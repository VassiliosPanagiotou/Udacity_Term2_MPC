#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //std::cout << sdata << endl; // not needed
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
          double delta0 = - double(j[1]["steering_angle"]) * deg2rad(25);
          double a0     = j[1]["throttle"];

          // Reference trajectory (in vehicle coordinate system)
          size_t dim = ptsy.size();
		  if (ptsx.size() < ptsy.size()) ptsx.size();
		  
          Eigen::VectorXd xvals(dim), yvals(dim);
		  
		  //
          for (size_t i=0; i<dim; i++)
          {
            xvals[i] = cos(-psi) * (ptsx[i] - px) - sin(-psi) * (ptsy[i] - py);
            yvals[i] = sin(-psi) * (ptsx[i] - px) + cos(-psi) * (ptsy[i] - py);
          }
		  
		  // Fit a third degree polynomial
          auto coeffs = polyfit(xvals, yvals, 3);

          // MPC in vehicle coordinate system
          const double latency = 0.1;
          double yawRate = v * delta0 / MPC::Lf;
          Eigen::VectorXd state(6);
		  
          if ( fabs(yawRate) < 0.0001 )
          {
            state[0] = v * latency;
            state[1] = 0.0;
            state[2] = 0.0;
          }
          else
          {
            state[0] = v / yawRate * ( sin(0.0 + yawRate * latency) - sin(0.0) );
            state[1] = v / yawRate * ( cos(0.0)                     - cos(0.0 + yawRate * latency) );
            state[2] = yawRate * latency;
          }
          state[3] = v + latency * a0;
          state[4] = polyeval(coeffs, state[0]) - state[1];
          state[5] = state[2] - atan(coeffs[1] + 2.0 * coeffs[2] * state[0] + 3.0 * coeffs[3] * state[0] * state[0]);
          double delta, a;

		  // Calculate MPC predicted trajectory 
		  vector<double> mpc_x_vals;
		  vector<double> mpc_y_vals;

          std::cout << "INPUT:     {x, y, psi, v, cte, epsi} = " << std::endl
                    << "           {"
                    << state[0] << ", " << state[1] << ", " << state[2] << ", " << state[3] << ", " << state[4] << "}" << std::endl;

          mpc.Solve(state, coeffs, delta, a, mpc_x_vals, mpc_y_vals);

          /*
          * Calculate steering angle and throttle using MPC.
          * Both are in between [-1, 1].
          */
          double steer_value    = -delta / deg2rad(25);
          double throttle_value = a;
		  if (steer_value < -1.0) steer_value = -1.0;
		  if (steer_value >  1.0) steer_value = 1.0;
		  if (throttle_value < -1.0) throttle_value = -1.0;
		  if (throttle_value >  1.0) throttle_value = 1.0;

		  std::cout << "OUTPUT:    {steer, a}                = " << std::endl
			  << "           {"
			  << steer_value << ", " << a << "}" << std::endl;
		  std::cout << "MPC:       {(x, y) -> (x, y)} = " << std::endl
			  << "           {("
			  << mpc_x_vals[0] << ", " << mpc_y_vals[0] << ") -> ("
			  << mpc_x_vals[1] << ", " << mpc_y_vals[1] << ") }" << std::endl << std::endl;

          json msgJson;
          
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;
		  //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          next_x_vals.resize(mpc_x_vals.size());
          next_y_vals.resize(mpc_y_vals.size());
		  
          for (size_t i=0; i<mpc_x_vals.size(); i++)
          {
            next_x_vals[i] = mpc_x_vals[i];
            next_y_vals[i] = polyeval(coeffs, next_x_vals[i]);
          }

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl; // not needed
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
