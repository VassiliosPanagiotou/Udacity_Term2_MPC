#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

//////////
// Init static const variables
//////////
// Set the timestep length and duration
const size_t N = 10;
const double dt = 0.1;
// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
// Reference velocity set to 40 mph
 const double ref_v = 40.0;

// State variables and actuator
size_t x_start     = 0;
size_t y_start     = x_start     + N;
size_t psi_start   = y_start     + N;
size_t v_start     = psi_start   + N;
size_t cte_start   = v_start     + N;
size_t epsi_start  = cte_start   + N;
size_t delta_start = epsi_start  + N;
size_t a_start     = delta_start + N - 1;

// Definition of class that computes objective and constraints
class FG_eval
{
public:
  // Definition of double vector in CppAD
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Implement MPC
    // Cost is the first element of `fg`.
    fg[0] = 0;
    // Cost based on the reference state.
    for (size_t t = 0; t < MPC::N; t++)
    {
      fg[0] += 1.0  * CppAD::pow(vars[cte_start  + t], 2);
      fg[0] += 1.0  * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += 0.01 * CppAD::pow(vars[v_start    + t] - ref_v, 2);
    }
    // Minimize the use of actuators.
    for (size_t t = 0; t < MPC::N - 1; t++)
    {
      fg[0] += 1.0  * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 0.25 * CppAD::pow(vars[a_start     + t], 2);
    }
    // Minimize the value gap between sequential actuations.
    for (size_t t = 0; t < N - 2; t++)
    {
      fg[0] += 4.0 * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += 1.0 * CppAD::pow(vars[a_start     + t + 1] - vars[a_start     + t], 2);
    }

    // Initial constraints
    fg[x_start + 1]    = vars[x_start];
    fg[y_start + 1]    = vars[y_start];
    fg[psi_start + 1]  = vars[psi_start];
    fg[v_start + 1]    = vars[v_start];
    fg[cte_start + 1]  = vars[cte_start];
    fg[epsi_start + 1] = vars[epsi_start];
    // Constraints
    for (size_t t = 1; t < MPC::N; t++)
    {
      // State at time t+1 .
      AD<double> x1    = vars[x_start    + t];
      AD<double> y1    = vars[y_start    + t];
      AD<double> psi1  = vars[psi_start  + t];
      AD<double> v1    = vars[v_start    + t];
      AD<double> cte1  = vars[cte_start  + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // State at time t.
      AD<double> x0    = vars[x_start    + t - 1];
      AD<double> y0    = vars[y_start    + t - 1];
      AD<double> psi0  = vars[psi_start  + t - 1];
      AD<double> v0    = vars[v_start    + t - 1];
      AD<double> cte0  = vars[cte_start  + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0     = vars[a_start     + t - 1];

      // Yaw rate at time t
      AD<double> yawRate0 = v0 * delta0 / MPC::Lf;

      // desired y0 value according to polynomial coefficients
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      // desired angle value according to polynomial coefficients
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2.0 * coeffs[2] * x0 + 3.0 * coeffs[3] * CppAD::pow(x0, 2));

      // State transitions
      if ( fabs(yawRate0) < 0.0001 )
      {
        fg[1 + x_start    + t] = x1 -  (x0               + v0 * CppAD::cos(psi0) * dt);
        fg[1 + y_start    + t] = y1 -  (y0               + v0 * CppAD::sin(psi0) * dt);
      }
      else
      {
        fg[1 + x_start    + t] = x1 -  (x0               + (v0 / yawRate0) * (CppAD::sin(psi0 + yawRate0 * MPC::dt) - CppAD::sin(psi0)) );
        fg[1 + y_start    + t] = y1 -  (y0               + (v0 / yawRate0) * (CppAD::cos(psi0)                      - CppAD::cos(psi0 + yawRate0 * MPC::dt)) );
      }
      fg[1 + psi_start  + t] = psi1 -  (psi0             + v0 * delta0 / MPC::Lf * MPC::dt);
      fg[1 + v_start    + t] = v1 -    (v0               + a0 * MPC::dt);
      fg[1 + cte_start  + t] = cte1 -  ((f0 - y0)        + (v0 * CppAD::sin(epsi0) * MPC::dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / MPC::Lf * MPC::dt);
    }
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

bool MPC::Solve(const Eigen::VectorXd& state, const Eigen::VectorXd& coeffs,	double& delta, double& a, vector<double>& trajectory_x, vector<double>& trajectory_y) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Set the number of model variables (includes both states and inputs).
  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2;
  // Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  // Initial variable values
  vars[x_start]    = x;
  vars[y_start]    = y;
  vars[psi_start]  = psi;
  vars[v_start]    = v;
  vars[cte_start]  = cte;
  vars[epsi_start] = epsi;
  
  // Lower and upper limits for state variables
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // Set lower and upper limits for variables.
  for (size_t i = 0; i < delta_start; i++)
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  for (size_t i = delta_start; i < a_start; i++)
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  for (size_t i = a_start; i < n_vars; i++)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  
  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  // Initial lower bounds
  constraints_lowerbound[x_start]    = x;
  constraints_lowerbound[y_start]    = y;
  constraints_lowerbound[psi_start]  = psi;
  constraints_lowerbound[v_start]    = v;
  constraints_lowerbound[cte_start]  = cte;
  constraints_lowerbound[epsi_start] = epsi;
  // Initial upper bounds
  constraints_upperbound[x_start]    = x;
  constraints_upperbound[y_start]    = y;
  constraints_upperbound[psi_start]  = psi;
  constraints_upperbound[v_start]    = v;
  constraints_upperbound[cte_start]  = cte;
  constraints_upperbound[epsi_start] = epsi;
  
  // Compute objective and constraints
  FG_eval fg_eval(coeffs);
  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  if ( ok )
  {
    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;
  }
  else
  {
    // Error message
    std::cout << "MPC::solve - Error occured!" << std::endl;
  }

  // Set the output values
  delta = solution.x[delta_start];
  a     = solution.x[a_start];
  trajectory_x.resize(N);
  trajectory_y.resize(N);
  for (size_t i=0; i<N; i++)
  {
    trajectory_x[i] = solution.x[x_start + i];
    trajectory_y[i] = solution.x[y_start + i];
  }
  return ok;
}
