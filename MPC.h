#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {

public:
	//C'tor - D'tor
	MPC();
	virtual ~MPC();

	// Solve the model given an initial state and polynomial coefficients.
	//@param  state          initial state {x, y, psi, v, cte, epsi}
	//@param  coeffs         polynomial coefficients of ref trajectory
	//@param  delta          steering angle
	//@param  a              acceleration
	//@param  trajectory_x   x component of the trajectory
	//@param  trajectory_y   y component of the trajectory
	bool Solve(const Eigen::VectorXd& state, const Eigen::VectorXd& coeffs,	double& delta, double& a, vector<double>& trajectory_x, vector<double>& trajectory_y);

	// Timestep length and duration
	const size_t N;
	const double dt;
	// Length from front to CoG
	const double Lf;
	// Reference velocity
	const double ref_v;
};

#endif /* MPC_H */
