#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// Reference velocity is used in the cost function to punish any solution with
// a lower speed than the reference velocity. This could also be adapted dynamically,
// e.g. according to speed limits or lower it in curves
double ref_v = 60;

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

// Defines the starting point of each variable set in the state vector
uint x_start = 0;
uint y_start = x_start + N;
uint psi_start = y_start + N;
uint v_start = psi_start + N;
uint cte_start = v_start + N;
uint epsi_start = cte_start + N;
uint delta_start = epsi_start + N;
uint a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    size_t i;

    // Cost influence factors for fine-tuning the driving behavior
    double factor_cte = 1000.0;
    double factor_epsi = 1000.0;
    double factor_velocity = 10.0;
    double factor_delta_usage = 50.0;
    double factor_a_usage = 100.0;
    double factor_delta_delta = 100.0;
    double factor_a_delta = 20.0;

    // Define the cost function fg[0]
    // Note: This is taken from the solution to the MPC line fitting quiz in the Udacity course
    fg[0] = 0;
    // The part of the cost based on the reference state.
    for (i = 0; i < N; i++) {
      fg[0] += factor_cte * CppAD::pow(vars[cte_start + i], 2);
      fg[0] += factor_epsi * CppAD::pow(vars[epsi_start + i], 2);
      fg[0] += factor_velocity * CppAD::pow(vars[v_start + i] - ref_v, 2);
    }
    // Minimizes the usage of actuators, i.e. use them as little as possible
    for (i = 0; i < N - 1; i++) {
      fg[0] += factor_delta_usage * CppAD::pow(vars[delta_start + i], 2);
      fg[0] += factor_a_usage * CppAD::pow(vars[a_start + i], 2);
    }
    // Minimizes changes of actuators between two timesteps to ensure a smoother riding experience
    for (i = 0; i < N - 2; i++) {
      fg[0] += factor_delta_delta * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += factor_a_delta * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    // Take over all state variables from the vars vector
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];   

    // The remaining constraints are defined with the state vector transition from timestep t to t+1
    for (i = 1; i < N; i++) {
      // The state at time t+1 .
      AD<double> x1 = vars[x_start + i];
      AD<double> y1 = vars[y_start + i];
      AD<double> psi1 = vars[psi_start + i];
      AD<double> v1 = vars[v_start + i];
      AD<double> cte1 = vars[cte_start + i];
      AD<double> epsi1 = vars[epsi_start + i];

      // The state at time t.
      AD<double> x0 = vars[x_start + i - 1];
      AD<double> y0 = vars[y_start + i - 1];
      AD<double> psi0 = vars[psi_start + i - 1];
      AD<double> v0 = vars[v_start + i - 1];
      AD<double> cte0 = vars[cte_start + i - 1];
      AD<double> epsi0 = vars[epsi_start + i - 1];
      AD<double> delta0 = vars[delta_start + i - 1];
      AD<double> a0 = vars[a_start + i - 1];

      // Helper calculations
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0; // Evaluate the desired y position at x0 to calculate CTE
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2* coeffs[2] * x0 + 3 * coeffs[3] * x0 * x0); // Desired steering angle to calculate e_psi. Note derivative of second and third order components

      fg[1 + x_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + i] = psi1 - (psi0 - (v0 / Lf) * delta0 * dt);
      fg[1 + v_start + i] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + i] = epsi1 - ((psi0 - psides0) - (v0 / Lf) * delta0 * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9

  // 6 state variables (x,y,psi,v,cte,Epsi) and 2 actuators (delta, a)
  size_t n_vars = (6 * N) + (2 * (N-1));
  // TODO: Set the number of constraints
  // Each state variable is used as a constraint (relation between variable in timestep t+1 and t = 0)
  size_t n_constraints = 6 * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set initial state from given state vector
  vars[x_start] = state[0];
  vars[y_start] = state[1];
  vars[psi_start] = state[2];
  vars[v_start] = state[3];
  vars[cte_start] = state[4];
  vars[epsi_start] = state[5];

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // State variables are not really limited (of course in practice velocity e.g. is limited, but this will be handled by the optimizer)
  // Therefore all variables are set to the numeric minimum and maximum limits of the double variable type
  for (i=0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
    //vars_lowerbound[i] = std::numeric_limits<double>::min();
    //vars_upperbound[i] = std::numeric_limits<double>::max();
  }

  // Steering angle is limited to -25 and +25 degrees
  for (i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332; // - 25° in radians
    vars_upperbound[i] = 0.436332; // + 25° in radians
  }

  // acceleration is limited to -1 to +1. Can be higher values, but may result in erratic driving behavior
  for (i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Set the initial state to exactly the given state (x=x0, ...)
  constraints_lowerbound[x_start] = state[0];
  constraints_lowerbound[y_start] = state[1];
  constraints_lowerbound[psi_start] = state[2];
  constraints_lowerbound[v_start] = state[3];
  constraints_lowerbound[cte_start] = state[4];
  constraints_lowerbound[epsi_start] = state[5];

  constraints_upperbound[x_start] = state[0];
  constraints_upperbound[y_start] = state[1];
  constraints_upperbound[psi_start] = state[2];
  constraints_upperbound[v_start] = state[3];
  constraints_upperbound[cte_start] = state[4];
  constraints_upperbound[epsi_start] = state[5];

  // object that computes objective and constraints
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

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  // Solutions vector consists of all (x,y) prediction pairs, 
  // followed by the one delta and a value for the next actuation
  std::vector<double> solutions;

  for (i=0; i< N; i++){
    solutions.push_back(solution.x[x_start + i]);
    solutions.push_back(solution.x[y_start + i]);
  }

  solutions.push_back(solution.x[delta_start]);
  solutions.push_back(solution.x[a_start]);

  return solutions;
}
