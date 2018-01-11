# Model Predictive Control Vehicle controller 
This project is part of the Udacity Self-Driving Car Engineer Nanodegree Program. The project base code can be found under [https://github.com/udacity/CarND-MPC-Project](https://github.com/udacity/CarND-MPC-Project). The program interacts with a simulator that can be downloaded at [https://github.com/udacity/self-driving-car-sim/releases](https://github.com/udacity/self-driving-car-sim/releases) -> Term 2 simulator.

---

## Project outline
In this project, the task was to implement a Controller based on a more complex vehicle model. The model includes a state with 6 variables: $(x,y)$ position of the vehicle, car heading $\psi$, velocity $v$, current cross track error $CTE$ and car heading error $e_{\psi}$. Additionally, two actuator variables steering angle $\delta$ and acceleration $a$ are used.

Solving the model is an optimization task in which a cost function is minimized under given constraints. The constraints are given as the current state of the vehicle in timestep $t=0$ as well as state transition functions from timestep $t$ to $t+1$. These are given by the equations
$$ x_{t+1} = x_t + v * cos(\psi) * dt $$
$$ y_{t+1} = y_t + v * sin(\psi) * dt $$
$$ \psi_{t+1} = \psi_t - v * dt * \delta / L_f $$
$$ v_{t+1} = v_t + a * dt $$

with $L_f$ a constant that is derived from the real steering behavior of the car (steering angle translates into a steering circle)


The cross track error at time 0 is given by the difference of the vehicle's position (defined in vehicle coordinates as 0) and the current waypoint path, which is the evaluation of the polynomial fitting at position 0, or coeffs[0]. 

The model evaluates N timesteps of each dt seconds into the future (total prediction timeframe: N*dt) to minimize the cost function that as best as possible matches the given waypoint. This ensures that, unlike a PID controller, not only the next command is optimized, but the next command is executed in a way so as to minimize the error N timesteps in the future. Because errors accumulate, I have chosen a relatively short period of 10 timesteps of each 100ms, or a total of 1s prediction into the future.

Regarding the waypoints, the points are given in map coordinates, but need to be handed to the simulator in vehicle coordinates to be displayed correctly, which also greatly simplifies the model calculations, because the CTE calculation from the vehicle perspective is much simpler: it only involves the deviation in $y$ direction, since the $x$ axis is right ahead of the vehicle. Therefore, I am transforming the waypoints into vehicle coordinates first.

In the time of the latency, the vehicle travels with velocity v further in x direction. Since $\psi$ is by definition $0$ (in vehicle coordinates), we remain on the imagined trajectory of the car, therefore $x = v * latency$, otherwise additional multiplication with $cos(\psi)$ is required. With the same logic, y remains zero. psi changes slightly over time according to the vehicle model equation $(\psi = v * latency * \frac{-steer\_value}{L_f}$. Steer value needs to be negative because in the simulator, a positive value implies a right turn, whereas the mathematial definition is the other way around. The velocity changes with the acceleration, $v = v + a * latency$. Finally, CTE needs to be slightly adapted if at timepoint $t=0$ the heading error $e_\psi$ is non-zero, which causes CTE to grow by $ v * sin( \e_\psi) * latency $.


## Dependencies and basic build instructions

See the Udacity starter code page [https://github.com/udacity/CarND-MPC-Project#dependencies](https://github.com/udacity/CarND-MPC-Project#dependencies)


