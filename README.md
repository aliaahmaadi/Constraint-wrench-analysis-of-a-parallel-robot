# This repository contains MATLAB scripts for simulating the inverse dynamics of a parallel robot using screw motion theory. The analysis focuses on constraint wrench evaluation and trajectory planning.

Features
Inverse Displacement Analysis: Computes the joint angles for given end-effector positions and orientations.
Trajectory Simulation: Implements predefined trajectories to test robot dynamics.
Constraint Wrench Evaluation: Uses screw theory to analyze forces and moments in the system.
Visualization: Plots of joint displacements, velocities, accelerations, and test trajectories.
Files and Their Functions
Finding_home.m
Determines the home posture and geometry parameters of the parallel robot.

Inverse_displacement.m
Solves the inverse displacement problem for the robot.

Inverse_displacement_method2.m
Alternative method for inverse displacement calculation.

Inverse_Displacement_with_trajectory.m
Simulates joint angles over time for a predefined trajectory.

Inverse_dynamics.m
Simulates the inverse dynamics using constraint wrench analysis.

CPM.m
Computes the cross-product matrix for a given vector.

Requirements
MATLAB (R2020b or later recommended)
Basic knowledge of robot kinematics and dynamics
Usage
Clone the repository:
bash
Copy code
git clone https://github.com/your-username/inverse-dynamics-parallel-robot.git
Open MATLAB and navigate to the cloned repository folder.
Run the desired script (.m file) to perform the analysis.
Output
Numerical results for joint angles, velocities, and accelerations.
Graphs for displacement, velocity, and acceleration over time.
Plots of the robot's trajectory in the workspace.
