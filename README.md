# This repository contains MATLAB scripts for simulating the inverse dynamics of a parallel robot using screw motion theory. The analysis focuses on constraint wrench evaluation and trajectory planning.

# Inverse Dynamics Simulation of a Parallel Robot

This repository contains MATLAB scripts for simulating the inverse dynamics of a parallel robot using screw motion theory. The analysis focuses on constraint wrench evaluation and trajectory planning.

## Features

- **Inverse Displacement Analysis:**  
  Computes the joint angles for given end-effector positions and orientations.
  
- **Trajectory Simulation:**  
  Implements predefined trajectories to test robot dynamics.
  
- **Constraint Wrench Evaluation:**  
  Uses screw theory to analyze forces and moments in the system.
  
- **Visualization:**  
  Plots of joint displacements, velocities, accelerations, and test trajectories.

## Files and Their Functions

1. **`Finding_home.m`**  
   Determines the home posture and geometry parameters of the parallel robot.

2. **`Inverse_displacement.m`**  
   Solves the inverse displacement problem for the robot.

3. **`Inverse_displacement_method2.m`**  
   Alternative method for inverse displacement calculation.

4. **`Inverse_Displacement_with_trajectory.m`**  
   Simulates joint angles over time for a predefined trajectory.

5. **`Inverse_dynamics.m`**  
   Simulates the inverse dynamics using constraint wrench analysis.

6. **`CPM.m`**  
   Computes the cross-product matrix for a given vector.

## Requirements

- **MATLAB**: R2020b or later (recommended).  
- **Knowledge**: Basic understanding of robot kinematics and dynamics.

## Usage

1. Clone the repository:  
   ```bash
   git clone https://github.com/your-username/inverse-dynamics-parallel-robot.git

