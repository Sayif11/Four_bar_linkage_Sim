# Four_bar_linkage_Sim

## Overview

Simulates the motion of a planar four-bar linkage consisting of a fixed ground link, an input crank, a coupler, and an output rocker. The mechanism configuration is obtained from the loop-closure equations, while the crank motion is computed by integrating the equations of motion using MATLAB's `ode45` solver.

* Geometric constraint equations determine the rocker angle from the crank angle.
* Numerical differentiation is used to evaluate the first and second derivatives of the rocker motion.
* The crank dynamics are integrated in time using the specified inertia and applied torque.
* Link positions are updated and animated at each timestep.

## Outputs

* Input crank angular velocity as a function of time.
* Animated visualization of the four-bar linkage motion.
* Time-varying positions of all links and joints.

## Assumptions

* Planar rigid-body motion.
* Ideal pin joints.
* Constant link lengths.
* No friction or damping.
* Ground link aligned with the global (y)-axis.
