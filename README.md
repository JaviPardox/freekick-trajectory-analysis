# Free kick trajectory analysis

This code simulates the motion of a ball starting at rest during a free kick.

In association football, or soccer, a free kick is a method of restarting a play.
It happens after an infringement of the rules by one of the players.

This code takes different parameters to create any type of free kick.

Parameters that can be changed:

	- Vo [x,y,z] - An array for the guess of initial velocities

	- Rf [x,y,z] - An array for the target position

	- time - The time it takes to reach the target position 

	- Wx, Wy and Wz - Angular velocities

### FK_t.py

This file simulates a standart free kick over a range of different times, from 1 second to 15 seconds.

Assumptions for a standart free kick: 

The player is right footed hence a positive spin
in the z direction, an average distance of 20 meters to the goal, with a wall of
players 9.1 meters away from the ball. The said spin is constant, drag and lift
coefficients also remain constant.

### FK_wx.py

This file simulates a standart free kick over a range of different possitive values of spin, only in the x direction.
Spin in the y and z direction are set to 0.

The values range from 0 rad/s to 4.5 rad/s.

### FK_wy.py

This file simulates a standart free kick over a range of different possitive values of spin, only in the y direction.
Spin in the x and z direction are set to 0.

The values range from 0 rad/s to 4.5 rad/s

### FK_wz.py

This file simulates a standart free kick over a range of different possitive values of spin, only in the z direction.
Spin in the y and x direction are set to 0.

The values range from 0 rad/s to 4.5 rad/s

## Historical free kicks

Using the model, animated simulation of two examples of sublime free kicks.

### MESSI_FK.py

Lionel Messi's free kick against Liverpool on May 2019

### ROBERTO_CARLOS_FK.py

Roberto Carlos' free kick against France on June 1997

