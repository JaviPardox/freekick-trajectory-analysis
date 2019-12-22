# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 03:38:04 2019

@author: Javier Pardo
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scp
from scipy import optimize
import matplotlib.animation as animation
from mpl_toolkits import mplot3d

#DRAG AND LIFT COEFFICIENTS
Cd = 0.25
Cl = 0.28

#MASS AND RADIUS OF THE BALL
m = 0.4
r = 0.1095

#ANGULAR VELOCITIES
Wx = 0
Wz = -1
Wy = 0

#CONSTANTS
g = 9.81
p = 1.2
A = 0.0375
B = (p*A)/2
I = (2/3)*m*r**2 #moment of inertia

plt.close('all')

#RIGHT HAND SIDE EQUATIONS
def rhs(t,state):
      
    x = state[0]
    Vx = state[1]
    
    y = state[2]
    Vy = state[3]
    
    z = state[4]
    Vz = state[5]
    
    fx = Vx
    fVx = B*np.sqrt(Vx**2+Vy**2+Vz**2)*(-(Cd/m)*Vx + (Cl/m)*(Wy*Vz-Wz*Vy))
    
    fy = Vy
    fVy = B*np.sqrt(Vx**2+Vy**2+Vz**2)*(-(Cd/m)*Vy + (Cl/m)*(Wx*Vz-Wz*Vx))
    
    fz = Vz
    fVz = B*np.sqrt(Vx**2+Vy**2+Vz**2)*(-(Cd/m)*Vz + (Cl/m)*(Wx*Vy-Wy*Vx)) - g
    
    return np.array([fx,fVx,fy,fVy,fz,fVz],float)


#TIME
#initial, final, number of points
time = np.linspace(0, 1.0, 100)

#INITIAL CONDITIONS
Vo = np.array([10,5,3]) #initial velocities
Ro = np.array([0,0,0]) #initial positions

#ROBERTO CARLOS TARGET
Rf = np.array([35,-3.66,1.5]) 

def dist_to_target(Vo):
      
    init = np.array([Ro[0],Vo[0],Ro[1],Vo[1],Ro[2],Vo[2]])
    sol = scp.solve_ivp(rhs,[0,10],init, t_eval=time)
    
    x = sol['y'][0,:]
    y = sol['y'][2,:]
    z = sol['y'][4,:]
    
    dist_x = x[-1] - Rf[0]
    dist_y = y[-1] - Rf[1]
    dist_z = z[-1] - Rf[2]
    
    return np.array([dist_x,dist_y,dist_z],float)

best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)
w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)
r_magnitude = np.sqrt(Rf[0]**2 + Rf[1]**2 + Rf[2]**2)
    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")

print("Velocity magnitude: ",v_magnitude,"m/s")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")

#using the velocities to find the position
init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

data = np.array([x,y,z])
N = 300

def update(num, data, line):
   line.set_data(data[:2, :num])
   line.set_3d_properties(data[2, :num])

#ROBERTO CARLOS FREE KICK
fig, cx = plt.subplots()

cx = plt.axes(projection='3d')
line, = cx.plot3D(x,y,z,'-r',label='soccer ball')

ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)

#SAVE ANIMATION
#Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#ani.save('Roberto_carlos.mp4', writer=writer)


#GOAL PLACEMENT
goal_z_p = np.linspace(2.44,2.44,50)
goal_x_p = np.linspace(Rf[0],Rf[0],50)
goal_y_p = np.linspace(Rf[1],-Rf[1],50)

cx.plot3D(goal_x_p,goal_y_p,goal_z_p,'-b',label='goal')

goal_z_pr = np.linspace(0,2.44,50)
goal_x_pr = np.linspace(Rf[0],Rf[0],50)
goal_y_pr = np.linspace(Rf[1],Rf[1],50)

cx.plot3D(goal_x_pr,goal_y_pr,goal_z_pr,'-b')

goal_z_pl = np.linspace(0,2.44,50)
goal_x_pl = np.linspace(Rf[0],Rf[0],50)
goal_y_pl = np.linspace(-Rf[1],-Rf[1],50)

cx.plot3D(goal_x_pl,goal_y_pl,goal_z_pl,'-b')

#WALL OF PLAYERS PLACEMENT
wall_z_top = np.linspace(2.0,2.0,50)
wall_x_top = np.linspace(9.1,9.1,50)
wall_y_top = np.linspace(0,-1.5,50)

cx.plot3D(wall_x_top,wall_y_top,wall_z_top,'-g',label='wall of players')

wall_z_r = np.linspace(2.0,0.0,50)
wall_x_r = np.linspace(9.1,9.1,50)
wall_y_r = np.linspace(-1.5,-1.5,50)

cx.plot3D(wall_x_r,wall_y_r,wall_z_r,'-g')

wall_z_l = np.linspace(2.0,0.0,50)
wall_x_l = np.linspace(9.1,9.1,50)
wall_y_l = np.linspace(0.0,0.0,50)

cx.plot3D(wall_x_l,wall_y_l,wall_z_l,'-g')

cx.set_title('Roberto Carlos Free Kick')

cx.set_xlabel('x (m)')
cx.set_ylabel('y (m)')
cx.set_zlabel('z (m)')

cx.set_xlim(0,35)
cx.set_ylim(-10,10)
cx.set_zlim(0,15)

cx.legend()
