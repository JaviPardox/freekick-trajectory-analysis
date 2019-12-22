# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 02:57:48 2019

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
Wz = 0
Wy = 0

#CONSTANTS
g = 9.81
p = 1.2
A = 0.0375
B = (p*A)/2
I = (2/3)*m*r**2 #moment of inertia

#LISTS
Vz = []
Vy = []
Vx = []
V = []
W = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5]

#plt.close('all')


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


#time: initial, final, number of points
time = np.linspace(0, 2.0, 200)

#INITIAL CONDITIONS
Vo = np.array([10,5,3]) #initial velocities
Ro = np.array([0,0,0]) #initial positions

#TARGET
Rf = np.array([20,3.66,2.4]) 

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
V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)
r_magnitude = np.sqrt(Rf[0]**2 + Rf[1]**2 + Rf[2]**2)
    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter
print("0.0:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (1s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

#using the velocities to find the position
init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

#UNCOMMENT SECTION FOR ANIMATION (CPU HEAVY)
#data = np.array([x,y,z])
#N = 500

#def update(num, data, line):
#   line.set_data(data[:2, :num])
#   line.set_3d_properties(data[2, :num])

#1S FREE KICK
fig, cx = plt.subplots()
cx = plt.axes(projection='3d')
line, = cx.plot3D(x,y,z,'-r',label='0.0')

fig, ax = plt.subplots()
ax.plot(x,z,label='0.0')

fig, dx = plt.subplots()
dx.plot(y,z,label='0.0')

#UNCOMMENT FOR ANIMATION
#ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)

  ################################
####### 0.5 Wx FREE KICK #########
  ################################

Wx = 0.5
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("0.5:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (0.5): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,'-b',label='0.5')
ax.plot(x,z,label='0.5')
dx.plot(y,z,label='0.5')

  ################################
####### 1.0 Wx FREE KICK #########
  ################################

Wx = 1.0
Vo = np.array([5,3,1])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]
v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("1.0:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (1.0): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,'-g',label='1.0')
ax.plot(x,z,label='1.0')
dx.plot(y,z,label='1.0')

  ################################
####### 1.5 Wx FREE KICK #########
  ################################

Wx = 1.5
Vo = np.array([2,1,1])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("1.5:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (1.5): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,color='orange',label='1.5')
ax.plot(x,z,label='1.5')
dx.plot(y,z,label='1.5')

  ################################
####### 2.0 Wx FREE KICK #########
  ################################

Wx = 2.0
Vo = np.array([2,1,1])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("2.0:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (2.0): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,color='violet',label='2.0')
ax.plot(x,z,label='2.0')
dx.plot(y,z,label='2.0')

  ################################
####### 2.5 Wx FREE KICK #########
  ################################

Wx = 2.5
Vo = np.array([1,1,1])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("2.5:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (2.5): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,color='grey',label='2.5')
ax.plot(x,z,label='2.5')
dx.plot(y,z,label='2.5')

  ################################
####### 3.0 Wx FREE KICK #########
  ################################

Wx = 3.0
Vo = np.array([0.01,0.01,0.01])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("3.0:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (3.0): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='3.0')
ax.plot(x,z,label='3.0')
dx.plot(y,z,label='3.0')

  ################################
####### 3.5 Wx FREE KICK #########
  ################################

Wx = 3.5
Vo = np.array([0.01,0.01,0.01])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("3.5")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (3.5): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='3.5')
ax.plot(x,z,label='3.5')
dx.plot(y,z,label='3.5')

  ################################
####### 4 Wx FREE KICK #########
  ################################

Wx = 4.0
Vo = np.array([0.01,0.01,0.01])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("4:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (4): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='4')
ax.plot(x,z,label='4')
dx.plot(y,z,label='4')

  ################################
####### 4.5 Wx FREE KICK #########
  ################################

Wx = 4.5
Vo = np.array([0.01,0.01,0.01])
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

print("4.5:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (4.5): ",v_magnitude,"m/s")
print("")


init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[0,10],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,color='black',label='4.5')
ax.plot(x,z,label='4.5')
dx.plot(y,z,label='4.5')

# 3D PLOT'S PARAMETERS
cx.set_title('Average Free Kick - Trajectories changing $\omega_x$')

cx.set_xlabel('x (m)')
cx.set_ylabel('y (m)')
cx.set_zlabel('z (m)')

#cx.set_xlim(0,25)
#cx.set_ylim(-10,10)
#cx.set_zlim(0,15)

cx.legend()

#V vs T plots
fig, bx = plt.subplots()
bx.plot(W,V,label='V')
bx.plot(W,Vx,label='Vx')
bx.plot(W,Vy,label='Vy')
bx.plot(W,Vz,label='Vz')

bx.set_title('Velocities vs Different $\omega_x$')

bx.set_xlabel('$\omega_x$ (rad/s)')
bx.set_ylabel('velocity (m/s)')
bx.legend()

#X vs Z plot
ax.set_title('X vs Z')
ax.set_xlabel('x (m)')
ax.set_ylabel('z (m)')
ax.legend()

#Y vs Z plot
dx.set_title('Y vs Z')
dx.set_xlabel('Y (m)')
dx.set_ylabel('z (m)')
dx.legend()