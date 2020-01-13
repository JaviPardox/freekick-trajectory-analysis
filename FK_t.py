# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 19:40:51 2019

@author: Javier Pardo 
         https://www.linkedin.com/in/javier-pardo-fernandez-87b565124/
         javiyupipa@gmail.com
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
Wz = -2
Wy = 0

#CONSTANTS
g = 9.81
p = 1.2
A = 0.0375
B = (p*A)/2
I = (2/3)*m*r**2 #moment of inertia

#LISTS for graphing purposes
Vz = []
Vy = []
Vx = []
V = []
E = []

#TIMES for graphing purposes
#number of elements = number of shots
#element value = time for the shot
T = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15] 

plt.close('all')

#RIGHT HAND SIDE EQUATIONS
#returns position and velocities at all times
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


Ro = np.array([0,0,0]) #initial positions
Rf = np.array([20,3.66,2.4]) #final position

 
def dist_to_target(Vo):
      
    init = np.array([Ro[0],Vo[0],Ro[1],Vo[1],Ro[2],Vo[2]])
    sol = scp.solve_ivp(rhs,[time[0],time[-1]],init, t_eval=time) #force the amount of points
    
    x = sol['y'][0,:]
    y = sol['y'][2,:]
    z = sol['y'][4,:]
    
    #shooting method
    #when dist_r = 0 -> arrived to target
    dist_x = x[-1] - Rf[0]
    dist_y = y[-1] - Rf[1]
    dist_z = z[-1] - Rf[2]
    
    return np.array([dist_x,dist_y,dist_z],float)

  ################################
####### 1 SECOND FREE KICK #########
  ################################
  
#TIME
#(initial, final, number of points)
time = np.linspace(0, 1.0, 200)

#INITIAL CONDITIONS
Vo = np.array([10,5,3]) #initial velocities

#returns the initial velocities needed to reach the desired target
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("1 SECOND:")
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

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

#UNCOMMENT FOR ANIMATION (CPU HEAVY)
#data = np.array([x,y,z])
#N = 500  #adjust speed of animation

#def update(num, data, line):
#   line.set_data(data[:2, :num])
#   line.set_3d_properties(data[2, :num])

cx = plt.axes(projection='3d')
line, = cx.plot3D(x,y,z,label='1s')

fig, ax = plt.subplots()
ax.plot(x,z,label='1s')

fig, dx = plt.subplots()
dx.plot(y,z,label='1s')

fig, fx = plt.subplots()
fx.plot(x,y,label='1s')
#UNCOMMENT FOR ANIMATION (CPU HEAVY)
#ani = animation.FuncAnimation(fig, update, N, fargs=(data, line), interval=10000/N, blit=False)

  ################################
####### 2 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 2.0, 200)
best_init_Vo = optimize.root(dist_to_target, Vo)

Vox = best_init_Vo['x'][0]
Voy = best_init_Vo['x'][1]
Voz = best_init_Vo['x'][2]

v_magnitude = np.sqrt(Vox**2 + Voy**2 + Voz**2)
w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)
r_magnitude = np.sqrt(Rf[0]**2 + Rf[1]**2 + Rf[2]**2)
    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

V.append(v_magnitude)
Vx.append(Vox)
Vy.append(Voy)
Vz.append(Voz)

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("2 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (2s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='2s')
ax.plot(x,z,label='2s')
dx.plot(y,z,label='2s')
fx.plot(x,y,label='2s')

  ################################
####### 3 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 3.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("3 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (3s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='3s')
ax.plot(x,z,label='3s')
dx.plot(y,z,label='3s')
fx.plot(x,y,label='3s')

  ################################
####### 4 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 4.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("4 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (4s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='4s')
ax.plot(x,z,label='4s')
dx.plot(y,z,label='4s')
fx.plot(x,y,label='4s')

  ################################
####### 5 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 5.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("5 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (5s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='5s')
ax.plot(x,z,label='5s')
dx.plot(y,z,label='5s')
fx.plot(x,y,label='5s')

  ################################
####### 6 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 6.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("6 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (6s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='6s')
ax.plot(x,z,label='6s')
dx.plot(y,z,label='6s')
fx.plot(x,y,label='6s')

  ################################
####### 7 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 7.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("7 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (7s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='7s')
ax.plot(x,z,label='7s')
dx.plot(y,z,label='7s')
fx.plot(x,y,label='7s')

  ################################
####### 8 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 8.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("8 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (8s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='8s')
ax.plot(x,z,label='8s')
dx.plot(y,z,label='8s')
fx.plot(x,y,label='8s')

  ################################
####### 9 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 9.0, 200)
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

w_magnitude = np.sqrt(Wx**2 + Wy**2 + Wz**2)    
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("9 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (9s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='9s')
ax.plot(x,z,label='9s')
dx.plot(y,z,label='9s')
fx.plot(x,y,label='9s')

  ################################
####### 10 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 10.0, 500)
Vo = np.array([0.0225,0.0225,0.0225])
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("10 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (10s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='10s')
ax.plot(x,z,label='10s')
dx.plot(y,z,label='10s')
fx.plot(x,y,label='10s')

  ################################
####### 11 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 11.0, 500)
Vo = np.array([0.0225,0.0225,0.0225])
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("11 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (11s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='11s')
ax.plot(x,z,label='11s')
dx.plot(y,z,label='11s')
fx.plot(x,y,label='11s')

  ################################
####### 12 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 12.0, 500)
Vo = np.array([0.0225,0.0225,0.0225])
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("12 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (12s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='12s')
ax.plot(x,z,label='12s')
dx.plot(y,z,label='12s')
fx.plot(x,y,label='12s')

  ################################
####### 13 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 13.0, 500)
Vo = np.array([0.0225,0.0225,0.0225])
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("13 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (13s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='13s')
ax.plot(x,z,label='13s')
dx.plot(y,z,label='13s')
fx.plot(x,y,label='13s')

  ################################
####### 14 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 14.0, 500)
Vo = np.array([0.0225,0.0225,0.0225])
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("14 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (14s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='14s')
ax.plot(x,z,label='14s')
dx.plot(y,z,label='14s')
fx.plot(x,y,label='14s')

  ################################
####### 15 SECONDS FREE KICK #########
  ################################

time = np.linspace(0, 15.0, 500)
Vo = np.array([0.0225,0.0225,0.0225])
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
Energy = (1/2)*m*(v_magnitude**2) + (1/2)*I*(w_magnitude**2)
Work_equivalent_on_mass = Energy/g*1  #mass moved against gravity for 1 meter

E.append(Energy)

print("15 SECONDS:")
print("Velocity necessary x-direction: ",Vox,"m/s")
print("Velocity necessary y-direction: ",Voy,"m/s")
print("Velocity necessary z-direction: ",Voz,"m/s")
print("Velocity magnitude (15s): ",v_magnitude,"m/s")
print("")
print("Energy used: ",Energy,"J")
print("Is equivalent to lift a mass of ",Work_equivalent_on_mass,"kg for a meter")
print("")

init_graph = np.array([Ro[0],Vox,Ro[1],Voy,Ro[2],Voz])

sol = scp.solve_ivp(rhs,[time[0],time[-1]],init_graph, t_eval=time) #force the amount of points
x = sol['y'][0,:]
y = sol['y'][2,:]
z = sol['y'][4,:]

line, = cx.plot3D(x,y,z,label='15s')
ax.plot(x,z,label='15s')
dx.plot(y,z,label='15s')
fx.plot(x,y,label='15s')

# 3D PLOT'S PARAMETERS
cx.set_title('Average Free Kick - Trajectories')

cx.set_xlabel('x (m)')
cx.set_ylabel('y (m)')
cx.set_zlabel('z (m)')

#set limits
#cx.set_xlim(0,25)
#cx.set_ylim(-10,10)
#cx.set_zlim(0,15)

cx.legend()

#V vs T plots
fig, bx = plt.subplots()
bx.plot(T,V,label='V')
bx.plot(T,Vx,label='Vx')
bx.plot(T,Vy,label='Vy')
bx.plot(T,Vz,label='Vz')

bx.set_title('Velocities vs Time of ball in the air')

bx.set_xlabel('time (s)')
bx.set_ylabel('velocity (m/s)')
bx.legend()

#X vs Z plot
ax.set_title('X vs Z')
ax.set_xlabel('x (m)')
ax.set_ylabel('z (m)')
ax.legend()

#Y vs Z plot
dx.set_title('Y vs Z')
dx.set_xlabel('y (m)')
dx.set_ylabel('z (m)')
dx.legend()

#X vs Y plot
fx.set_title('X vs Y')
fx.set_xlabel('x (m)')
fx.set_ylabel('y (m)')
fx.legend()

#Energy vs times plot
fig, ex = plt.subplots()
ex.plot(T,E)
ex.set_title('Energy vs Times')
ex.set_xlabel('Times (s)')
ex.set_ylabel('Energy (J)')




