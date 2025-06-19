# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:57:05 2025

@author: phan

Equation from this reference:
G. Akoun, and J.-P. Yonnet, "3D analytical calculation of the forces exerted between two cuboidal magnets," IEEE Trans. Magn., vol. MAG-20, no.5, Sept. 1984.

"""
import numpy as np
import matplotlib.pyplot as plt
# Define constants
mu0 = 4*np.pi*1e-7
J1 = 0.38 #in Tesla
J2 = 0.38 #in Tesla

# Define sizes of magnets
# Size of magnet 1 -- 2a, 2b, 2c and sizes of magnet 2 -- 2A, 2B, 2C
# Displacements of magnet 2 relative to magnet 1 are: alpha, beta, gamma

a = 20.0/2*1e-3 #in meter
b = 12.0/2*1e-3
c = 6.0/2*1e-3

A = 12.0/2*1e-3
B = 20.0/2*1e-3
C = 6.0/2*1e-3

alpha = -4*1e-3
beta = -4*1e-3
gamma = 8*1e-3
def MagneticForceCalc(a,b,c,J1, A,B,C,J2, alpha, beta, gamma):
    # Initializing the interaction energy W
    W = 0
    Fx = 0
    Fy = 0
    Fz = 0
    count = 0

    for i in range(2):
        for j in range(2):
            uij = alpha + (-1)**j*A - (-1)**i*a
            for k in range(2):
                for l in range(2):
                    vkl = beta + (-1)**l*B - (-1)**k*b
                    for p in range(2):
                        for q in range (2):
                            wpq = gamma + (-1)**q*C - (-1)**p*c
                            r = np.sqrt(uij**2 + vkl**2 + wpq**2)
                            psi_uij_vkl_wpq_r = 1.0/2*uij*(vkl**2 - wpq**2)*np.log(r - uij) \
                                                + 1.0/2*vkl*(uij**2 - wpq**2)*np.log(r - vkl) \
                                                + uij*vkl*wpq*np.arctan(uij*vkl/r/wpq) \
                                                + r/6*(uij**2 + vkl**2 -2*wpq**2)
                            W += J1*J2/(4*np.pi*mu0)*(-1)**(i+j + k+l + p+q)*psi_uij_vkl_wpq_r
                            phi_x_uij_vkl_wpq_r = 1.0/2*(vkl**2 - wpq**2)*np.log(r - uij) \
                                                + uij*vkl*np.log(r - vkl) \
                                                + vkl*wpq*np.arctan(uij*vkl/r/wpq) \
                                                + 1.0/2*r*uij
                            Fx += J1*J2/(4*np.pi*mu0)*(-1)**(i+j + k+l + p+q)*phi_x_uij_vkl_wpq_r
                            phi_y_uij_vkl_wpq_r = 1.0/2*(uij**2 - wpq**2)*np.log(r - vkl) \
                                                + uij*vkl*np.log(r - uij) \
                                                + uij*wpq*np.arctan(uij*vkl/r/wpq) \
                                                + 1.0/2*r*vkl
                            Fy += J1*J2/(4*np.pi*mu0)*(-1)**(i+j + k+l + p+q)*phi_y_uij_vkl_wpq_r
                            phi_z_uij_vkl_wpq_r = -uij*wpq*np.log(r - uij) \
                                                - vkl*wpq*np.log(r - vkl) \
                                                + uij*vkl*np.arctan(uij*vkl/r/wpq) \
                                                - r*wpq
                            Fz += J1*J2/(4*np.pi*mu0)*(-1)**(i+j + k+l + p+q)*phi_z_uij_vkl_wpq_r   
                            count += 1
    return [W, Fx, Fy, Fz, count]


## testing purpose
#d=0
#Output_list = MagneticForceCalc(a,b,c,J1, A,B,C,J2, alpha + d, beta, gamma)
#print(Output_list)
Output_list = []

for d in range(0, 32, 2):
    Output_list.append(MagneticForceCalc(a,b,c,J1, A,B,C,J2, alpha + d*1e-3, beta, gamma))
    
    
reshaped_array = np.reshape(Output_list, (16, 5))

plt.plot(range(0, 32, 2), reshaped_array[:,1], label = "Force_X", linestyle = "--")
plt.plot(range(0, 32, 2), reshaped_array[:,2], label = "Force_Y", linestyle = "-.")
plt.plot(range(0, 32, 2), reshaped_array[:,3], label = "Force_Z")
plt.xlabel("Displacement in +X direction [mm]")
plt.ylabel("EM force [N]")
plt.legend()
plt.grid(True)
    

    
                    
                                    