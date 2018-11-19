#list of constants

import numpy as np
k_0 = 20.
hbar = 1.
m = 1.
L = 20.
N_x = int(2.5*k_0*L) + 1

del_x = L /(N_x-1)

x_s = 5
v_g = hbar * k_0/m

sigma = 1.0

#for sigmax in [0.1, 0.5, 1.0, 2.0]: #Problem 2

rho = 0.1    #scaling dt woop woop
# c = 0 for problem 1, c = 0 for 2, c = 0.5 for 3, c = 0.0 1.5 og 50 for 4
c = 0.5
E = hbar**2*k_0**2/2*m
w = E/hbar
dt = rho  *  2*m*hbar*del_x**2 / (2*m*c*E*del_x**2 + hbar**2)

x = np.linspace(0.0 , L, N_x)