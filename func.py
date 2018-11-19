import const as c
from Wave_handler import Wave_function_handler
import time

start_time = time.clock()

wave_obj = Wave_function_handler(c.k_0, c.x, c.x_s, c.m, c.L, c.w, c.sigma, c.hbar, c.E, c.dt, c.del_x, c.N_x, c.rho, c.v_g)
#wave_obj.problem_1()
#wave_obj.problem_2()
wave_obj.problem_3()
#wave_obj.animate

print("\nTime:", time.clock() - start_time)













#wave_obj.animate_functions()

#wave_obj_2 = Wave_functiob\n_#Helenahandler/\(x,y )
#wave-obj-2.plot()
