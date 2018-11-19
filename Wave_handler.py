import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math as ma
import pylab


class Wave_function_handler:
    def __init__(self, k_0, x, x_s, m, L, w, sigma, hbar, E, dt, del_x, N_x, rho, v_g):
        self.__k_0 = k_0
        self.__x = x
        self.__x_s = x_s
        self.__m = m
        self.__L = L
        self.__w = w
        self.__sigma = sigma
        self.__hbar = hbar
        self.__E = E
        self.__dt = dt
        self.__del_x = del_x
        self.__N_x = N_x
        self.__rho = rho
        self.__v_g = v_g

    def get_del_x(self):
        return self.__del_x

    def get_k_0(self):
        return self.__k_0

    def get_E(self):
        return self.__E

    def get_x(self):
        return self.__x

    def get_hbar(self):
        return self.__hbar

    def get_x_s(self):
        return self.__x_s

    def get_m(self):
        return self.__m

    def get_L(self):
        return self.__L

    def get_w(self):
        return self.__w

    def get_sigma(self):
        return self.__sigma

    def get_dt(self):
        return self.__dt

    def get_N_x(self):
        return self.__N_x

    def get_rho(self):
        return self.__rho
    def get_v_g(self):
        return self.__v_g



    def set_sigma(self, sigma):
        self.__sigma = sigma

    def set_del_x(self, del_x):
        self.__del_x = del_x

    def set_k_o(self, k_o):
        self.__k_o = k_o

    def set_w(self, w):
        self.__w = w

    def set_x_s(self, x_s):
        self.__x_s = x_s

    def set_hbar(self, hbar):
        self.__hbar = hbar

    def set_x(self, x):
        self.__x = x

    def set_e(self, E):
        self.__E = E

    def set_m(self, m):
        self.__m = m

    def set_N_x(self, N_x):
        self.__N_x = N_x


    def set_rho(self, rho):
        self.__rho = rho

    def set_v_g(self, v_g):
        self.__v_g = v_g

    #Constant C for norm

    def c_factor(self):

        return 1 / ((np.pi * self.get_sigma() ** 2) ** (1 / 4))



    # Funksjoner ============================================================================================================#

    #Potensiale
    def V_func(self, x, c, b):
        l = int(b / self.get_del_x())

        V_list = [0.0] * int(self.get_N_x() / 2 + self.get_N_x() % 2 - l / 2 - l % 2)

        V_list += [c * self.get_E()] * (l)

        V_list += [0.0] * int(self.get_N_x() / 2 - l / 2)

        return np.array(V_list)



    #gauss
    def gaussian(self):
        return self.c_factor() * np.exp(-(self.get_x() - self.get_x_s()) ** (2) / (2*self.get_sigma() ** (2)))


    #Real and Imag
    def calc_psi_r(self):

        return self.gaussian() * np.cos(self.get_k_0() * self.get_x())

    def calc_psi_i(self):

        return self.gaussian() * np.sin(self.get_k_0() * self.get_x() - self.get_w() * self.get_dt() / 2)


    def calc_matrix(self):

        new_psi_r = self.calc_psi_r()
        new_psi_i = self.calc_psi_i()

        new_psi_r[0] = 0
        new_psi_r[-1] = 0
        new_psi_i[0] = 0
        new_psi_i[-1] = 0

        return new_psi_r + 1j * new_psi_i



    def Psi_plot(self, x, Psi_R, Psi_I, plot_squared=False):
        matplotlib.rcParams.update({'font.size': 22})

        plt.figure(figsize=(15, 10))

        plt.plot(x, Psi_R,'m', label='$\Psi_R$')

        plt.plot(x, Psi_I,'lightgreen', label='$\Psi_I$')

        plt.legend()

        plt.show()

        if plot_squared:
            Psi = Psi_R ** 2 + Psi_I ** 2

            plt.figure(figsize=(15, 10))

            plt.plot(x, Psi, label='$|\Psi(x,t=T)|^2$')

            plt.legend()

            plt.show()



    #Udef
    def Psi_propagate(self, c, b, plot=False, plot_squared=False):
        # Define the last remaining variables

        T = self.get_L() / (2 * self.get_v_g())

        Nt = int(T / self.get_dt())

        # Defining arrays

        V = self.V_func(self.get_x(), c, b)

        Psi = self.calc_matrix()
        PsiR0 = Psi.real
        PsiI0 = Psi.imag

        #To be iterated over
        PsiR = Psi.real
        PsiI = Psi.imag
        # Defining matrices


        A = -2 * np.diagflat([0.0] + [1.0] * (self.get_N_x() - 2) + [0.0])

        A += np.diagflat([0.0] + [1.0] * (self.get_N_x() - 2), 1)

        A += np.diagflat([1.0] * (self.get_N_x() - 2) + [0.0], -1)

        A = A * (self.get_hbar() * self.get_dt() / (2 * self.get_m() * self.get_del_x() ** 2))

        P = np.diagflat(V)  # I will assume that V is always zero at the endpoints.

        P = P * (self.get_dt() / self.get_hbar())

        Matrix = - A + P

        # Propagating wave function

        for t in range(Nt):

            PsiR += Matrix.dot(PsiI)

            PsiI += - Matrix.dot(PsiR)

            if t % 100 == 0:
                C = np.sqrt(sum(PsiR ** 2 + PsiI ** 2) * self.get_del_x())

                PsiR = PsiR/C #normalization preservation

                PsiI = PsiI/C

        if plot:
            self.Psi_plot(self.get_x(), PsiR0, PsiI0, plot_squared=plot_squared)

            self.Psi_plot(self.get_x(), PsiR, PsiI, plot_squared=plot_squared)

        return PsiR, PsiI

    def probabilities_ref_tra(self, PsiR, PsiI):

        Psi_squared = PsiR ** 2 + PsiI ** 2

        if self.get_N_x() % 2 == 1:

            p_ref = (sum(Psi_squared[:int(self.get_N_x() / 2)]) + 0.5 * Psi_squared[int(self.get_N_x() / 2)]) * self.get_del_x()

            p_tra = (sum(Psi_squared[int(self.get_N_x() / 2 + 1):]) + 0.5 * Psi_squared[int(self.get_N_x() / 2)]) * self.get_del_x()

        else:

            p_ref = sum(Psi_squared[:self.get_N_x() / 2]) * self.get_del_x()

            p_tra = sum(Psi_squared[self.get_N_x() / 2:]) * self.get_del_x()

        return p_ref, p_tra





    def problem_1(self):
        #sigma = 1.0
        self.Psi_propagate(0, 0, plot=True, plot_squared=True)

    def problem_2(self):
            self.Psi_propagate(0, 0, plot=True, plot_squared=False)

    def problem_3(self):

        c = 0.5

        width = self.get_L() / 50

        PsiR, PsiI = self.Psi_propagate(c, width, plot=False, plot_squared=False)

        x = np.linspace(0.0, self.get_L(), self.get_N_x())

        V = self.V_func(x, c, width)

        plt.figure(figsize=(15, 10))

        plt.plot(x, PsiR, label='$\Psi_R$')

        plt.plot(x, PsiI, label='$\Psi_I$')

        plt.plot(x, V / self.get_E(), label='$V$', linewidth=2)

        plt.legend()

        plt.show()

        p_ref, p_tra = self.probabilities_ref_tra(PsiR, PsiI)

        print("Probabilities of ")

        print("reflection: \t", p_ref)

        print("transmission: \t", p_tra)

        print("total: \t\t", p_ref + p_tra)