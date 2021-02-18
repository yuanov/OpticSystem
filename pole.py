import numpy as np
import matplotlib.pyplot as plt
from cmath import sqrt, exp
from math import log10
from dictionary import pc_cells
from tqdm import tqdm

ReW_min = 5.25
ReW_max = 7
ImW_min = -0.2
ImW_max = 0.9

graph_t = 3 * 10 ** -3
smooth_out = False


class Pole:
    def __init__(self, thickness, epsilon, tls, rest):
        self.thickness = thickness
        self.epsilon = epsilon
        self.tls = tls
        self.rest = rest

        self.w = None
        self.alpha = 0.998

        self.material_tags = ['vac', 'pas', 'act']

        self.surface_in, self.surface_out, self.medium = {}, {}, {}
        for material_tag in self.material_tags:
            self.surface_in[material_tag] = np.ndarray((2, 2), dtype=np.complex128)
            self.surface_out[material_tag] = np.ndarray((2, 2), dtype=np.complex128)
            self.medium[material_tag] = np.eye(2, dtype=np.complex128)

        self.size_re = 150
        self.size_im = 150

        increment_ReW = (ReW_max - ReW_min) / self.size_re  # grid spacing
        increment_ImW = (ImW_max - ImW_min) / self.size_im

        self.window = {'ReW_min': ReW_min,
                       'ReW_max': ReW_max,
                       'ImW_min': ImW_min,
                       'ImW_max': ImW_max,
                       'increment_ReW': increment_ReW,
                       'increment_ImW': increment_ImW}

        self.Re = np.zeros(self.size_re)
        re_w = self.window['ReW_min']
        for count1 in range(self.size_re):
            self.Re[count1] = re_w
            re_w += self.window['increment_ReW']

        self.Im = np.zeros(self.size_im)
        im_w = self.window['ImW_min']
        for count1 in range(self.size_im):
            self.Im[count1] = im_w
            im_w += self.window['increment_ImW']

        self.t = np.zeros(self.size_re * self.size_im)
        self.t = self.t.reshape(self.size_im, self.size_re)

    def calculate_eps(self):
        w = self.w
        epsilon = self.epsilon
        alpha = self.alpha
        omega_tls = self.tls['omega']
        gamma_tls = self.tls['gamma']

        if abs(w - omega_tls + 1j * gamma_tls) < pow(10, -6):
            return 1
        else:
            return epsilon['material'] + alpha * gamma_tls / (w - omega_tls + 1j * gamma_tls)

    def refresh_surface(self):
        for material_tag in self.material_tags:
            if material_tag == 'vac':
                continue
            self.surface_in[material_tag][0][0] = (1 / sqrt(self.epsilon[material_tag]) + 1) / 2
            self.surface_in[material_tag][0][1] = (-1 / sqrt(self.epsilon[material_tag]) + 1) / 2
            self.surface_in[material_tag][1][0] = (-1 / sqrt(self.epsilon[material_tag]) + 1) / 2
            self.surface_in[material_tag][1][1] = (1 / sqrt(self.epsilon[material_tag]) + 1) / 2

            self.surface_out[material_tag][0][0] = (sqrt(self.epsilon[material_tag]) + 1) / 2
            self.surface_out[material_tag][0][1] = (-sqrt(self.epsilon[material_tag]) + 1) / 2
            self.surface_out[material_tag][1][0] = (-sqrt(self.epsilon[material_tag]) + 1) / 2
            self.surface_out[material_tag][1][1] = (sqrt(self.epsilon[material_tag]) + 1) / 2

    def refresh_medium(self):
        w = self.w
        c = self.rest['c']

        for material_tag in self.material_tags:
            self.medium[material_tag][0][0] = exp(
                w * 1j * sqrt(self.epsilon[material_tag]) * self.thickness[material_tag] / c)
            self.medium[material_tag][1][1] = exp(
                -w * 1j * sqrt(self.epsilon[material_tag]) * self.thickness[material_tag] / c)

    def transition_pc(self, matrix):
        for i in range(pc_cells):
            matrix = self.transition(matrix, 'act')
            matrix = self.transition(matrix, 'pas')
        return matrix

    def transition(self, matrix, material_tag):
        if material_tag == 'vac':
            matrix = np.dot(self.medium[material_tag], matrix)
            return matrix

        matrix = np.dot(self.surface_in[material_tag], matrix)
        matrix = np.dot(self.medium[material_tag], matrix)
        matrix = np.dot(self.surface_out[material_tag], matrix)
        return matrix

    def config(self):
        T = np.eye(2)
        self.refresh_surface()
        self.refresh_medium()
        T = self.transition(T, 'vac')
        T = self.transition_pc(T)
        T = self.transition(T, 'vac')

        return T

    @staticmethod
    def transmission_coefficient(T):
        return min(abs(T[0][0] - T[0][1] * T[1][0] / T[1][1]), graph_t)

    def plot_distribution(self):
        window = self.window
        transmission_coefficient = self.transmission_coefficient

        self.w = window['ReW_min'] + 1j * window['ImW_min']
        for count1 in tqdm(range(self.size_re)):
            self.w = self.w.real + 1j * self.window['ImW_min']
            for count2 in range(self.size_im):
                self.epsilon['act'] = self.calculate_eps()
                T = self.config()
                if smooth_out:
                    self.t[count2][count1] = log10(1 + transmission_coefficient(T))
                else:
                    self.t[count2][count1] = transmission_coefficient(T)
                self.w += 1j * window['increment_ImW']
            self.w += window['increment_ReW']

        plt.figure(figsize=(18, 8))
        font = {'size': 14}
        plt.rc('font', **font)
        cc = plt.contourf(self.Re, self.Im, self.t, cmap="coolwarm", levels=10)

        plt.colorbar(cc)
        plt.title('t(ω)', color='black')
        plt.xlabel('Re ω')
        plt.ylabel('Im ω')
        plt.grid(True)
        plt.show()
        plt.close()
