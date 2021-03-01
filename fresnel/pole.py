import numpy as np
import time
import matplotlib.pyplot as plt
from cmath import sqrt, exp
from math import log10
from fresnel.constants import Constants
from tqdm import tqdm


class Pole(Constants):
    def __init__(self):
        super().__init__()

        self.alpha = None
        self.w = None

        self.surface_in, self.surface_out, self.medium = {}, {}, {}
        for material_tag in self.material_tags:
            self.surface_in[material_tag] = np.ndarray((2, 2), dtype=np.complex128)
            self.surface_out[material_tag] = np.ndarray((2, 2), dtype=np.complex128)
            self.medium[material_tag] = np.eye(2, dtype=np.complex128)

    def calculate_eps(self):
        w = self.w
        epsilon = self.epsilon
        alpha = self.alpha
        omega_tls = self.tls['omega']
        gamma_tls = self.tls['gamma']

        if abs(w - omega_tls + 1j * gamma_tls) < 0.2:
            return 1
        else:
            return epsilon['mat'] + alpha * gamma_tls / (w - omega_tls + 1j * gamma_tls)

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
        c = self.c

        for material_tag in self.material_tags:
            self.medium[material_tag][0][0] = exp(
                w * 1j * sqrt(self.epsilon[material_tag]) * self.thickness[material_tag] / c)
            self.medium[material_tag][1][1] = exp(
                -w * 1j * sqrt(self.epsilon[material_tag]) * self.thickness[material_tag] / c)

    def transition(self, matrix, material_tag):
        if material_tag == 'vac' or material_tag == 'vac_eq_pc':
            matrix = np.dot(self.medium[material_tag], matrix)
            return matrix
        if material_tag == 'pc':
            for i in range(self.pc_cells):
                matrix = self.transition(matrix, 'pas')
                matrix = self.transition(matrix, 'act')
                return matrix

        matrix = np.dot(self.surface_in[material_tag], matrix)
        matrix = np.dot(self.medium[material_tag], matrix)
        matrix = np.dot(self.surface_out[material_tag], matrix)
        return matrix

    def T_matrix(self):
        T = np.eye(2)
        self.epsilon['act'] = self.calculate_eps()
        self.epsilon['res1'] = self.epsilon['act']
        self.epsilon['res2'] = self.epsilon['act']
        self.refresh_surface()
        self.refresh_medium()

        T = self.transition(T, 'neg')
        T = self.transition(T, 'res1')
        T = self.transition(T, 'neg')
        T = self.transition(T, 'res2')
        T = self.transition(T, 'neg')

        return T

    @staticmethod
    def transmission_coefficient(T, log_scale, clip_t):
        t = abs(T[0][0] - T[0][1] * T[1][0] / T[1][1])
        if log_scale:
            return min(log10(1 + t), clip_t)
        else:
            return min(t, clip_t)

    def plot(self, re_w, im_w, num_re=150, num_im=75, log_scale=False, clip_t=float('inf'), alpha=0):
        t = np.zeros((num_im, num_re), dtype=float)
        re_w, im_w = np.linspace(re_w[0], re_w[1], num=num_re), np.linspace(im_w[0], im_w[1], num=num_im)
        for i, re in tqdm(enumerate(re_w), total=len(re_w), ncols=100):
            for j, im in enumerate(im_w):
                self.w = re + 1j * im
                T = self.T_matrix()
                t[j][i] = self.transmission_coefficient(T, clip_t=clip_t, log_scale=log_scale)

        self._plot(re_w, im_w, t)

    @staticmethod
    def _plot(re_w, im_w, t):
        plt.figure(figsize=(9.6, 5))
        font = {'size': 14}
        plt.rc('font', **font)
        cc = plt.contourf(re_w, im_w, t, cmap="coolwarm", levels=10)

        plt.colorbar(cc)
        plt.title('t(ω)', color='black')
        plt.xlabel('Re ω')
        plt.ylabel('Im ω')
        plt.grid(True)
        plt.tight_layout()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.show()
        plt.close()

    @staticmethod
    def _show(re_w, im_w, t):
        plt.rc('font', **{'size': 14})

        plt.title('t(ω)', color='black')
        plt.xlabel('Re ω')
        plt.ylabel('Im ω')
        cc = plt.contourf(re_w, im_w, t, cmap="coolwarm", levels=10)
        plt.colorbar(cc)
        plt.grid(True)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        plt.savefig("pictures/" + str(time.time_ns()) + ".png")
        if plt.waitforbuttonpress(timeout=0.2):
            plt.waitforbuttonpress()
        plt.clf()

    def inverted_pc_experiment(self, re_w, im_w, alpha, num_re=150, num_im=75, log_scale=False, clip_t=float('inf')):
        alphas = np.arange(alpha[0], alpha[1], 0.1) / self.epsilon['mat']

        self.epsilon['pas'] /= self.epsilon['mat']
        self.epsilon['mat'] = 1

        plt.ion()
        plt.figure(figsize=(9.6, 5))
        plt.draw()

        re_w, im_w = np.linspace(re_w[0], re_w[1], num=num_re), np.linspace(im_w[0], im_w[1], num=num_im)
        for self.alpha in alphas:
            t = np.zeros((num_im, num_re), dtype=float)
            for i, re in enumerate(re_w):
                for j, im in enumerate(im_w):
                    self.w = re + 1j * im
                    T = self.T_matrix()
                    t[j][i] = self.transmission_coefficient(T, clip_t=clip_t, log_scale=log_scale)
            self._show(re_w, im_w, t)
