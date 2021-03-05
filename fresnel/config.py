import numpy as np
from cmath import sqrt, exp
from math import log10


class Material:
    def __init__(self, thickness, epsilon, active=False):
        self.thickness = thickness
        self.epsilon = epsilon
        self.active = active

        if active:
            self.omega_tls = 4.5
            self.gamma_tls = 0.1 * self.omega_tls

        self.T = np.eye(2, dtype=np.complex128)

    def _epsilon_dispersion(self, w, alpha):
        omega_tls, gamma_tls = self.omega_tls, self.gamma_tls

        if abs(w - omega_tls + 1j * gamma_tls) < 0.2:
            return 0
        else:
            return alpha * gamma_tls / (w - omega_tls + 1j * gamma_tls)

    @staticmethod
    def _surface(epsilon):
        surface_in = np.eye(2, dtype=np.complex128)
        surface_out = np.eye(2, dtype=np.complex128)

        sq = sqrt(epsilon)
        expr1 = (1 / sq + 1) / 2
        expr2 = (-1 / sq + 1) / 2

        surface_in[0][0] = expr1
        surface_in[0][1] = expr2
        surface_in[1][0] = expr2
        surface_in[1][1] = expr1

        expr1 = (sq + 1) / 2
        expr2 = (-sq + 1) / 2

        surface_out[0][0] = expr1
        surface_out[0][1] = expr2
        surface_out[1][0] = expr2
        surface_out[1][1] = expr1

        return surface_in, surface_out

    @staticmethod
    def _medium(w, epsilon, thickness):
        medium = np.eye(2, dtype=np.complex128)

        sq = sqrt(epsilon)

        medium[0][0] = exp(w * 1j * sq * thickness)
        medium[1][1] = exp(-w * 1j * sq * thickness)

        return medium

    def refresh(self, w, alpha):
        epsilon = self.epsilon
        if self.active:
            epsilon += self._epsilon_dispersion(w, alpha)

        surface_in, surface_out = self._surface(epsilon)
        medium = self._medium(w, epsilon, self.thickness)

        T = np.eye(2)
        T = np.dot(surface_in, T)
        T = np.dot(medium, T)
        T = np.dot(surface_out, T)

        self.T = T


class Config:
    def _sequence(self):
        pass

    def _refresh(self, w, alpha):
        for tag in self.__dict__:
            material = getattr(self, tag)
            material.refresh(w=w, alpha=alpha)

    def _T(self, w, alpha):
        self._refresh(w, alpha)

        T = np.eye(2)  # todo check complexity
        for cell in self._sequence():
            T = np.dot(cell.T, T)
        return T

    def t(self, w, alpha, log_scale, clip_t):
        """
        transmission_coefficient
        :param w:
        :param alpha:
        :param log_scale:
        :param clip_t:
        :return:
        """
        return self._t(self._T(w, alpha), log_scale, clip_t)

    @staticmethod
    def _t(T, log_scale, clip_t):
        t = abs(T[0][0] - T[0][1] * T[1][0] / T[1][1])
        if log_scale:
            return min(log10(1 + t), clip_t)
        else:
            return min(t, clip_t)


class NegActPC(Config):
    def __init__(self):
        self.negative = Material(thickness=0.1, epsilon=-7)
        self.active = Material(thickness=0.8, epsilon=2, active=True)

    def _sequence(self):
        sequence = []
        num_cells = 8
        for _ in range(num_cells):
            sequence.append(self.negative)
            sequence.append(self.active)
        return sequence


class Neg2Res(Config):
    def __init__(self):
        self.negative = Material(thickness=0.1, epsilon=-7)
        self.res1 = Material(thickness=0.8, epsilon=10)
        self.res2 = Material(thickness=0.9, epsilon=10)

    def _sequence(self):
        sequence = [self.negative, self.res1, self.negative, self.res2, self.negative]

        return sequence
