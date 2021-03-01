import numpy as np
import matplotlib.pyplot as plt
from cmath import sqrt, exp
from fresnel.constants import Constants


class Field(Constants):
    def __init__(self, w, alpha=0):
        super().__init__()

        self.alpha = alpha
        self.w = w

        self.surface_in, self.surface_out, self.medium = {}, {}, {}
        for material_tag in self.material_tags:
            self.surface_in[material_tag] = np.ndarray((2, 2), dtype=np.complex128)
            self.surface_out[material_tag] = np.ndarray((2, 2), dtype=np.complex128)
            self.medium[material_tag] = np.eye(2, dtype=np.complex128)

        self.x = 0
        self.increment_x = 3 * 10 ** -5

        self.E = np.ndarray((2, 1), dtype=np.complex128)
        self.E[0] = 0
        self.E[1] = 1

        self.boxes = []

        self.x_graph = []
        self.y_graph = []

    def calculate_eps(self):
        w = self.w
        epsilon = self.epsilon
        alpha = self.alpha
        omega_tls = self.tls['omega']
        gamma_tls = self.tls['gamma']

        if abs(w - omega_tls + 1j * gamma_tls) < pow(10, -6):
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
                w * 1j * sqrt(self.epsilon[material_tag]) * self.increment_x / c)
            self.medium[material_tag][1][1] = exp(
                -w * 1j * sqrt(self.epsilon[material_tag]) * self.increment_x / c)

    @staticmethod
    def get_color(material_tag):
        if material_tag == 'vac':
            return 'black'
        elif material_tag == 'act':
            return 'red'
        elif material_tag == 'pas':
            return 'blue'

        raise BaseException('color not found')

    def _plot(self):
        x = np.array(self.x_graph)
        y = np.array(self.y_graph)

        dx = self.increment_x
        y /= np.max(y)
        x -= np.max(x)
        x = -x

        plt.figure(figsize=(18, 8))
        font = {'size': 18}
        plt.rc('font', **font)
        # plt.title('|E|^2 (x)', color='black')
        plt.xlabel('x')
        for x1, x2, material_tag in self.boxes:
            material_color = self.get_color(material_tag)
            plt.plot([x1 + dx, x1 + dx], [0, 1], color='gray')
            plt.plot([x1 + dx, x2 - dx], [1, 1], color=material_color)
            plt.plot([x2 - dx, x2 - dx], [0, 1], color='gray')

        plt.plot(x, y, '-b')

        plt.show()

    def y_quantity(self, E, material_tag):
        return sqrt(self.epsilon[material_tag]).real * pow(abs(E[0] + E[1]), 2)

    def transition(self, material_tag, thickness=None, add_boxes=True):
        E = self.E
        x = self.x
        increment_x = self.increment_x
        medium = self.medium[material_tag]
        surface_in, surface_out = self.surface_in[material_tag], self.surface_out[material_tag]
        if thickness is None:
            thickness = self.thickness[material_tag]

        if add_boxes:
            self.boxes.append((x, x + thickness, material_tag))
        if not material_tag == 'vac':
            E = np.dot(surface_in, E)

        x_current = thickness
        while x_current > 0:
            E = np.dot(medium, E)

            self.x_graph.append(x)
            self.y_graph.append(self.y_quantity(E, material_tag))

            x_current = x_current - increment_x
            x = x + increment_x
        if not material_tag == 'vac':
            E = np.dot(surface_out, E)

        self.E = E
        self.x = x

    def plot(self):
        self.epsilon['act'] = self.calculate_eps()
        self.refresh_medium()
        self.refresh_surface()

        # backward pass!
        self.transition('vac')
        for _ in range(self.pc_cells):
            self.transition('act')
            self.transition('pas')
        self.transition('vac')

        self._plot()

