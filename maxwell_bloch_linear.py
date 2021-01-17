import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la
import pickle
from math import exp, log10, log
from random import uniform, seed
from tqdm import tqdm

num_a = 50
num_sigma_d = num_a

gamma_sigma = 0.1


class MBLSystem:
    def __init__(self):
        self.delta = np.ndarray(shape=(num_a,), dtype=np.complex128)
        self.gamma = np.ndarray(shape=(num_a,), dtype=np.complex128)

        self.omega = np.ndarray(shape=(num_a, num_sigma_d), dtype=np.complex128)

        self.quantities = {}

    def generate_consts(self, alpha_average=1):
        omega_0 = 1 * pow(10, -3)
        alpha_dispersion = 1

        self.omega = np.zeros((num_a, num_a))
        for i in range(num_a):
            alpha = uniform(alpha_average - alpha_dispersion, alpha_average + alpha_dispersion)
            for j in range(num_sigma_d):
                self.omega[i][j] = omega_0 * exp(-abs(i - j) * alpha)
                self.omega[j][i] = omega_0 * exp(-abs(i - j) * alpha)

        for i in range(num_a):
            self.delta[i] = uniform(- gamma_sigma, gamma_sigma)
            self.gamma[i] = gamma_sigma * uniform(0.5, 0.8)

        pickle.dump((self.omega, self.delta, self.gamma), open("consts_MBL.p", "wb"))

    def load_consts(self):
        self.omega, self.delta, self.gamma = pickle.load(open("consts_MBL.p", "rb"))

    def update_omega(self, alpha_cur, alpha_prev):
        for i in range(num_a):
            for j in range(num_sigma_d):
                self.omega[i][j] = self.omega[i][j] * exp(-abs(i - j) * (alpha_cur - alpha_prev))
                self.omega[j][i] = self.omega[i][j]

    def build_matrix(self, D0):
        num_dim = num_a + num_sigma_d
        Nat = 10 ** 5
        omega = self.omega
        delta = self.delta
        gamma = self.gamma

        system_matrix = np.zeros((num_dim, num_dim), dtype=np.complex128)

        for i in range(num_a):
            system_matrix[i, i] = -1j * delta[i] - gamma[i]
            j = i + num_a
            system_matrix[j, j] = - gamma_sigma

        system_matrix[:num_a, num_a:] = -1j * Nat * omega
        system_matrix[num_a:, :num_a] = 1j * D0 * np.conjugate(np.transpose(omega))

        return system_matrix

    @staticmethod
    def eigen(matrix):
        freq, vec = la.eig(matrix)
        vec = np.transpose(vec)
        return freq, vec

    @staticmethod
    def sort_by_decay(freq, vec):
        idx = freq.real.argsort()[::-1]
        freq = freq[idx]
        vec = vec[idx, :]
        return freq, vec

    @staticmethod
    def sort_by_frequency(freq, vec):
        idx = freq.imag.argsort()[::-1]
        freq = freq[idx]
        vec = vec[idx, :]
        return freq, vec

    @staticmethod
    def mode_volume(vec):
        vol = np.ndarray(shape=(len(vec),))
        for i in range(len(vec)):
            vol[i] = np.sum(abs(vec[i, :num_a]) ** 2) / np.max(abs(vec[i, :num_a]) ** 2)
        return vol

    @staticmethod
    def correlation_function(freq, div=5):
        """
        hardly depends on delta_w
        """
        freq = sorted(freq.imag)
        delta_w = gamma_sigma / div

        num = len(freq)
        num_cells = int((np.max(freq) - np.min(freq)) / delta_w)

        photons = np.zeros(num_cells, dtype=int)
        w = np.min(freq) + delta_w / 2
        k = 0
        for i in range(num_cells):
            for j in range(k, num):
                if abs(freq[j] - w) <= delta_w / 2 + 10 ** -10:
                    photons[i] += 1
                else:
                    k = j
                    break
            w += delta_w

        probabilities = np.zeros(num)
        for i in range(num_cells):
            probabilities[photons[i]] += 1

        probabilities /= np.sum(probabilities)

        sum_up = 0
        for n in range(len(probabilities)):
            sum_up += probabilities[n] * n * (n - 1)

        sum_down = 0
        for n in range(len(probabilities)):
            sum_down += probabilities[n] * n

        return sum_up / (sum_down ** 2)

    @staticmethod
    def probabilities(freq, delta, incr):
        freq = sorted(freq.imag)
        probab = 0
        for i in range(len(freq) - 1):
            if incr > freq[i + 1] - freq[i] - delta > 0:
                probab += 1 / (len(freq) - 1)
        return probab

    @staticmethod
    def update_probabilities(avg, new, num):
        return avg * num / (num + 1) + new / (num + 1)

    def plot_graph(self, x, xlabel=None):
        for tag in self.quantities:
            quantity = self.quantities[tag]
            y = np.transpose(quantity)
            if np.shape(y) == (len(y),):
                y = y[np.newaxis, :]

            plt.figure(figsize=(18, 8))
            font = {'size': 18}
            plt.rc('font', **font)
            plt.title(tag, color='black')
            plt.xlabel(xlabel)
            for line in y:
                plt.plot(x, line, 'k', linestyle='-', linewidth=0.8)
            plt.grid(True)
            plt.show()

    def append_quantity(self, tag, quantity):
        try:
            self.quantities[tag].append(quantity)
        except KeyError:
            self.quantities[tag] = [quantity]

    def append_correlation_function(self, freq):
        self.append_quantity(tag='correlation function', quantity=self.correlation_function(freq[:num_a], div=8))

    def append_average_mode_volume(self, vec):
        self.append_quantity(tag='<V>', quantity=np.average(self.mode_volume(vec)))

    def append_max_mode_volume(self, vec):
        self.append_quantity(tag='max V', quantity=np.max(self.mode_volume(vec)))

    def append_decay(self, freq):
        self.append_quantity(tag='γ', quantity=sorted(freq.imag))

    def run_d(self):
        D0_start = 10 ** -3
        D0_stop = 1
        D0_div = 1.01
        D0_steps = int(log(D0_stop / D0_start, D0_div))
        D0_stamps = [D0_start * D0_div ** n for n in range(D0_steps)]

        for D0 in tqdm(D0_stamps, ncols=100):
            system_matrix = self.build_matrix(D0=D0)
            freq, vec = self.eigen(matrix=system_matrix)
            freq, vec = self.sort_by_decay(freq, vec)
            # self.append_correlation_function(freq=freq)
            self.append_average_mode_volume(vec=vec[:num_a])
            self.append_max_mode_volume(vec=vec[:num_a])
            self.append_decay(freq=freq)

        D0_log_stamps = np.log10(np.array(D0_stamps)).tolist()
        self.plot_graph(x=D0_log_stamps, xlabel='lg D')

    def run_alpha(self):
        alpha_start = 1
        alpha_stop = 3
        alpha_stamps = np.linspace(alpha_start, alpha_stop, num=1000)

        for alpha in tqdm(alpha_stamps, ncols=100):
            self.generate_consts(alpha_average=alpha)
            system_matrix = self.build_matrix(D0=0.01)
            freq, vec = self.eigen(matrix=system_matrix)
            freq, vec = self.sort_by_decay(freq, vec)
            # self.append_correlation_function(freq=freq)
            self.append_average_mode_volume(vec=vec[:num_a])
            # self.append_max_mode_volume(vec=vec[:num_a])
            # self.append_decay(freq=freq)

        self.plot_graph(x=alpha_stamps, xlabel='<α>')

    def average_probabilities(self):
        self.generate_consts()
        self.init_matrix()
        self.sort_vectors_by_frequency()
        incr_delta = 0.0005
        deltas = np.arange(0, 0.0125, incr_delta)
        i = 0

        probabilities, number_of_averages = pickle.load(open("probab0.p", "rb"))
        # probabilities = np.array([self.probabilities(self.freq[:num_a], delta, incr_delta) for delta in deltas])
        # number_of_averages = 1

        while True:
            self.generate_consts()
            self.init_matrix()
            self.sort_vectors_by_frequency()

            new_probabilities = np.array([self.probabilities(self.freq[:num_a], delta, incr_delta) for delta in deltas])
            probabilities = self.update_probabilities(probabilities, new_probabilities, number_of_averages)
            pickle.dump((probabilities, number_of_averages + 1), open("probab" + str(i) + ".p", "wb"))
            i += 1
            i %= 2
            number_of_averages += 1
            print(number_of_averages)
            break
        probabilities1, number_of_averages = pickle.load(open("probab3.p", "rb"))
        self.quantities.append(
            Plot(graph_x=deltas, graph_y=np.array([np.log(probabilities), np.log(probabilities1)]).T))
        self.plot_graph(xlabel='Δ')

    def animate_re_im(self):
        plt.ion()
        plt.figure(figsize=(19, 10))
        font = {'size': 18}
        plt.draw()

        D0_start = 10 ** -3
        D0_stop = 1
        D0_div = 2
        D0_steps = int(log(D0_stop / D0_start, D0_div))
        D0_stamps = [D0_start * D0_div ** n for n in range(D0_steps)]
        for D0 in D0_stamps:
            system_matrix = self.build_matrix(D0=D0)
            freq, vec = self.eigen(matrix=system_matrix)

            graph_x = freq.imag
            graph_y = freq.real

            plt.rc('font', **font)
            plt.xlabel('Re ω')
            plt.ylabel('Im ω')
            # plt.title(round(self.D0, 3))
            plt.scatter(graph_x, graph_y, c='k')
            if plt.waitforbuttonpress(timeout=1):
                plt.waitforbuttonpress()
            plt.clf()
        plt.close()

    def heatmap_logD_alpha_V(self):
        alpha_start = 1
        alpha_stop = 3
        alpha_stamps = np.linspace(alpha_start, alpha_stop, num=100)

        D0_start = 10 ** -3
        D0_stop = 1
        D0_div = 1.1
        D0_steps = int(log(D0_stop / D0_start, D0_div))
        D0_stamps = [D0_start * D0_div ** n for n in range(D0_steps)]
        D0_log_stamps = np.log10(np.array(D0_stamps)).tolist()

        graph_z = np.zeros((len(alpha_stamps), len(D0_log_stamps)))
        graph_z = np.transpose(graph_z)

        self.generate_consts(alpha_average=1)

        alpha_previous = alpha_start
        for i, alpha in enumerate(tqdm(alpha_stamps)):
            self.update_omega(alpha_cur=alpha, alpha_prev=alpha_previous)
            for j, D0 in enumerate(D0_stamps):
                system_matrix = self.build_matrix(D0=D0)
                freq, vec = self.eigen(matrix=system_matrix)
                freq, vec = self.sort_by_decay(freq, vec)

                graph_z[j][i] = np.average(self.mode_volume(vec[:num_a]))
            alpha_previous = alpha

        plt.figure(figsize=(18, 8))
        font = {'size': 14}
        plt.rc('font', **font)
        cc = plt.contourf(alpha_stamps, D0_log_stamps, graph_z, cmap="coolwarm", levels=10)
        plt.colorbar(cc)
        plt.grid(True)
        plt.xlabel('<α>')
        plt.ylabel('lg D')
        plt.show()
