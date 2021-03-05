import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from tqdm import tqdm


class Plot:
    @staticmethod
    def start():
        plt.ion()
        plt.figure(figsize=(9.6, 5))
        plt.draw()

    @staticmethod
    def show(x, y, xlabel=None, ylabel=None, title=None, save_as=None):
        plt.rc('font', **{'size': 14})

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(x, y)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.makedirs("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        if plt.waitforbuttonpress(timeout=0.4):
            plt.waitforbuttonpress()
        plt.clf()

    @staticmethod
    def plot(x, y, xlabel=None, ylabel=None, title=None, save_as=None):
        plt.figure(figsize=(9.6, 5))
        plt.rc('font', **{'size': 14})

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(x, y)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.makedirs("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        plt.show()


class Scatter:
    @staticmethod
    def start():
        plt.ion()
        plt.figure(figsize=(9.6, 5))
        plt.draw()

    @staticmethod
    def show(x, y, xlabel=None, ylabel=None, title=None, save_as=None):
        plt.rc('font', **{'size': 14})

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.scatter(x, y)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.makedirs("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        if plt.waitforbuttonpress(timeout=0.4):
            plt.waitforbuttonpress()
        plt.clf()

    @staticmethod
    def plot(x, y, xlabel=None, ylabel=None, title=None, save_as=None):
        plt.figure(figsize=(9.6, 5))
        plt.rc('font', **{'size': 14})

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.scatter(x, y)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.makedirs("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        plt.show()


class Heatmap:
    @staticmethod
    def start():
        plt.ion()
        plt.figure(figsize=(9.6, 5))
        plt.draw()

    @staticmethod
    def show(x, y, z, xlabel=None, ylabel=None, title=None, clip_t=float('inf'), save_as=None):
        plt.rc('font', **{'size': 14})

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if clip_t == float('inf'):
            cc = plt.contourf(x, y, z, cmap="coolwarm", levels=100)
        else:
            levels = np.linspace(0, clip_t, 101)
            cc = plt.contourf(x, y, z, cmap="coolwarm", levels=levels)
        plt.colorbar(cc)
        plt.grid(True)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.mkdir("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        if plt.waitforbuttonpress(timeout=0.2):
            plt.waitforbuttonpress()
        plt.clf()

    @staticmethod
    def plot(x, y, z, xlabel=None, ylabel=None, title=None, clip_t=float('inf'), save_as=None):
        plt.figure(figsize=(9.6, 5))
        font = {'size': 14}
        plt.rc('font', **font)
        if clip_t == float('inf'):
            cc = plt.contourf(x, y, z, cmap="coolwarm", levels=100)
        else:
            levels = np.linspace(0, clip_t, 101)
            cc = plt.contourf(x, y, z, cmap="coolwarm", levels=levels)

        plt.colorbar(cc)
        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.mkdir("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        plt.show()
        plt.close()


class Surface:
    @staticmethod
    def start():
        plt.ion()
        plt.figure(figsize=(9.6, 5))
        plt.draw()

    @staticmethod
    def show(x, y, z, xlabel=None, ylabel=None, title=None, clip_t=float('inf'), save_as=None):
        plt.rc('font', **{'size': 14})

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if clip_t == float('inf'):
            cc = plt.contourf(x, y, z, cmap="coolwarm", levels=100)
        else:
            levels = np.linspace(0, clip_t, 101)
            cc = plt.contourf(x, y, z, cmap="coolwarm", levels=levels)
        plt.colorbar(cc)
        plt.grid(True)
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.mkdir("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        if plt.waitforbuttonpress(timeout=0.2):
            plt.waitforbuttonpress()
        plt.clf()

    @staticmethod
    def plot(x, y, z, xlabel=None, ylabel=None, title=None, clip_t=float('inf'), save_as=None):
        x, y = np.meshgrid(x, y)

        # plt.figure(figsize=(9.6, 5))
        font = {'size': 11}
        plt.rc('font', **font)
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.plot_surface(x, y, z, rstride=1, cstride=1,
                        cmap='coolwarm', edgecolor='none')

        plt.title(title, color='black')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        figManager = plt.get_current_fig_manager()
        figManager.window.showMaximized()
        if save_as:
            try:
                os.mkdir("pictures/" + str(save_as) + '/')
            except FileExistsError:
                pass
            plt.savefig("pictures/" + str(save_as) + '/' + str(time.time_ns()) + ".png")
        plt.show()
        plt.close()


class ReImT:
    @staticmethod
    def _loop(re_w, im_w, alpha, config, log_scale=False, clip_t=float('inf')):
        t = np.zeros((len(im_w), len(re_w)), dtype=float)
        for i, re in enumerate(re_w):
            for j, im in enumerate(im_w):
                w = re + 1j * im
                t[j][i] = config.t(w=w, alpha=alpha, clip_t=clip_t, log_scale=log_scale)
        return t

    def _animate(self, re_w, im_w, alpha, num_re=150, num_im=75, num_alpha=15, clip_t=float('inf'), save_as=None,
                 *args, **kwargs):
        alphas = np.linspace(alpha[0], alpha[1], num_alpha)
        re_w = np.linspace(re_w[0], re_w[1], num_re)
        im_w = np.linspace(im_w[0], im_w[1], num_im)

        Heatmap.start()

        for alpha in alphas:
            t = self._loop(re_w, im_w, alpha, clip_t=clip_t, *args, **kwargs)
            Heatmap.show(re_w, im_w, t, xlabel='Re ω', ylabel='Im ω', title='t(ω), α =' + str(round(alpha, 2)),
                         clip_t=clip_t, save_as=save_as)

    def _plot(self, re_w, im_w, alpha, num_re=150, num_im=75, clip_t=float('inf'), save_as=None, *args, **kwargs):
        re_w = np.linspace(re_w[0], re_w[1], num_re)
        im_w = np.linspace(im_w[0], im_w[1], num_im)
        t = self._loop(re_w, im_w, alpha, clip_t=clip_t, *args, **kwargs)

        Heatmap.plot(re_w, im_w, t, xlabel='Re ω', ylabel='Im ω', title='t(ω), α =' + str(round(alpha, 2)),
                     save_as=save_as)

    def plot(self, alpha, *args, **kwargs):
        if type(alpha) is not list or len(alpha) == 1:
            animate = False
        elif len(alpha) == 2:
            animate = True
        else:
            raise ValueError('alpha is a value or a pair of boundary values')

        if animate:
            self._animate(alpha=alpha, *args, **kwargs)
        else:
            self._plot(alpha=alpha, *args, **kwargs)

    def animate_eps(self, re_w, im_w, epsilon, config, num_re=150, num_im=75, num_eps=73, clip_t=float('inf'),
                    save_as=None, *args, **kwargs):
        epsilon = np.linspace(epsilon[0], epsilon[1], num_eps)
        re_w = np.linspace(re_w[0], re_w[1], num_re)
        im_w = np.linspace(im_w[0], im_w[1], num_im)

        Heatmap.start()

        for eps in epsilon:
            config.negative.epsilon = eps
            t = self._loop(re_w, im_w, alpha=0, clip_t=clip_t, config=config, *args, **kwargs)
            Heatmap.show(re_w, im_w, t, xlabel='Re ω', ylabel='Im ω', title='t(ω), eps =' + str(round(eps, 1)),
                         clip_t=clip_t, save_as=save_as)


class ReT:
    @staticmethod
    def _loop(re_w, alpha, config, log_scale=False, clip_t=float('inf')):
        t = np.zeros((len(re_w),), dtype=float)
        for i, w in enumerate(re_w):
            t[i] = config.t(w=w, alpha=alpha, clip_t=clip_t, log_scale=log_scale)
        return t

    def _animate(self, re_w, alpha, config, num_re=150, num_alpha=15, save_as=None, *args, **kwargs):
        alphas = np.linspace(alpha[0], alpha[1], num_alpha)
        re_w = np.linspace(re_w[0], re_w[1], num_re)

        Plot.start()

        for alpha in alphas:
            t = self._loop(re_w, alpha, config, *args, **kwargs)
            Plot.show(re_w, t, xlabel='ω', ylabel='t(ω)', title='α =' + str(round(alpha, 4)), save_as=save_as)

    def _plot(self, re_w, alpha, config, num_re=150, save_as=None, *args, **kwargs):
        re_w = np.linspace(re_w[0], re_w[1], num_re)
        t = self._loop(re_w, alpha, config, *args, **kwargs)

        Plot.plot(re_w, t, xlabel='ω', ylabel='t(ω)', title='α =' + str(round(alpha, 2)), save_as=save_as)

    def plot(self, alpha, *args, **kwargs):
        if type(alpha) is not list or len(alpha) == 1:
            print(type(alpha), len(alpha))
            animate = False
        elif len(alpha) == 2:
            animate = True
        else:
            raise ValueError('alpha is a value or a pair of boundary values')

        if animate:
            self._animate(alpha=alpha, *args, **kwargs)
        else:
            self._plot(alpha=alpha, *args, **kwargs)


class ReAlphaT:
    @staticmethod
    def _loop(re_w, alpha, config, log_scale=False, clip_t=float('inf')):
        t = np.zeros((len(alpha), len(re_w)), dtype=float)
        for i, w in enumerate(re_w):
            for j, alp in enumerate(alpha):
                t[j][i] = config.t(w=w, alpha=alp, clip_t=clip_t, log_scale=log_scale)
        return t

    def plot(self, re_w, alpha, num_re=150, num_alpha=75, clip_t=float('inf'), graph='heatmap', save_as=None, *args,
             **kwargs):
        re_w = np.linspace(re_w[0], re_w[1], num_re)
        alpha = np.linspace(alpha[0], alpha[1], num_alpha)
        t = self._loop(re_w, alpha, clip_t=clip_t, *args, **kwargs)

        if graph == 'heatmap':
            Heatmap.plot(re_w, alpha, t, xlabel='ω', ylabel='α', title='t(ω,α)', save_as=save_as)
        elif graph == '3d':
            Surface.plot(re_w, alpha, t, xlabel='ω', ylabel='α', title='t(ω,α)', save_as=save_as)


class EpsT:
    @staticmethod
    def _optimize(t, height=0):
        peaks, _ = find_peaks(t, height=height)
        return t[peaks]

    def _loop(self, re_w, epsilon, config, log_scale=False, clip_t=float('inf')):
        t_peaks = []
        for eps in tqdm(epsilon, ncols=100):
            t = np.zeros((len(re_w, )))
            config.negative.epsilon = eps
            for j, w in enumerate(re_w):
                t[j] = (config.t(w=w, alpha=0, clip_t=clip_t, log_scale=log_scale))

            for t_peak in self._optimize(t):
                t_peaks.append([eps, t_peak])
        return np.array(t_peaks)

    def plot(self, re_w, epsilon, config, save_as=None, num_re=150, num_epsilon=73, *args, **kwargs):
        epsilon = np.linspace(epsilon[0], epsilon[1], num_epsilon)
        re_w = np.linspace(re_w[0], re_w[1], num_re)
        t = self._loop(re_w, epsilon, config, *args, **kwargs)
        Scatter.plot(t[:, 0], t[:, 1], xlabel='ω', title='t(ω)', save_as=save_as)


re_im_t = ReImT()
re_t = ReT()
re_alpha_t = ReAlphaT()
eps_t = EpsT()
