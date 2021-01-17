from dictionary import Dictionary
from pole import Pole
from maxwell_bloch import MBSystem
from maxwell_bloch_linear import MBLSystem

# mbl
mbl = MBLSystem()

mbl.generate_consts()
# mbl.load_consts()

mbl.heatmap_logD_alpha_V()
# solve.average_probabilities()
# mbl.run_d()
# mbl.run_alpha()
# mbl.animate_re_im()


# mb
# mb = MBSystem()

# mb.generate_consts()
# mb.read_consts()


# mb.solve_equations()
# mb.plot_dynamic()
# mb.plot_fourier()

# Pole
# dict = Dictionary()
#
# pole = Pole(thickness=dict.thickness, epsilon=dict.epsilon, tls=dict.tls, rest=dict.rest)
# pole.plot_distribution()
