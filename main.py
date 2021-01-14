from dictionary import Dictionary
from pole import Pole
from maxwell_bloch import MBSystem

# mb
mb = MBSystem()

mb.generate_consts()
# solve.read_consts()


mb.solve_equations()
mb.plot_dynamic()
mb.plot_fourier()

# Pole
# dict = Dictionary()
#
# pole = Pole(thickness=dict.thickness, epsilon=dict.epsilon, tls=dict.tls, rest=dict.rest)
# pole.plot_distribution()
