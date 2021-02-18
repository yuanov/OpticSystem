from fresnel.pole import Pole
from fresnel.field import Field
from maxwell_bloch.non_linear import MBSystem
from maxwell_bloch.linear import MBLSystem

# mbl
# mbl = MBLSystem()

# mbl.generate_consts()
# mbl.load_consts()

# mbl.heatmap_logD_alpha_V()
# mbl.average_probabilities(plot=True, start_over=True)
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
pole = Pole(w=4.6, alpha=0.5)
pole.plot_distribution()

# Field
# field = Field(w=4.6, alpha=0.2)
# field.plot()

