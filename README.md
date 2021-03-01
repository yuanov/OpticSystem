# OpticSystem

#mbl
```
from maxwell_bloch.linear import MBLSystem

mbl.generate_consts()
mbl.load_consts()

mbl.heatmap_logD_alpha_V()

mbl.average_probabilities(plot=True, start_over=True)

mbl.run_d()

mbl.run_alpha()

mbl.animate_re_im()
```

# mb
```
from maxwell_bloch.non_linear import MBSystem

mb = MBSystem()

mb.generate_consts()
mb.read_consts()


mb.solve_equations()
mb.plot_dynamic()
mb.plot_fourier()
```

# Pole
```
from fresnel.pole import Pole

pole = Pole(alpha=0.5, log_scale=True, clip_t=10)
pole.plot(re_w=[3, 6], im_w=[-0.5, 0.5])
pole.inverted_pc_experiment(re_w=[3, 6], im_w=[-0.5, 0.5], alpha=[0, 1])
```

# Field
```
from fresnel.field import Field

field = Field(w=4.6, alpha=0.2)
field.plot()
```


