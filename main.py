from fresnel.pole import Pole

pole = Pole()
pole.inverted_pc_experiment(re_w=[3, 6], im_w=[-0.5, 0.5], alpha=[0, 1], log_scale=True, clip_t=0.5)
