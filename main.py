import fresnel.config as config
from fresnel.dependency import re_im_t, re_t, re_alpha_t, eps_t

# configuration = config.NegActPC()
# re_im_t.plot(re_w=[3, 6], im_w=[-0.5, 0.5], alpha=[0.1, 0.3],
#              config=configuration,
#              log_scale=True,
#              clip_t=0.5,
#              save_as='T(re, im)')

# re_t.plot(re_w=[3, 6], alpha=[0, 1],
#           num_re=1000,
#           num_alpha=50,
#           config=configuration,
#           log_scale=True,
#           clip_t=0.5,
#           save_as='T(w)/4.03'
#           )
# re_alpha_t.plot(re_w=[4, 5.2], alpha=[0, 0.01],
#                 num_alpha=20,
#                 num_re=200,
#                 config=configuration,
#                 log_scale=True,
#                 clip_t=1,
#                 graph='3d',
# #                save_as='T(w, alpha)/4.03'
#                 )
configuration = config.Neg2Res()
# re_im_t.animate_eps(re_w=[3, 6], im_w=[-0.5, 0.5], epsilon=[-10, 11],
#                     config=configuration,
#                     log_scale=True,
#                     clip_t=0.5,
#                     )
eps_t.plot(re_w=[3.48, 5.01], epsilon=[-10, 11],
           num_epsilon=7314,
           config=configuration)
