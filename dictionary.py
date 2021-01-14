epsilon_material = 1.2  # + 0.9j
epsilon_pas = -6.6  # - 0.9 * 2.34j  # non-active material, border
epsilon_vac = 1

pc_cells = 16
thickness_pas = 0.0595
thickness_act = 0.778
thickness_vac = 10  # brd

omega_tls = 4.82  # 4.75  # position of TLS point on diagram is (omega_tls, -gamma_tls)
gamma_tls = 0.1 * omega_tls  # poles should be at the top


class Dictionary:
    def __init__(self):
        self.thickness = {'pas': thickness_pas,
                          'vac': thickness_vac,
                          'act': thickness_act
                          }

        self.epsilon = {'material': epsilon_material,
                        'pas': epsilon_pas,
                        'vac': epsilon_vac
                        }

        self.tls = {'omega': omega_tls,
                    'gamma': gamma_tls
                    }

        self.rest = {'c': 1}
