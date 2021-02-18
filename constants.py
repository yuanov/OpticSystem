class Constants:
    def __init__(self):
        self.material_tags = ['vac', 'pas', 'act', 'res1', 'res2', 'vac_eq_pc']

        self.thickness = {'pas': 0.1,
                          'vac': 0,
                          'act': 0.8,
                          'res1': 5,
                          'res2': 7,
                          'vac_eq_pc': (0.1 + 0.8) * 8
                          }

        self.epsilon = {'mat': 2,  # background material for active medium
                        'pas': -7 + 5j,
                        'vac': 1,
                        'res1': 3,
                        'res2': 4,
                        'vac_eq_pc': 1
                        }

        omega_tls = 4.5
        self.tls = {'omega': omega_tls,
                    'gamma': 0.1 * omega_tls
                    }

        self.c = 1  # speed of light
        self.pc_cells = 8  # number for cells in photon crystal
        self.alpha = 0
