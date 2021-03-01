class Constants:
    def __init__(self):
        self.thickness = {'pas': 0.1,
                          'vac': 0,
                          'act': 0.8,
                          'res1': 5,
                          'res2': 7,
                          'neg': 0.1
                          }
        self.thickness['vac_eq_pc'] = (self.thickness['pas'] + self.thickness['act']) * 8

        self.epsilon = {'mat': 2,  # background material for active medium
                        'pas': -7,  # + 5j
                        'vac': 1,
                        'res1': 3,
                        'res2': 4,
                        'vac_eq_pc': 1,
                        'neg': -7
                        }

        self.material_tags = self.thickness.keys()

        omega_tls = 4.5
        self.tls = {'omega': omega_tls,
                    'gamma': 0.1 * omega_tls
                    }

        self.c = 1  # speed of light
        self.pc_cells = 8  # number for cells in photon crystal
        self.alpha = 0
