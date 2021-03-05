class Constants:
    def __init__(self):
        self.thickness = {'pas': 0.1,
                          'vac': 0,
                          'act': 0.8,
                          'res1': 5,
                          'res2': 7,
                          'neg': 0.1,
                          'pc1': 0.8,
                          'pc2': 0.1,
                          }

        self.pc_cells = 8  # number for cells in photon crystal
        self.thickness['vac_eq_pc'] = (self.thickness['pc1'] + self.thickness['pc2']) * self.pc_cells

        self.epsilon = {'mat': 2,  # background material for active medium
                        'pas': -7,  # + 5j
                        'vac': 1,
                        'res1': 3,
                        'res2': 4,
                        'vac_eq_pc': 1,
                        'neg': -7,
                        'pc1': 2,
                        'pc2': -7,
                        'act': 2
                        }

        self.material_tags = self.thickness.keys()

        omega_tls = 4.5
        self.tls = {'omega': omega_tls,
                    'gamma': 0.1 * omega_tls
                    }

        self.c = 1  # speed of light
