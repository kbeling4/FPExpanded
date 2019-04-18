import numpy as np

class Particle:
    def __init__(self):
        self.name   = "Proton"
        self.mass   = 938.27231
        self.Z      = 1
        self.qmax   = 0
        self.beta2  = 0
        self.gamma2 = 0
        
    def get_qmax( self, energy ):
        self.gamma2 = np.power( ( energy + self.mass ) / self.mass, 2 )
        self.beta2  = 1 - ( 1 / self.gamma2  )
        self.qmax =  1.022 * self.beta2 * self.gamma2
