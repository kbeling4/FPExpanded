import particle as prt
import weights as wgt
import material as mat
import spectrum as spec
import grid as gd

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg


def main():
    #------- Setup Problem ------------------------
    Emin =  100
    Emax = 1000
    Ne = 100
    
    Xmin = 0
    Xmax = 30
    Nx = 1000

    particle = prt.Particle()
    material = mat.Material()

    #------ Generate A matrix --------------------
    grid = gd.Grid( Ne, Emin, Emax, Nx, Xmin, Xmax )
    Enodes = grid.get_EnodesChebyshev()
    Ewts = wgt.Weights2( Enodes )
    A1 = Ewts.get_A1()
    A2 = Ewts.get_A2()

    A = grid.get_Agrid( A1, A2, particle, material )

    # ----- Initial Spectrum -------------------------------
    idx1 = grid.find_Enode( Emax )
    idx2 = grid.find_Enode( Emin )

    Spec = spec.Spectrum( Enodes, idx1, idx2 )  

    B = -1*np.array( Spec.gaussian( 700, 500 ) )
    B1 = -1*B
    
    # ----- Solver ---------------------------------------------
    for i in range( 0, Nx ):
        x = linalg.solve( A, B )
        grid.phi[:,i] = x
        B = (-1)*x
        if i % 100 == 0:
            print( 'step: ', i )

    idx4 = grid.find_Xnode( 25.0 )

    np.savetxt('output.txt', np.column_stack((grid.Enodes, grid.phi[:,idx4]) ), fmt="%1.6e", delimiter=' ')

    plt.figure( 1 )
    plt.plot( grid.Enodes, B1, 'r' ) 
    plt.plot( grid.Enodes, grid.phi[:,idx4], 'g' )
    plt.plot( grid.Enodes, grid.phi[:,idx4], 'o' )

    plt.show()
            
    
if __name__ == "__main__": main() 
