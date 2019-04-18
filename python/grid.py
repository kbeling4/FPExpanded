import numpy as np
from scipy import special
import physics as phy

class Grid:
    def __init__(self, Ne, Emin, Emax, Nx, Xmin, Xmax):
        self.Ne     = Ne
        self.Nx     = Nx
        self.Emin   = Emin
        self.Emax   = Emax
        self.Xmin   = Xmin
        self.Xmax   = Xmax
        self.Enodes = np.zeros( Ne )
        self.Nnodes = np.zeros( Ne )
        self.b      = np.zeros( Ne )
        self.Xnodes = np.linspace( Xmin, Xmax, self.Nx, endpoint=True)
        self.Svec   = np.zeros( self.Ne )
        self.dSvec   = np.zeros( self.Ne )
        self.Tvec   = np.zeros( self.Ne )
        self.dTvec   = np.zeros( self.Ne )
        self.ddTvec   = np.zeros( self.Ne )
        self.Agrid  = np.array( np.zeros( ( self.Ne, self.Ne ) ) )
        self.phi    = np.array( np.zeros( ( self.Ne, self.Nx ) ) )

    def get_EnodesChebyshev( self ):
        self.Enodes[0]         = self.Emin
        self.Enodes[self.Ne-1] = self.Emax
        for i in range( 2, self.Ne ):
            self.Enodes[i-1] = self.Emin + 0.5*( 1 - np.cos( (i - 1)/(self.Ne - 1)*np.pi ) )*(self.Emax - self.Emin)
        return self.Enodes

    def get_EnodesUniform( self ):
        self.Enodes[0]         = self.Emin
        self.Enodes[self.Ne-1] = self.Emax
        for i in range( 2, self.Ne ):
            self.Enodes[i-1] = self.Emin + ((i - 1)/( self.Ne - 1))*( self.Emax - self.Emin )
        return self.Enodes

    def get_EnodesGauss( self ):
        x, wi = special.roots_legendre( self.Ne, mu=False )
        for i in range( 0, self.Ne ):
            self.Enodes[i] = x[i]*( self.Emax - self.Emin )/2 + ( self.Emin + self.Emax )/2
        return self.Enodes
            
    def get_Normal( self ):
        for i in range( 0, self.Ne ):
            self.Nnodes[i] = ( self.Enodes[i] - self.Emin ) / ( self.Emax - self.Emin )
        return self.Nnodes
        
    def find_Enode( self, value):
        return (np.abs(self.Enodes - value)).argmin()

    def find_Xnode( self, value):
        return (np.abs(self.Xnodes - value)).argmin()

    def get_Svec( self, particle, material ):
        physics = phy.Physics
        for i in range( 0, self.Ne ):
            self.Svec[i] = physics.S( self.Enodes[i], particle, material )
        return self.Svec

    def get_dSvec( self, particle, material ):
        physics = phy.Physics
        for i in range( 0, self.Ne ):
            self.dSvec[i] = physics.dS( self.Enodes[i], particle, material )
        return self.dSvec

    def get_Tvec( self, particle, material ):
        physics = phy.Physics
        for i in range( 0, self.Ne ):
            self.Tvec[i] = physics.T( self.Enodes[i], particle, material )
        return self.Tvec

    def get_dTvec( self, particle, material ):
        physics = phy.Physics
        for i in range( 0, self.Ne ):
            self.dTvec[i] = physics.dT( self.Enodes[i], particle, material )
        return self.dTvec

    def get_ddTvec( self, particle, material ):
        physics = phy.Physics
        for i in range( 0, self.Ne ):
            self.ddTvec[i] = physics.ddT( self.Enodes[i], particle, material )
        return self.ddTvec

    def get_Agrid( self, A1, A2, particle, material ):
        self.get_Svec( particle, material )
        self.get_dSvec( particle, material )
        self.get_Tvec( particle, material )
        self.get_dTvec( particle, material )
        self.get_ddTvec( particle, material )
        for i in range(0, self.Ne ):
            for j in range(0, self.Ne ):
                del_x = self.Xnodes[i+1] - self.Xnodes[i]
                if j != i:
                    self.Agrid[i,j] = del_x*( (self.Svec[j] + self.dTvec[j])*A1[i,j] + 0.5*self.Tvec[j]*A2[i,j] )
                    
                if i == j:
                    self.Agrid[i,i] = del_x*( (self.dSvec[i] + 0.5*self.ddTvec[i])
                                              + (self.Svec[i] + self.dTvec[i])*A1[i,i]
                                              + 0.5*self.Tvec[i]*A2[i,i] ) - 1
        return self.Agrid

