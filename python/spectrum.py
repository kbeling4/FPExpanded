import numpy as np

class Spectrum:
    def __init__( self, Egrid, idx1, idx2 ):
        self.Egrid = Egrid
        self.idx1 = idx1
        self.idx2 = idx2
        self.norm = np.zeros( len( Egrid ) )
        self.spec = np.zeros( len( Egrid ) )
        
    def normalizer( self ):
        diff = self.Egrid[self.idx1] - self.Egrid[self.idx2]
        for i in range( self.idx2, self.idx1+1 ):
            self.norm[i] = ( self.Egrid[i] - self.Egrid[self.idx2] ) / diff
        return self.norm

    def exponential( self, lam ):
        for i in range( self.idx2, self.idx1+1 ):
            self.spec[i] = (np.exp( -lam*self.norm[i] ) - np.exp( -lam ) ) / ( 1 - np.exp( -lam ) )
        return self.spec

    def gaussian( self, mean, var ):
        const = 1/np.sqrt(2*np.pi*var)
        for i in range( self.idx2, self.idx1+1 ):
            self.spec[i] = const*np.exp( -1*( self.Egrid[i] - mean )**2 / ( 2 * var ) )
        return self.spec
