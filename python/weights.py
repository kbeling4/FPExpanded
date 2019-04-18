import numpy as np

class Weights2:
    def __init__( self, nodes ):
        self.nodes = nodes
        self.N = len( self.nodes )
        self.A1 = np.array( np.zeros( ( self.N, self.N ) ) )
        self.A2 = np.array( np.zeros( ( self.N, self.N ) ) )

    def get_A1( self ):
        for i in range( 0, self.N ):
            for j in range( 0, self.N ):
                sum_i = 1
                sum_j = 1
                if j != i:
                    for k in range( 0, self.N ):
                        if k != i:
                            sum_i *= ( self.nodes[i] - self.nodes[k] )
                        if k != j: 
                            sum_j *= ( self.nodes[j] - self.nodes[k] )
                    self.A1[i,j] = sum_i / ( ( self.nodes[i] - self.nodes[j] ) * sum_j )
                if j == self.N - 1:
                    self.A1[i,i] = -1*sum( self.A1[i,:] )
        return self.A1

    def get_A2( self ):
        for i in range( 0, self.N ):
            for j in range( 0, self.N ):
                if j != i:
                    self.A2[i,j] = 2 * self.A1[i,j] * ( self.A1[i,j] - 1 / ( self.nodes[i] - self.nodes[j] ) )
                if j == self.N - 1:
                    self.A2[i,i] = -1*sum( self.A2[i,:] )
        return self.A2
