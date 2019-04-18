import numpy as np

class Physics:
    def S( E, particle, material ):
        particle.get_qmax( E )
        amp = 0.1536*(particle.Z)**2*material.Z*material.rho/(material.M)
        return amp*( np.log( particle.qmax/material.qmin )*( E + particle.mass )**2/(E*(E+2*particle.mass))
                     + material.qmin/particle.qmax - 1 )
    
    def dS( E, particle, material ):
        particle.get_qmax( E )
        amp = (-1)*0.1536*(particle.Z)**2*material.Z*material.rho/(material.M)
        return amp*np.log( particle.qmax/material.qmin ) * ( 2*(particle.mass)**2 * (particle.mass + E )
                                                             / ( E**2 * ( 2*particle.mass + E )**2 ) )
    
    def T( E, particle, material ):
        particle.get_qmax( E )
        amp  = 0.1536*(particle.Z)**2*material.Z*material.rho/(material.M)
        return amp*( particle.qmax*( ( E + particle.mass)**2 / ( E*(2*particle.mass + E) ) ) - 0.5*particle.qmax
                     - material.qmin*( ( E + particle.mass)**2 / ( E*(2*particle.mass + E) ) )
                     + 0.5*material.qmin**2/particle.qmax )

    def dT( E, particle, material ):
        particle.get_qmax(E)
        amp = 0.1536*(particle.Z)**2*material.Z*material.rho/(material.M)
        return amp*( 2*particle.mass**2 * (particle.mass + E ) / ( E**2 *( 2*particle.mass + E )**2 )
                     * ( material.qmin - particle.qmax ) )

    def ddT( E, particle, material ):
        particle.get_qmax(E)
        amp = 0.1536*(particle.Z)**2*material.Z*material.rho/(material.M)
        return amp*( ( 2*particle.mass**2 ) * ( 2*particle.mass - 3*E ) / ( E**3 * ( 2*particle.mass + E )**2 )
                     * ( particle.qmax - material.qmin ) )


    
    # def TT( E, particle, material ):
    #     particle.get_qmax( E )
    #     amp = 0.1536*(particle.Z)**2*material.Z*material.rho/(material.M)
        
