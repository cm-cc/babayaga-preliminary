        dpsi = -(1D0*(4D0 + 2D0*Log(lambda2/ame**2) + 
     &      Log(mu2dim/ame**2) + deltauv))
        dphi = -(2D0*(Log(lambda2/ampi**2) - 
     &      1D0*Log(mu2dim/ampi**2) - 1D0*deltauv)
     &    )
        dmpi = -(1D0*(7D0 + 3D0*Log(mu2dim/ampi**2) + 
     &      3D0*deltauv))
        dme = -(1D0*(4D0 + 3D0*Log(mu2dim/ame**2) + 
     &      3D0*deltauv))
        da = -(1.33333333333333D0*
     &    (Log(mu2dim/ame**2) + deltauv))
        de = 6.66666666666667D-1*
     &  (Log(mu2dim/ame**2) + deltauv)
                  
