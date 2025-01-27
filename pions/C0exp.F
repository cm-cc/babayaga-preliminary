        VirtPi = (1.05542899627435D-3*Sqrt(s*(s - ampi**2*4D0))*
     &    (3.94784176043574D1 + 
     &      (2.07944154167984D0 + 1.2D1*Log(ampi) - 
     &         6D0*Log(lambda2) - 6D0*Log(mu2dim**(-2)) - 
     &         6D0*Log(mu2dim**2) - 
     &         3D0*Log(ampi**2/(s - ampi**2*4D0)**2) - 
     &         3D0*Log(s - ampi**2*2D0 + 
     &            Sqrt(s*(s - ampi**2*4D0))))*
     &       (Log(5.D-1/ampi**2) + 
     &         Log(s - ampi**2*2D0 - 
     &           1D0*Sqrt(s*(s - ampi**2*4D0)))) + 
     &      1.2D1*ddilog(
     &        (5.D-1*(s - ampi**2*2D0 - 
     &             1D0*Sqrt(s*(s - ampi**2*4D0))))/ampi**2)))/
     &  (s*(-(s*1D0) + ampi**2*4D0))
                  
