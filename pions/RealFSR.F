      RealFSR = (2D0*(-(p4q**3*(p1p2 + ame**2)*ampi**2*1D0) + 
     &      p3q**2*ampi**2*
     &       (-(1D0*(p1p3 - p1p4*1D0 - p1q*1D0)*
     &            (p2p3 - p2p4*1D0 - p2q*1D0)) + 
     &         p1p2*(ampi**2 - p3p4*1D0 - p3q*1D0) + 
     &         ame**2*(ampi**2 - p3p4*1D0 - p3q*1D0)) - 
     &      p3q*p4q*1D0*
     &       (p3p4*(p1q*p2q - 
     &            1D0*(p1p3 - p1p4*1D0)*
     &             (p2p3 - p2p4*1D0) + 
     &            p1p2*(ampi**2 - p3p4*1D0) + 
     &            ame**2*(ampi**2 - p3p4*1D0))*2D0 + 
     &         p3q*(p2p4*(p1p3 - p1q*1D0)*2D0 + 
     &            p1p4*2D0*
     &             (p2p3 - p2q*1D0 - p2p4*2D0) + 
     &            p1p2*(ampi**2 - p3p4*2D0) + 
     &            ame**2*(ampi**2 - p3p4*2D0))) + 
     &      p4q**2*(ampi**2*
     &          (-(1D0*(p1p3 + p1q - p1p4*1D0)*
     &               (p2p3 + p2q - p2p4*1D0)) + 
     &            p1p2*(ampi**2 - p3p4*1D0) + 
     &            ame**2*(ampi**2 - p3p4*1D0)) + 
     &         p3q*(p2p3*(p1q - p1p4*1D0)*2D0 + 
     &            p1p3*2D0*
     &             (p2q - p2p4*1D0 + p2p3*2D0) - 
     &            p1p2*1D0*(ampi**2 - p3p4*2D0) - 
     &            ame**2*1D0*(ampi**2 - p3p4*2D0)) + 
     &         p3q**2*ame**2*4D0)))/
     &  (p3q**2*p4q**2*
     &    (p3p4 + p3q + p4q + ampi**2)**2)
