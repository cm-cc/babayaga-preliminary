        elmat2 = 8D0*((2D0*(((p3q + ampi**2)*2.5D-1*
     &            (-((p1p3 + p1q - p1p4*1D0)*
     &                 (p2p3 + p2q - p2p4*1D0)*2D0) + 
     &              ame**2*
     &               (-(p3p4*2D0) + p3q*2D0 - p4q*2D0 + 
     &                 ampi**2*2D0)))/p3q**2 + 
     &         ((p4q + ampi**2)*2.5D-1*
     &            ((p1p3 - p1p4*1D0 - p1q*1D0)*
     &               (p2p4 + p2q - p2p3*1D0)*2D0 + 
     &              ame**2*
     &               (-(p3p4*2D0) - p3q*2D0 + p4q*2D0 + 
     &                 ampi**2*2D0)))/p4q**2 + 
     &         p1p2*(2D0 + 
     &            ((p3q + ampi**2)*2.5D-1*
     &               (-(p3p4*2D0) + p3q*2D0 - p4q*2D0 + 
     &                 ampi**2*2D0))/p3q**2 + 
     &            ((p4q + ampi**2)*2.5D-1*
     &               (-(p3p4*2D0) - p3q*2D0 + p4q*2D0 + 
     &                 ampi**2*2D0))/p4q**2) + ame**2*4D0))/
     &     (p3p4*2D0 + p3q*2D0 + p4q*2D0 + ampi**2*2D0)**
     &      2 + ((2.5D-1*
     &          (p1p4*2D0*
     &             (p2p3*ame**2*2D0 - p2p4*ame**2*2D0 - 
     &               (p2q + ame**2)*(p3q - p4q*1D0)*2D0
     &               ) + 
     &            (-(p3p4*2D0) + ampi**2*2D0)*
     &             (p1p2*ame**2*2D0 + 
     &               (ame**4 - p1q*ame**2*1D0 - 
     &                  p2q*(p1q + ame**2)*1D0)*2D0) + 
     &            p1p3*
     &             (-(p2p3*ame**2*4D0) + 
     &               p2p4*ame**2*4D0 + 
     &               (p2q + ame**2)*(p3q - p4q*1D0)*4D0
     &               )))/p2q**2 + 
     &       (-(p1p4*2D0*
     &             (p2p3**2*4D0 + p2p4*ampi**2*4D0 + 
     &               p2p3*(-(p2p4*4D0) - p3p4*4D0))) + 
     &          p1p3*2D0*
     &           (-(p2p4**2*4D0) + p2p4*p3p4*4D0 + 
     &             p2p3*(p2p4*4D0 - ampi**2*4D0)) + 
     &          (p1p2 + ame**2*2D0)*
     &           ((-(p3p4*2D0) + ampi**2*2D0)*
     &              (p3p4*2D0 + ampi**2*2D0) + 
     &             p2p3**2*4D0 + p2p4**2*4D0 - 
     &             p2p3*p2p4*8D0))/
     &        (-(p2p3*2D0) - p2p4*2D0 + p3p4*2D0 + 
     &           ampi**2*2D0)**2)/(p3p4*2D0 + ampi**2*2D0)**2)
                  
