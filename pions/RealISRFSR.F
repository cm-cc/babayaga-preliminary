        RealISRFSR = -((1D0*(p2p4**3*p3q*
     &         (p1q - p1p3*2D0 + p1p4*2D0) + 
     &        p2q**2*(-(p1q*p2p3*(p3q + p4q)*1D0) - 
     &           p1p4*1D0*
     &            (p3q*(p2p3 + p3p4 + ampi**2) + 
     &              p4q*(p3p4 + ampi**2 - p2p3*1D0)) + 
     &           p1p3*(p3q*(p3p4 + ampi**2) + 
     &              p4q*(p3p4 + ampi**2 - p2p3*2D0))) - 
     &        p4q*1D0*(p3p4 + ampi**2 - p2p3*1D0)*
     &         (p3q*ame**2*
     &            (p1p3 + p3q - p1p4*1D0 - p4q*1D0)*
     &            2D0 - p2p3**2*1D0*
     &            (p1q + p1p3*2D0 - p1p4*2D0) + 
     &           p2p3*((p1p3 - p1p4*1D0)*
     &               (p3q - p4q*1D0) + 
     &              p1q*(p3p4 - ampi**2*1D0) + 
     &              p1p2*
     &               (p3q - p4q*1D0 - p3p4*2D0 + 
     &                 ampi**2*2D0) + 
     &              ame**2*
     &               (-(p3q*1D0) - p4q*1D0 - p3p4*2D0 + 
     &                 ampi**2*2D0))) + 
     &        p2p4**2*(p3q**2*
     &            (p1p2 + p1p4 + ame**2 - p1p3*1D0) + 
     &           p4q*(-(p1p4*p2q*1D0) + 
     &              p2p3*(p1q + p1p3*2D0 - p1p4*2D0))
     &            + p3q*
     &            (p4q*(ame**2 - p1p4*1D0) + 
     &              p1q*p2q*2D0 - p1p4*ampi**2*2D0 - 
     &              ame**2*ampi**2*2D0 + 
     &              p3p4*(ame**2 - p1p4*1D0 - p1q*1D0)*
     &               2D0 - 
     &              p1p2*1D0*
     &               (p4q - p3p4*2D0 + ampi**2*2D0) + 
     &              p1p3*
     &               (p4q + (p3p4 + ampi**2)*2D0 - 
     &                 p2q*3D0) + p1p4*p2q*4D0)) + 
     &        p2p4*(p3q**2*
     &            (-(p1q*p2q*1D0) + 
     &              (p1p3 - p1p4*1D0)*
     &               (p3p4 + ampi**2 - p2p3*1D0) + 
     &              p1p2*
     &               (p2p3 + p2q - p3p4*1D0 - 
     &                 ampi**2*1D0) + 
     &              ame**2*
     &               (p2p3 + p2q - p3p4*1D0 - 
     &                 ampi**2*1D0 + p4q*2D0)) + 
     &           p4q*(p2p3*
     &               (-(p4q*1D0*
     &                    (p1p3 + ame**2 - p1p4*1D0)) + 
     &                 (p1p4*(p3p4 + ampi**2) - 
     &                    p1p3*1D0*
     &                     (p3p4 + ampi**2 - p2q*1D0) + 
     &                    ame**2*(ampi**2 - p3p4*1D0) + 
     &                    p1q*(p2q - ampi**2*1D0))*2D0 + 
     &                 p1p2*
     &                  (-(p4q*1D0) - p3p4*2D0 + 
     &                    ampi**2*2D0)) + 
     &              p2q*
     &               (p1p3*
     &                  (p2q - 
     &                    (p3p4 + p4q + ampi**2)*2D0) + 
     &                 p1q*(p2q + p3p4 - ampi**2*3D0) + 
     &                 p1p4*ampi**2*4D0)) + 
     &           p3q*(-(p1p3*p2q**2*1D0) + 
     &              p1p3*p2p3**2*2D0 - 
     &              p1p4*p2p3**2*2D0 - 
     &              p1p4*p2p3*p2q*2D0 + 
     &              p1p4*p2q**2*2D0 - 
     &              p1p3*p2p3*p3p4*2D0 + 
     &              p1p4*p2p3*p3p4*2D0 + 
     &              p2p3*p3p4*ame**2*2D0 + 
     &              p2q*p3p4*ame**2*2D0 - 
     &              p3p4**2*ame**2*2D0 - 
     &              p4q**2*ame**2*2D0 - 
     &              p1p3*p2p3*ampi**2*2D0 + 
     &              p1p4*p2p3*ampi**2*2D0 + 
     &              p1p3*p2q*ampi**2*2D0 - 
     &              p1p4*p2q*ampi**2*2D0 - 
     &              p2p3*ame**2*ampi**2*2D0 - 
     &              p2q*ame**2*ampi**2*2D0 + 
     &              ame**2*ampi**4*2D0 + 
     &              p1p2*
     &               (p4q*(p2q + p3p4 + ampi**2) + 
     &                 (p3p4 + ampi**2 - p2p3*1D0 - 
     &                    p2q*1D0)*(ampi**2 - p3p4*1D0)*2D0
     &                 ) + 
     &              p1q*
     &               (p2q**2 + p3p4**2 - p2p3**2*1D0 - 
     &                 ampi**4*1D0 + 
     &                 p2p3*(ampi**2 - p2q*1D0)*2D0 + 
     &                 p2q*(ampi**2 - p4q*1D0 - p3p4*3D0)
     &                 ) + p1p3*p2q*p3p4*4D0 - 
     &              p1p4*p2q*p3p4*4D0 - 
     &              p4q*1D0*
     &               (p1p3*
     &                  (p3p4 + ampi**2 - p2p3*2D0 - 
     &                    ame**2*2D0) + 
     &                 p1p4*
     &                  (-(p3p4*1D0) - ampi**2*1D0 + 
     &                    p2p3*2D0 + p2q*2D0 + ame**2*2D0
     &                    ) + 
     &                 ame**2*(p3p4 + ampi**2 - p2q*5D0)))
     &           ) + p2q*
     &         (-(p3q**2*1D0*
     &              (p1p2*(p3p4 + ampi**2) + 
     &                ame**2*(p3p4 + ampi**2) + 
     &                p1p4*(p3p4 + ampi**2 - p2p3*2D0)))
     &            + p4q*
     &            (p1q*ampi**4 - p1q*p3p4**2*1D0 - 
     &              p1q*p2p3*ampi**2*1D0 + 
     &              p4q*
     &               (p1p3*(p3p4 + ampi**2) + 
     &                 ame**2*(p3p4 + ampi**2) + 
     &                 p2p3*(p1q - ame**2*1D0)) - 
     &              p1q*p2p3**2*2D0 - 
     &              p2p3*p3p4*ame**2*2D0 + 
     &              p3p4**2*ame**2*2D0 + 
     &              p1p3*p2p3*ampi**2*2D0 + 
     &              p2p3*ame**2*ampi**2*2D0 - 
     &              ame**2*ampi**4*2D0 - 
     &              p1p2*1D0*(p3p4 + ampi**2 - p2p3*1D0)*
     &               (-(p4q*1D0) - p3p4*2D0 + ampi**2*2D0)
     &               + p1q*p2p3*p3p4*3D0 - 
     &              p1p4*1D0*(p3p4 + ampi**2 - p2p3*1D0)*
     &               (ampi**2 - p3p4*1D0 + p2p3*3D0) - 
     &              p1p3*p2p3**2*4D0 + 
     &              p1p3*p2p3*p3p4*4D0) + 
     &           p3q*(p1q*p3p4**2 - 
     &              p1q*p2p3*p3p4*1D0 - 
     &              p1q*ampi**4*1D0 + 
     &              p1p4*p2p3*p3p4*2D0 - 
     &              p3p4**2*ame**2*2D0 + 
     &              p1p4*p2p3*ampi**2*2D0 + 
     &              ame**2*ampi**4*2D0 + 
     &              p1p2*
     &               (-(p2p3*p4q*1D0) - p3p4**2*2D0 + 
     &                 ampi**4*2D0) + 
     &              p1q*p2p3*ampi**2*3D0 + 
     &              p1p3*
     &               (p2p3**2 + 
     &                 p2p3*2D0*(p4q - ampi**2*2D0) + 
     &                 (p3p4 + ampi**2)*
     &                  (ampi**2 - p3p4*1D0 - p4q*3D0)) + 
     &              p4q*
     &               (p1p4*(p3p4 + ampi**2)*3D0 + 
     &                 p2p3*(p1q - ame**2*5D0))))))/
     &    (p2q*p3q*p4q*(p3p4 + ampi**2)*
     &      (p3p4 + p3q + p4q + ampi**2)*
     &      (p3p4 + ampi**2 - p2p3*1D0 - p2p4*1D0)))

   