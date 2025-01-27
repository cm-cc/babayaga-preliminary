      double precision mxx2min,mxx2max,mxx,mxx2,mxxg,mxxgmin
      double precision thghgfmax
      double precision qhardest(0:3),qinside(40,0:3),qfake(0:3)
      logical thereisoneinside
      common/cutsvarsstrong2020/
!     CMD
     .     thavgmin,thavgmax,dfmax,ppmmin,csimax,pp,pm,thm,thp,z,thavg,
     .     dphi,
!     KLOE[I,II] tagged and untagged
     .     pzmin,pperpmin,mxx2min,mxx2max,egammamin,thgtlow,thgthigh,
     .     thgulow,thguhigh,pzm,pzp,pperpm,pperpp,mxx2,mxx,qhardest,
     .     qinside,qfake,thgtilda,
!     BES3
     .     cthpmmin,cthpmmax,egammamin1,egammamin2,
     .     acthgmax1,acthgmax2,acthgmin2,cm,cp,
!     B
     .     mxxgmin,mxxg,thghgfmax,thghgf,thmincutB,thmaxcutB,
     .     thereisoneinside
