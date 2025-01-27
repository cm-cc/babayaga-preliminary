      implicit real*8 (a-h,o-z)
      common/ecms/ecms,ecmsnom
      common/epssoft/eps
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      dimension p3(0:3),p4(0:3),qph(40,0:3)

#if CMD
      parameter (nbin   = 600)
      parameter (ndistr =  10)
#elif KLOEI
      parameter (nbin   = 600)
      parameter (ndistr =  18)
#elif KLOEII
      parameter (nbin   = 600)
      parameter (ndistr =  17)
#elif BES3
      parameter (nbin   = 600)
      parameter (ndistr =  14)
#elif B
      parameter (nbin   = 600)
      parameter (ndistr =  13)
#else
! default
      parameter (nbin   = 100)
      parameter (ndistr =  18)
#endif      

      integer lung(ndistr)
