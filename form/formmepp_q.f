*
* quadruple precision version for handling isolated points with elmat2 < 0
*     
*      cfsr = (-im)*(2.d0*1.d0/(qp+qm+k)^2)            !c5 phase: -i
*      cfsrp= (+im)*1.d0/((qp+k)^2-mpi2)/((qp+qm+k)^2) !c3 phase: +i
*      cfsrm= (+im)*1.d0/((qm+k)^2-mpi2)/((qp+qm+k)^2) !c4 phase: +i
*      cisrp= (+im)*1.d0/((pp-k)^2-me2)/(qp+qm)^2      !c1 phase: +i
*      cisrm= (+im)*1.d0/((pm-k)^2-me2)/(qp+qm)^2      !c2 phase: +i
*
      im= (0.d0,1.d0)
      cfsr = (-im)*2.q0*1.q0/qppqmpk2                 !c5 phase: -i
      cfsrp= (+im)*1.q0/qppk2mmpi2/qppqmpk2           !c3 phase: +i
      cfsrm= (+im)*1.q0/qmpk2mmpi2/qppqmpk2           !c4 phase: +i
      cisrp= (-im)*1.q0/ppmk2mme2/qppqm2              !c1 phase: -i
      cisrm= (-im)*1.q0/pmmk2mme2/qppqm2              !c2 phase: -i

      cfsrd = conjg(cfsr)
      cfsrpd= conjg(cfsrp)
      cfsrmd= conjg(cfsrm)
      cisrpd= conjg(cisrp)
      cisrmd= conjg(cisrm)

c$$$      cfsr2 = abs((-im)*2.d0*1.d0/qppqmpk2)                 !c5 phase: -i
c$$$      cfsrp2= abs((+im)*1.d0/qppk2mmpi2/qppqmpk2)           !c3 phase: +i
c$$$      cfsrm2= abs((+im)*1.d0/qmpk2mmpi2/qppqmpk2)           !c4 phase: +i
c$$$      cisrp2= abs((+im)*1.d0/ppmk2mme2/qppqm2)              !c1 phase: +i
c$$$      cisrm2= abs((+im)*1.d0/pmmk2mme2/qppqm2)              !c2 phase: +i

      interf =
     &  + cfsr * ( 16*cfsrd*me2*vp122 + 8*cfsrpd*me2*mpi2*vp122 + 8*
     &    cfsrmd*me2*mpi2*vp122 - 8*pm_pp*qm_qp*cfsrpd*vp122 - 8*pm_pp*
     &    qm_qp*cfsrmd*vp122 - 4*pm_pp*qm_k*cfsrpd*vp122 + 12*pm_pp*
     &    qm_k*cfsrmd*vp122 - 8*pm_pp*qm_k*cisrpd*vp12*vp34 + 8*pm_pp*
     &    qm_k*cisrmd*vp12*vp34 + 12*pm_pp*qp_k*cfsrpd*vp122 - 4*pm_pp*
     &    qp_k*cfsrmd*vp122 + 8*pm_pp*qp_k*cisrpd*vp12*vp34 - 8*pm_pp*
     &    qp_k*cisrmd*vp12*vp34 + 8*pm_pp*cfsrd*vp122 + 8*pm_pp*cfsrpd*
     &    mpi2*vp122 + 8*pm_pp*cfsrmd*mpi2*vp122 - 16*pm_qm*pp_qm*
     &    cfsrmd*vp122 + 8*pm_qm*pp_qp*cfsrpd*vp122 + 8*pm_qm*pp_qp*
     &    cfsrmd*vp122 + 4*pm_qm*pp_k*cfsrpd*vp122 - 12*pm_qm*pp_k*
     &    cfsrmd*vp122 - 8*pm_qm*pp_k*cisrpd*vp12*vp34 - 8*pm_qm*pp_k*
     &    cisrmd*vp12*vp34 - 8*pm_qm*cisrpd*me2*vp12*vp34 - 8*pm_qm*
     &    cisrmd*me2*vp12*vp34 + 8*pm_qp*pp_qm*cfsrpd*vp122 + 8*pm_qp*
     &    pp_qm*cfsrmd*vp122 - 16*pm_qp*pp_qp*cfsrpd*vp122 - 12*pm_qp*
     &    pp_k*cfsrpd*vp122 + 4*pm_qp*pp_k*cfsrmd*vp122 + 8*pm_qp*pp_k*
     &    cisrpd*vp12*vp34 )
      interf = interf + cfsr * ( 8*pm_qp*pp_k*cisrmd*vp12*vp34 + 8*
     &    pm_qp*cisrpd*me2*vp12*vp34 + 8*pm_qp*cisrmd*me2*vp12*vp34 + 4
     &    *pm_k*pp_qm*cfsrpd*vp122 - 12*pm_k*pp_qm*cfsrmd*vp122 + 8*
     &    pm_k*pp_qm*cisrpd*vp12*vp34 + 8*pm_k*pp_qm*cisrmd*vp12*vp34
     &     - 12*pm_k*pp_qp*cfsrpd*vp122 + 4*pm_k*pp_qp*cfsrmd*vp122 - 8
     &    *pm_k*pp_qp*cisrpd*vp12*vp34 - 8*pm_k*pp_qp*cisrmd*vp12*vp34
     &     - 8*pm_k*pp_k*cfsrpd*vp122 - 8*pm_k*pp_k*cfsrmd*vp122 + 8*
     &    pp_qm*cisrpd*me2*vp12*vp34 + 8*pp_qm*cisrmd*me2*vp12*vp34 - 8
     &    *pp_qp*cisrpd*me2*vp12*vp34 - 8*pp_qp*cisrmd*me2*vp12*vp34 - 
     &    8*qm_qp*cfsrpd*me2*vp122 - 8*qm_qp*cfsrmd*me2*vp122 - 4*qm_k*
     &    cfsrpd*me2*vp122 + 12*qm_k*cfsrmd*me2*vp122 - 16*qm_k*cisrpd*
     &    me2*vp12*vp34 + 16*qm_k*cisrmd*me2*vp12*vp34 + 12*qp_k*cfsrpd
     &    *me2*vp122 - 4*qp_k*cfsrmd*me2*vp122 + 16*qp_k*cisrpd*me2*
     &    vp12*vp34 - 16*qp_k*cisrmd*me2*vp12*vp34 )
      interf = interf + cfsrp * ( 8*cfsrd*me2*mpi2*vp122 + 32*cfsrpd*
     &    me2*mpi4*vp122 - 32*pm_pp*pm_qp*qm_qp*cisrmd*vp12*vp34 - 16*
     &    pm_pp*pm_qp*qm_k*cisrmd*vp12*vp34 + 16*pm_pp*pm_qp*qp_k*
     &    cisrmd*vp12*vp34 + 32*pm_pp*pm_qp*cisrmd*mpi2*vp12*vp34 - 16*
     &    pm_pp*pm_k*qm_qp*cisrmd*vp12*vp34 - 8*pm_pp*pm_k*qm_k*cisrmd*
     &    vp12*vp34 + 8*pm_pp*pm_k*qp_k*cisrmd*vp12*vp34 + 16*pm_pp*
     &    pm_k*cisrmd*mpi2*vp12*vp34 + 32*pm_pp*pp_qp*qm_qp*cisrpd*vp12
     &    *vp34 + 16*pm_pp*pp_qp*qm_k*cisrpd*vp12*vp34 - 16*pm_pp*pp_qp
     &    *qp_k*cisrpd*vp12*vp34 - 32*pm_pp*pp_qp*cisrpd*mpi2*vp12*vp34
     &     + 16*pm_pp*pp_k*qm_qp*cisrpd*vp12*vp34 + 8*pm_pp*pp_k*qm_k*
     &    cisrpd*vp12*vp34 - 8*pm_pp*pp_k*qp_k*cisrpd*vp12*vp34 - 16*
     &    pm_pp*pp_k*cisrpd*mpi2*vp12*vp34 + 16*pm_pp*qm_qp*qm_k*cfsrmd
     &    *vp122 - 32*pm_pp*qm_qp*qp_k*cfsrpd*vp122 + 16*pm_pp*qm_qp*
     &    qp_k*cfsrmd*vp122 - 16*pm_pp*qm_qp*qp_k*cisrpd*vp12*vp34 + 16
     &    *pm_pp*qm_qp*qp_k*cisrmd*vp12*vp34 - 8*pm_pp*qm_qp*cfsrd*
     &    vp122 )
      interf = interf + cfsrp * (  - 32*pm_pp*qm_qp*cfsrpd*mpi2*vp122
     &     - 32*pm_pp*qm_qp*cfsrmd*mpi2*vp122 + 32*pm_pp*qm_qp2*
     &    cfsrmd*vp122 - 32*pm_pp*qm_k*qp_k*cfsrpd*vp122 - 16*pm_pp*
     &    qm_k*qp_k*cisrpd*vp12*vp34 + 16*pm_pp*qm_k*qp_k*cisrmd*vp12*
     &    vp34 - 4*pm_pp*qm_k*cfsrd*vp122 - 32*pm_pp*qm_k*cfsrpd*mpi2*
     &    vp122 - 16*pm_pp*qm_k*cfsrmd*mpi2*vp122 + 12*pm_pp*qp_k*cfsrd
     &    *vp122 + 64*pm_pp*qp_k*cfsrpd*mpi2*vp122 - 16*pm_pp*qp_k*
     &    cfsrmd*mpi2*vp122 + 16*pm_pp*qp_k*cisrpd*mpi2*vp12*vp34 - 16*
     &    pm_pp*qp_k*cisrmd*mpi2*vp12*vp34 + 32*pm_pp*qp_k2*cfsrpd*
     &    vp122 + 16*pm_pp*qp_k2*cisrpd*vp12*vp34 - 16*pm_pp*qp_k2*
     &    cisrmd*vp12*vp34 + 8*pm_pp*cfsrd*mpi2*vp122 + 32*pm_pp*cfsrpd
     &    *mpi4*vp122 - 32*pm_qm*pm_qp*pp_qm*cisrmd*vp12*vp34 + 32*
     &    pm_qm*pm_qp*pp_qp*cisrmd*vp12*vp34 + 16*pm_qm*pm_qp*pp_k*
     &    cisrmd*vp12*vp34 - 16*pm_qm*pm_k*pp_qm*cisrmd*vp12*vp34 + 16*
     &    pm_qm*pm_k*pp_qp*cisrmd*vp12*vp34 + 8*pm_qm*pm_k*pp_k*cisrmd*
     &    vp12*vp34 )
      interf = interf + cfsrp * ( 32*pm_qm*pp_qm*pp_qp*cisrpd*vp12*vp34
     &     + 16*pm_qm*pp_qm*pp_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qm*
     &    qm_qp*cfsrmd*vp122 + 16*pm_qm*pp_qm*qm_k*cfsrmd*vp122 - 32*
     &    pm_qm*pp_qm*qp_k*cfsrpd*vp122 + 16*pm_qm*pp_qm*qp_k*cfsrmd*
     &    vp122 - 16*pm_qm*pp_qm*qp_k*cisrpd*vp12*vp34 + 16*pm_qm*pp_qm
     &    *qp_k*cisrmd*vp12*vp34 - 32*pm_qm*pp_qm*cfsrpd*mpi2*vp122 - 
     &    32*pm_qm*pp_qp*pp_k*cisrpd*vp12*vp34 - 32*pm_qm*pp_qp*qm_qp*
     &    cfsrmd*vp122 - 16*pm_qm*pp_qp*qm_k*cfsrmd*vp122 - 16*pm_qm*
     &    pp_qp*qm_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qp*qp_k*cfsrpd*
     &    vp122 - 16*pm_qm*pp_qp*qp_k*cfsrmd*vp122 + 32*pm_qm*pp_qp*
     &    qp_k*cisrpd*vp12*vp34 - 16*pm_qm*pp_qp*qp_k*cisrmd*vp12*vp34
     &     + 8*pm_qm*pp_qp*cfsrd*vp122 + 32*pm_qm*pp_qp*cfsrpd*mpi2*
     &    vp122 - 32*pm_qm*pp_qp2*cisrpd*vp12*vp34 + 16*pm_qm*pp_k*
     &    qm_qp*cisrpd*vp12*vp34 + 32*pm_qm*pp_k*qp_k*cfsrpd*vp122 - 16
     &    *pm_qm*pp_k*qp_k*cisrmd*vp12*vp34 + 4*pm_qm*pp_k*cfsrd*vp122
     &     + 32*pm_qm*pp_k*cfsrpd*mpi2*vp122 )
      interf = interf + cfsrp * (  - 16*pm_qm*pp_k*cisrpd*mpi2*vp12*
     &    vp34 - 8*pm_qm*pp_k2*cisrpd*vp12*vp34 + 32*pm_qp*pm_k*pp_qm
     &    *cisrmd*vp12*vp34 - 32*pm_qp*pm_k*pp_qp*cisrmd*vp12*vp34 - 8*
     &    pm_qp*pm_k*pp_k*cisrmd*vp12*vp34 - 32*pm_qp*pp_qm*pp_qp*
     &    cisrpd*vp12*vp34 - 16*pm_qp*pp_qm*pp_k*cisrpd*vp12*vp34 - 32*
     &    pm_qp*pp_qm*qm_qp*cfsrmd*vp122 - 16*pm_qp*pp_qm*qm_k*cfsrmd*
     &    vp122 + 16*pm_qp*pp_qm*qm_k*cisrmd*vp12*vp34 + 32*pm_qp*pp_qm
     &    *qp_k*cfsrpd*vp122 - 16*pm_qp*pp_qm*qp_k*cfsrmd*vp122 + 16*
     &    pm_qp*pp_qm*qp_k*cisrpd*vp12*vp34 - 32*pm_qp*pp_qm*qp_k*
     &    cisrmd*vp12*vp34 + 8*pm_qp*pp_qm*cfsrd*vp122 + 32*pm_qp*pp_qm
     &    *cfsrpd*mpi2*vp122 + 32*pm_qp*pp_qp*pp_k*cisrpd*vp12*vp34 + 
     &    32*pm_qp*pp_qp*qm_qp*cfsrmd*vp122 + 16*pm_qp*pp_qp*qm_k*
     &    cfsrmd*vp122 + 16*pm_qp*pp_qp*qm_k*cisrpd*vp12*vp34 - 16*
     &    pm_qp*pp_qp*qm_k*cisrmd*vp12*vp34 - 32*pm_qp*pp_qp*qp_k*
     &    cfsrpd*vp122 + 16*pm_qp*pp_qp*qp_k*cfsrmd*vp122 - 32*pm_qp*
     &    pp_qp*qp_k*cisrpd*vp12*vp34 )
      interf = interf + cfsrp * ( 32*pm_qp*pp_qp*qp_k*cisrmd*vp12*vp34
     &     - 16*pm_qp*pp_qp*cfsrd*vp122 - 32*pm_qp*pp_qp*cfsrpd*mpi2*
     &    vp122 + 32*pm_qp*pp_qp2*cisrpd*vp12*vp34 + 16*pm_qp*pp_k*
     &    qm_qp*cisrmd*vp12*vp34 + 16*pm_qp*pp_k*qm_k*cisrpd*vp12*vp34
     &     - 32*pm_qp*pp_k*qp_k*cfsrpd*vp122 - 16*pm_qp*pp_k*qp_k*
     &    cisrpd*vp12*vp34 + 16*pm_qp*pp_k*qp_k*cisrmd*vp12*vp34 - 12*
     &    pm_qp*pp_k*cfsrd*vp122 - 32*pm_qp*pp_k*cfsrpd*mpi2*vp122 - 16
     &    *pm_qp*pp_k*cisrmd*mpi2*vp12*vp34 + 8*pm_qp*pp_k2*cisrpd*
     &    vp12*vp34 - 32*pm_qp*qm_qp*cisrmd*me2*vp12*vp34 - 16*pm_qp*
     &    qm_k*cisrmd*me2*vp12*vp34 + 16*pm_qp*qp_k*cisrmd*me2*vp12*
     &    vp34 + 32*pm_qp*cisrmd*me2*mpi2*vp12*vp34 + 32*pm_qp2*pp_qm
     &    *cisrmd*vp12*vp34 - 32*pm_qp2*pp_qp*cisrmd*vp12*vp34 - 16*
     &    pm_qp2*pp_k*cisrmd*vp12*vp34 - 16*pm_k*pp_qm*pp_qp*cisrpd*
     &    vp12*vp34 - 8*pm_k*pp_qm*pp_k*cisrpd*vp12*vp34 - 16*pm_k*
     &    pp_qm*qm_qp*cisrmd*vp12*vp34 + 32*pm_k*pp_qm*qp_k*cfsrpd*
     &    vp122 )
      interf = interf + cfsrp * ( 16*pm_k*pp_qm*qp_k*cisrpd*vp12*vp34
     &     + 4*pm_k*pp_qm*cfsrd*vp122 + 32*pm_k*pp_qm*cfsrpd*mpi2*vp122
     &     + 16*pm_k*pp_qm*cisrmd*mpi2*vp12*vp34 + 8*pm_k*pp_qp*pp_k*
     &    cisrpd*vp12*vp34 - 16*pm_k*pp_qp*qm_qp*cisrpd*vp12*vp34 - 16*
     &    pm_k*pp_qp*qm_k*cisrmd*vp12*vp34 - 32*pm_k*pp_qp*qp_k*cfsrpd*
     &    vp122 - 16*pm_k*pp_qp*qp_k*cisrpd*vp12*vp34 + 16*pm_k*pp_qp*
     &    qp_k*cisrmd*vp12*vp34 - 12*pm_k*pp_qp*cfsrd*vp122 - 32*pm_k*
     &    pp_qp*cfsrpd*mpi2*vp122 + 16*pm_k*pp_qp*cisrpd*mpi2*vp12*vp34
     &     + 16*pm_k*pp_qp2*cisrpd*vp12*vp34 - 32*pm_k*pp_k*qm_qp*
     &    cfsrmd*vp122 - 16*pm_k*pp_k*qm_qp*cisrpd*vp12*vp34 + 16*pm_k*
     &    pp_k*qm_qp*cisrmd*vp12*vp34 - 16*pm_k*pp_k*qm_k*cfsrmd*vp122
     &     - 32*pm_k*pp_k*qp_k*cfsrpd*vp122 - 16*pm_k*pp_k*qp_k*cfsrmd*
     &    vp122 - 8*pm_k*pp_k*cfsrd*vp122 - 32*pm_k*pp_k*cfsrpd*mpi2*
     &    vp122 + 16*pm_k*pp_k*cisrpd*mpi2*vp12*vp34 - 16*pm_k*pp_k*
     &    cisrmd*mpi2*vp12*vp34 - 16*pm_k*qm_qp*cisrmd*me2*vp12*vp34 - 
     &    8*pm_k*qm_k*cisrmd*me2*vp12*vp34 )
      interf = interf + cfsrp * ( 8*pm_k*qp_k*cisrmd*me2*vp12*vp34 + 16
     &    *pm_k*cisrmd*me2*mpi2*vp12*vp34 + 8*pm_k2*pp_qm*cisrmd*vp12
     &    *vp34 - 8*pm_k2*pp_qp*cisrmd*vp12*vp34 + 32*pp_qp*qm_qp*
     &    cisrpd*me2*vp12*vp34 + 16*pp_qp*qm_k*cisrpd*me2*vp12*vp34 - 
     &    16*pp_qp*qp_k*cisrpd*me2*vp12*vp34 - 32*pp_qp*cisrpd*me2*mpi2
     &    *vp12*vp34 + 16*pp_k*qm_qp*cisrpd*me2*vp12*vp34 + 8*pp_k*qm_k
     &    *cisrpd*me2*vp12*vp34 - 8*pp_k*qp_k*cisrpd*me2*vp12*vp34 - 16
     &    *pp_k*cisrpd*me2*mpi2*vp12*vp34 + 16*qm_qp*qm_k*cfsrmd*me2*
     &    vp122 - 32*qm_qp*qp_k*cfsrpd*me2*vp122 + 16*qm_qp*qp_k*cfsrmd
     &    *me2*vp122 - 16*qm_qp*qp_k*cisrpd*me2*vp12*vp34 + 16*qm_qp*
     &    qp_k*cisrmd*me2*vp12*vp34 - 8*qm_qp*cfsrd*me2*vp122 - 32*
     &    qm_qp*cfsrpd*me2*mpi2*vp122 - 32*qm_qp*cfsrmd*me2*mpi2*vp122
     &     + 32*qm_qp2*cfsrmd*me2*vp122 - 32*qm_k*qp_k*cfsrpd*me2*
     &    vp122 - 16*qm_k*qp_k*cisrpd*me2*vp12*vp34 + 16*qm_k*qp_k*
     &    cisrmd*me2*vp12*vp34 - 4*qm_k*cfsrd*me2*vp122 - 32*qm_k*
     &    cfsrpd*me2*mpi2*vp122 )
      interf = interf + cfsrp * (  - 16*qm_k*cfsrmd*me2*mpi2*vp122 + 12
     &    *qp_k*cfsrd*me2*vp122 + 64*qp_k*cfsrpd*me2*mpi2*vp122 - 16*
     &    qp_k*cfsrmd*me2*mpi2*vp122 + 16*qp_k*cisrpd*me2*mpi2*vp12*
     &    vp34 - 16*qp_k*cisrmd*me2*mpi2*vp12*vp34 + 32*qp_k2*cfsrpd*
     &    me2*vp122 + 16*qp_k2*cisrpd*me2*vp12*vp34 - 16*qp_k2*
     &    cisrmd*me2*vp12*vp34 )
      interf = interf + cfsrm * ( 8*cfsrd*me2*mpi2*vp122 + 32*cfsrmd*
     &    me2*mpi4*vp122 + 32*pm_pp*pm_qm*qm_qp*cisrmd*vp12*vp34 - 16*
     &    pm_pp*pm_qm*qm_k*cisrmd*vp12*vp34 + 16*pm_pp*pm_qm*qp_k*
     &    cisrmd*vp12*vp34 - 32*pm_pp*pm_qm*cisrmd*mpi2*vp12*vp34 + 16*
     &    pm_pp*pm_k*qm_qp*cisrmd*vp12*vp34 - 8*pm_pp*pm_k*qm_k*cisrmd*
     &    vp12*vp34 + 8*pm_pp*pm_k*qp_k*cisrmd*vp12*vp34 - 16*pm_pp*
     &    pm_k*cisrmd*mpi2*vp12*vp34 - 32*pm_pp*pp_qm*qm_qp*cisrpd*vp12
     &    *vp34 + 16*pm_pp*pp_qm*qm_k*cisrpd*vp12*vp34 - 16*pm_pp*pp_qm
     &    *qp_k*cisrpd*vp12*vp34 + 32*pm_pp*pp_qm*cisrpd*mpi2*vp12*vp34
     &     - 16*pm_pp*pp_k*qm_qp*cisrpd*vp12*vp34 + 8*pm_pp*pp_k*qm_k*
     &    cisrpd*vp12*vp34 - 8*pm_pp*pp_k*qp_k*cisrpd*vp12*vp34 + 16*
     &    pm_pp*pp_k*cisrpd*mpi2*vp12*vp34 + 16*pm_pp*qm_qp*qm_k*cfsrpd
     &    *vp122 - 32*pm_pp*qm_qp*qm_k*cfsrmd*vp122 + 16*pm_pp*qm_qp*
     &    qm_k*cisrpd*vp12*vp34 - 16*pm_pp*qm_qp*qm_k*cisrmd*vp12*vp34
     &     + 16*pm_pp*qm_qp*qp_k*cfsrpd*vp122 - 8*pm_pp*qm_qp*cfsrd*
     &    vp122 )
      interf = interf + cfsrm * (  - 32*pm_pp*qm_qp*cfsrpd*mpi2*vp122
     &     - 32*pm_pp*qm_qp*cfsrmd*mpi2*vp122 + 32*pm_pp*qm_qp2*
     &    cfsrpd*vp122 - 32*pm_pp*qm_k*qp_k*cfsrmd*vp122 + 16*pm_pp*
     &    qm_k*qp_k*cisrpd*vp12*vp34 - 16*pm_pp*qm_k*qp_k*cisrmd*vp12*
     &    vp34 + 12*pm_pp*qm_k*cfsrd*vp122 - 16*pm_pp*qm_k*cfsrpd*mpi2*
     &    vp122 + 64*pm_pp*qm_k*cfsrmd*mpi2*vp122 - 16*pm_pp*qm_k*
     &    cisrpd*mpi2*vp12*vp34 + 16*pm_pp*qm_k*cisrmd*mpi2*vp12*vp34
     &     + 32*pm_pp*qm_k2*cfsrmd*vp122 - 16*pm_pp*qm_k2*cisrpd*
     &    vp12*vp34 + 16*pm_pp*qm_k2*cisrmd*vp12*vp34 - 4*pm_pp*qp_k*
     &    cfsrd*vp122 - 16*pm_pp*qp_k*cfsrpd*mpi2*vp122 - 32*pm_pp*qp_k
     &    *cfsrmd*mpi2*vp122 + 8*pm_pp*cfsrd*mpi2*vp122 + 32*pm_pp*
     &    cfsrmd*mpi4*vp122 - 32*pm_qm*pm_qp*pp_qm*cisrmd*vp12*vp34 + 
     &    32*pm_qm*pm_qp*pp_qp*cisrmd*vp12*vp34 - 16*pm_qm*pm_qp*pp_k*
     &    cisrmd*vp12*vp34 + 32*pm_qm*pm_k*pp_qm*cisrmd*vp12*vp34 - 32*
     &    pm_qm*pm_k*pp_qp*cisrmd*vp12*vp34 + 8*pm_qm*pm_k*pp_k*cisrmd*
     &    vp12*vp34 )
      interf = interf + cfsrm * ( 32*pm_qm*pp_qm*pp_qp*cisrpd*vp12*vp34
     &     - 32*pm_qm*pp_qm*pp_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qm*
     &    qm_qp*cfsrpd*vp122 + 16*pm_qm*pp_qm*qm_k*cfsrpd*vp122 - 32*
     &    pm_qm*pp_qm*qm_k*cfsrmd*vp122 + 32*pm_qm*pp_qm*qm_k*cisrpd*
     &    vp12*vp34 - 32*pm_qm*pp_qm*qm_k*cisrmd*vp12*vp34 + 16*pm_qm*
     &    pp_qm*qp_k*cfsrpd*vp122 - 16*pm_qm*pp_qm*qp_k*cisrpd*vp12*
     &    vp34 + 16*pm_qm*pp_qm*qp_k*cisrmd*vp12*vp34 - 16*pm_qm*pp_qm*
     &    cfsrd*vp122 - 32*pm_qm*pp_qm*cfsrmd*mpi2*vp122 - 32*pm_qm*
     &    pp_qm2*cisrpd*vp12*vp34 + 16*pm_qm*pp_qp*pp_k*cisrpd*vp12*
     &    vp34 - 32*pm_qm*pp_qp*qm_qp*cfsrpd*vp122 - 16*pm_qm*pp_qp*
     &    qm_k*cfsrpd*vp122 + 32*pm_qm*pp_qp*qm_k*cfsrmd*vp122 - 16*
     &    pm_qm*pp_qp*qm_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qp*qm_k*
     &    cisrmd*vp12*vp34 - 16*pm_qm*pp_qp*qp_k*cfsrpd*vp122 - 16*
     &    pm_qm*pp_qp*qp_k*cisrmd*vp12*vp34 + 8*pm_qm*pp_qp*cfsrd*vp122
     &     + 32*pm_qm*pp_qp*cfsrmd*mpi2*vp122 - 16*pm_qm*pp_k*qm_qp*
     &    cisrmd*vp12*vp34 )
      interf = interf + cfsrm * (  - 32*pm_qm*pp_k*qm_k*cfsrmd*vp122 + 
     &    16*pm_qm*pp_k*qm_k*cisrpd*vp12*vp34 - 16*pm_qm*pp_k*qm_k*
     &    cisrmd*vp12*vp34 - 16*pm_qm*pp_k*qp_k*cisrpd*vp12*vp34 - 12*
     &    pm_qm*pp_k*cfsrd*vp122 - 32*pm_qm*pp_k*cfsrmd*mpi2*vp122 + 16
     &    *pm_qm*pp_k*cisrmd*mpi2*vp12*vp34 - 8*pm_qm*pp_k2*cisrpd*
     &    vp12*vp34 + 32*pm_qm*qm_qp*cisrmd*me2*vp12*vp34 - 16*pm_qm*
     &    qm_k*cisrmd*me2*vp12*vp34 + 16*pm_qm*qp_k*cisrmd*me2*vp12*
     &    vp34 - 32*pm_qm*cisrmd*me2*mpi2*vp12*vp34 + 32*pm_qm2*pp_qm
     &    *cisrmd*vp12*vp34 - 32*pm_qm2*pp_qp*cisrmd*vp12*vp34 + 16*
     &    pm_qm2*pp_k*cisrmd*vp12*vp34 - 16*pm_qp*pm_k*pp_qm*cisrmd*
     &    vp12*vp34 + 16*pm_qp*pm_k*pp_qp*cisrmd*vp12*vp34 - 8*pm_qp*
     &    pm_k*pp_k*cisrmd*vp12*vp34 - 32*pm_qp*pp_qm*pp_qp*cisrpd*vp12
     &    *vp34 + 32*pm_qp*pp_qm*pp_k*cisrpd*vp12*vp34 - 32*pm_qp*pp_qm
     &    *qm_qp*cfsrpd*vp122 - 16*pm_qp*pp_qm*qm_k*cfsrpd*vp122 + 32*
     &    pm_qp*pp_qm*qm_k*cfsrmd*vp122 - 32*pm_qp*pp_qm*qm_k*cisrpd*
     &    vp12*vp34 )
      interf = interf + cfsrm * ( 16*pm_qp*pp_qm*qm_k*cisrmd*vp12*vp34
     &     - 16*pm_qp*pp_qm*qp_k*cfsrpd*vp122 + 16*pm_qp*pp_qm*qp_k*
     &    cisrpd*vp12*vp34 + 8*pm_qp*pp_qm*cfsrd*vp122 + 32*pm_qp*pp_qm
     &    *cfsrmd*mpi2*vp122 + 32*pm_qp*pp_qm2*cisrpd*vp12*vp34 - 16*
     &    pm_qp*pp_qp*pp_k*cisrpd*vp12*vp34 + 32*pm_qp*pp_qp*qm_qp*
     &    cfsrpd*vp122 + 16*pm_qp*pp_qp*qm_k*cfsrpd*vp122 - 32*pm_qp*
     &    pp_qp*qm_k*cfsrmd*vp122 + 16*pm_qp*pp_qp*qm_k*cisrpd*vp12*
     &    vp34 - 16*pm_qp*pp_qp*qm_k*cisrmd*vp12*vp34 + 16*pm_qp*pp_qp*
     &    qp_k*cfsrpd*vp122 - 32*pm_qp*pp_qp*cfsrmd*mpi2*vp122 - 16*
     &    pm_qp*pp_k*qm_qp*cisrpd*vp12*vp34 + 32*pm_qp*pp_k*qm_k*cfsrmd
     &    *vp122 + 16*pm_qp*pp_k*qm_k*cisrmd*vp12*vp34 + 4*pm_qp*pp_k*
     &    cfsrd*vp122 + 32*pm_qp*pp_k*cfsrmd*mpi2*vp122 + 16*pm_qp*pp_k
     &    *cisrpd*mpi2*vp12*vp34 + 8*pm_qp*pp_k2*cisrpd*vp12*vp34 + 
     &    16*pm_k*pp_qm*pp_qp*cisrpd*vp12*vp34 - 8*pm_k*pp_qm*pp_k*
     &    cisrpd*vp12*vp34 + 16*pm_k*pp_qm*qm_qp*cisrpd*vp12*vp34 - 32*
     &    pm_k*pp_qm*qm_k*cfsrmd*vp122 )
      interf = interf + cfsrm * ( 16*pm_k*pp_qm*qm_k*cisrpd*vp12*vp34
     &     - 16*pm_k*pp_qm*qm_k*cisrmd*vp12*vp34 + 16*pm_k*pp_qm*qp_k*
     &    cisrmd*vp12*vp34 - 12*pm_k*pp_qm*cfsrd*vp122 - 32*pm_k*pp_qm*
     &    cfsrmd*mpi2*vp122 - 16*pm_k*pp_qm*cisrpd*mpi2*vp12*vp34 - 16*
     &    pm_k*pp_qm2*cisrpd*vp12*vp34 + 8*pm_k*pp_qp*pp_k*cisrpd*
     &    vp12*vp34 + 16*pm_k*pp_qp*qm_qp*cisrmd*vp12*vp34 + 32*pm_k*
     &    pp_qp*qm_k*cfsrmd*vp122 - 16*pm_k*pp_qp*qm_k*cisrpd*vp12*vp34
     &     + 4*pm_k*pp_qp*cfsrd*vp122 + 32*pm_k*pp_qp*cfsrmd*mpi2*vp122
     &     - 16*pm_k*pp_qp*cisrmd*mpi2*vp12*vp34 - 32*pm_k*pp_k*qm_qp*
     &    cfsrpd*vp122 + 16*pm_k*pp_k*qm_qp*cisrpd*vp12*vp34 - 16*pm_k*
     &    pp_k*qm_qp*cisrmd*vp12*vp34 - 16*pm_k*pp_k*qm_k*cfsrpd*vp122
     &     - 32*pm_k*pp_k*qm_k*cfsrmd*vp122 - 16*pm_k*pp_k*qp_k*cfsrpd*
     &    vp122 - 8*pm_k*pp_k*cfsrd*vp122 - 32*pm_k*pp_k*cfsrmd*mpi2*
     &    vp122 - 16*pm_k*pp_k*cisrpd*mpi2*vp12*vp34 + 16*pm_k*pp_k*
     &    cisrmd*mpi2*vp12*vp34 + 16*pm_k*qm_qp*cisrmd*me2*vp12*vp34 - 
     &    8*pm_k*qm_k*cisrmd*me2*vp12*vp34 )
      interf = interf + cfsrm * ( 8*pm_k*qp_k*cisrmd*me2*vp12*vp34 - 16
     &    *pm_k*cisrmd*me2*mpi2*vp12*vp34 + 8*pm_k2*pp_qm*cisrmd*vp12
     &    *vp34 - 8*pm_k2*pp_qp*cisrmd*vp12*vp34 - 32*pp_qm*qm_qp*
     &    cisrpd*me2*vp12*vp34 + 16*pp_qm*qm_k*cisrpd*me2*vp12*vp34 - 
     &    16*pp_qm*qp_k*cisrpd*me2*vp12*vp34 + 32*pp_qm*cisrpd*me2*mpi2
     &    *vp12*vp34 - 16*pp_k*qm_qp*cisrpd*me2*vp12*vp34 + 8*pp_k*qm_k
     &    *cisrpd*me2*vp12*vp34 - 8*pp_k*qp_k*cisrpd*me2*vp12*vp34 + 16
     &    *pp_k*cisrpd*me2*mpi2*vp12*vp34 + 16*qm_qp*qm_k*cfsrpd*me2*
     &    vp122 - 32*qm_qp*qm_k*cfsrmd*me2*vp122 + 16*qm_qp*qm_k*cisrpd
     &    *me2*vp12*vp34 - 16*qm_qp*qm_k*cisrmd*me2*vp12*vp34 + 16*
     &    qm_qp*qp_k*cfsrpd*me2*vp122 - 8*qm_qp*cfsrd*me2*vp122 - 32*
     &    qm_qp*cfsrpd*me2*mpi2*vp122 - 32*qm_qp*cfsrmd*me2*mpi2*vp122
     &     + 32*qm_qp2*cfsrpd*me2*vp122 - 32*qm_k*qp_k*cfsrmd*me2*
     &    vp122 + 16*qm_k*qp_k*cisrpd*me2*vp12*vp34 - 16*qm_k*qp_k*
     &    cisrmd*me2*vp12*vp34 + 12*qm_k*cfsrd*me2*vp122 - 16*qm_k*
     &    cfsrpd*me2*mpi2*vp122 )
      interf = interf + cfsrm * ( 64*qm_k*cfsrmd*me2*mpi2*vp122 - 16*
     &    qm_k*cisrpd*me2*mpi2*vp12*vp34 + 16*qm_k*cisrmd*me2*mpi2*vp12
     &    *vp34 + 32*qm_k2*cfsrmd*me2*vp122 - 16*qm_k2*cisrpd*me2*
     &    vp12*vp34 + 16*qm_k2*cisrmd*me2*vp12*vp34 - 4*qp_k*cfsrd*
     &    me2*vp122 - 16*qp_k*cfsrpd*me2*mpi2*vp122 - 32*qp_k*cfsrmd*
     &    me2*mpi2*vp122 )
      interf = interf + cisrp * ( 32*cisrpd*me4*mpi2*vp342 + 32*pm_pp*
     &    pm_qm*pp_qm*cisrmd*vp342 - 32*pm_pp*pm_qm*pp_qp*cisrmd*vp342
     &     - 16*pm_pp*pm_qm*qm_k*cisrmd*vp342 + 16*pm_pp*pm_qm*qp_k*
     &    cisrmd*vp342 - 32*pm_pp*pm_qp*pp_qm*cisrmd*vp342 + 32*pm_pp*
     &    pm_qp*pp_qp*cisrmd*vp342 + 16*pm_pp*pm_qp*qm_k*cisrmd*vp342
     &     - 16*pm_pp*pm_qp*qp_k*cisrmd*vp342 - 32*pm_pp*pm_k*qm_qp*
     &    cisrmd*vp342 + 32*pm_pp*pm_k*cisrmd*mpi2*vp342 - 32*pm_pp*
     &    pp_qm*qm_qp*cfsrmd*vp12*vp34 + 16*pm_pp*pp_qm*qm_k*cfsrmd*
     &    vp12*vp34 - 16*pm_pp*pp_qm*qm_k*cisrmd*vp342 - 16*pm_pp*pp_qm
     &    *qp_k*cfsrmd*vp12*vp34 + 16*pm_pp*pp_qm*qp_k*cisrmd*vp342 + 
     &    32*pm_pp*pp_qm*cfsrmd*mpi2*vp12*vp34 + 32*pm_pp*pp_qp*qm_qp*
     &    cfsrpd*vp12*vp34 + 16*pm_pp*pp_qp*qm_k*cfsrpd*vp12*vp34 + 16*
     &    pm_pp*pp_qp*qm_k*cisrmd*vp342 - 16*pm_pp*pp_qp*qp_k*cfsrpd*
     &    vp12*vp34 - 16*pm_pp*pp_qp*qp_k*cisrmd*vp342 - 32*pm_pp*pp_qp
     &    *cfsrpd*mpi2*vp12*vp34 + 16*pm_pp*pp_k*qm_qp*cfsrpd*vp12*vp34
     &     - 16*pm_pp*pp_k*qm_qp*cfsrmd*vp12*vp34 )
      interf = interf + cisrp * (  - 32*pm_pp*pp_k*qm_qp*cisrmd*vp342
     &     + 8*pm_pp*pp_k*qm_k*cfsrpd*vp12*vp34 + 8*pm_pp*pp_k*qm_k*
     &    cfsrmd*vp12*vp34 - 8*pm_pp*pp_k*qp_k*cfsrpd*vp12*vp34 - 8*
     &    pm_pp*pp_k*qp_k*cfsrmd*vp12*vp34 - 16*pm_pp*pp_k*cfsrpd*mpi2*
     &    vp12*vp34 + 16*pm_pp*pp_k*cfsrmd*mpi2*vp12*vp34 + 32*pm_pp*
     &    pp_k*cisrmd*mpi2*vp342 + 16*pm_pp*qm_qp*qm_k*cfsrmd*vp12*vp34
     &     - 16*pm_pp*qm_qp*qp_k*cfsrpd*vp12*vp34 - 32*pm_pp*qm_qp*
     &    cisrpd*me2*vp342 + 32*pm_pp*qm_qp*cisrmd*me2*vp342 - 16*pm_pp
     &    *qm_k*qp_k*cfsrpd*vp12*vp34 + 16*pm_pp*qm_k*qp_k*cfsrmd*vp12*
     &    vp34 - 8*pm_pp*qm_k*cfsrd*vp12*vp34 - 16*pm_pp*qm_k*cfsrmd*
     &    mpi2*vp12*vp34 - 16*pm_pp*qm_k2*cfsrmd*vp12*vp34 + 8*pm_pp*
     &    qp_k*cfsrd*vp12*vp34 + 16*pm_pp*qp_k*cfsrpd*mpi2*vp12*vp34 + 
     &    16*pm_pp*qp_k2*cfsrpd*vp12*vp34 + 32*pm_pp*cisrpd*me2*mpi2*
     &    vp342 - 32*pm_pp*cisrmd*me2*mpi2*vp342 + 32*pm_pp2*qm_qp*
     &    cisrmd*vp342 - 32*pm_pp2*cisrmd*mpi2*vp342 - 32*pm_qm*pm_qp
     &    *pp_k*cisrmd*vp342 )
      interf = interf + cisrp * (  - 16*pm_qm*pm_k*pp_qm*cisrmd*vp342
     &     + 16*pm_qm*pm_k*pp_qp*cisrmd*vp342 + 32*pm_qm*pp_qm*pp_qp*
     &    cfsrpd*vp12*vp34 + 32*pm_qm*pp_qm*pp_qp*cfsrmd*vp12*vp34 + 16
     &    *pm_qm*pp_qm*pp_k*cfsrpd*vp12*vp34 - 32*pm_qm*pp_qm*pp_k*
     &    cfsrmd*vp12*vp34 - 16*pm_qm*pp_qm*pp_k*cisrmd*vp342 + 32*
     &    pm_qm*pp_qm*qm_k*cfsrmd*vp12*vp34 - 16*pm_qm*pp_qm*qp_k*
     &    cfsrpd*vp12*vp34 - 16*pm_qm*pp_qm*qp_k*cfsrmd*vp12*vp34 - 32*
     &    pm_qm*pp_qm*cisrpd*me2*vp342 - 32*pm_qm*pp_qm2*cfsrmd*vp12*
     &    vp34 - 32*pm_qm*pp_qp*pp_k*cfsrpd*vp12*vp34 + 16*pm_qm*pp_qp*
     &    pp_k*cfsrmd*vp12*vp34 + 16*pm_qm*pp_qp*pp_k*cisrmd*vp342 - 16
     &    *pm_qm*pp_qp*qm_k*cfsrpd*vp12*vp34 - 16*pm_qm*pp_qp*qm_k*
     &    cfsrmd*vp12*vp34 + 32*pm_qm*pp_qp*qp_k*cfsrpd*vp12*vp34 + 32*
     &    pm_qm*pp_qp*cisrpd*me2*vp342 - 32*pm_qm*pp_qp2*cfsrpd*vp12*
     &    vp34 + 16*pm_qm*pp_k*qm_qp*cfsrpd*vp12*vp34 + 16*pm_qm*pp_k*
     &    qm_k*cfsrmd*vp12*vp34 + 32*pm_qm*pp_k*qm_k*cisrpd*vp342 - 16*
     &    pm_qm*pp_k*qp_k*cfsrmd*vp12*vp34 )
      interf = interf + cisrp * (  - 32*pm_qm*pp_k*qp_k*cisrpd*vp342 - 
     &    8*pm_qm*pp_k*cfsrd*vp12*vp34 - 16*pm_qm*pp_k*cfsrpd*mpi2*vp12
     &    *vp34 - 8*pm_qm*pp_k2*cfsrpd*vp12*vp34 - 8*pm_qm*pp_k2*
     &    cfsrmd*vp12*vp34 + 32*pm_qm*qm_k*cisrpd*me2*vp342 - 32*pm_qm*
     &    qp_k*cisrpd*me2*vp342 - 8*pm_qm*cfsrd*me2*vp12*vp34 + 16*
     &    pm_qm2*pp_k*cisrmd*vp342 + 16*pm_qp*pm_k*pp_qm*cisrmd*vp342
     &     - 16*pm_qp*pm_k*pp_qp*cisrmd*vp342 - 32*pm_qp*pp_qm*pp_qp*
     &    cfsrpd*vp12*vp34 - 32*pm_qp*pp_qm*pp_qp*cfsrmd*vp12*vp34 - 16
     &    *pm_qp*pp_qm*pp_k*cfsrpd*vp12*vp34 + 32*pm_qp*pp_qm*pp_k*
     &    cfsrmd*vp12*vp34 + 16*pm_qp*pp_qm*pp_k*cisrmd*vp342 - 32*
     &    pm_qp*pp_qm*qm_k*cfsrmd*vp12*vp34 + 16*pm_qp*pp_qm*qp_k*
     &    cfsrpd*vp12*vp34 + 16*pm_qp*pp_qm*qp_k*cfsrmd*vp12*vp34 + 32*
     &    pm_qp*pp_qm*cisrpd*me2*vp342 + 32*pm_qp*pp_qm2*cfsrmd*vp12*
     &    vp34 + 32*pm_qp*pp_qp*pp_k*cfsrpd*vp12*vp34 - 16*pm_qp*pp_qp*
     &    pp_k*cfsrmd*vp12*vp34 - 16*pm_qp*pp_qp*pp_k*cisrmd*vp342 + 16
     &    *pm_qp*pp_qp*qm_k*cfsrpd*vp12*vp34 )
      interf = interf + cisrp * ( 16*pm_qp*pp_qp*qm_k*cfsrmd*vp12*vp34
     &     - 32*pm_qp*pp_qp*qp_k*cfsrpd*vp12*vp34 - 32*pm_qp*pp_qp*
     &    cisrpd*me2*vp342 + 32*pm_qp*pp_qp2*cfsrpd*vp12*vp34 - 16*
     &    pm_qp*pp_k*qm_qp*cfsrmd*vp12*vp34 + 16*pm_qp*pp_k*qm_k*cfsrpd
     &    *vp12*vp34 - 32*pm_qp*pp_k*qm_k*cisrpd*vp342 - 16*pm_qp*pp_k*
     &    qp_k*cfsrpd*vp12*vp34 + 32*pm_qp*pp_k*qp_k*cisrpd*vp342 + 8*
     &    pm_qp*pp_k*cfsrd*vp12*vp34 + 16*pm_qp*pp_k*cfsrmd*mpi2*vp12*
     &    vp34 + 8*pm_qp*pp_k2*cfsrpd*vp12*vp34 + 8*pm_qp*pp_k2*
     &    cfsrmd*vp12*vp34 - 32*pm_qp*qm_k*cisrpd*me2*vp342 + 32*pm_qp*
     &    qp_k*cisrpd*me2*vp342 + 8*pm_qp*cfsrd*me2*vp12*vp34 + 16*
     &    pm_qp2*pp_k*cisrmd*vp342 - 16*pm_k*pp_qm*pp_qp*cfsrpd*vp12*
     &    vp34 + 16*pm_k*pp_qm*pp_qp*cfsrmd*vp12*vp34 - 32*pm_k*pp_qm*
     &    pp_qp*cisrmd*vp342 - 8*pm_k*pp_qm*pp_k*cfsrpd*vp12*vp34 - 8*
     &    pm_k*pp_qm*pp_k*cfsrmd*vp12*vp34 + 16*pm_k*pp_qm*qm_qp*cfsrmd
     &    *vp12*vp34 + 16*pm_k*pp_qm*qm_k*cfsrmd*vp12*vp34 + 16*pm_k*
     &    pp_qm*qp_k*cfsrpd*vp12*vp34 )
      interf = interf + cisrp * ( 8*pm_k*pp_qm*cfsrd*vp12*vp34 - 16*
     &    pm_k*pp_qm*cfsrmd*mpi2*vp12*vp34 - 16*pm_k*pp_qm2*cfsrmd*
     &    vp12*vp34 + 16*pm_k*pp_qm2*cisrmd*vp342 + 8*pm_k*pp_qp*pp_k
     &    *cfsrpd*vp12*vp34 + 8*pm_k*pp_qp*pp_k*cfsrmd*vp12*vp34 - 16*
     &    pm_k*pp_qp*qm_qp*cfsrpd*vp12*vp34 - 16*pm_k*pp_qp*qm_k*cfsrmd
     &    *vp12*vp34 - 16*pm_k*pp_qp*qp_k*cfsrpd*vp12*vp34 - 8*pm_k*
     &    pp_qp*cfsrd*vp12*vp34 + 16*pm_k*pp_qp*cfsrpd*mpi2*vp12*vp34
     &     + 16*pm_k*pp_qp2*cfsrpd*vp12*vp34 + 16*pm_k*pp_qp2*
     &    cisrmd*vp342 - 16*pm_k*pp_k*qm_qp*cfsrpd*vp12*vp34 + 16*pm_k*
     &    pp_k*qm_qp*cfsrmd*vp12*vp34 + 32*pm_k*pp_k*qm_qp*cisrpd*vp342
     &     + 16*pm_k*pp_k*cfsrpd*mpi2*vp12*vp34 - 16*pm_k*pp_k*cfsrmd*
     &    mpi2*vp12*vp34 - 32*pm_k*pp_k*cisrpd*mpi2*vp342 + 32*pm_k*
     &    qm_qp*cisrpd*me2*vp342 - 32*pm_k*cisrpd*me2*mpi2*vp342 - 32*
     &    pp_qm*qm_qp*cfsrmd*me2*vp12*vp34 + 16*pp_qm*qm_k*cfsrmd*me2*
     &    vp12*vp34 - 16*pp_qm*qp_k*cfsrmd*me2*vp12*vp34 + 8*pp_qm*
     &    cfsrd*me2*vp12*vp34 )
      interf = interf + cisrp * ( 32*pp_qm*cfsrmd*me2*mpi2*vp12*vp34 + 
     &    32*pp_qp*qm_qp*cfsrpd*me2*vp12*vp34 + 16*pp_qp*qm_k*cfsrpd*
     &    me2*vp12*vp34 - 16*pp_qp*qp_k*cfsrpd*me2*vp12*vp34 - 8*pp_qp*
     &    cfsrd*me2*vp12*vp34 - 32*pp_qp*cfsrpd*me2*mpi2*vp12*vp34 + 16
     &    *pp_k*qm_qp*cfsrpd*me2*vp12*vp34 - 16*pp_k*qm_qp*cfsrmd*me2*
     &    vp12*vp34 + 32*pp_k*qm_qp*cisrpd*me2*vp342 + 8*pp_k*qm_k*
     &    cfsrpd*me2*vp12*vp34 + 8*pp_k*qm_k*cfsrmd*me2*vp12*vp34 - 8*
     &    pp_k*qp_k*cfsrpd*me2*vp12*vp34 - 8*pp_k*qp_k*cfsrmd*me2*vp12*
     &    vp34 - 16*pp_k*cfsrpd*me2*mpi2*vp12*vp34 + 16*pp_k*cfsrmd*me2
     &    *mpi2*vp12*vp34 - 32*pp_k*cisrpd*me2*mpi2*vp342 + 16*qm_qp*
     &    qm_k*cfsrmd*me2*vp12*vp34 - 16*qm_qp*qp_k*cfsrpd*me2*vp12*
     &    vp34 - 32*qm_qp*cisrpd*me4*vp342 - 16*qm_k*qp_k*cfsrpd*me2*
     &    vp12*vp34 + 16*qm_k*qp_k*cfsrmd*me2*vp12*vp34 + 32*qm_k*qp_k*
     &    cisrmd*me2*vp342 - 16*qm_k*cfsrd*me2*vp12*vp34 - 16*qm_k*
     &    cfsrmd*me2*mpi2*vp12*vp34 - 16*qm_k2*cfsrmd*me2*vp12*vp34
     &     - 16*qm_k2*cisrmd*me2*vp342 )
      interf = interf + cisrp * ( 16*qp_k*cfsrd*me2*vp12*vp34 + 16*qp_k
     &    *cfsrpd*me2*mpi2*vp12*vp34 + 16*qp_k2*cfsrpd*me2*vp12*vp34
     &     - 16*qp_k2*cisrmd*me2*vp342 )
      interf = interf + cisrm * ( 32*cisrmd*me4*mpi2*vp342 + 32*pm_pp*
     &    pm_qm*pp_qm*cisrpd*vp342 - 32*pm_pp*pm_qm*pp_qp*cisrpd*vp342
     &     + 32*pm_pp*pm_qm*qm_qp*cfsrmd*vp12*vp34 - 16*pm_pp*pm_qm*
     &    qm_k*cfsrmd*vp12*vp34 - 16*pm_pp*pm_qm*qm_k*cisrpd*vp342 + 16
     &    *pm_pp*pm_qm*qp_k*cfsrmd*vp12*vp34 + 16*pm_pp*pm_qm*qp_k*
     &    cisrpd*vp342 - 32*pm_pp*pm_qm*cfsrmd*mpi2*vp12*vp34 - 32*
     &    pm_pp*pm_qp*pp_qm*cisrpd*vp342 + 32*pm_pp*pm_qp*pp_qp*cisrpd*
     &    vp342 - 32*pm_pp*pm_qp*qm_qp*cfsrpd*vp12*vp34 - 16*pm_pp*
     &    pm_qp*qm_k*cfsrpd*vp12*vp34 + 16*pm_pp*pm_qp*qm_k*cisrpd*
     &    vp342 + 16*pm_pp*pm_qp*qp_k*cfsrpd*vp12*vp34 - 16*pm_pp*pm_qp
     &    *qp_k*cisrpd*vp342 + 32*pm_pp*pm_qp*cfsrpd*mpi2*vp12*vp34 - 
     &    16*pm_pp*pm_k*qm_qp*cfsrpd*vp12*vp34 + 16*pm_pp*pm_k*qm_qp*
     &    cfsrmd*vp12*vp34 - 32*pm_pp*pm_k*qm_qp*cisrpd*vp342 - 8*pm_pp
     &    *pm_k*qm_k*cfsrpd*vp12*vp34 - 8*pm_pp*pm_k*qm_k*cfsrmd*vp12*
     &    vp34 + 8*pm_pp*pm_k*qp_k*cfsrpd*vp12*vp34 + 8*pm_pp*pm_k*qp_k
     &    *cfsrmd*vp12*vp34 )
      interf = interf + cisrm * ( 16*pm_pp*pm_k*cfsrpd*mpi2*vp12*vp34
     &     - 16*pm_pp*pm_k*cfsrmd*mpi2*vp12*vp34 + 32*pm_pp*pm_k*cisrpd
     &    *mpi2*vp342 - 16*pm_pp*pp_qm*qm_k*cisrpd*vp342 + 16*pm_pp*
     &    pp_qm*qp_k*cisrpd*vp342 + 16*pm_pp*pp_qp*qm_k*cisrpd*vp342 - 
     &    16*pm_pp*pp_qp*qp_k*cisrpd*vp342 - 32*pm_pp*pp_k*qm_qp*cisrpd
     &    *vp342 + 32*pm_pp*pp_k*cisrpd*mpi2*vp342 - 16*pm_pp*qm_qp*
     &    qm_k*cfsrmd*vp12*vp34 + 16*pm_pp*qm_qp*qp_k*cfsrpd*vp12*vp34
     &     + 32*pm_pp*qm_qp*cisrpd*me2*vp342 - 32*pm_pp*qm_qp*cisrmd*
     &    me2*vp342 + 16*pm_pp*qm_k*qp_k*cfsrpd*vp12*vp34 - 16*pm_pp*
     &    qm_k*qp_k*cfsrmd*vp12*vp34 + 8*pm_pp*qm_k*cfsrd*vp12*vp34 + 
     &    16*pm_pp*qm_k*cfsrmd*mpi2*vp12*vp34 + 16*pm_pp*qm_k2*cfsrmd
     &    *vp12*vp34 - 8*pm_pp*qp_k*cfsrd*vp12*vp34 - 16*pm_pp*qp_k*
     &    cfsrpd*mpi2*vp12*vp34 - 16*pm_pp*qp_k2*cfsrpd*vp12*vp34 - 
     &    32*pm_pp*cisrpd*me2*mpi2*vp342 + 32*pm_pp*cisrmd*me2*mpi2*
     &    vp342 + 32*pm_pp2*qm_qp*cisrpd*vp342 - 32*pm_pp2*cisrpd*
     &    mpi2*vp342 )
      interf = interf + cisrm * (  - 32*pm_qm*pm_qp*pp_qm*cfsrpd*vp12*
     &    vp34 - 32*pm_qm*pm_qp*pp_qm*cfsrmd*vp12*vp34 + 32*pm_qm*pm_qp
     &    *pp_qp*cfsrpd*vp12*vp34 + 32*pm_qm*pm_qp*pp_qp*cfsrmd*vp12*
     &    vp34 + 16*pm_qm*pm_qp*pp_k*cfsrpd*vp12*vp34 - 16*pm_qm*pm_qp*
     &    pp_k*cfsrmd*vp12*vp34 - 32*pm_qm*pm_qp*pp_k*cisrpd*vp342 - 16
     &    *pm_qm*pm_k*pp_qm*cfsrpd*vp12*vp34 + 32*pm_qm*pm_k*pp_qm*
     &    cfsrmd*vp12*vp34 - 16*pm_qm*pm_k*pp_qm*cisrpd*vp342 + 16*
     &    pm_qm*pm_k*pp_qp*cfsrpd*vp12*vp34 - 32*pm_qm*pm_k*pp_qp*
     &    cfsrmd*vp12*vp34 + 16*pm_qm*pm_k*pp_qp*cisrpd*vp342 + 8*pm_qm
     &    *pm_k*pp_k*cfsrpd*vp12*vp34 + 8*pm_qm*pm_k*pp_k*cfsrmd*vp12*
     &    vp34 - 16*pm_qm*pp_qm*pp_k*cisrpd*vp342 - 32*pm_qm*pp_qm*qm_k
     &    *cfsrmd*vp12*vp34 + 16*pm_qm*pp_qm*qp_k*cfsrpd*vp12*vp34 + 16
     &    *pm_qm*pp_qm*qp_k*cfsrmd*vp12*vp34 - 32*pm_qm*pp_qm*cisrmd*
     &    me2*vp342 + 16*pm_qm*pp_qp*pp_k*cisrpd*vp342 + 32*pm_qm*pp_qp
     &    *qm_k*cfsrmd*vp12*vp34 - 16*pm_qm*pp_qp*qp_k*cfsrpd*vp12*vp34
     &     - 16*pm_qm*pp_qp*qp_k*cfsrmd*vp12*vp34 )
      interf = interf + cisrm * ( 32*pm_qm*pp_qp*cisrmd*me2*vp342 - 16*
     &    pm_qm*pp_k*qm_qp*cfsrmd*vp12*vp34 - 16*pm_qm*pp_k*qm_k*cfsrmd
     &    *vp12*vp34 - 16*pm_qm*pp_k*qp_k*cfsrpd*vp12*vp34 - 8*pm_qm*
     &    pp_k*cfsrd*vp12*vp34 + 16*pm_qm*pp_k*cfsrmd*mpi2*vp12*vp34 + 
     &    32*pm_qm*qm_qp*cfsrmd*me2*vp12*vp34 - 16*pm_qm*qm_k*cfsrmd*
     &    me2*vp12*vp34 + 16*pm_qm*qp_k*cfsrmd*me2*vp12*vp34 - 8*pm_qm*
     &    cfsrd*me2*vp12*vp34 - 32*pm_qm*cfsrmd*me2*mpi2*vp12*vp34 + 32
     &    *pm_qm2*pp_qm*cfsrmd*vp12*vp34 - 32*pm_qm2*pp_qp*cfsrmd*
     &    vp12*vp34 + 16*pm_qm2*pp_k*cfsrmd*vp12*vp34 + 16*pm_qm2*
     &    pp_k*cisrpd*vp342 + 32*pm_qp*pm_k*pp_qm*cfsrpd*vp12*vp34 - 16
     &    *pm_qp*pm_k*pp_qm*cfsrmd*vp12*vp34 + 16*pm_qp*pm_k*pp_qm*
     &    cisrpd*vp342 - 32*pm_qp*pm_k*pp_qp*cfsrpd*vp12*vp34 + 16*
     &    pm_qp*pm_k*pp_qp*cfsrmd*vp12*vp34 - 16*pm_qp*pm_k*pp_qp*
     &    cisrpd*vp342 - 8*pm_qp*pm_k*pp_k*cfsrpd*vp12*vp34 - 8*pm_qp*
     &    pm_k*pp_k*cfsrmd*vp12*vp34 + 16*pm_qp*pp_qm*pp_k*cisrpd*vp342
     &     + 16*pm_qp*pp_qm*qm_k*cfsrpd*vp12*vp34 )
      interf = interf + cisrm * ( 16*pm_qp*pp_qm*qm_k*cfsrmd*vp12*vp34
     &     - 32*pm_qp*pp_qm*qp_k*cfsrpd*vp12*vp34 + 32*pm_qp*pp_qm*
     &    cisrmd*me2*vp342 - 16*pm_qp*pp_qp*pp_k*cisrpd*vp342 - 16*
     &    pm_qp*pp_qp*qm_k*cfsrpd*vp12*vp34 - 16*pm_qp*pp_qp*qm_k*
     &    cfsrmd*vp12*vp34 + 32*pm_qp*pp_qp*qp_k*cfsrpd*vp12*vp34 - 32*
     &    pm_qp*pp_qp*cisrmd*me2*vp342 + 16*pm_qp*pp_k*qm_qp*cfsrpd*
     &    vp12*vp34 + 16*pm_qp*pp_k*qm_k*cfsrmd*vp12*vp34 + 16*pm_qp*
     &    pp_k*qp_k*cfsrpd*vp12*vp34 + 8*pm_qp*pp_k*cfsrd*vp12*vp34 - 
     &    16*pm_qp*pp_k*cfsrpd*mpi2*vp12*vp34 - 32*pm_qp*qm_qp*cfsrpd*
     &    me2*vp12*vp34 - 16*pm_qp*qm_k*cfsrpd*me2*vp12*vp34 + 16*pm_qp
     &    *qp_k*cfsrpd*me2*vp12*vp34 + 8*pm_qp*cfsrd*me2*vp12*vp34 + 32
     &    *pm_qp*cfsrpd*me2*mpi2*vp12*vp34 + 32*pm_qp2*pp_qm*cfsrpd*
     &    vp12*vp34 - 32*pm_qp2*pp_qp*cfsrpd*vp12*vp34 - 16*pm_qp2*
     &    pp_k*cfsrpd*vp12*vp34 + 16*pm_qp2*pp_k*cisrpd*vp342 - 32*
     &    pm_k*pp_qm*pp_qp*cisrpd*vp342 - 16*pm_k*pp_qm*qm_qp*cfsrpd*
     &    vp12*vp34 )
      interf = interf + cisrm * (  - 16*pm_k*pp_qm*qm_k*cfsrmd*vp12*
     &    vp34 + 32*pm_k*pp_qm*qm_k*cisrmd*vp342 + 16*pm_k*pp_qm*qp_k*
     &    cfsrmd*vp12*vp34 - 32*pm_k*pp_qm*qp_k*cisrmd*vp342 + 8*pm_k*
     &    pp_qm*cfsrd*vp12*vp34 + 16*pm_k*pp_qm*cfsrpd*mpi2*vp12*vp34
     &     + 16*pm_k*pp_qm2*cisrpd*vp342 + 16*pm_k*pp_qp*qm_qp*cfsrmd
     &    *vp12*vp34 - 16*pm_k*pp_qp*qm_k*cfsrpd*vp12*vp34 - 32*pm_k*
     &    pp_qp*qm_k*cisrmd*vp342 + 16*pm_k*pp_qp*qp_k*cfsrpd*vp12*vp34
     &     + 32*pm_k*pp_qp*qp_k*cisrmd*vp342 - 8*pm_k*pp_qp*cfsrd*vp12*
     &    vp34 - 16*pm_k*pp_qp*cfsrmd*mpi2*vp12*vp34 + 16*pm_k*pp_qp2
     &    *cisrpd*vp342 + 16*pm_k*pp_k*qm_qp*cfsrpd*vp12*vp34 - 16*pm_k
     &    *pp_k*qm_qp*cfsrmd*vp12*vp34 + 32*pm_k*pp_k*qm_qp*cisrmd*
     &    vp342 - 16*pm_k*pp_k*cfsrpd*mpi2*vp12*vp34 + 16*pm_k*pp_k*
     &    cfsrmd*mpi2*vp12*vp34 - 32*pm_k*pp_k*cisrmd*mpi2*vp342 - 16*
     &    pm_k*qm_qp*cfsrpd*me2*vp12*vp34 + 16*pm_k*qm_qp*cfsrmd*me2*
     &    vp12*vp34 + 32*pm_k*qm_qp*cisrmd*me2*vp342 - 8*pm_k*qm_k*
     &    cfsrpd*me2*vp12*vp34 )
      interf = interf + cisrm * (  - 8*pm_k*qm_k*cfsrmd*me2*vp12*vp34
     &     + 8*pm_k*qp_k*cfsrpd*me2*vp12*vp34 + 8*pm_k*qp_k*cfsrmd*me2*
     &    vp12*vp34 + 16*pm_k*cfsrpd*me2*mpi2*vp12*vp34 - 16*pm_k*
     &    cfsrmd*me2*mpi2*vp12*vp34 - 32*pm_k*cisrmd*me2*mpi2*vp342 + 8
     &    *pm_k2*pp_qm*cfsrpd*vp12*vp34 + 8*pm_k2*pp_qm*cfsrmd*vp12
     &    *vp34 - 8*pm_k2*pp_qp*cfsrpd*vp12*vp34 - 8*pm_k2*pp_qp*
     &    cfsrmd*vp12*vp34 + 32*pp_qm*qm_k*cisrmd*me2*vp342 - 32*pp_qm*
     &    qp_k*cisrmd*me2*vp342 + 8*pp_qm*cfsrd*me2*vp12*vp34 - 32*
     &    pp_qp*qm_k*cisrmd*me2*vp342 + 32*pp_qp*qp_k*cisrmd*me2*vp342
     &     - 8*pp_qp*cfsrd*me2*vp12*vp34 + 32*pp_k*qm_qp*cisrmd*me2*
     &    vp342 - 32*pp_k*cisrmd*me2*mpi2*vp342 - 16*qm_qp*qm_k*cfsrmd*
     &    me2*vp12*vp34 + 16*qm_qp*qp_k*cfsrpd*me2*vp12*vp34 - 32*qm_qp
     &    *cisrmd*me4*vp342 + 16*qm_k*qp_k*cfsrpd*me2*vp12*vp34 - 16*
     &    qm_k*qp_k*cfsrmd*me2*vp12*vp34 + 32*qm_k*qp_k*cisrpd*me2*
     &    vp342 + 16*qm_k*cfsrd*me2*vp12*vp34 + 16*qm_k*cfsrmd*me2*mpi2
     &    *vp12*vp34 )
      interf = interf + cisrm * ( 16*qm_k2*cfsrmd*me2*vp12*vp34 - 16*
     &    qm_k2*cisrpd*me2*vp342 - 16*qp_k*cfsrd*me2*vp12*vp34 - 16*
     &    qp_k*cfsrpd*me2*mpi2*vp12*vp34 - 16*qp_k2*cfsrpd*me2*vp12*
     &    vp34 - 16*qp_k2*cisrpd*me2*vp342 )
