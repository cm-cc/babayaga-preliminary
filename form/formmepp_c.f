*      
*      cfsr = (-im)*(2.d0*1.d0/(qp+qm+k)^2)            !c5 phase: -i
*      cfsrp= (+im)*1.d0/((qp+k)^2-mpi2)/((qp+qm+k)^2) !c3 phase: +i
*      cfsrm= (+im)*1.d0/((qm+k)^2-mpi2)/((qp+qm+k)^2) !c4 phase: +i
*      cisrp= (+im)*1.d0/((pp-k)^2-me2)/(qp+qm)^2      !c1 phase: +i
*      cisrm= (+im)*1.d0/((pm-k)^2-me2)/(qp+qm)^2      !c2 phase: +i
*
c      im= (0.d0,1.d0)
      cfsr = (-im)*2.d0*1.d0/qppqmpk2                 !c5 phase: -i
      cfsrp= (+im)*1.d0/qppk2mmpi2/qppqmpk2           !c3 phase: +i
      cfsrm= (+im)*1.d0/qmpk2mmpi2/qppqmpk2           !c4 phase: +i
      cisrp= (-im)*1.d0/ppmk2mme2/qppqm2              !c1 phase: -i
      cisrm= (-im)*1.d0/pmmk2mme2/qppqm2              !c2 phase: -i
*
* inclusion of the form factor (called with q2 of the virtual photon)
*
      cfsr = cfsr * ffpi12
      cfsrp= cfsrp * ffpi12
      cfsrm= cfsrm * ffpi12
      cisrp= cisrp * ffpi34
      cisrm= cisrm * ffpi34
*
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

c$$$      interf =
c$$$     &  + cfsr * ( 16*cfsrd*me2*vp122 + 8*cfsrpd*me2*mpi2*vp122 + 8*
c$$$     &    cfsrmd*me2*mpi2*vp122 - 8*pm_pp*qm_qp*cfsrpd*vp122 - 8*pm_pp*
c$$$     &    qm_qp*cfsrmd*vp122 - 4*pm_pp*qm_k*cfsrpd*vp122 + 12*pm_pp*
c$$$     &    qm_k*cfsrmd*vp122 - 8*pm_pp*qm_k*cisrpd*vp12*vp34 + 8*pm_pp*
c$$$     &    qm_k*cisrmd*vp12*vp34 + 12*pm_pp*qp_k*cfsrpd*vp122 - 4*pm_pp*
c$$$     &    qp_k*cfsrmd*vp122 + 8*pm_pp*qp_k*cisrpd*vp12*vp34 - 8*pm_pp*
c$$$     &    qp_k*cisrmd*vp12*vp34 + 8*pm_pp*cfsrd*vp122 + 8*pm_pp*cfsrpd*
c$$$     &    mpi2*vp122 + 8*pm_pp*cfsrmd*mpi2*vp122 - 16*pm_qm*pp_qm*
c$$$     &    cfsrmd*vp122 + 8*pm_qm*pp_qp*cfsrpd*vp122 + 8*pm_qm*pp_qp*
c$$$     &    cfsrmd*vp122 + 4*pm_qm*pp_k*cfsrpd*vp122 - 12*pm_qm*pp_k*
c$$$     &    cfsrmd*vp122 - 8*pm_qm*pp_k*cisrpd*vp12*vp34 - 8*pm_qm*pp_k*
c$$$     &    cisrmd*vp12*vp34 - 8*pm_qm*cisrpd*me2*vp12*vp34 - 8*pm_qm*
c$$$     &    cisrmd*me2*vp12*vp34 + 8*pm_qp*pp_qm*cfsrpd*vp122 + 8*pm_qp*
c$$$     &    pp_qm*cfsrmd*vp122 - 16*pm_qp*pp_qp*cfsrpd*vp122 - 12*pm_qp*
c$$$     &    pp_k*cfsrpd*vp122 + 4*pm_qp*pp_k*cfsrmd*vp122 + 8*pm_qp*pp_k*
c$$$     &    cisrpd*vp12*vp34 )
c$$$      interf = interf + cfsr * ( 8*pm_qp*pp_k*cisrmd*vp12*vp34 + 8*
c$$$     &    pm_qp*cisrpd*me2*vp12*vp34 + 8*pm_qp*cisrmd*me2*vp12*vp34 + 4
c$$$     &    *pm_k*pp_qm*cfsrpd*vp122 - 12*pm_k*pp_qm*cfsrmd*vp122 + 8*
c$$$     &    pm_k*pp_qm*cisrpd*vp12*vp34 + 8*pm_k*pp_qm*cisrmd*vp12*vp34
c$$$     &     - 12*pm_k*pp_qp*cfsrpd*vp122 + 4*pm_k*pp_qp*cfsrmd*vp122 - 8
c$$$     &    *pm_k*pp_qp*cisrpd*vp12*vp34 - 8*pm_k*pp_qp*cisrmd*vp12*vp34
c$$$     &     - 8*pm_k*pp_k*cfsrpd*vp122 - 8*pm_k*pp_k*cfsrmd*vp122 + 8*
c$$$     &    pp_qm*cisrpd*me2*vp12*vp34 + 8*pp_qm*cisrmd*me2*vp12*vp34 - 8
c$$$     &    *pp_qp*cisrpd*me2*vp12*vp34 - 8*pp_qp*cisrmd*me2*vp12*vp34 - 
c$$$     &    8*qm_qp*cfsrpd*me2*vp122 - 8*qm_qp*cfsrmd*me2*vp122 - 4*qm_k*
c$$$     &    cfsrpd*me2*vp122 + 12*qm_k*cfsrmd*me2*vp122 - 16*qm_k*cisrpd*
c$$$     &    me2*vp12*vp34 + 16*qm_k*cisrmd*me2*vp12*vp34 + 12*qp_k*cfsrpd
c$$$     &    *me2*vp122 - 4*qp_k*cfsrmd*me2*vp122 + 16*qp_k*cisrpd*me2*
c$$$     &    vp12*vp34 - 16*qp_k*cisrmd*me2*vp12*vp34 )
c$$$      interf = interf + cfsrp * ( 8*cfsrd*me2*mpi2*vp122 + 32*cfsrpd*
c$$$     &    me2*mpi4*vp122 - 32*pm_pp*pm_qp*qm_qp*cisrmd*vp12*vp34 - 16*
c$$$     &    pm_pp*pm_qp*qm_k*cisrmd*vp12*vp34 + 16*pm_pp*pm_qp*qp_k*
c$$$     &    cisrmd*vp12*vp34 + 32*pm_pp*pm_qp*cisrmd*mpi2*vp12*vp34 - 16*
c$$$     &    pm_pp*pm_k*qm_qp*cisrmd*vp12*vp34 - 8*pm_pp*pm_k*qm_k*cisrmd*
c$$$     &    vp12*vp34 + 8*pm_pp*pm_k*qp_k*cisrmd*vp12*vp34 + 16*pm_pp*
c$$$     &    pm_k*cisrmd*mpi2*vp12*vp34 + 32*pm_pp*pp_qp*qm_qp*cisrpd*vp12
c$$$     &    *vp34 + 16*pm_pp*pp_qp*qm_k*cisrpd*vp12*vp34 - 16*pm_pp*pp_qp
c$$$     &    *qp_k*cisrpd*vp12*vp34 - 32*pm_pp*pp_qp*cisrpd*mpi2*vp12*vp34
c$$$     &     + 16*pm_pp*pp_k*qm_qp*cisrpd*vp12*vp34 + 8*pm_pp*pp_k*qm_k*
c$$$     &    cisrpd*vp12*vp34 - 8*pm_pp*pp_k*qp_k*cisrpd*vp12*vp34 - 16*
c$$$     &    pm_pp*pp_k*cisrpd*mpi2*vp12*vp34 + 16*pm_pp*qm_qp*qm_k*cfsrmd
c$$$     &    *vp122 - 32*pm_pp*qm_qp*qp_k*cfsrpd*vp122 + 16*pm_pp*qm_qp*
c$$$     &    qp_k*cfsrmd*vp122 - 16*pm_pp*qm_qp*qp_k*cisrpd*vp12*vp34 + 16
c$$$     &    *pm_pp*qm_qp*qp_k*cisrmd*vp12*vp34 - 8*pm_pp*qm_qp*cfsrd*
c$$$     &    vp122 )
c$$$      interf = interf + cfsrp * (  - 32*pm_pp*qm_qp*cfsrpd*mpi2*vp122
c$$$     &     - 32*pm_pp*qm_qp*cfsrmd*mpi2*vp122 + 32*pm_pp*qm_qp**2*
c$$$     &    cfsrmd*vp122 - 32*pm_pp*qm_k*qp_k*cfsrpd*vp122 - 16*pm_pp*
c$$$     &    qm_k*qp_k*cisrpd*vp12*vp34 + 16*pm_pp*qm_k*qp_k*cisrmd*vp12*
c$$$     &    vp34 - 4*pm_pp*qm_k*cfsrd*vp122 - 32*pm_pp*qm_k*cfsrpd*mpi2*
c$$$     &    vp122 - 16*pm_pp*qm_k*cfsrmd*mpi2*vp122 + 12*pm_pp*qp_k*cfsrd
c$$$     &    *vp122 + 64*pm_pp*qp_k*cfsrpd*mpi2*vp122 - 16*pm_pp*qp_k*
c$$$     &    cfsrmd*mpi2*vp122 + 16*pm_pp*qp_k*cisrpd*mpi2*vp12*vp34 - 16*
c$$$     &    pm_pp*qp_k*cisrmd*mpi2*vp12*vp34 + 32*pm_pp*qp_k**2*cfsrpd*
c$$$     &    vp122 + 16*pm_pp*qp_k**2*cisrpd*vp12*vp34 - 16*pm_pp*qp_k**2*
c$$$     &    cisrmd*vp12*vp34 + 8*pm_pp*cfsrd*mpi2*vp122 + 32*pm_pp*cfsrpd
c$$$     &    *mpi4*vp122 - 32*pm_qm*pm_qp*pp_qm*cisrmd*vp12*vp34 + 32*
c$$$     &    pm_qm*pm_qp*pp_qp*cisrmd*vp12*vp34 + 16*pm_qm*pm_qp*pp_k*
c$$$     &    cisrmd*vp12*vp34 - 16*pm_qm*pm_k*pp_qm*cisrmd*vp12*vp34 + 16*
c$$$     &    pm_qm*pm_k*pp_qp*cisrmd*vp12*vp34 + 8*pm_qm*pm_k*pp_k*cisrmd*
c$$$     &    vp12*vp34 )
c$$$      interf = interf + cfsrp * ( 32*pm_qm*pp_qm*pp_qp*cisrpd*vp12*vp34
c$$$     &     + 16*pm_qm*pp_qm*pp_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qm*
c$$$     &    qm_qp*cfsrmd*vp122 + 16*pm_qm*pp_qm*qm_k*cfsrmd*vp122 - 32*
c$$$     &    pm_qm*pp_qm*qp_k*cfsrpd*vp122 + 16*pm_qm*pp_qm*qp_k*cfsrmd*
c$$$     &    vp122 - 16*pm_qm*pp_qm*qp_k*cisrpd*vp12*vp34 + 16*pm_qm*pp_qm
c$$$     &    *qp_k*cisrmd*vp12*vp34 - 32*pm_qm*pp_qm*cfsrpd*mpi2*vp122 - 
c$$$     &    32*pm_qm*pp_qp*pp_k*cisrpd*vp12*vp34 - 32*pm_qm*pp_qp*qm_qp*
c$$$     &    cfsrmd*vp122 - 16*pm_qm*pp_qp*qm_k*cfsrmd*vp122 - 16*pm_qm*
c$$$     &    pp_qp*qm_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qp*qp_k*cfsrpd*
c$$$     &    vp122 - 16*pm_qm*pp_qp*qp_k*cfsrmd*vp122 + 32*pm_qm*pp_qp*
c$$$     &    qp_k*cisrpd*vp12*vp34 - 16*pm_qm*pp_qp*qp_k*cisrmd*vp12*vp34
c$$$     &     + 8*pm_qm*pp_qp*cfsrd*vp122 + 32*pm_qm*pp_qp*cfsrpd*mpi2*
c$$$     &    vp122 - 32*pm_qm*pp_qp**2*cisrpd*vp12*vp34 + 16*pm_qm*pp_k*
c$$$     &    qm_qp*cisrpd*vp12*vp34 + 32*pm_qm*pp_k*qp_k*cfsrpd*vp122 - 16
c$$$     &    *pm_qm*pp_k*qp_k*cisrmd*vp12*vp34 + 4*pm_qm*pp_k*cfsrd*vp122
c$$$     &     + 32*pm_qm*pp_k*cfsrpd*mpi2*vp122 )
c$$$      interf = interf + cfsrp * (  - 16*pm_qm*pp_k*cisrpd*mpi2*vp12*
c$$$     &    vp34 - 8*pm_qm*pp_k**2*cisrpd*vp12*vp34 + 32*pm_qp*pm_k*pp_qm
c$$$     &    *cisrmd*vp12*vp34 - 32*pm_qp*pm_k*pp_qp*cisrmd*vp12*vp34 - 8*
c$$$     &    pm_qp*pm_k*pp_k*cisrmd*vp12*vp34 - 32*pm_qp*pp_qm*pp_qp*
c$$$     &    cisrpd*vp12*vp34 - 16*pm_qp*pp_qm*pp_k*cisrpd*vp12*vp34 - 32*
c$$$     &    pm_qp*pp_qm*qm_qp*cfsrmd*vp122 - 16*pm_qp*pp_qm*qm_k*cfsrmd*
c$$$     &    vp122 + 16*pm_qp*pp_qm*qm_k*cisrmd*vp12*vp34 + 32*pm_qp*pp_qm
c$$$     &    *qp_k*cfsrpd*vp122 - 16*pm_qp*pp_qm*qp_k*cfsrmd*vp122 + 16*
c$$$     &    pm_qp*pp_qm*qp_k*cisrpd*vp12*vp34 - 32*pm_qp*pp_qm*qp_k*
c$$$     &    cisrmd*vp12*vp34 + 8*pm_qp*pp_qm*cfsrd*vp122 + 32*pm_qp*pp_qm
c$$$     &    *cfsrpd*mpi2*vp122 + 32*pm_qp*pp_qp*pp_k*cisrpd*vp12*vp34 + 
c$$$     &    32*pm_qp*pp_qp*qm_qp*cfsrmd*vp122 + 16*pm_qp*pp_qp*qm_k*
c$$$     &    cfsrmd*vp122 + 16*pm_qp*pp_qp*qm_k*cisrpd*vp12*vp34 - 16*
c$$$     &    pm_qp*pp_qp*qm_k*cisrmd*vp12*vp34 - 32*pm_qp*pp_qp*qp_k*
c$$$     &    cfsrpd*vp122 + 16*pm_qp*pp_qp*qp_k*cfsrmd*vp122 - 32*pm_qp*
c$$$     &    pp_qp*qp_k*cisrpd*vp12*vp34 )
c$$$      interf = interf + cfsrp * ( 32*pm_qp*pp_qp*qp_k*cisrmd*vp12*vp34
c$$$     &     - 16*pm_qp*pp_qp*cfsrd*vp122 - 32*pm_qp*pp_qp*cfsrpd*mpi2*
c$$$     &    vp122 + 32*pm_qp*pp_qp**2*cisrpd*vp12*vp34 + 16*pm_qp*pp_k*
c$$$     &    qm_qp*cisrmd*vp12*vp34 + 16*pm_qp*pp_k*qm_k*cisrpd*vp12*vp34
c$$$     &     - 32*pm_qp*pp_k*qp_k*cfsrpd*vp122 - 16*pm_qp*pp_k*qp_k*
c$$$     &    cisrpd*vp12*vp34 + 16*pm_qp*pp_k*qp_k*cisrmd*vp12*vp34 - 12*
c$$$     &    pm_qp*pp_k*cfsrd*vp122 - 32*pm_qp*pp_k*cfsrpd*mpi2*vp122 - 16
c$$$     &    *pm_qp*pp_k*cisrmd*mpi2*vp12*vp34 + 8*pm_qp*pp_k**2*cisrpd*
c$$$     &    vp12*vp34 - 32*pm_qp*qm_qp*cisrmd*me2*vp12*vp34 - 16*pm_qp*
c$$$     &    qm_k*cisrmd*me2*vp12*vp34 + 16*pm_qp*qp_k*cisrmd*me2*vp12*
c$$$     &    vp34 + 32*pm_qp*cisrmd*me2*mpi2*vp12*vp34 + 32*pm_qp**2*pp_qm
c$$$     &    *cisrmd*vp12*vp34 - 32*pm_qp**2*pp_qp*cisrmd*vp12*vp34 - 16*
c$$$     &    pm_qp**2*pp_k*cisrmd*vp12*vp34 - 16*pm_k*pp_qm*pp_qp*cisrpd*
c$$$     &    vp12*vp34 - 8*pm_k*pp_qm*pp_k*cisrpd*vp12*vp34 - 16*pm_k*
c$$$     &    pp_qm*qm_qp*cisrmd*vp12*vp34 + 32*pm_k*pp_qm*qp_k*cfsrpd*
c$$$     &    vp122 )
c$$$      interf = interf + cfsrp * ( 16*pm_k*pp_qm*qp_k*cisrpd*vp12*vp34
c$$$     &     + 4*pm_k*pp_qm*cfsrd*vp122 + 32*pm_k*pp_qm*cfsrpd*mpi2*vp122
c$$$     &     + 16*pm_k*pp_qm*cisrmd*mpi2*vp12*vp34 + 8*pm_k*pp_qp*pp_k*
c$$$     &    cisrpd*vp12*vp34 - 16*pm_k*pp_qp*qm_qp*cisrpd*vp12*vp34 - 16*
c$$$     &    pm_k*pp_qp*qm_k*cisrmd*vp12*vp34 - 32*pm_k*pp_qp*qp_k*cfsrpd*
c$$$     &    vp122 - 16*pm_k*pp_qp*qp_k*cisrpd*vp12*vp34 + 16*pm_k*pp_qp*
c$$$     &    qp_k*cisrmd*vp12*vp34 - 12*pm_k*pp_qp*cfsrd*vp122 - 32*pm_k*
c$$$     &    pp_qp*cfsrpd*mpi2*vp122 + 16*pm_k*pp_qp*cisrpd*mpi2*vp12*vp34
c$$$     &     + 16*pm_k*pp_qp**2*cisrpd*vp12*vp34 - 32*pm_k*pp_k*qm_qp*
c$$$     &    cfsrmd*vp122 - 16*pm_k*pp_k*qm_qp*cisrpd*vp12*vp34 + 16*pm_k*
c$$$     &    pp_k*qm_qp*cisrmd*vp12*vp34 - 16*pm_k*pp_k*qm_k*cfsrmd*vp122
c$$$     &     - 32*pm_k*pp_k*qp_k*cfsrpd*vp122 - 16*pm_k*pp_k*qp_k*cfsrmd*
c$$$     &    vp122 - 8*pm_k*pp_k*cfsrd*vp122 - 32*pm_k*pp_k*cfsrpd*mpi2*
c$$$     &    vp122 + 16*pm_k*pp_k*cisrpd*mpi2*vp12*vp34 - 16*pm_k*pp_k*
c$$$     &    cisrmd*mpi2*vp12*vp34 - 16*pm_k*qm_qp*cisrmd*me2*vp12*vp34 - 
c$$$     &    8*pm_k*qm_k*cisrmd*me2*vp12*vp34 )
c$$$      interf = interf + cfsrp * ( 8*pm_k*qp_k*cisrmd*me2*vp12*vp34 + 16
c$$$     &    *pm_k*cisrmd*me2*mpi2*vp12*vp34 + 8*pm_k**2*pp_qm*cisrmd*vp12
c$$$     &    *vp34 - 8*pm_k**2*pp_qp*cisrmd*vp12*vp34 + 32*pp_qp*qm_qp*
c$$$     &    cisrpd*me2*vp12*vp34 + 16*pp_qp*qm_k*cisrpd*me2*vp12*vp34 - 
c$$$     &    16*pp_qp*qp_k*cisrpd*me2*vp12*vp34 - 32*pp_qp*cisrpd*me2*mpi2
c$$$     &    *vp12*vp34 + 16*pp_k*qm_qp*cisrpd*me2*vp12*vp34 + 8*pp_k*qm_k
c$$$     &    *cisrpd*me2*vp12*vp34 - 8*pp_k*qp_k*cisrpd*me2*vp12*vp34 - 16
c$$$     &    *pp_k*cisrpd*me2*mpi2*vp12*vp34 + 16*qm_qp*qm_k*cfsrmd*me2*
c$$$     &    vp122 - 32*qm_qp*qp_k*cfsrpd*me2*vp122 + 16*qm_qp*qp_k*cfsrmd
c$$$     &    *me2*vp122 - 16*qm_qp*qp_k*cisrpd*me2*vp12*vp34 + 16*qm_qp*
c$$$     &    qp_k*cisrmd*me2*vp12*vp34 - 8*qm_qp*cfsrd*me2*vp122 - 32*
c$$$     &    qm_qp*cfsrpd*me2*mpi2*vp122 - 32*qm_qp*cfsrmd*me2*mpi2*vp122
c$$$     &     + 32*qm_qp**2*cfsrmd*me2*vp122 - 32*qm_k*qp_k*cfsrpd*me2*
c$$$     &    vp122 - 16*qm_k*qp_k*cisrpd*me2*vp12*vp34 + 16*qm_k*qp_k*
c$$$     &    cisrmd*me2*vp12*vp34 - 4*qm_k*cfsrd*me2*vp122 - 32*qm_k*
c$$$     &    cfsrpd*me2*mpi2*vp122 )
c$$$      interf = interf + cfsrp * (  - 16*qm_k*cfsrmd*me2*mpi2*vp122 + 12
c$$$     &    *qp_k*cfsrd*me2*vp122 + 64*qp_k*cfsrpd*me2*mpi2*vp122 - 16*
c$$$     &    qp_k*cfsrmd*me2*mpi2*vp122 + 16*qp_k*cisrpd*me2*mpi2*vp12*
c$$$     &    vp34 - 16*qp_k*cisrmd*me2*mpi2*vp12*vp34 + 32*qp_k**2*cfsrpd*
c$$$     &    me2*vp122 + 16*qp_k**2*cisrpd*me2*vp12*vp34 - 16*qp_k**2*
c$$$     &    cisrmd*me2*vp12*vp34 )
c$$$      interf = interf + cfsrm * ( 8*cfsrd*me2*mpi2*vp122 + 32*cfsrmd*
c$$$     &    me2*mpi4*vp122 + 32*pm_pp*pm_qm*qm_qp*cisrmd*vp12*vp34 - 16*
c$$$     &    pm_pp*pm_qm*qm_k*cisrmd*vp12*vp34 + 16*pm_pp*pm_qm*qp_k*
c$$$     &    cisrmd*vp12*vp34 - 32*pm_pp*pm_qm*cisrmd*mpi2*vp12*vp34 + 16*
c$$$     &    pm_pp*pm_k*qm_qp*cisrmd*vp12*vp34 - 8*pm_pp*pm_k*qm_k*cisrmd*
c$$$     &    vp12*vp34 + 8*pm_pp*pm_k*qp_k*cisrmd*vp12*vp34 - 16*pm_pp*
c$$$     &    pm_k*cisrmd*mpi2*vp12*vp34 - 32*pm_pp*pp_qm*qm_qp*cisrpd*vp12
c$$$     &    *vp34 + 16*pm_pp*pp_qm*qm_k*cisrpd*vp12*vp34 - 16*pm_pp*pp_qm
c$$$     &    *qp_k*cisrpd*vp12*vp34 + 32*pm_pp*pp_qm*cisrpd*mpi2*vp12*vp34
c$$$     &     - 16*pm_pp*pp_k*qm_qp*cisrpd*vp12*vp34 + 8*pm_pp*pp_k*qm_k*
c$$$     &    cisrpd*vp12*vp34 - 8*pm_pp*pp_k*qp_k*cisrpd*vp12*vp34 + 16*
c$$$     &    pm_pp*pp_k*cisrpd*mpi2*vp12*vp34 + 16*pm_pp*qm_qp*qm_k*cfsrpd
c$$$     &    *vp122 - 32*pm_pp*qm_qp*qm_k*cfsrmd*vp122 + 16*pm_pp*qm_qp*
c$$$     &    qm_k*cisrpd*vp12*vp34 - 16*pm_pp*qm_qp*qm_k*cisrmd*vp12*vp34
c$$$     &     + 16*pm_pp*qm_qp*qp_k*cfsrpd*vp122 - 8*pm_pp*qm_qp*cfsrd*
c$$$     &    vp122 )
c$$$      interf = interf + cfsrm * (  - 32*pm_pp*qm_qp*cfsrpd*mpi2*vp122
c$$$     &     - 32*pm_pp*qm_qp*cfsrmd*mpi2*vp122 + 32*pm_pp*qm_qp**2*
c$$$     &    cfsrpd*vp122 - 32*pm_pp*qm_k*qp_k*cfsrmd*vp122 + 16*pm_pp*
c$$$     &    qm_k*qp_k*cisrpd*vp12*vp34 - 16*pm_pp*qm_k*qp_k*cisrmd*vp12*
c$$$     &    vp34 + 12*pm_pp*qm_k*cfsrd*vp122 - 16*pm_pp*qm_k*cfsrpd*mpi2*
c$$$     &    vp122 + 64*pm_pp*qm_k*cfsrmd*mpi2*vp122 - 16*pm_pp*qm_k*
c$$$     &    cisrpd*mpi2*vp12*vp34 + 16*pm_pp*qm_k*cisrmd*mpi2*vp12*vp34
c$$$     &     + 32*pm_pp*qm_k**2*cfsrmd*vp122 - 16*pm_pp*qm_k**2*cisrpd*
c$$$     &    vp12*vp34 + 16*pm_pp*qm_k**2*cisrmd*vp12*vp34 - 4*pm_pp*qp_k*
c$$$     &    cfsrd*vp122 - 16*pm_pp*qp_k*cfsrpd*mpi2*vp122 - 32*pm_pp*qp_k
c$$$     &    *cfsrmd*mpi2*vp122 + 8*pm_pp*cfsrd*mpi2*vp122 + 32*pm_pp*
c$$$     &    cfsrmd*mpi4*vp122 - 32*pm_qm*pm_qp*pp_qm*cisrmd*vp12*vp34 + 
c$$$     &    32*pm_qm*pm_qp*pp_qp*cisrmd*vp12*vp34 - 16*pm_qm*pm_qp*pp_k*
c$$$     &    cisrmd*vp12*vp34 + 32*pm_qm*pm_k*pp_qm*cisrmd*vp12*vp34 - 32*
c$$$     &    pm_qm*pm_k*pp_qp*cisrmd*vp12*vp34 + 8*pm_qm*pm_k*pp_k*cisrmd*
c$$$     &    vp12*vp34 )
c$$$      interf = interf + cfsrm * ( 32*pm_qm*pp_qm*pp_qp*cisrpd*vp12*vp34
c$$$     &     - 32*pm_qm*pp_qm*pp_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qm*
c$$$     &    qm_qp*cfsrpd*vp122 + 16*pm_qm*pp_qm*qm_k*cfsrpd*vp122 - 32*
c$$$     &    pm_qm*pp_qm*qm_k*cfsrmd*vp122 + 32*pm_qm*pp_qm*qm_k*cisrpd*
c$$$     &    vp12*vp34 - 32*pm_qm*pp_qm*qm_k*cisrmd*vp12*vp34 + 16*pm_qm*
c$$$     &    pp_qm*qp_k*cfsrpd*vp122 - 16*pm_qm*pp_qm*qp_k*cisrpd*vp12*
c$$$     &    vp34 + 16*pm_qm*pp_qm*qp_k*cisrmd*vp12*vp34 - 16*pm_qm*pp_qm*
c$$$     &    cfsrd*vp122 - 32*pm_qm*pp_qm*cfsrmd*mpi2*vp122 - 32*pm_qm*
c$$$     &    pp_qm**2*cisrpd*vp12*vp34 + 16*pm_qm*pp_qp*pp_k*cisrpd*vp12*
c$$$     &    vp34 - 32*pm_qm*pp_qp*qm_qp*cfsrpd*vp122 - 16*pm_qm*pp_qp*
c$$$     &    qm_k*cfsrpd*vp122 + 32*pm_qm*pp_qp*qm_k*cfsrmd*vp122 - 16*
c$$$     &    pm_qm*pp_qp*qm_k*cisrpd*vp12*vp34 + 32*pm_qm*pp_qp*qm_k*
c$$$     &    cisrmd*vp12*vp34 - 16*pm_qm*pp_qp*qp_k*cfsrpd*vp122 - 16*
c$$$     &    pm_qm*pp_qp*qp_k*cisrmd*vp12*vp34 + 8*pm_qm*pp_qp*cfsrd*vp122
c$$$     &     + 32*pm_qm*pp_qp*cfsrmd*mpi2*vp122 - 16*pm_qm*pp_k*qm_qp*
c$$$     &    cisrmd*vp12*vp34 )
c$$$      interf = interf + cfsrm * (  - 32*pm_qm*pp_k*qm_k*cfsrmd*vp122 + 
c$$$     &    16*pm_qm*pp_k*qm_k*cisrpd*vp12*vp34 - 16*pm_qm*pp_k*qm_k*
c$$$     &    cisrmd*vp12*vp34 - 16*pm_qm*pp_k*qp_k*cisrpd*vp12*vp34 - 12*
c$$$     &    pm_qm*pp_k*cfsrd*vp122 - 32*pm_qm*pp_k*cfsrmd*mpi2*vp122 + 16
c$$$     &    *pm_qm*pp_k*cisrmd*mpi2*vp12*vp34 - 8*pm_qm*pp_k**2*cisrpd*
c$$$     &    vp12*vp34 + 32*pm_qm*qm_qp*cisrmd*me2*vp12*vp34 - 16*pm_qm*
c$$$     &    qm_k*cisrmd*me2*vp12*vp34 + 16*pm_qm*qp_k*cisrmd*me2*vp12*
c$$$     &    vp34 - 32*pm_qm*cisrmd*me2*mpi2*vp12*vp34 + 32*pm_qm**2*pp_qm
c$$$     &    *cisrmd*vp12*vp34 - 32*pm_qm**2*pp_qp*cisrmd*vp12*vp34 + 16*
c$$$     &    pm_qm**2*pp_k*cisrmd*vp12*vp34 - 16*pm_qp*pm_k*pp_qm*cisrmd*
c$$$     &    vp12*vp34 + 16*pm_qp*pm_k*pp_qp*cisrmd*vp12*vp34 - 8*pm_qp*
c$$$     &    pm_k*pp_k*cisrmd*vp12*vp34 - 32*pm_qp*pp_qm*pp_qp*cisrpd*vp12
c$$$     &    *vp34 + 32*pm_qp*pp_qm*pp_k*cisrpd*vp12*vp34 - 32*pm_qp*pp_qm
c$$$     &    *qm_qp*cfsrpd*vp122 - 16*pm_qp*pp_qm*qm_k*cfsrpd*vp122 + 32*
c$$$     &    pm_qp*pp_qm*qm_k*cfsrmd*vp122 - 32*pm_qp*pp_qm*qm_k*cisrpd*
c$$$     &    vp12*vp34 )
c$$$      interf = interf + cfsrm * ( 16*pm_qp*pp_qm*qm_k*cisrmd*vp12*vp34
c$$$     &     - 16*pm_qp*pp_qm*qp_k*cfsrpd*vp122 + 16*pm_qp*pp_qm*qp_k*
c$$$     &    cisrpd*vp12*vp34 + 8*pm_qp*pp_qm*cfsrd*vp122 + 32*pm_qp*pp_qm
c$$$     &    *cfsrmd*mpi2*vp122 + 32*pm_qp*pp_qm**2*cisrpd*vp12*vp34 - 16*
c$$$     &    pm_qp*pp_qp*pp_k*cisrpd*vp12*vp34 + 32*pm_qp*pp_qp*qm_qp*
c$$$     &    cfsrpd*vp122 + 16*pm_qp*pp_qp*qm_k*cfsrpd*vp122 - 32*pm_qp*
c$$$     &    pp_qp*qm_k*cfsrmd*vp122 + 16*pm_qp*pp_qp*qm_k*cisrpd*vp12*
c$$$     &    vp34 - 16*pm_qp*pp_qp*qm_k*cisrmd*vp12*vp34 + 16*pm_qp*pp_qp*
c$$$     &    qp_k*cfsrpd*vp122 - 32*pm_qp*pp_qp*cfsrmd*mpi2*vp122 - 16*
c$$$     &    pm_qp*pp_k*qm_qp*cisrpd*vp12*vp34 + 32*pm_qp*pp_k*qm_k*cfsrmd
c$$$     &    *vp122 + 16*pm_qp*pp_k*qm_k*cisrmd*vp12*vp34 + 4*pm_qp*pp_k*
c$$$     &    cfsrd*vp122 + 32*pm_qp*pp_k*cfsrmd*mpi2*vp122 + 16*pm_qp*pp_k
c$$$     &    *cisrpd*mpi2*vp12*vp34 + 8*pm_qp*pp_k**2*cisrpd*vp12*vp34 + 
c$$$     &    16*pm_k*pp_qm*pp_qp*cisrpd*vp12*vp34 - 8*pm_k*pp_qm*pp_k*
c$$$     &    cisrpd*vp12*vp34 + 16*pm_k*pp_qm*qm_qp*cisrpd*vp12*vp34 - 32*
c$$$     &    pm_k*pp_qm*qm_k*cfsrmd*vp122 )
c$$$      interf = interf + cfsrm * ( 16*pm_k*pp_qm*qm_k*cisrpd*vp12*vp34
c$$$     &     - 16*pm_k*pp_qm*qm_k*cisrmd*vp12*vp34 + 16*pm_k*pp_qm*qp_k*
c$$$     &    cisrmd*vp12*vp34 - 12*pm_k*pp_qm*cfsrd*vp122 - 32*pm_k*pp_qm*
c$$$     &    cfsrmd*mpi2*vp122 - 16*pm_k*pp_qm*cisrpd*mpi2*vp12*vp34 - 16*
c$$$     &    pm_k*pp_qm**2*cisrpd*vp12*vp34 + 8*pm_k*pp_qp*pp_k*cisrpd*
c$$$     &    vp12*vp34 + 16*pm_k*pp_qp*qm_qp*cisrmd*vp12*vp34 + 32*pm_k*
c$$$     &    pp_qp*qm_k*cfsrmd*vp122 - 16*pm_k*pp_qp*qm_k*cisrpd*vp12*vp34
c$$$     &     + 4*pm_k*pp_qp*cfsrd*vp122 + 32*pm_k*pp_qp*cfsrmd*mpi2*vp122
c$$$     &     - 16*pm_k*pp_qp*cisrmd*mpi2*vp12*vp34 - 32*pm_k*pp_k*qm_qp*
c$$$     &    cfsrpd*vp122 + 16*pm_k*pp_k*qm_qp*cisrpd*vp12*vp34 - 16*pm_k*
c$$$     &    pp_k*qm_qp*cisrmd*vp12*vp34 - 16*pm_k*pp_k*qm_k*cfsrpd*vp122
c$$$     &     - 32*pm_k*pp_k*qm_k*cfsrmd*vp122 - 16*pm_k*pp_k*qp_k*cfsrpd*
c$$$     &    vp122 - 8*pm_k*pp_k*cfsrd*vp122 - 32*pm_k*pp_k*cfsrmd*mpi2*
c$$$     &    vp122 - 16*pm_k*pp_k*cisrpd*mpi2*vp12*vp34 + 16*pm_k*pp_k*
c$$$     &    cisrmd*mpi2*vp12*vp34 + 16*pm_k*qm_qp*cisrmd*me2*vp12*vp34 - 
c$$$     &    8*pm_k*qm_k*cisrmd*me2*vp12*vp34 )
c$$$      interf = interf + cfsrm * ( 8*pm_k*qp_k*cisrmd*me2*vp12*vp34 - 16
c$$$     &    *pm_k*cisrmd*me2*mpi2*vp12*vp34 + 8*pm_k**2*pp_qm*cisrmd*vp12
c$$$     &    *vp34 - 8*pm_k**2*pp_qp*cisrmd*vp12*vp34 - 32*pp_qm*qm_qp*
c$$$     &    cisrpd*me2*vp12*vp34 + 16*pp_qm*qm_k*cisrpd*me2*vp12*vp34 - 
c$$$     &    16*pp_qm*qp_k*cisrpd*me2*vp12*vp34 + 32*pp_qm*cisrpd*me2*mpi2
c$$$     &    *vp12*vp34 - 16*pp_k*qm_qp*cisrpd*me2*vp12*vp34 + 8*pp_k*qm_k
c$$$     &    *cisrpd*me2*vp12*vp34 - 8*pp_k*qp_k*cisrpd*me2*vp12*vp34 + 16
c$$$     &    *pp_k*cisrpd*me2*mpi2*vp12*vp34 + 16*qm_qp*qm_k*cfsrpd*me2*
c$$$     &    vp122 - 32*qm_qp*qm_k*cfsrmd*me2*vp122 + 16*qm_qp*qm_k*cisrpd
c$$$     &    *me2*vp12*vp34 - 16*qm_qp*qm_k*cisrmd*me2*vp12*vp34 + 16*
c$$$     &    qm_qp*qp_k*cfsrpd*me2*vp122 - 8*qm_qp*cfsrd*me2*vp122 - 32*
c$$$     &    qm_qp*cfsrpd*me2*mpi2*vp122 - 32*qm_qp*cfsrmd*me2*mpi2*vp122
c$$$     &     + 32*qm_qp**2*cfsrpd*me2*vp122 - 32*qm_k*qp_k*cfsrmd*me2*
c$$$     &    vp122 + 16*qm_k*qp_k*cisrpd*me2*vp12*vp34 - 16*qm_k*qp_k*
c$$$     &    cisrmd*me2*vp12*vp34 + 12*qm_k*cfsrd*me2*vp122 - 16*qm_k*
c$$$     &    cfsrpd*me2*mpi2*vp122 )
c$$$      interf = interf + cfsrm * ( 64*qm_k*cfsrmd*me2*mpi2*vp122 - 16*
c$$$     &    qm_k*cisrpd*me2*mpi2*vp12*vp34 + 16*qm_k*cisrmd*me2*mpi2*vp12
c$$$     &    *vp34 + 32*qm_k**2*cfsrmd*me2*vp122 - 16*qm_k**2*cisrpd*me2*
c$$$     &    vp12*vp34 + 16*qm_k**2*cisrmd*me2*vp12*vp34 - 4*qp_k*cfsrd*
c$$$     &    me2*vp122 - 16*qp_k*cfsrpd*me2*mpi2*vp122 - 32*qp_k*cfsrmd*
c$$$     &    me2*mpi2*vp122 )
c$$$      interf = interf + cisrp * ( 32*cisrpd*me4*mpi2*vp342 + 32*pm_pp*
c$$$     &    pm_qm*pp_qm*cisrmd*vp342 - 32*pm_pp*pm_qm*pp_qp*cisrmd*vp342
c$$$     &     - 16*pm_pp*pm_qm*qm_k*cisrmd*vp342 + 16*pm_pp*pm_qm*qp_k*
c$$$     &    cisrmd*vp342 - 32*pm_pp*pm_qp*pp_qm*cisrmd*vp342 + 32*pm_pp*
c$$$     &    pm_qp*pp_qp*cisrmd*vp342 + 16*pm_pp*pm_qp*qm_k*cisrmd*vp342
c$$$     &     - 16*pm_pp*pm_qp*qp_k*cisrmd*vp342 - 32*pm_pp*pm_k*qm_qp*
c$$$     &    cisrmd*vp342 + 32*pm_pp*pm_k*cisrmd*mpi2*vp342 - 32*pm_pp*
c$$$     &    pp_qm*qm_qp*cfsrmd*vp12*vp34 + 16*pm_pp*pp_qm*qm_k*cfsrmd*
c$$$     &    vp12*vp34 - 16*pm_pp*pp_qm*qm_k*cisrmd*vp342 - 16*pm_pp*pp_qm
c$$$     &    *qp_k*cfsrmd*vp12*vp34 + 16*pm_pp*pp_qm*qp_k*cisrmd*vp342 + 
c$$$     &    32*pm_pp*pp_qm*cfsrmd*mpi2*vp12*vp34 + 32*pm_pp*pp_qp*qm_qp*
c$$$     &    cfsrpd*vp12*vp34 + 16*pm_pp*pp_qp*qm_k*cfsrpd*vp12*vp34 + 16*
c$$$     &    pm_pp*pp_qp*qm_k*cisrmd*vp342 - 16*pm_pp*pp_qp*qp_k*cfsrpd*
c$$$     &    vp12*vp34 - 16*pm_pp*pp_qp*qp_k*cisrmd*vp342 - 32*pm_pp*pp_qp
c$$$     &    *cfsrpd*mpi2*vp12*vp34 + 16*pm_pp*pp_k*qm_qp*cfsrpd*vp12*vp34
c$$$     &     - 16*pm_pp*pp_k*qm_qp*cfsrmd*vp12*vp34 )
c$$$      interf = interf + cisrp * (  - 32*pm_pp*pp_k*qm_qp*cisrmd*vp342
c$$$     &     + 8*pm_pp*pp_k*qm_k*cfsrpd*vp12*vp34 + 8*pm_pp*pp_k*qm_k*
c$$$     &    cfsrmd*vp12*vp34 - 8*pm_pp*pp_k*qp_k*cfsrpd*vp12*vp34 - 8*
c$$$     &    pm_pp*pp_k*qp_k*cfsrmd*vp12*vp34 - 16*pm_pp*pp_k*cfsrpd*mpi2*
c$$$     &    vp12*vp34 + 16*pm_pp*pp_k*cfsrmd*mpi2*vp12*vp34 + 32*pm_pp*
c$$$     &    pp_k*cisrmd*mpi2*vp342 + 16*pm_pp*qm_qp*qm_k*cfsrmd*vp12*vp34
c$$$     &     - 16*pm_pp*qm_qp*qp_k*cfsrpd*vp12*vp34 - 32*pm_pp*qm_qp*
c$$$     &    cisrpd*me2*vp342 + 32*pm_pp*qm_qp*cisrmd*me2*vp342 - 16*pm_pp
c$$$     &    *qm_k*qp_k*cfsrpd*vp12*vp34 + 16*pm_pp*qm_k*qp_k*cfsrmd*vp12*
c$$$     &    vp34 - 8*pm_pp*qm_k*cfsrd*vp12*vp34 - 16*pm_pp*qm_k*cfsrmd*
c$$$     &    mpi2*vp12*vp34 - 16*pm_pp*qm_k**2*cfsrmd*vp12*vp34 + 8*pm_pp*
c$$$     &    qp_k*cfsrd*vp12*vp34 + 16*pm_pp*qp_k*cfsrpd*mpi2*vp12*vp34 + 
c$$$     &    16*pm_pp*qp_k**2*cfsrpd*vp12*vp34 + 32*pm_pp*cisrpd*me2*mpi2*
c$$$     &    vp342 - 32*pm_pp*cisrmd*me2*mpi2*vp342 + 32*pm_pp**2*qm_qp*
c$$$     &    cisrmd*vp342 - 32*pm_pp**2*cisrmd*mpi2*vp342 - 32*pm_qm*pm_qp
c$$$     &    *pp_k*cisrmd*vp342 )
c$$$      interf = interf + cisrp * (  - 16*pm_qm*pm_k*pp_qm*cisrmd*vp342
c$$$     &     + 16*pm_qm*pm_k*pp_qp*cisrmd*vp342 + 32*pm_qm*pp_qm*pp_qp*
c$$$     &    cfsrpd*vp12*vp34 + 32*pm_qm*pp_qm*pp_qp*cfsrmd*vp12*vp34 + 16
c$$$     &    *pm_qm*pp_qm*pp_k*cfsrpd*vp12*vp34 - 32*pm_qm*pp_qm*pp_k*
c$$$     &    cfsrmd*vp12*vp34 - 16*pm_qm*pp_qm*pp_k*cisrmd*vp342 + 32*
c$$$     &    pm_qm*pp_qm*qm_k*cfsrmd*vp12*vp34 - 16*pm_qm*pp_qm*qp_k*
c$$$     &    cfsrpd*vp12*vp34 - 16*pm_qm*pp_qm*qp_k*cfsrmd*vp12*vp34 - 32*
c$$$     &    pm_qm*pp_qm*cisrpd*me2*vp342 - 32*pm_qm*pp_qm**2*cfsrmd*vp12*
c$$$     &    vp34 - 32*pm_qm*pp_qp*pp_k*cfsrpd*vp12*vp34 + 16*pm_qm*pp_qp*
c$$$     &    pp_k*cfsrmd*vp12*vp34 + 16*pm_qm*pp_qp*pp_k*cisrmd*vp342 - 16
c$$$     &    *pm_qm*pp_qp*qm_k*cfsrpd*vp12*vp34 - 16*pm_qm*pp_qp*qm_k*
c$$$     &    cfsrmd*vp12*vp34 + 32*pm_qm*pp_qp*qp_k*cfsrpd*vp12*vp34 + 32*
c$$$     &    pm_qm*pp_qp*cisrpd*me2*vp342 - 32*pm_qm*pp_qp**2*cfsrpd*vp12*
c$$$     &    vp34 + 16*pm_qm*pp_k*qm_qp*cfsrpd*vp12*vp34 + 16*pm_qm*pp_k*
c$$$     &    qm_k*cfsrmd*vp12*vp34 + 32*pm_qm*pp_k*qm_k*cisrpd*vp342 - 16*
c$$$     &    pm_qm*pp_k*qp_k*cfsrmd*vp12*vp34 )
c$$$      interf = interf + cisrp * (  - 32*pm_qm*pp_k*qp_k*cisrpd*vp342 - 
c$$$     &    8*pm_qm*pp_k*cfsrd*vp12*vp34 - 16*pm_qm*pp_k*cfsrpd*mpi2*vp12
c$$$     &    *vp34 - 8*pm_qm*pp_k**2*cfsrpd*vp12*vp34 - 8*pm_qm*pp_k**2*
c$$$     &    cfsrmd*vp12*vp34 + 32*pm_qm*qm_k*cisrpd*me2*vp342 - 32*pm_qm*
c$$$     &    qp_k*cisrpd*me2*vp342 - 8*pm_qm*cfsrd*me2*vp12*vp34 + 16*
c$$$     &    pm_qm**2*pp_k*cisrmd*vp342 + 16*pm_qp*pm_k*pp_qm*cisrmd*vp342
c$$$     &     - 16*pm_qp*pm_k*pp_qp*cisrmd*vp342 - 32*pm_qp*pp_qm*pp_qp*
c$$$     &    cfsrpd*vp12*vp34 - 32*pm_qp*pp_qm*pp_qp*cfsrmd*vp12*vp34 - 16
c$$$     &    *pm_qp*pp_qm*pp_k*cfsrpd*vp12*vp34 + 32*pm_qp*pp_qm*pp_k*
c$$$     &    cfsrmd*vp12*vp34 + 16*pm_qp*pp_qm*pp_k*cisrmd*vp342 - 32*
c$$$     &    pm_qp*pp_qm*qm_k*cfsrmd*vp12*vp34 + 16*pm_qp*pp_qm*qp_k*
c$$$     &    cfsrpd*vp12*vp34 + 16*pm_qp*pp_qm*qp_k*cfsrmd*vp12*vp34 + 32*
c$$$     &    pm_qp*pp_qm*cisrpd*me2*vp342 + 32*pm_qp*pp_qm**2*cfsrmd*vp12*
c$$$     &    vp34 + 32*pm_qp*pp_qp*pp_k*cfsrpd*vp12*vp34 - 16*pm_qp*pp_qp*
c$$$     &    pp_k*cfsrmd*vp12*vp34 - 16*pm_qp*pp_qp*pp_k*cisrmd*vp342 + 16
c$$$     &    *pm_qp*pp_qp*qm_k*cfsrpd*vp12*vp34 )
c$$$      interf = interf + cisrp * ( 16*pm_qp*pp_qp*qm_k*cfsrmd*vp12*vp34
c$$$     &     - 32*pm_qp*pp_qp*qp_k*cfsrpd*vp12*vp34 - 32*pm_qp*pp_qp*
c$$$     &    cisrpd*me2*vp342 + 32*pm_qp*pp_qp**2*cfsrpd*vp12*vp34 - 16*
c$$$     &    pm_qp*pp_k*qm_qp*cfsrmd*vp12*vp34 + 16*pm_qp*pp_k*qm_k*cfsrpd
c$$$     &    *vp12*vp34 - 32*pm_qp*pp_k*qm_k*cisrpd*vp342 - 16*pm_qp*pp_k*
c$$$     &    qp_k*cfsrpd*vp12*vp34 + 32*pm_qp*pp_k*qp_k*cisrpd*vp342 + 8*
c$$$     &    pm_qp*pp_k*cfsrd*vp12*vp34 + 16*pm_qp*pp_k*cfsrmd*mpi2*vp12*
c$$$     &    vp34 + 8*pm_qp*pp_k**2*cfsrpd*vp12*vp34 + 8*pm_qp*pp_k**2*
c$$$     &    cfsrmd*vp12*vp34 - 32*pm_qp*qm_k*cisrpd*me2*vp342 + 32*pm_qp*
c$$$     &    qp_k*cisrpd*me2*vp342 + 8*pm_qp*cfsrd*me2*vp12*vp34 + 16*
c$$$     &    pm_qp**2*pp_k*cisrmd*vp342 - 16*pm_k*pp_qm*pp_qp*cfsrpd*vp12*
c$$$     &    vp34 + 16*pm_k*pp_qm*pp_qp*cfsrmd*vp12*vp34 - 32*pm_k*pp_qm*
c$$$     &    pp_qp*cisrmd*vp342 - 8*pm_k*pp_qm*pp_k*cfsrpd*vp12*vp34 - 8*
c$$$     &    pm_k*pp_qm*pp_k*cfsrmd*vp12*vp34 + 16*pm_k*pp_qm*qm_qp*cfsrmd
c$$$     &    *vp12*vp34 + 16*pm_k*pp_qm*qm_k*cfsrmd*vp12*vp34 + 16*pm_k*
c$$$     &    pp_qm*qp_k*cfsrpd*vp12*vp34 )
c$$$      interf = interf + cisrp * ( 8*pm_k*pp_qm*cfsrd*vp12*vp34 - 16*
c$$$     &    pm_k*pp_qm*cfsrmd*mpi2*vp12*vp34 - 16*pm_k*pp_qm**2*cfsrmd*
c$$$     &    vp12*vp34 + 16*pm_k*pp_qm**2*cisrmd*vp342 + 8*pm_k*pp_qp*pp_k
c$$$     &    *cfsrpd*vp12*vp34 + 8*pm_k*pp_qp*pp_k*cfsrmd*vp12*vp34 - 16*
c$$$     &    pm_k*pp_qp*qm_qp*cfsrpd*vp12*vp34 - 16*pm_k*pp_qp*qm_k*cfsrmd
c$$$     &    *vp12*vp34 - 16*pm_k*pp_qp*qp_k*cfsrpd*vp12*vp34 - 8*pm_k*
c$$$     &    pp_qp*cfsrd*vp12*vp34 + 16*pm_k*pp_qp*cfsrpd*mpi2*vp12*vp34
c$$$     &     + 16*pm_k*pp_qp**2*cfsrpd*vp12*vp34 + 16*pm_k*pp_qp**2*
c$$$     &    cisrmd*vp342 - 16*pm_k*pp_k*qm_qp*cfsrpd*vp12*vp34 + 16*pm_k*
c$$$     &    pp_k*qm_qp*cfsrmd*vp12*vp34 + 32*pm_k*pp_k*qm_qp*cisrpd*vp342
c$$$     &     + 16*pm_k*pp_k*cfsrpd*mpi2*vp12*vp34 - 16*pm_k*pp_k*cfsrmd*
c$$$     &    mpi2*vp12*vp34 - 32*pm_k*pp_k*cisrpd*mpi2*vp342 + 32*pm_k*
c$$$     &    qm_qp*cisrpd*me2*vp342 - 32*pm_k*cisrpd*me2*mpi2*vp342 - 32*
c$$$     &    pp_qm*qm_qp*cfsrmd*me2*vp12*vp34 + 16*pp_qm*qm_k*cfsrmd*me2*
c$$$     &    vp12*vp34 - 16*pp_qm*qp_k*cfsrmd*me2*vp12*vp34 + 8*pp_qm*
c$$$     &    cfsrd*me2*vp12*vp34 )
c$$$      interf = interf + cisrp * ( 32*pp_qm*cfsrmd*me2*mpi2*vp12*vp34 + 
c$$$     &    32*pp_qp*qm_qp*cfsrpd*me2*vp12*vp34 + 16*pp_qp*qm_k*cfsrpd*
c$$$     &    me2*vp12*vp34 - 16*pp_qp*qp_k*cfsrpd*me2*vp12*vp34 - 8*pp_qp*
c$$$     &    cfsrd*me2*vp12*vp34 - 32*pp_qp*cfsrpd*me2*mpi2*vp12*vp34 + 16
c$$$     &    *pp_k*qm_qp*cfsrpd*me2*vp12*vp34 - 16*pp_k*qm_qp*cfsrmd*me2*
c$$$     &    vp12*vp34 + 32*pp_k*qm_qp*cisrpd*me2*vp342 + 8*pp_k*qm_k*
c$$$     &    cfsrpd*me2*vp12*vp34 + 8*pp_k*qm_k*cfsrmd*me2*vp12*vp34 - 8*
c$$$     &    pp_k*qp_k*cfsrpd*me2*vp12*vp34 - 8*pp_k*qp_k*cfsrmd*me2*vp12*
c$$$     &    vp34 - 16*pp_k*cfsrpd*me2*mpi2*vp12*vp34 + 16*pp_k*cfsrmd*me2
c$$$     &    *mpi2*vp12*vp34 - 32*pp_k*cisrpd*me2*mpi2*vp342 + 16*qm_qp*
c$$$     &    qm_k*cfsrmd*me2*vp12*vp34 - 16*qm_qp*qp_k*cfsrpd*me2*vp12*
c$$$     &    vp34 - 32*qm_qp*cisrpd*me4*vp342 - 16*qm_k*qp_k*cfsrpd*me2*
c$$$     &    vp12*vp34 + 16*qm_k*qp_k*cfsrmd*me2*vp12*vp34 + 32*qm_k*qp_k*
c$$$     &    cisrmd*me2*vp342 - 16*qm_k*cfsrd*me2*vp12*vp34 - 16*qm_k*
c$$$     &    cfsrmd*me2*mpi2*vp12*vp34 - 16*qm_k**2*cfsrmd*me2*vp12*vp34
c$$$     &     - 16*qm_k**2*cisrmd*me2*vp342 )
c$$$      interf = interf + cisrp * ( 16*qp_k*cfsrd*me2*vp12*vp34 + 16*qp_k
c$$$     &    *cfsrpd*me2*mpi2*vp12*vp34 + 16*qp_k**2*cfsrpd*me2*vp12*vp34
c$$$     &     - 16*qp_k**2*cisrmd*me2*vp342 )
c$$$      interf = interf + cisrm * ( 32*cisrmd*me4*mpi2*vp342 + 32*pm_pp*
c$$$     &    pm_qm*pp_qm*cisrpd*vp342 - 32*pm_pp*pm_qm*pp_qp*cisrpd*vp342
c$$$     &     + 32*pm_pp*pm_qm*qm_qp*cfsrmd*vp12*vp34 - 16*pm_pp*pm_qm*
c$$$     &    qm_k*cfsrmd*vp12*vp34 - 16*pm_pp*pm_qm*qm_k*cisrpd*vp342 + 16
c$$$     &    *pm_pp*pm_qm*qp_k*cfsrmd*vp12*vp34 + 16*pm_pp*pm_qm*qp_k*
c$$$     &    cisrpd*vp342 - 32*pm_pp*pm_qm*cfsrmd*mpi2*vp12*vp34 - 32*
c$$$     &    pm_pp*pm_qp*pp_qm*cisrpd*vp342 + 32*pm_pp*pm_qp*pp_qp*cisrpd*
c$$$     &    vp342 - 32*pm_pp*pm_qp*qm_qp*cfsrpd*vp12*vp34 - 16*pm_pp*
c$$$     &    pm_qp*qm_k*cfsrpd*vp12*vp34 + 16*pm_pp*pm_qp*qm_k*cisrpd*
c$$$     &    vp342 + 16*pm_pp*pm_qp*qp_k*cfsrpd*vp12*vp34 - 16*pm_pp*pm_qp
c$$$     &    *qp_k*cisrpd*vp342 + 32*pm_pp*pm_qp*cfsrpd*mpi2*vp12*vp34 - 
c$$$     &    16*pm_pp*pm_k*qm_qp*cfsrpd*vp12*vp34 + 16*pm_pp*pm_k*qm_qp*
c$$$     &    cfsrmd*vp12*vp34 - 32*pm_pp*pm_k*qm_qp*cisrpd*vp342 - 8*pm_pp
c$$$     &    *pm_k*qm_k*cfsrpd*vp12*vp34 - 8*pm_pp*pm_k*qm_k*cfsrmd*vp12*
c$$$     &    vp34 + 8*pm_pp*pm_k*qp_k*cfsrpd*vp12*vp34 + 8*pm_pp*pm_k*qp_k
c$$$     &    *cfsrmd*vp12*vp34 )
c$$$      interf = interf + cisrm * ( 16*pm_pp*pm_k*cfsrpd*mpi2*vp12*vp34
c$$$     &     - 16*pm_pp*pm_k*cfsrmd*mpi2*vp12*vp34 + 32*pm_pp*pm_k*cisrpd
c$$$     &    *mpi2*vp342 - 16*pm_pp*pp_qm*qm_k*cisrpd*vp342 + 16*pm_pp*
c$$$     &    pp_qm*qp_k*cisrpd*vp342 + 16*pm_pp*pp_qp*qm_k*cisrpd*vp342 - 
c$$$     &    16*pm_pp*pp_qp*qp_k*cisrpd*vp342 - 32*pm_pp*pp_k*qm_qp*cisrpd
c$$$     &    *vp342 + 32*pm_pp*pp_k*cisrpd*mpi2*vp342 - 16*pm_pp*qm_qp*
c$$$     &    qm_k*cfsrmd*vp12*vp34 + 16*pm_pp*qm_qp*qp_k*cfsrpd*vp12*vp34
c$$$     &     + 32*pm_pp*qm_qp*cisrpd*me2*vp342 - 32*pm_pp*qm_qp*cisrmd*
c$$$     &    me2*vp342 + 16*pm_pp*qm_k*qp_k*cfsrpd*vp12*vp34 - 16*pm_pp*
c$$$     &    qm_k*qp_k*cfsrmd*vp12*vp34 + 8*pm_pp*qm_k*cfsrd*vp12*vp34 + 
c$$$     &    16*pm_pp*qm_k*cfsrmd*mpi2*vp12*vp34 + 16*pm_pp*qm_k**2*cfsrmd
c$$$     &    *vp12*vp34 - 8*pm_pp*qp_k*cfsrd*vp12*vp34 - 16*pm_pp*qp_k*
c$$$     &    cfsrpd*mpi2*vp12*vp34 - 16*pm_pp*qp_k**2*cfsrpd*vp12*vp34 - 
c$$$     &    32*pm_pp*cisrpd*me2*mpi2*vp342 + 32*pm_pp*cisrmd*me2*mpi2*
c$$$     &    vp342 + 32*pm_pp**2*qm_qp*cisrpd*vp342 - 32*pm_pp**2*cisrpd*
c$$$     &    mpi2*vp342 )
c$$$      interf = interf + cisrm * (  - 32*pm_qm*pm_qp*pp_qm*cfsrpd*vp12*
c$$$     &    vp34 - 32*pm_qm*pm_qp*pp_qm*cfsrmd*vp12*vp34 + 32*pm_qm*pm_qp
c$$$     &    *pp_qp*cfsrpd*vp12*vp34 + 32*pm_qm*pm_qp*pp_qp*cfsrmd*vp12*
c$$$     &    vp34 + 16*pm_qm*pm_qp*pp_k*cfsrpd*vp12*vp34 - 16*pm_qm*pm_qp*
c$$$     &    pp_k*cfsrmd*vp12*vp34 - 32*pm_qm*pm_qp*pp_k*cisrpd*vp342 - 16
c$$$     &    *pm_qm*pm_k*pp_qm*cfsrpd*vp12*vp34 + 32*pm_qm*pm_k*pp_qm*
c$$$     &    cfsrmd*vp12*vp34 - 16*pm_qm*pm_k*pp_qm*cisrpd*vp342 + 16*
c$$$     &    pm_qm*pm_k*pp_qp*cfsrpd*vp12*vp34 - 32*pm_qm*pm_k*pp_qp*
c$$$     &    cfsrmd*vp12*vp34 + 16*pm_qm*pm_k*pp_qp*cisrpd*vp342 + 8*pm_qm
c$$$     &    *pm_k*pp_k*cfsrpd*vp12*vp34 + 8*pm_qm*pm_k*pp_k*cfsrmd*vp12*
c$$$     &    vp34 - 16*pm_qm*pp_qm*pp_k*cisrpd*vp342 - 32*pm_qm*pp_qm*qm_k
c$$$     &    *cfsrmd*vp12*vp34 + 16*pm_qm*pp_qm*qp_k*cfsrpd*vp12*vp34 + 16
c$$$     &    *pm_qm*pp_qm*qp_k*cfsrmd*vp12*vp34 - 32*pm_qm*pp_qm*cisrmd*
c$$$     &    me2*vp342 + 16*pm_qm*pp_qp*pp_k*cisrpd*vp342 + 32*pm_qm*pp_qp
c$$$     &    *qm_k*cfsrmd*vp12*vp34 - 16*pm_qm*pp_qp*qp_k*cfsrpd*vp12*vp34
c$$$     &     - 16*pm_qm*pp_qp*qp_k*cfsrmd*vp12*vp34 )
c$$$      interf = interf + cisrm * ( 32*pm_qm*pp_qp*cisrmd*me2*vp342 - 16*
c$$$     &    pm_qm*pp_k*qm_qp*cfsrmd*vp12*vp34 - 16*pm_qm*pp_k*qm_k*cfsrmd
c$$$     &    *vp12*vp34 - 16*pm_qm*pp_k*qp_k*cfsrpd*vp12*vp34 - 8*pm_qm*
c$$$     &    pp_k*cfsrd*vp12*vp34 + 16*pm_qm*pp_k*cfsrmd*mpi2*vp12*vp34 + 
c$$$     &    32*pm_qm*qm_qp*cfsrmd*me2*vp12*vp34 - 16*pm_qm*qm_k*cfsrmd*
c$$$     &    me2*vp12*vp34 + 16*pm_qm*qp_k*cfsrmd*me2*vp12*vp34 - 8*pm_qm*
c$$$     &    cfsrd*me2*vp12*vp34 - 32*pm_qm*cfsrmd*me2*mpi2*vp12*vp34 + 32
c$$$     &    *pm_qm**2*pp_qm*cfsrmd*vp12*vp34 - 32*pm_qm**2*pp_qp*cfsrmd*
c$$$     &    vp12*vp34 + 16*pm_qm**2*pp_k*cfsrmd*vp12*vp34 + 16*pm_qm**2*
c$$$     &    pp_k*cisrpd*vp342 + 32*pm_qp*pm_k*pp_qm*cfsrpd*vp12*vp34 - 16
c$$$     &    *pm_qp*pm_k*pp_qm*cfsrmd*vp12*vp34 + 16*pm_qp*pm_k*pp_qm*
c$$$     &    cisrpd*vp342 - 32*pm_qp*pm_k*pp_qp*cfsrpd*vp12*vp34 + 16*
c$$$     &    pm_qp*pm_k*pp_qp*cfsrmd*vp12*vp34 - 16*pm_qp*pm_k*pp_qp*
c$$$     &    cisrpd*vp342 - 8*pm_qp*pm_k*pp_k*cfsrpd*vp12*vp34 - 8*pm_qp*
c$$$     &    pm_k*pp_k*cfsrmd*vp12*vp34 + 16*pm_qp*pp_qm*pp_k*cisrpd*vp342
c$$$     &     + 16*pm_qp*pp_qm*qm_k*cfsrpd*vp12*vp34 )
c$$$      interf = interf + cisrm * ( 16*pm_qp*pp_qm*qm_k*cfsrmd*vp12*vp34
c$$$     &     - 32*pm_qp*pp_qm*qp_k*cfsrpd*vp12*vp34 + 32*pm_qp*pp_qm*
c$$$     &    cisrmd*me2*vp342 - 16*pm_qp*pp_qp*pp_k*cisrpd*vp342 - 16*
c$$$     &    pm_qp*pp_qp*qm_k*cfsrpd*vp12*vp34 - 16*pm_qp*pp_qp*qm_k*
c$$$     &    cfsrmd*vp12*vp34 + 32*pm_qp*pp_qp*qp_k*cfsrpd*vp12*vp34 - 32*
c$$$     &    pm_qp*pp_qp*cisrmd*me2*vp342 + 16*pm_qp*pp_k*qm_qp*cfsrpd*
c$$$     &    vp12*vp34 + 16*pm_qp*pp_k*qm_k*cfsrmd*vp12*vp34 + 16*pm_qp*
c$$$     &    pp_k*qp_k*cfsrpd*vp12*vp34 + 8*pm_qp*pp_k*cfsrd*vp12*vp34 - 
c$$$     &    16*pm_qp*pp_k*cfsrpd*mpi2*vp12*vp34 - 32*pm_qp*qm_qp*cfsrpd*
c$$$     &    me2*vp12*vp34 - 16*pm_qp*qm_k*cfsrpd*me2*vp12*vp34 + 16*pm_qp
c$$$     &    *qp_k*cfsrpd*me2*vp12*vp34 + 8*pm_qp*cfsrd*me2*vp12*vp34 + 32
c$$$     &    *pm_qp*cfsrpd*me2*mpi2*vp12*vp34 + 32*pm_qp**2*pp_qm*cfsrpd*
c$$$     &    vp12*vp34 - 32*pm_qp**2*pp_qp*cfsrpd*vp12*vp34 - 16*pm_qp**2*
c$$$     &    pp_k*cfsrpd*vp12*vp34 + 16*pm_qp**2*pp_k*cisrpd*vp342 - 32*
c$$$     &    pm_k*pp_qm*pp_qp*cisrpd*vp342 - 16*pm_k*pp_qm*qm_qp*cfsrpd*
c$$$     &    vp12*vp34 )
c$$$      interf = interf + cisrm * (  - 16*pm_k*pp_qm*qm_k*cfsrmd*vp12*
c$$$     &    vp34 + 32*pm_k*pp_qm*qm_k*cisrmd*vp342 + 16*pm_k*pp_qm*qp_k*
c$$$     &    cfsrmd*vp12*vp34 - 32*pm_k*pp_qm*qp_k*cisrmd*vp342 + 8*pm_k*
c$$$     &    pp_qm*cfsrd*vp12*vp34 + 16*pm_k*pp_qm*cfsrpd*mpi2*vp12*vp34
c$$$     &     + 16*pm_k*pp_qm**2*cisrpd*vp342 + 16*pm_k*pp_qp*qm_qp*cfsrmd
c$$$     &    *vp12*vp34 - 16*pm_k*pp_qp*qm_k*cfsrpd*vp12*vp34 - 32*pm_k*
c$$$     &    pp_qp*qm_k*cisrmd*vp342 + 16*pm_k*pp_qp*qp_k*cfsrpd*vp12*vp34
c$$$     &     + 32*pm_k*pp_qp*qp_k*cisrmd*vp342 - 8*pm_k*pp_qp*cfsrd*vp12*
c$$$     &    vp34 - 16*pm_k*pp_qp*cfsrmd*mpi2*vp12*vp34 + 16*pm_k*pp_qp**2
c$$$     &    *cisrpd*vp342 + 16*pm_k*pp_k*qm_qp*cfsrpd*vp12*vp34 - 16*pm_k
c$$$     &    *pp_k*qm_qp*cfsrmd*vp12*vp34 + 32*pm_k*pp_k*qm_qp*cisrmd*
c$$$     &    vp342 - 16*pm_k*pp_k*cfsrpd*mpi2*vp12*vp34 + 16*pm_k*pp_k*
c$$$     &    cfsrmd*mpi2*vp12*vp34 - 32*pm_k*pp_k*cisrmd*mpi2*vp342 - 16*
c$$$     &    pm_k*qm_qp*cfsrpd*me2*vp12*vp34 + 16*pm_k*qm_qp*cfsrmd*me2*
c$$$     &    vp12*vp34 + 32*pm_k*qm_qp*cisrmd*me2*vp342 - 8*pm_k*qm_k*
c$$$     &    cfsrpd*me2*vp12*vp34 )
c$$$      interf = interf + cisrm * (  - 8*pm_k*qm_k*cfsrmd*me2*vp12*vp34
c$$$     &     + 8*pm_k*qp_k*cfsrpd*me2*vp12*vp34 + 8*pm_k*qp_k*cfsrmd*me2*
c$$$     &    vp12*vp34 + 16*pm_k*cfsrpd*me2*mpi2*vp12*vp34 - 16*pm_k*
c$$$     &    cfsrmd*me2*mpi2*vp12*vp34 - 32*pm_k*cisrmd*me2*mpi2*vp342 + 8
c$$$     &    *pm_k**2*pp_qm*cfsrpd*vp12*vp34 + 8*pm_k**2*pp_qm*cfsrmd*vp12
c$$$     &    *vp34 - 8*pm_k**2*pp_qp*cfsrpd*vp12*vp34 - 8*pm_k**2*pp_qp*
c$$$     &    cfsrmd*vp12*vp34 + 32*pp_qm*qm_k*cisrmd*me2*vp342 - 32*pp_qm*
c$$$     &    qp_k*cisrmd*me2*vp342 + 8*pp_qm*cfsrd*me2*vp12*vp34 - 32*
c$$$     &    pp_qp*qm_k*cisrmd*me2*vp342 + 32*pp_qp*qp_k*cisrmd*me2*vp342
c$$$     &     - 8*pp_qp*cfsrd*me2*vp12*vp34 + 32*pp_k*qm_qp*cisrmd*me2*
c$$$     &    vp342 - 32*pp_k*cisrmd*me2*mpi2*vp342 - 16*qm_qp*qm_k*cfsrmd*
c$$$     &    me2*vp12*vp34 + 16*qm_qp*qp_k*cfsrpd*me2*vp12*vp34 - 32*qm_qp
c$$$     &    *cisrmd*me4*vp342 + 16*qm_k*qp_k*cfsrpd*me2*vp12*vp34 - 16*
c$$$     &    qm_k*qp_k*cfsrmd*me2*vp12*vp34 + 32*qm_k*qp_k*cisrpd*me2*
c$$$     &    vp342 + 16*qm_k*cfsrd*me2*vp12*vp34 + 16*qm_k*cfsrmd*me2*mpi2
c$$$     &    *vp12*vp34 )
c$$$      interf = interf + cisrm * ( 16*qm_k**2*cfsrmd*me2*vp12*vp34 - 16*
c$$$     &    qm_k**2*cisrpd*me2*vp342 - 16*qp_k*cfsrd*me2*vp12*vp34 - 16*
c$$$     &    qp_k*cfsrpd*me2*mpi2*vp12*vp34 - 16*qp_k**2*cfsrpd*me2*vp12*
c$$$     &    vp34 - 16*qp_k**2*cisrpd*me2*vp342 )
c$$$

      Z1_=pm_pp
      Z2_=qm_qp
      Z3_=qm_k
      Z4_=qp_k
      Z5_=pm_qm
      Z6_=pp_qp
      Z7_=pp_k
      Z8_=pm_qp
      Z9_=pp_qm
      Z10_=pm_k
      Z11_=Z1_ + 2*me2
      Z11_=8*Z11_
      Z12_=3*Z7_
      Z13_=2*Z6_
      Z14_=Z13_ - Z9_
      Z15_=Z12_ + 2*Z14_
      Z15_=Z15_*Z8_
      Z16_=me2 + Z1_
      Z17_=Z2_ - mpi2
      Z18_=2*Z17_
      Z19_=Z18_ + Z3_
      Z20_= - Z19_ + 3*Z4_
      Z20_=Z16_*Z20_
      Z21_=2*Z7_
      Z22_=3*Z6_ + Z21_ - Z9_
      Z22_=Z22_*Z10_
      Z23_=Z13_ + Z7_
      Z23_=Z23_*Z5_
      Z15_= - Z15_ - Z22_ + Z23_ + Z20_
      Z15_=4*Z15_
      Z20_=2*Z9_
      Z22_=Z20_ - Z6_
      Z12_=Z12_ + 2*Z22_
      Z12_=Z12_*Z5_
      Z23_= - Z4_ - Z18_ + 3*Z3_
      Z23_=Z16_*Z23_
      Z21_=3*Z9_ + Z21_ - Z6_
      Z21_=Z21_*Z10_
      Z24_=Z20_ + Z7_
      Z24_=Z24_*Z8_
      Z12_= - Z12_ - Z21_ + Z24_ + Z23_
      Z12_=4*Z12_
      Z21_=Z8_ - Z5_
      Z23_=Z6_ - Z9_
      Z24_=Z23_ - Z21_
      Z25_=2*Z3_
      Z26_= - Z25_ + 2*Z4_
      Z27_=Z24_ - Z26_
      Z27_=Z27_*me2
      Z28_=Z8_*Z7_
      Z29_=Z5_*Z7_
      Z30_=Z28_ - Z29_
      Z31_=Z23_*Z10_
      Z32_= - Z31_ + Z30_
      Z33_=Z4_ - Z3_
      Z34_=Z33_*Z1_
      Z27_=Z27_ - Z32_ - Z34_
      Z27_=8*Z27_
      Z24_=Z24_ + Z26_
      Z24_=Z24_*me2
      Z24_=Z24_ - Z32_ + Z34_
      Z24_=8*Z24_
      Z26_= - Z2_ + 2*mpi2
      Z32_=Z26_ + Z33_
      Z32_=Z32_*Z4_
      Z34_=Z3_*mpi2
      Z35_=Z2_*mpi2
      Z35_=Z35_ - mpi4
      Z32_= - Z32_ + Z34_ + Z35_
      Z32_= - Z32_*Z16_
      Z34_=Z4_ + mpi2
      Z34_=Z34_*Z7_
      Z36_= - mpi2*Z23_
      Z37_=Z23_*Z4_
      Z34_= - Z34_ + Z36_ - Z37_
      Z36_=Z10_ + Z21_
      Z34_=Z34_*Z36_
      Z32_=Z34_ + Z32_
      Z32_=32*Z32_
      Z34_=Z16_*Z17_
      Z36_=Z10_*Z7_
      Z34_= - Z36_ + Z34_
      Z36_=Z3_ + 2*Z2_
      Z38_=Z36_ + Z4_
      Z34_=Z38_*Z34_
      Z36_= - Z36_*Z23_
      Z36_=Z36_ - Z37_
      Z36_=Z21_*Z36_
      Z34_= - Z36_ + Z34_
      Z34_=16*Z34_
      Z36_=Z17_ - Z23_
      Z36_=Z36_*Z6_
      Z36_=Z36_ + Z37_
      Z38_=Z23_ - Z18_
      Z38_=Z38_*Z7_
      Z36_= - Z38_ + 2*Z36_
      Z36_=Z36_*Z10_
      Z38_=Z17_ - Z14_
      Z38_= - Z7_ + 2*Z38_
      Z38_=Z38_*Z7_
      Z13_=Z13_ - Z20_
      Z39_=Z13_ + Z3_
      Z39_=Z39_*Z6_
      Z40_=Z14_*Z4_
      Z39_=Z39_ - Z40_
      Z39_=2*Z39_
      Z38_=Z38_ - Z39_
      Z38_=Z38_*Z5_
      Z14_=Z14_ - Z33_
      Z14_=Z7_ + 2*Z14_
      Z14_=Z14_*Z7_
      Z14_=Z14_ + Z39_
      Z14_=Z14_*Z8_
      Z39_=Z33_ - Z17_
      Z40_=Z39_ - Z6_
      Z40_=Z40_*Z4_
      Z19_=Z19_*Z6_
      Z19_=Z40_ + Z19_
      Z40_=Z18_ - Z33_
      Z41_=Z40_*Z7_
      Z19_=Z41_ + 2*Z19_
      Z19_=Z16_*Z19_
      Z14_=Z19_ - Z36_ + Z38_ + Z14_
      Z14_=8*Z14_
      Z19_=Z9_*Z3_
      Z36_=Z6_*Z3_
      Z38_=Z19_ - Z36_
      Z41_=Z17_ + Z4_
      Z41_=Z41_*Z7_
      Z42_=2*Z23_
      Z43_=Z42_ + Z7_
      Z44_=Z43_*Z8_
      Z43_=Z43_*Z5_
      Z41_= - Z38_ - Z41_ - 2*Z37_ + Z44_ - Z43_
      Z41_=Z41_*Z8_
      Z44_=Z7_*Z4_
      Z44_=Z44_ + Z37_
      Z44_=Z44_*Z5_
      Z41_=Z41_ + Z44_
      Z44_=Z4_*Z6_
      Z45_=Z17_*Z9_
      Z46_=Z17_*Z7_
      Z36_=Z36_ - Z44_ + Z45_ - Z46_
      Z44_=4*Z23_
      Z45_=Z44_ + Z7_
      Z45_=Z45_*Z8_
      Z36_=Z45_ - Z43_ + Z31_ + 2*Z36_
      Z36_=Z36_*Z10_
      Z43_=Z40_*Z8_
      Z39_=Z39_*Z4_
      Z39_=Z43_ + Z39_
      Z40_=Z40_*Z10_
      Z39_=Z40_ + 2*Z39_
      Z39_=Z16_*Z39_
      Z36_=2*Z41_ + Z36_ + Z39_
      Z36_=8*Z36_
      Z26_=Z26_ + Z3_
      Z26_=Z26_*Z3_
      Z39_=Z3_ + mpi2
      Z40_=Z39_*Z4_
      Z26_= - Z26_ + Z40_ + Z35_
      Z26_= - Z26_*Z16_
      Z35_= - Z23_ + Z7_
      Z40_= - Z10_ + Z21_
      Z35_=Z40_*Z39_*Z35_
      Z26_=Z35_ + Z26_
      Z26_=32*Z26_
      Z35_=Z9_ - Z3_
      Z39_=Z17_ - Z35_
      Z39_=Z39_*Z9_
      Z40_=Z35_*Z6_
      Z39_=Z39_ + Z40_
      Z40_=Z23_ + Z18_
      Z40_=Z40_*Z7_
      Z39_=Z40_ + 2*Z39_
      Z39_=Z39_*Z10_
      Z40_=Z35_*Z20_
      Z20_=Z20_ - Z3_
      Z20_=Z20_*Z6_
      Z41_=Z4_*Z9_
      Z20_=Z41_ + Z40_ - Z20_
      Z20_=2*Z20_
      Z40_=Z17_ - Z22_
      Z40_= - Z7_ + 2*Z40_
      Z40_=Z40_*Z7_
      Z40_=Z40_ - Z20_
      Z40_=Z40_*Z8_
      Z43_=Z22_ + Z33_
      Z43_=Z7_ + 2*Z43_
      Z43_=Z43_*Z7_
      Z20_=Z43_ + Z20_
      Z20_=Z20_*Z5_
      Z43_=Z18_ - Z3_
      Z43_=Z43_*Z9_
      Z35_=Z35_*Z4_
      Z45_=Z17_ - Z3_
      Z45_=Z45_*Z3_
      Z35_=Z43_ + Z35_ - Z45_
      Z43_=Z18_ + Z33_
      Z47_=Z43_*Z7_
      Z35_=Z47_ + 2*Z35_
      Z35_=Z16_*Z35_
      Z20_=Z35_ - Z39_ + Z40_ + Z20_
      Z20_=8*Z20_
      Z35_=Z17_*Z6_
      Z19_= - Z35_ - Z41_ + Z46_ + Z19_
      Z35_=Z44_ - Z7_
      Z35_=Z35_*Z5_
      Z39_=Z42_ - Z7_
      Z40_=Z39_*Z8_
      Z19_=2*Z19_ + Z35_ - Z40_ + Z31_
      Z19_=Z19_*Z10_
      Z31_=Z17_ + Z3_
      Z31_=Z31_*Z7_
      Z35_=Z39_*Z5_
      Z31_=Z31_ + 2*Z38_ + Z35_ + Z37_
      Z31_=Z31_*Z5_
      Z39_=Z7_*Z3_
      Z35_=Z39_ + Z35_ + Z38_
      Z35_=Z35_*Z8_
      Z31_=Z31_ - Z35_
      Z35_=Z43_*Z5_
      Z39_=Z4_*Z3_
      Z35_=Z45_ + Z39_ - Z35_
      Z39_=Z43_*Z10_
      Z35_= - Z39_ + 2*Z35_
      Z16_=Z16_*Z35_
      Z16_=Z16_ + Z19_ + 2*Z31_
      Z16_=8*Z16_
      Z19_=Z1_ - Z10_
      Z19_=Z19_*Z17_
      Z31_=Z23_ - Z33_
      Z31_= - Z31_*Z21_
      Z31_=Z46_ - Z19_ + Z31_
      Z31_=me2*Z31_
      Z30_=Z33_*Z30_
      Z17_=Z17_*me4
      Z35_=Z10_*Z46_
      Z30_=Z31_ + Z35_ + Z30_ - Z17_
      Z30_=32*Z30_
      Z31_=Z10_ + Z7_
      Z31_=Z31_*Z18_
      Z13_=Z13_ - Z33_
      Z13_=Z21_*Z13_
      Z33_=Z38_ + Z37_
      Z18_=Z18_*Z1_
      Z13_= - Z31_ + Z18_ + Z13_ - Z33_
      Z13_=Z13_*Z1_
      Z25_=Z25_ - Z4_
      Z25_=Z25_*Z4_
      Z31_=Z3_**2
      Z18_=Z25_ - Z31_ + Z18_
      Z18_=Z18_*me2
      Z25_=Z23_*Z7_
      Z28_=2*Z29_ + Z25_ - Z28_
      Z28_=Z28_*Z8_
      Z21_= - Z23_*Z21_
      Z22_=Z22_*Z6_
      Z23_=Z9_**2
      Z22_= - Z22_ + Z23_ + Z21_
      Z22_=Z22_*Z10_
      Z23_=Z25_ + Z29_
      Z23_=Z23_*Z5_
      Z13_= - Z23_ - Z13_ - Z18_ + Z28_ - Z22_
      Z13_=16*Z13_
      Z18_=Z33_ + Z46_
      Z19_=Z21_ - Z19_ + Z18_
      Z19_=me2*Z19_
      Z18_=Z10_*Z18_
      Z17_=Z19_ + Z18_ - Z17_
      Z17_=32*Z17_
      interf=cfsr*cfsrd*vp122*Z11_ + cfsr*cfsrpd*vp122*Z15_ + cfsr*
     & cfsrmd*vp122*Z12_ - cfsr*cisrpd*vp12*vp34*Z27_ - cfsr*cisrmd*
     & vp12*vp34*Z24_ + cfsrp*cfsrd*vp122*Z15_ + cfsrp*cfsrpd*vp122*
     & Z32_ + cfsrp*cfsrmd*vp122*Z34_ + cfsrp*cisrpd*vp12*vp34*Z14_ - 
     & cfsrp*cisrmd*vp12*vp34*Z36_ + cfsrm*cfsrd*vp122*Z12_ + cfsrm*
     & cfsrpd*vp122*Z34_ + cfsrm*cfsrmd*vp122*Z26_ - cfsrm*cisrpd*vp12*
     & vp34*Z20_ - cfsrm*cisrmd*vp12*vp34*Z16_ - cisrp*cfsrd*vp12*vp34*
     & Z27_ + cisrp*cfsrpd*vp12*vp34*Z14_ - cisrp*cfsrmd*vp12*vp34*Z20_
     &  + cisrp*cisrpd*vp342*Z30_ - cisrp*cisrmd*vp342*Z13_ - cisrm*
     & cfsrd*vp12*vp34*Z24_ - cisrm*cfsrpd*vp12*vp34*Z36_ - cisrm*
     & cfsrmd*vp12*vp34*Z16_ - cisrm*cisrpd*vp342*Z13_ + cisrm*cisrmd*
     & vp342*Z17_

*** STATS: original	0P 2222M 757A : 2979
*** STATS: optimized 0P 165M 172A : 337
      
