*      
*      cfsr = (-im)*(2.d0*1.d0/(qp+qm+k)^2)            !c5 phase: -i
*      cfsrp= (+im)*1.d0/((qp+k)^2-mpi2)/((qp+qm+k)^2) !c3 phase: +i
*      cfsrm= (+im)*1.d0/((qm+k)^2-mpi2)/((qp+qm+k)^2) !c4 phase: +i
*      cisrp= (+im)*1.d0/((pp-k)^2-me2)/(qp+qm)^2      !c1 phase: +i
*      cisrm= (+im)*1.d0/((pm-k)^2-me2)/(qp+qm)^2      !c2 phase: +i
*
      im= (0.d0,1.d0)
      cfsr = (-im)*2.d0*1.d0/qppqmpk2                 !c5 phase: -i
      cfsrp= (+im)*1.d0/qppk2mmpi2/qppqmpk2           !c3 phase: +i
      cfsrm= (+im)*1.d0/qmpk2mmpi2/qppqmpk2           !c4 phase: +i
      cisrp= (+im)*1.d0/ppmk2mme2/qppqm2              !c1 phase: +i
      cisrm= (+im)*1.d0/pmmk2mme2/qppqm2              !c2 phase: +i

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
     &  + cisrp * ( 32*cisrpd*me4*mpi2*vp342 + 32*pm_pp*pm_qm*pp_qm*
     &    cisrmd*vp342 - 32*pm_pp*pm_qm*pp_qp*cisrmd*vp342 - 16*pm_pp*
     &    pm_qm*qm_k*cisrmd*vp342 + 16*pm_pp*pm_qm*qp_k*cisrmd*vp342 - 
     &    32*pm_pp*pm_qp*pp_qm*cisrmd*vp342 + 32*pm_pp*pm_qp*pp_qp*
     &    cisrmd*vp342 + 16*pm_pp*pm_qp*qm_k*cisrmd*vp342 - 16*pm_pp*
     &    pm_qp*qp_k*cisrmd*vp342 - 32*pm_pp*pm_k*qm_qp*cisrmd*vp342 + 
     &    32*pm_pp*pm_k*cisrmd*mpi2*vp342 - 16*pm_pp*pp_qm*qm_k*cisrmd*
     &    vp342 + 16*pm_pp*pp_qm*qp_k*cisrmd*vp342 + 16*pm_pp*pp_qp*
     &    qm_k*cisrmd*vp342 - 16*pm_pp*pp_qp*qp_k*cisrmd*vp342 - 32*
     &    pm_pp*pp_k*qm_qp*cisrmd*vp342 + 32*pm_pp*pp_k*cisrmd*mpi2*
     &    vp342 - 32*pm_pp*qm_qp*cisrpd*me2*vp342 + 32*pm_pp*qm_qp*
     &    cisrmd*me2*vp342 + 32*pm_pp*cisrpd*me2*mpi2*vp342 - 32*pm_pp*
     &    cisrmd*me2*mpi2*vp342 + 32*pm_pp**2*qm_qp*cisrmd*vp342 - 32*
     &    pm_pp**2*cisrmd*mpi2*vp342 - 32*pm_qm*pm_qp*pp_k*cisrmd*vp342
     &     - 16*pm_qm*pm_k*pp_qm*cisrmd*vp342 + 16*pm_qm*pm_k*pp_qp*
     &    cisrmd*vp342 )
      interf = interf + cisrp * (  - 16*pm_qm*pp_qm*pp_k*cisrmd*vp342
     &     - 32*pm_qm*pp_qm*cisrpd*me2*vp342 + 16*pm_qm*pp_qp*pp_k*
     &    cisrmd*vp342 + 32*pm_qm*pp_qp*cisrpd*me2*vp342 + 32*pm_qm*
     &    pp_k*qm_k*cisrpd*vp342 - 32*pm_qm*pp_k*qp_k*cisrpd*vp342 + 32
     &    *pm_qm*qm_k*cisrpd*me2*vp342 - 32*pm_qm*qp_k*cisrpd*me2*vp342
     &     + 16*pm_qm**2*pp_k*cisrmd*vp342 + 16*pm_qp*pm_k*pp_qm*cisrmd
     &    *vp342 - 16*pm_qp*pm_k*pp_qp*cisrmd*vp342 + 16*pm_qp*pp_qm*
     &    pp_k*cisrmd*vp342 + 32*pm_qp*pp_qm*cisrpd*me2*vp342 - 16*
     &    pm_qp*pp_qp*pp_k*cisrmd*vp342 - 32*pm_qp*pp_qp*cisrpd*me2*
     &    vp342 - 32*pm_qp*pp_k*qm_k*cisrpd*vp342 + 32*pm_qp*pp_k*qp_k*
     &    cisrpd*vp342 - 32*pm_qp*qm_k*cisrpd*me2*vp342 + 32*pm_qp*qp_k
     &    *cisrpd*me2*vp342 + 16*pm_qp**2*pp_k*cisrmd*vp342 - 32*pm_k*
     &    pp_qm*pp_qp*cisrmd*vp342 + 16*pm_k*pp_qm**2*cisrmd*vp342 + 16
     &    *pm_k*pp_qp**2*cisrmd*vp342 + 32*pm_k*pp_k*qm_qp*cisrpd*vp342
     &     - 32*pm_k*pp_k*cisrpd*mpi2*vp342 + 32*pm_k*qm_qp*cisrpd*me2*
     &    vp342 )
      interf = interf + cisrp * (  - 32*pm_k*cisrpd*me2*mpi2*vp342 + 32
     &    *pp_k*qm_qp*cisrpd*me2*vp342 - 32*pp_k*cisrpd*me2*mpi2*vp342
     &     - 32*qm_qp*cisrpd*me4*vp342 + 32*qm_k*qp_k*cisrmd*me2*vp342
     &     - 16*qm_k**2*cisrmd*me2*vp342 - 16*qp_k**2*cisrmd*me2*vp342
     &     )
      interf = interf + cisrm * ( 32*cisrmd*me4*mpi2*vp342 + 32*pm_pp*
     &    pm_qm*pp_qm*cisrpd*vp342 - 32*pm_pp*pm_qm*pp_qp*cisrpd*vp342
     &     - 16*pm_pp*pm_qm*qm_k*cisrpd*vp342 + 16*pm_pp*pm_qm*qp_k*
     &    cisrpd*vp342 - 32*pm_pp*pm_qp*pp_qm*cisrpd*vp342 + 32*pm_pp*
     &    pm_qp*pp_qp*cisrpd*vp342 + 16*pm_pp*pm_qp*qm_k*cisrpd*vp342
     &     - 16*pm_pp*pm_qp*qp_k*cisrpd*vp342 - 32*pm_pp*pm_k*qm_qp*
     &    cisrpd*vp342 + 32*pm_pp*pm_k*cisrpd*mpi2*vp342 - 16*pm_pp*
     &    pp_qm*qm_k*cisrpd*vp342 + 16*pm_pp*pp_qm*qp_k*cisrpd*vp342 + 
     &    16*pm_pp*pp_qp*qm_k*cisrpd*vp342 - 16*pm_pp*pp_qp*qp_k*cisrpd
     &    *vp342 - 32*pm_pp*pp_k*qm_qp*cisrpd*vp342 + 32*pm_pp*pp_k*
     &    cisrpd*mpi2*vp342 + 32*pm_pp*qm_qp*cisrpd*me2*vp342 - 32*
     &    pm_pp*qm_qp*cisrmd*me2*vp342 - 32*pm_pp*cisrpd*me2*mpi2*vp342
     &     + 32*pm_pp*cisrmd*me2*mpi2*vp342 + 32*pm_pp**2*qm_qp*cisrpd*
     &    vp342 - 32*pm_pp**2*cisrpd*mpi2*vp342 - 32*pm_qm*pm_qp*pp_k*
     &    cisrpd*vp342 - 16*pm_qm*pm_k*pp_qm*cisrpd*vp342 + 16*pm_qm*
     &    pm_k*pp_qp*cisrpd*vp342 )
      interf = interf + cisrm * (  - 16*pm_qm*pp_qm*pp_k*cisrpd*vp342
     &     - 32*pm_qm*pp_qm*cisrmd*me2*vp342 + 16*pm_qm*pp_qp*pp_k*
     &    cisrpd*vp342 + 32*pm_qm*pp_qp*cisrmd*me2*vp342 + 16*pm_qm**2*
     &    pp_k*cisrpd*vp342 + 16*pm_qp*pm_k*pp_qm*cisrpd*vp342 - 16*
     &    pm_qp*pm_k*pp_qp*cisrpd*vp342 + 16*pm_qp*pp_qm*pp_k*cisrpd*
     &    vp342 + 32*pm_qp*pp_qm*cisrmd*me2*vp342 - 16*pm_qp*pp_qp*pp_k
     &    *cisrpd*vp342 - 32*pm_qp*pp_qp*cisrmd*me2*vp342 + 16*pm_qp**2
     &    *pp_k*cisrpd*vp342 - 32*pm_k*pp_qm*pp_qp*cisrpd*vp342 + 32*
     &    pm_k*pp_qm*qm_k*cisrmd*vp342 - 32*pm_k*pp_qm*qp_k*cisrmd*
     &    vp342 + 16*pm_k*pp_qm**2*cisrpd*vp342 - 32*pm_k*pp_qp*qm_k*
     &    cisrmd*vp342 + 32*pm_k*pp_qp*qp_k*cisrmd*vp342 + 16*pm_k*
     &    pp_qp**2*cisrpd*vp342 + 32*pm_k*pp_k*qm_qp*cisrmd*vp342 - 32*
     &    pm_k*pp_k*cisrmd*mpi2*vp342 + 32*pm_k*qm_qp*cisrmd*me2*vp342
     &     - 32*pm_k*cisrmd*me2*mpi2*vp342 + 32*pp_qm*qm_k*cisrmd*me2*
     &    vp342 - 32*pp_qm*qp_k*cisrmd*me2*vp342 - 32*pp_qp*qm_k*cisrmd
     &    *me2*vp342 )
      interf = interf + cisrm * ( 32*pp_qp*qp_k*cisrmd*me2*vp342 + 32*
     &    pp_k*qm_qp*cisrmd*me2*vp342 - 32*pp_k*cisrmd*me2*mpi2*vp342
     &     - 32*qm_qp*cisrmd*me4*vp342 + 32*qm_k*qp_k*cisrpd*me2*vp342
     &     - 16*qm_k**2*cisrpd*me2*vp342 - 16*qp_k**2*cisrpd*me2*vp342
     &     )

