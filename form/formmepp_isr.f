*      
*      cisrp= (+im)*1.d0/((pp-k)^2-me2)/(qp+qm)^2      !c1 phase: +i
*      cisrm= (+im)*1.d0/((pm-k)^2-me2)/(qp+qm)^2      !c2 phase: +i
*
      cisrp= 1.d0/ppmk2mme2/qppqm2              !c1 phase: +i
      cisrm= 1.d0/pmmk2mme2/qppqm2              !c2 phase: +i

c$$$      cisrp2= abs((+im)*1.d0/ppmk2mme2/qppqm2)              !c1 phase: +i
c$$$      cisrm2= abs((+im)*1.d0/pmmk2mme2/qppqm2)              !c2 phase: +i
      cisrp2= cisrp*cisrp
      cisrm2= cisrm*cisrm
      
c$$$      interf =
c$$$     &  + cisrp2*vp342 * ( 32*me4*mpi2 - 32*pm_pp*qm_qp*me2 + 32*pm_pp*
c$$$     &    me2*mpi2 - 32*pm_qm*pp_qm*me2 + 32*pm_qm*pp_qp*me2 + 32*pm_qm
c$$$     &    *pp_k*qm_k - 32*pm_qm*pp_k*qp_k + 32*pm_qm*qm_k*me2 - 32*
c$$$     &    pm_qm*qp_k*me2 + 32*pm_qp*pp_qm*me2 - 32*pm_qp*pp_qp*me2 - 32
c$$$     &    *pm_qp*pp_k*qm_k + 32*pm_qp*pp_k*qp_k - 32*pm_qp*qm_k*me2 + 
c$$$     &    32*pm_qp*qp_k*me2 + 32*pm_k*pp_k*qm_qp - 32*pm_k*pp_k*mpi2 + 
c$$$     &    32*pm_k*qm_qp*me2 - 32*pm_k*me2*mpi2 + 32*pp_k*qm_qp*me2 - 32
c$$$     &    *pp_k*me2*mpi2 - 32*qm_qp*me4 )
c$$$      interf = interf + cisrm2*vp342 * ( 32*me4*mpi2 - 32*pm_pp*qm_qp*
c$$$     &    me2 + 32*pm_pp*me2*mpi2 - 32*pm_qm*pp_qm*me2 + 32*pm_qm*pp_qp
c$$$     &    *me2 + 32*pm_qp*pp_qm*me2 - 32*pm_qp*pp_qp*me2 + 32*pm_k*
c$$$     &    pp_qm*qm_k - 32*pm_k*pp_qm*qp_k - 32*pm_k*pp_qp*qm_k + 32*
c$$$     &    pm_k*pp_qp*qp_k + 32*pm_k*pp_k*qm_qp - 32*pm_k*pp_k*mpi2 + 32
c$$$     &    *pm_k*qm_qp*me2 - 32*pm_k*me2*mpi2 + 32*pp_qm*qm_k*me2 - 32*
c$$$     &    pp_qm*qp_k*me2 - 32*pp_qp*qm_k*me2 + 32*pp_qp*qp_k*me2 + 32*
c$$$     &    pp_k*qm_qp*me2 - 32*pp_k*me2*mpi2 - 32*qm_qp*me4 )
c$$$      interf = interf + cisrm*cisrp*vp342 * ( 64*pm_pp*pm_qm*pp_qm - 64
c$$$     &    *pm_pp*pm_qm*pp_qp - 32*pm_pp*pm_qm*qm_k + 32*pm_pp*pm_qm*
c$$$     &    qp_k - 64*pm_pp*pm_qp*pp_qm + 64*pm_pp*pm_qp*pp_qp + 32*pm_pp
c$$$     &    *pm_qp*qm_k - 32*pm_pp*pm_qp*qp_k - 64*pm_pp*pm_k*qm_qp + 64*
c$$$     &    pm_pp*pm_k*mpi2 - 32*pm_pp*pp_qm*qm_k + 32*pm_pp*pp_qm*qp_k
c$$$     &     + 32*pm_pp*pp_qp*qm_k - 32*pm_pp*pp_qp*qp_k - 64*pm_pp*pp_k*
c$$$     &    qm_qp + 64*pm_pp*pp_k*mpi2 + 64*pm_pp*qm_qp*me2 - 64*pm_pp*
c$$$     &    me2*mpi2 + 64*pm_pp2*qm_qp - 64*pm_pp2*mpi2 - 64*pm_qm*
c$$$     &    pm_qp*pp_k - 32*pm_qm*pm_k*pp_qm + 32*pm_qm*pm_k*pp_qp - 32*
c$$$     &    pm_qm*pp_qm*pp_k + 32*pm_qm*pp_qp*pp_k + 32*pm_qm2*pp_k + 
c$$$     &    32*pm_qp*pm_k*pp_qm - 32*pm_qp*pm_k*pp_qp + 32*pm_qp*pp_qm*
c$$$     &    pp_k - 32*pm_qp*pp_qp*pp_k + 32*pm_qp2*pp_k - 64*pm_k*pp_qm
c$$$     &    *pp_qp + 32*pm_k*pp_qm2 + 32*pm_k*pp_qp2 + 64*qm_k*qp_k*
c$$$     &    me2 - 32*qm_k2*me2 - 32*qp_k2*me2 )

      Z1_=pm_pp
      Z2_=qm_qp
      Z3_=pm_qm
      Z4_=pp_qm
      Z5_=pp_qp
      Z6_=pp_k
      Z7_=qm_k
      Z8_=qp_k
      Z9_=pm_qp
      Z10_=pm_k
      Z11_= - Z1_ + Z10_ + Z6_
      Z12_=Z2_ - mpi2
      Z11_=Z11_*Z12_
      Z13_=Z7_ - Z8_
      Z14_=Z5_ - Z4_
      Z15_=Z13_ + Z14_
      Z16_=Z3_ - Z9_
      Z15_=Z16_*Z15_
      Z15_=Z11_ + Z15_
      Z15_=me2*Z15_
      Z17_=Z16_*Z13_
      Z18_=Z10_*Z12_
      Z17_=Z18_ + Z17_
      Z17_=Z6_*Z17_
      Z18_=Z12_*me4
      Z15_=Z15_ - Z18_ + Z17_
      Z15_=32*Z15_
      Z17_=Z16_ - Z13_
      Z17_=Z14_*Z17_
      Z11_=Z17_ + Z11_
      Z11_=me2*Z11_
      Z17_= - Z14_*Z13_
      Z19_=Z6_*Z12_
      Z17_=Z19_ + Z17_
      Z17_=Z10_*Z17_
      Z11_=Z11_ - Z18_ + Z17_
      Z11_=32*Z11_
      Z17_=Z16_ - Z14_
      Z13_= - Z17_*Z13_
      Z16_=Z16_*Z14_
      Z17_=2*Z1_
      Z12_=Z17_*Z12_
      Z13_=Z12_ - 2*Z16_ + Z13_
      Z13_=Z1_*Z13_
      Z17_=Z4_**2
      Z18_= - 2*Z4_ + Z5_
      Z18_=Z5_*Z18_
      Z16_= - Z12_ + Z17_ + Z18_ + Z16_
      Z16_=Z10_*Z16_
      Z17_=Z9_ - Z14_
      Z17_=Z9_*Z17_
      Z14_=Z3_ - 2*Z9_ + Z14_
      Z14_=Z3_*Z14_
      Z14_= - Z12_ + Z17_ + Z14_
      Z14_=Z6_*Z14_
      Z17_=Z8_**2
      Z18_=2*Z8_ - Z7_
      Z18_=Z7_*Z18_
      Z12_=Z18_ - Z17_ + Z12_
      Z12_=me2*Z12_
      Z12_=Z12_ + Z16_ + Z14_ + Z13_
      Z12_=32*Z12_
      interf=cisrp2*vp342*Z15_ + cisrm2*vp342*Z11_ + cisrm*cisrp*vp342*
     & Z12_

*** STATS: original	0P 239M 78A : 317
*** STATS: optimized 0P 33M 34A : 67
