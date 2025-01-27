      Z1_=Dtend(0,0,0,0)
      Z2_=Dtenc(0,0,0,0)
      Z3_=dottencmom(Dtend,-2,2,p1c*cu,cmucmn*cu,c5e*cu)
      Z4_=dottencmom(Dtend,-2,2,p2c*cu,c5mu*cu,cecmn*cu)
      Z5_=dottencmom(Dtend,-1,2,zerop*cu,zerop*cu,zerop*cu)
      Z6_=dottencmom(Dtend,1,2,p1c*cu,zerop*cu,zerop*cu)
      Z7_=dottencmom(Dtend,1,2,p2c*cu,zerop*cu,zerop*cu)
      Z8_=dottencmom(Dtend,1,2,cmucmn*cu,zerop*cu,zerop*cu)
      Z9_=dottencmom(Dtend,1,2,cecmn*cu,zerop*cu,zerop*cu)
      Z10_=dottencmom(Dtend,2,2,c5e*cu,c5mu*cu,zerop*cu)
      Z11_=dottencmom(Dtend,2,2,cecmn*cu,cmucmn*cu,zerop*cu)
      Z12_=dottencmom(Dtenc,-2,2,p1c*cu,cmucmn*cu,c5e*cu)
      Z13_=dottencmom(Dtenc,-2,2,p4c*cu,c5mu*cu,cecmn*cu)
      Z14_=dottencmom(Dtenc,-1,2,zerop*cu,zerop*cu,zerop*cu)
      Z15_=dottencmom(Dtenc,1,2,p1c*cu,zerop*cu,zerop*cu)
      Z16_=dottencmom(Dtenc,1,2,p4c*cu,zerop*cu,zerop*cu)
      Z17_=dottencmom(Dtenc,1,2,cmucmn*cu,zerop*cu,zerop*cu)
      Z18_=dottencmom(Dtenc,1,2,cecmn*cu,zerop*cu,zerop*cu)
      Z19_=dottencmom(Dtenc,2,2,c5e*cu,c5mu*cu,zerop*cu)
      Z20_=dottencmom(Dtenc,2,2,cecmn*cu,cmucmn*cu,zerop*cu)
      Z21_=p1p2*Z1_
      Z22_=Z2_*p1p4
      Z21_=Z21_ + Z22_
      Z21_=Z15_ + Z7_ - Z6_ - Z5_ + Z16_ + Z14_ + 2*Z21_
      Z21_=cmuce*Z21_
      Z22_= - Z14_ - Z5_
      Z22_=c5muc5e*Z22_
      Z23_=Z17_ - Z8_
      Z23_=p1ce*Z23_
      Z24_=cImu*mm
      Z25_= - Z24_ + p4cmu
      Z25_=Z18_*Z25_
      Z24_=Z24_ + p2cmu
      Z24_=Z9_*Z24_
      Z26_= - Z17_ - Z8_
      Z26_=me*cIe*Z26_
      boxes=Z21_ + Z26_ + Z24_ + Z25_ + Z13_ + Z23_ + Z12_ + Z22_ + 
     & Z20_ + Z19_ - Z11_ + Z10_ + Z3_ + Z4_

*** STATS: original	0P 26M 25A : 51
*** STATS: optimized 0P 11M 25A : 36


