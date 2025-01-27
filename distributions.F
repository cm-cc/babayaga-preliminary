*********************************************************************
      subroutine setfilenames(programma,corrections)
#include "shared.F"
      character*100  distname(ndistr),tmp(ndistr)
      common/filenames/distname
      character*(*) programma,corrections
      character*2   ec,ac
      character*7   bw

      do k = 1,7
         bw(k:k) = ' '
      enddo
      do k = 1,ndistr
         do i = 1,100
            tmp(k)(i:i) = ' '
         enddo
      enddo

      ebeam = ecmsnom/2.d0

      bw='.txt'

******************
#if CMD
      tmp(1)  = 'p+_'
      tmp(2)  = 'p-_'
      tmp(3)  = 'p+c_'
      tmp(4)  = 'p-c_'
      tmp(5)  = 'thav_'
      tmp(6)  = 'dphi_'
      tmp(7)  = 'xi_'
      tmp(8)  = 'mxx_'
      tmp(9)  = 'cth+_'
      tmp(10) = 'cth-_'
#elif KLOEI
      tmp(1)  = 'lth+_'
      tmp(2)  = 'lth-_'
      tmp(3)  = 'lp+_'
      tmp(4)  = 'lp-_'
      tmp(5)  = 'lpz+_'
      tmp(6)  = 'lpz-_'
      tmp(7)  = 'lpp+_'
      tmp(8)  = 'lpp-_'
      tmp(9)  = 'lmxx_'
      tmp(10) = 'lxi_'
      tmp(11) = 'ldphi_'
      tmp(12) = 'lthav_'
      tmp(13) = 'lthg_'
      tmp(14) = 'leg_'
      tmp(15) = 'lmxxg_'
      tmp(16) = 'lmtrk_'
      tmp(17) = 'lmxxc_'
      tmp(18) = 'lmtrkc_'
#elif KLOEII
      tmp(1)  = 'sth+_'
      tmp(2)  = 'sth-_'
      tmp(3)  = 'sp+_'
      tmp(4)  = 'sp-_'
      tmp(5)  = 'spz+_'
      tmp(6)  = 'spz-_'
      tmp(7)  = 'spp+_'
      tmp(8)  = 'spp-_'
      tmp(9)  = 'smxx_'
      tmp(10) = 'smtrk_'
      tmp(11) = 'sxi_'
      tmp(12) = 'sdphi_'
      tmp(13) = 'sthav_'
      tmp(14) = 'sthtgl_'
      tmp(15) = 'sthtgs_'
      tmp(16) = 'smxxc_'
      tmp(17) = 'smtrkc_'
#elif BES3
      tmp(1)  = 'cth+_'
      tmp(2)  = 'cth-_'
      tmp(3)  = 'pp+_'
      tmp(4)  = 'pp-_'
      tmp(5)  = 'p+_'
      tmp(6)  = 'p-_'
      tmp(7)  = 'eg_'
      tmp(8)  = 'dphi_'
      tmp(9)  = 'xi_'
      tmp(10) = 'mxx_'
      tmp(11) = 'mxxg_'
      tmp(12) = 'cthg_'
      tmp(13) = 'thav_'
      tmp(14) = 'mxxrw_'
#elif B
      tmp(1)  = 'p+_'
      tmp(2)  = 'th+_'
      tmp(3)  = 'p-_'
      tmp(4)  = 'th-_'
      tmp(5)  = 'Eg_'
      tmp(6)  = 'thg_'
      tmp(7)  = 'thav_'
      tmp(8)  = 'dphi_'
      tmp(9)  = 'xi_'
      tmp(10) = 'mxx_'
      tmp(11) = 'mxxrw_'
      tmp(12) = 'mxxg_'
      tmp(13) = 'thgg_'
#else
! default
      tmp(1)  = 'el_th_'
      tmp(2)  = 'el_en_'
      tmp(3)  = 'acoll_'
      tmp(4)  = 'g1_en_'
      tmp(5)  = 'g1_th_'
      tmp(6)  = 'ep_th_'
      tmp(7)  = 'ep_en_'
      tmp(8)  = 'e+e-m_'
      tmp(9)  = 'g2_en_'
      tmp(10) = 'g2_th_'
      tmp(11) = 'e1r_th_'
      tmp(12) = 'e1r_en_'
      tmp(13) = 'e2r_th_'
      tmp(14) = 'e2r_en_'
      tmp(15) = 'aco2_'
      tmp(16) = 't_'
      tmp(17) = 'cthe_'      
      tmp(18) = 'cthp_'
#endif
******************
      do k=1,ndistr
         lung(k) = len_trim(tmp(k))
      enddo

      if (corrections.eq.'oal') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'NLO'//bw
         enddo
      endif
      if (corrections.eq.'oal2') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'NNLO'//bw
         enddo
      endif

      if (corrections.eq.'exp') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'EXP'//bw
         enddo
      endif

      if (corrections.eq.'born') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'LO'//bw
         enddo
      endif

      if (corrections.eq.'struct') then
         do k=1,ndistr
            tmp(k)(lung(k)+1:) = 'stru_'//bw
         enddo
      endif

      do k = 1,len(programma)
         if (programma(k:k).ne.' ') lpath=k
      enddo

      do k = 1,ndistr
         distname(k)(1:lpath)  = programma
         distname(k)(lpath+1:) = tmp(k)
      enddo

      return
      end
*********************************************************************
      subroutine distributions(sd,nc,p3,p4,qph)
#include "shared.F"
      parameter (ibin = 3 * nbin)      
      parameter (iperj = ndistr*ibin)
      dimension pin1(0:3),pin2(0:3)
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b(0:3),p2b(0:3),iiiref
      character*2 fs
      common/finalstate/fs
      common/idebugging/idebug      
      common/distr/distr,s_distr,bmi,bma 
      common/ionlyfirsttime/ifirst,icount
      data ifirst,icount /1,0/ 

      common/tmintmax/tmin,tmax,tcur

      common/beamspreadsigmas/esig1,esig2

      common/parameters/ame,ammu,convfac,alpha,pi_i
      common/hparam/ampi
            
      common/par/pi,radtodeg,ebeam
      common/idarkon/idarkon
      dimension s_distr(ndistr,ibin)
      dimension distr(ndistr,ibin),bmi(ndistr),bma(ndistr),fun(ndistr)

      data ((distr(i,j),i=1,ndistr),j=1,ibin) / iperj * 0.d0/
      data ((s_distr(i,j),i=1,ndistr),j=1,ibin) / iperj * 0.d0/

      integer*8 nc
      double precision csi(1)
      
      double precision qa(0:3),qb(0:3),px(0:3),py(0:3)
      double precision pluto(0:3)
      common/massainv/amassainvmin,amassainvmax

      double precision qphtot(40,0:3)
      double precision p1inside(0:3),p2inside(0:3)
      double precision p1rad(0:3),p2rad(0:3)
      integer inside(40),icutgg
      common/cutsggcmn/qphtot,p1inside,p2inside,p1rad,p2rad,
     .     nph,npairs,inside,n_inside,icutgg
#include "strong2020common.F"

      ebeam = ecmsnom/2.d0      
#if CMD      
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi

         bmi(1)  =  0.4d0
         bma(1)  =  1.d0         
         bmi(2)  = bmi(1)
         bma(2)  = bma(1)         
         bmi(3)  = bmi(1)
         bma(3)  = bma(1)         
         bmi(4)  = bmi(1)
         bma(4)  = bma(1)         
         bmi(5)  =  1.d0
         bma(5)  =  pi - 1.d0
         bmi(6)  =  0.d0
         bma(6)  =  0.15d0
         bmi(7)  =  0.d0
         bma(7)  =  0.25d0
         bmi(8)  =  200.d0
         bma(8)  =  700.d0         
         bmi(9)  = -1.d0
         bma(9)  = -bmi(9)
         bmi(10) =  bmi(9)
         bma(10) =  bma(9)

! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo
         ifirst = 0
      endif
******************
      fun(1)  = 2.d0*pp/ecms
      fun(2)  = 2.d0*pm/ecms
      fun(3)  = bmi(3) - 10.d0
      fun(4)  = bmi(4) - 10.d0
      if (fun(2).gt.0.45d0.and.fun(2).lt.0.6d0) fun(3) = fun(1)
      if (fun(1).gt.0.45d0.and.fun(1).lt.0.6d0) fun(4) = fun(2)      
      fun(5)  = thavg
      fun(6)  = dphi
      fun(7)  = z
      fun(8)  = sqrt(dot(p3+p4,p3+p4))*1d3      
      fun(9)  = cos(thp)
      fun(10) = cos(thm)      
#elif KLOEI
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi

         bmi(1)  =  30d0
         bma(1)  = 130d0
         bmi(2)  =  bmi(1)          
         bma(2)  =  bma(1)
         bmi(3)  =   0d0
         bma(3)  = 600d0
         bmi(4)  = bmi(3)
         bma(4)  = bma(3)
         bmi(5)  = bmi(3)
         bma(5)  = bma(3)
         bmi(6)  = bmi(3)
         bma(6)  = bma(3)
         bmi(7)  =   7.d0
         bma(7)  = 557.d0
         bmi(8)  = bmi(7)
         bma(8)  = bma(7)
         bmi(9)  = 316d0
         bma(9)  = 922d0
         bmi(10) = 0.d0
         bma(10) = 80.d0
         bmi(11) = 0.d0
         bma(11) = 180d0
         bmi(12) = 30d0
         bma(12) = 130d0
         bmi(13) = 30d0
         bma(13) = 130d0
         bmi(14) = 0d0
         bma(14) = 500d0
         bmi(15) = 0d0
         bma(15) = 2000d0
         bmi(16) = 0.d0
         bma(16) = 500d0
         bmi(17) = 316d0
         bma(17) = 922d0
         bmi(18) = 0.d0
         bma(18) = 500d0
! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo
         ifirst = 0
      endif
******************
      fun(1) = thp*180d0/pi
      fun(2) = thm*180d0/pi
      fun(3) = pp*1d3
      fun(4) = pm*1d3
      fun(5) = pzp*1d3
      fun(6) = pzm*1d3
      fun(7) = pperpp*1d3
      fun(8) = pperpm*1d3
      fun(9) = mxx*1d3
      fun(10) = abs(thp+thm - pi)*180d0/pi

      phip = atan2(p3(2),p3(1))
      phim = atan2(p4(2),p4(1))      
      dphi = abs(abs(phip-phim) - pi)
      fun(11) = dphi*180d0/pi

      thav = (thm-thp+pi)*0.5d0
      fun(12) = thav * 180d0/pi

      fun(13:18) = bmi(13:18) -10.d0

      call mostenergeticphoton(qinside,pluto)
      fun(13) = acos(pluto(3)/pluto(0))*180d0/pi
      fun(14) = pluto(0)*1d3
      fun(15) = sqrt(dot(p3+p4+pluto,p3+p4+pluto))*1d3
      
      
      if (fs.eq.'ee') then
         trkm    = ame*1d3
         if (qph(2,0).gt.0.d0) then
            call getmasstrack(2,ecms,p3,p4,ame*ame,trkm2)
            trkm = sqrt(trkm2)*1d3
         endif
         if (trkm.gt.130d0.and.trkm.lt.220d0) then
            fun(17) = fun(9)
         endif
      elseif(fs.eq.'mm') then
         trkm    = ammu*1d3
         if (qph(2,0).gt.0.d0) then
            call getmasstrack(2,ecms,p3,p4,ammu*ammu,trkm2)
            trkm = sqrt(trkm2)*1d3
         endif         
!         if (trkm.gt.98d0.and.trkm.lt.166d0) then !!! WRONG
         if (trkm.gt.80d0.and.trkm.lt.115d0) then
            fun(17) = fun(9)
         endif
      elseif(fs.eq.'pp'.or.fs.eq.'pi') then
         trkm    = ampi*1d3
         if (qph(2,0).gt.0.d0) then
            call getmasstrack(2,ecms,p3,p4,ampi*ampi,trkm2)
            trkm = sqrt(trkm2)*1d3
         endif
         if (trkm.gt.130d0.and.trkm.lt.220d0) then
            fun(17) = fun(9)
         endif            
      endif
      fun(16) = trkm      
      if (mxx2.gt.0.5d0.and.mxx2.lt.0.7d0) fun(18) = fun(16)
#elif KLOEII
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi

         bmi(1)  =  30d0
         bma(1)  = 130d0
         bmi(2)  =  bmi(1)          
         bma(2)  =  bma(1)
         bmi(3)  = 215d0
         bma(3)  = 515d0
         bmi(4)  = bmi(3)
         bma(4)  = bma(3)
         bmi(5)  = 0.d0
         bma(5)  = 375.d0
         bma(5)  = 500.d0 ! aaargh
         bmi(6)  = bmi(5)
         bma(6)  = bma(5)
         bmi(7)  = 200.d0
         bma(7)  = 500.d0
         bmi(8)  = bmi(7)
         bma(8)  = bma(7)
         bmi(9)  = 591d0
         bma(9)  = 975d0
         bmi(10) = 0.d0
         bma(10) = 500.d0
         bmi(11) = 0.d0
         bma(11) = 60.d0
         bmi(12) = 0.d0
         bma(12) = 25d0
         bmi(13) = 30d0
         bma(13) = 130d0
         bmi(14) = 165d0
         bma(14) = 180d0
         bmi(15) = 0d0
         bma(15) = 15d0
         bmi(16) = 592d0
         bmi(16) = 591d0 ! aaargh
         bma(16) = 975d0
         bmi(17) = 0.d0
         bma(17) = 500d0         
! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo
         ifirst = 0
      endif
******************
      fun(1) = thp*180d0/pi
      fun(2) = thm*180d0/pi
      fun(3) = pp*1d3
      fun(4) = pm*1d3
      fun(5) = pzp*1d3
      fun(6) = pzm*1d3
      fun(7) = pperpp*1d3
      fun(8) = pperpm*1d3
      fun(9) = mxx*1d3

      fun(11) = abs(thp+thm - pi)*180d0/pi

      phip = atan2(p3(2),p3(1))
      phim = atan2(p4(2),p4(1))      
      dphi = abs(abs(phip-phim) - pi)
      fun(12) = dphi*180d0/pi

      thav = (thm-thp+pi)*0.5d0
      fun(13) = thav * 180d0/pi

      fun(14) = thgtilda*180.d0/pi
      fun(15) = fun(14)
      
      fun(16:17) = bmi(16:17) -10.d0
      if (fs.eq.'ee') then
         trkm    = ame*1d3
         if (qph(2,0).gt.0.d0) then
            call getmasstrack(2,ecms,p3,p4,ame*ame,trkm2)
            trkm = sqrt(trkm2)*1d3
         endif
         if (trkm.gt.130d0.and.trkm.lt.220d0) then
            fun(16) = fun(9)
         endif
      elseif(fs.eq.'mm') then
         trkm    = ammu*1d3
         if (qph(2,0).gt.0.d0) then
            call getmasstrack(2,ecms,p3,p4,ammu*ammu,trkm2)
            trkm = sqrt(trkm2)*1d3
         endif         
         if (trkm.gt.80d0.and.trkm.lt.115d0) then
            fun(16) = fun(9)
         endif
      elseif(fs.eq.'pp') then
         trkm    = ampi*1d3
         if (qph(2,0).gt.0.d0) then
            call getmasstrack(2,ecms,p3,p4,ampi*ampi,trkm2)
            trkm = sqrt(trkm2)*1d3
         endif
         if (trkm.gt.130d0.and.trkm.lt.220d0) then
            fun(16) = fun(9)
         endif            
      endif
      if (mxx2.gt.0.5d0.and.mxx2.lt.0.7d0) fun(17) = trkm
      fun(10) = trkm
#elif BES3
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi

         bmi(1)  =  0.d0
         bma(1)  =  1.8d0
         bmi(2)  =  bmi(1)          
         bma(2)  =  bma(1)
         bmi(3)  =    0d0
         bma(3)  = 2000d0
         bmi(4)  = bmi(3)
         bma(4)  = bma(3)
         bmi(5)  = bmi(3)
         bma(5)  = bma(3)
         bmi(6)  = bmi(3)
         bma(6)  = bma(3)
         bmi(7)  =    0.d0
         bma(7)  = 2500.d0
         bmi(8)  = 0.d0
         bma(8)  = 1.d0
         bmi(9)  = 0.d0
         bma(9)  = 2.4d0

         bmi(10) = 0.d0
         bma(10) = 4000.d0
         if (fs.eq.'mm') then
            bmi(10) = 80.d0
            bma(10) = 4080.d0
         elseif (fs.eq.'pp'.or.fs.eq.'pi') then
            bmi(10) = 80.d0
            bma(10) = 4080.d0
         endif

         bmi(11) = 0.d0
         bma(11) = 4000d0
         bmi(12) = 0d0
         bma(12) = 1.d0
         bmi(13) = 0.3d0
         bma(13) = 2.8d0

         bmi(14) = 0.d0
         bma(14) = 2400.d0
         
! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo
         ifirst = 0
      endif
******************
      fun(1) = abs(cp)
      fun(2) = abs(cm)
      fun(3) = pperpp*1d3
      fun(4) = pperpm*1d3
      fun(5) = pp*1d3
      fun(6) = pm*1d3

      call mostenergeticphoton(qinside,pluto)      
      fun(7) = pluto(0)*1d3

c      phip = atan2(p3(2),p3(1))
c      phim = atan2(p4(2),p4(1))      
c      dphi = abs(abs(phip-phim) - pi)  ! calculated in cuts, in common
      fun(8) = dphi/pi
      
      fun(9)  = z ! in common
      fun(10) = 1d3*sqrt(dot(p3+p4,p3+p4))
      fun(11) = 1d3*sqrt(dot(p3+p4+pluto,p3+p4+pluto))
      fun(12) = abs(pluto(3)/pluto(0))
      fun(13) = (thm-thp+pi)*0.5d0
      fun(14) = fun(10)
#elif B
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi

         bmi(1)  =  1000d0
         bma(1)  =  5000d0
         bmi(2)  =  0.25d0
         bma(2)  =  2.75d0
         bmi(3)  =  1000d0
         bma(3)  =  5000d0
         bmi(4)  =  0.25d0
         bma(4)  =  2.75d0
         bmi(5)  = 3000d0
         bma(5)  = 5000d0
         bmi(6)  = 0.2d0
         bma(6)  = 2.7d0
         bmi(7)  = 0.52d0
         bma(7)  = 2.63d0
         bmi(8)  = 0.d0
         bma(8)  = 1.d0
         bmi(9)  = 0.d0
         bma(9)  = 1.d0
         bmi(10) = 0.d0
         bma(10) = 8000.d0
         bmi(11) = 0.d0
         bma(11) = 2400.d0
         bmi(12) = 8000.d0
         bma(12) = 10000.d0
         bmi(13) = 0.d0
         bma(13) = 0.3d0
! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo
         ifirst = 0
      endif
******************
      fun(1) = pp*1d3
      fun(2) = thp !*180.d0/pi
      fun(3) = pm*1d3
      fun(4) = thm !* 180.d0/pi
      fun(5) = qhardest(0)*1d3
      fun(6) = acos(qhardest(3)/qhardest(0))
      fun(7) = (thm-thp+pi)*0.5d0

c done in cuts      phip    = atan2(p3(2),p3(1))
c      phim    = atan2(p4(2),p4(1))      
c      dphi    = abs(abs(phip-phim) - pi)            
      fun(8 ) = dphi/pi
      fun(9 ) = z/pi
      fun(10) = 1d3*sqrt(dot(p3+p4,p3+p4))      
      fun(11) = fun(10)
      fun(12) = mxxg*1d3
      fun(13) = thghgf
#else
! default
      if (ifirst.eq.1) then
         pi = 4.d0 * atan(1.d0)
         radtodeg = 180.d0/pi
         bmi(1) = thmine * radtodeg - 0.1d0
         bma(1) = thmaxe * radtodeg + 0.1d0

         bmi(1) = (thmine) * radtodeg - 1.d0
         bma(1) = (thmaxe) * radtodeg + 1.d0
         
         if (fs.eq.'gg') then ! for any fs!
            bmi(1) = 0.d0 * radtodeg - 0.1d0
            bma(1) = pi * radtodeg + 0.1d0
         endif
         bmi(2) = emin  - 0.01d0
         bma(2) = ebeam + 0.01d0!max(esig1,esig2)*2.d0 + 0.000001d0
         
         bmi(3) =  -1d0 !-.01d0
         bma(3) =  11d0!zmax*radtodeg + .01d0 
         bmi(4) =  -0.01d0!0.0d0*ebeam
         bma(4) =  0.7d0*ebeam!bma(2)
         bmi(5) =  -1d0!0.d0 * radtodeg - 0.01d0
         bma(5) =   181d0! pi * radtodeg + 0.01d0

         bmi(6) =  bmi(1)
         bma(6) =  bma(1)
         
         bmi(7) =  bmi(2)
         bma(7) =  bma(2)

         
         bmi(8) = 2.d0*0.511d-3  - 0.00001d0

         bmin = 1.d0
         cccc = cos(pi-zmax)
         bmi(8) = sqrt(2.d0*emin*emin*(1.d0-bmin*bmin*cccc))
         uuu = ame
         if (fs.eq.'mm') uuu = ammu
         if (fs.eq.'pp') uuu = ampi
         if (fs.eq.'pi') uuu = ampi
         bmi(8) = 2.d0*max(emin,uuu)
         bma(8) = ecms

         bmi(9)  =  bmi(4)
         bma(9)  =  bma(4)
         bmi(10) =  bmi(5)
         bma(10) =  bma(5)

         bmi(11) = bmi(1)
         bma(11) = bma(1)
         bmi(12) = bmi(2)
         bma(12) = bma(2)
         bmi(13) = bmi(6)
         bma(13) = bma(6)
         bmi(14) = bmi(7)
         bma(14) = bma(7)

         bmi(15) = bmi(3)
         bma(15) = bma(3)

         bmi(16) = -ecmsnom**2*0.5d0*(1.d0-cos(thmaxe))
         bma(16) = -ecmsnom**2*0.5d0*(1.d0-cos(thmine))

         bmi(17) = -1d0
         bma(17) =  1d0
         bmi(18) = -1d0
         bma(18) =  1d0
         
! filling distr. components when entering the first time
         do k = 1,ndistr
            d = (bma(k) - bmi(k))/nbin
            y = bmi(k)
            do i = 1,nbin
               y = y + d
               distr(k,i*3-2)   = y - d
               s_distr(k,i*3-2) = y - d
            enddo
         enddo
         ifirst = 0
      endif
******************
      if (fs.ne.'gg') then
         px = p3
         py = p4
         call duefotonipiuenergetici(qph,qa,qb)
      else
         px = p1inside
         py = p2inside
         qa = p1rad
         qb = p2rad
      endif
*** electron angle
      pe      = sqrt(px(1)**2+px(2)**2+px(3)**2)
      fun(17) = px(3)/pe
      fun(1)  = radtodeg*acos(fun(17))
**** electron energy
      fun(2) = px(0)
*** positron angle                                     
      pp      = sqrt(py(1)**2+py(2)**2+py(3)**2)
      fun(18) = py(3)/pp
      fun(6)  = radtodeg*acos(fun(18))
**** positron energy
      fun(7) = py(0)
         
**** acollinearity
      fun(3)  = abs(180.d0 - fun(1)- fun(6))
      fun(15) = abs(180.d0 - radtodeg*acos(tridot(px,py)/pe/pp))

      fun(16) = dot(pin1-px,pin1-px)
***   m.e. photon energy
      fun(4)  = bma(4) + 1.d0
      fun(5)  = bma(5) + 1.d0
      fun(9)  = bma(9) + 1.d0
      fun(10) = bma(10) + 1.d0

      soglia = 0.02d0 * ebeam
c      soglia = eps * ebeam
      if (qa(0).ge.soglia) then 
         fun(4) = qa(0)
         fun(5) = radtodeg*acos(qa(3)/qa(0))
      endif
      if (qb(0).ge.soglia) then 
         fun(9)  = qb(0)
         fun(10) = radtodeg*acos(qb(3)/qb(0))
      endif

      fun(8)  = sqrt(abs(dot(px+py,px+py)))

*** randomized p1,p2      
      call getrnd(csi,1)
      if (csi(1).lt.0.5d0) then
         fun(11) = fun(1)
         fun(12) = fun(2)
         fun(13) = fun(6)
         fun(14) = fun(7)
      else
         fun(13) = fun(1)
         fun(14) = fun(2)
         fun(11) = fun(6)
         fun(12) = fun(7)
      endif
#endif
***********************************************************************
      DO k = 1,ndistr
         bmax = bma(k)
         d = (bma(k) - bmi(k))/nbin
         x = bmi(k)
** NEW ***
         if (fun(k).le.bma(k).and.fun(k).ge.bmi(k)) then
            i  = int((fun(k) - bmi(k))/d)
c            if ((i+1).gt.nbin) print*,'WARNING!!',
c     .                            k,i+1,d,fun(k),bmi(k),bma(k)
            i  = min(i+1,nbin)*3

            distr(k,i-1) = distr(k,i-1) + sd
            distr(k,i)   = distr(k,i)   + sd**2                       
            su  = distr(k,i-1)
            su2 = distr(k,i)               
            s_distr(k,i-1) = su / nc
            argument = abs((su2/nc-s_distr(k,i-1)**2)/nc)
            s_distr(k,i) = sqrt(argument)
         endif         
** NEW ***
      ENDDO          
      return
      end
***      
      subroutine writedistributions
#include "shared.F"
      common/ifirstwritedistributions/ifirst
      
      parameter (ibin = 3 * nbin)      
      parameter (iperj = ndistr*ibin)
      
      common/filenames/distname
      common/distr/distr,s_distr,bmi,bma
      
      dimension s_distr(ndistr,ibin),bmi(ndistr),bma(ndistr)
      dimension distr(ndistr,ibin)
      
      character*100 distname(ndistr)

      data ifirst /1/

      do k = 1,ndistr
         open(12,file=distname(k),status='unknown')
         do i=1, nbin
            sezd  = s_distr(k,i*3-1)/(bma(k)-bmi(k))*nbin
            esezd = s_distr(k,i*3)  /(bma(k)-bmi(k))*nbin
            write(12,*)s_distr(k,i*3-2),sezd,esezd
         enddo
         close(12)
      enddo

      ifirst = 0
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mostenergeticphoton(q,q1)
      implicit double precision (a-h,o-z)
      dimension q(40,0:3),q1(0:3)
*  LEADING ENERGETIC PHOTON IS EXTRACTED
      q1(0) = 0.d0
      q1(1) = 0.d0
      q1(2) = 0.d0
      q1(3) = 0.d0

ccc      if (q(1,0).lt.1.d-11) return
C for backward compatibility with old ps-based BABAYAGA
      if (q(1,0).lt.1.d-12.and.q(11,0).lt.1.d-12.and.
     .     q(21,0).lt.1.d-12.and.q(31,0).lt.1.d-12) return
      
      ENPHOT = Q(1,0)
      ENHARD = ENPHOT
      J=1
*     
      DO I = 1,40
         if (q(i,0).gt.0.d0) then
            ENPHOT = Q(I,0)
            IF (ENHARD.GE.ENPHOT) THEN
               ENHARD = ENHARD
            ELSE
               J = I
               ENHARD = ENPHOT
            ENDIF
         endif
      ENDDO 
* 
      ENHARDL = ENHARD
      JL = J
*    
      if (j.gt.0) then
         DO I = 0,3
            Q1(I) = Q(J,I)
         ENDDO
      endif
      return
      end
*****************************************
      function get_dphi(p1,p2)
      implicit double precision (a-h,o-z)
      double precision p1(0:3),p2(0:3)
! from ALPGEN
      dphi = (p1(1)*p2(1)+p1(2)*p2(2))/
     .     sqrt(p1(1)**2+p1(2)**2)/sqrt(p2(1)**2+p2(2)**2)
      if (abs(dphi).gt.1.d0) then
         dphi = dphi/abs(dphi)
      endif
      get_dphi = acos(dphi)
      return
      end
***********************

      function getphi_here(ppp)
      implicit double precision (a-h,o-z)
      dimension ppp(0:3)
      parameter (pi = 3.141592653589793238462643383279502884197169399d0)
      pm = sqrt(ppp(1)**2+ppp(2)**2+ppp(3)**2)
      c  = ppp(3)/pm
      s  = sqrt(1.d0-c**2)
      if (s.eq.0.d0) then
         getphi_here = 0.d0
         return
      else   
         arg = ppp(1)/pm/s
*  avoiding numerical problems......
         if (abs(arg).ge.1.d0) then
            iarg = arg
            arg  = iarg
         endif   
         if (ppp(2).ge.0.d0) getphi_here = acos(arg)
         if (ppp(2).lt.0.d0) getphi_here = 2.d0*pi-acos(arg)
      endif             
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine duefotonipiuenergetici(q,q1,q2)
      implicit real*8 (a-h,o-z)
      dimension q(40,0:3),q1(0:3),q2(0:3),ql(40,0:3)
* 
*  LEADING ENERGETIC PHOTON IS EXTRACTED
*
      q1 = 0.d0
      q2 = 0.d0

      n = 0
      do k = 1,40
         if (q(k,0).gt.0.d0) then
            n = n+1
            ql(n,:) = q(k,:)
         endif
      enddo
      if (n.eq.0) return
      
      ENPHOT = QL(1,0)
      ENHARD = ENPHOT
      J=1
*     
      DO I = 1,N
         ENPHOT = QL(I,0)
         IF (enphot.ge.enhard) THEN
            J = I
            ENHARD = ENPHOT
         ENDIF
      ENDDO 
*    
      Q1 = QL(J,:)
      ENHARDL = ENHARD
      JL = J
*     
      j1 = j
*     
*     NEXT-TO-LEADING ENERGETIC PHOTON IS EXTRACTED
*     
      ENHARD = 0.D0
      J = 0
*     
      DO I = 1,N
         IF (I.NE.JL) THEN
            ENPHOT = QL(I,0)
            IF (ENPHOT.GE.ENHARD) THEN
               J = I
               ENHARD = ENPHOT
            ENDIF
         ENDIF
      ENDDO 
c
      if (j.gt.0) then
         Q2 = QL(J,:)
      endif
      j2=j

      return
      end      
************************************************
      subroutine duefotonipiuenergetici_bck(q,q1,q2)
      implicit real*8 (a-h,o-z)
      dimension q(40,0:3),q1(0:3),q2(0:3)
* 
*  LEADING ENERGETIC PHOTON IS EXTRACTED
*
        do i = 0,3
           q1(i) = 0.d0
           q2(i) = 0.d0
        enddo

        sum = 0.d0
	do k=1,40
          sum = sum + q(k,0)
        enddo

        if (sum.lt.1.d-11) return

        ENPHOT = Q(1,0)
        ENHARD = ENPHOT
        J=1
* 
        DO I = 1,40
           ENPHOT = Q(I,0)
           IF (ENHARD.GE.ENPHOT) THEN
              ENHARD = ENHARD
           ELSE
              J = I
              ENHARD = ENPHOT
           ENDIF
        ENDDO 
*   
        ENHARDL = ENHARD
        JL = J
*
        if (j.gt.0) then
           DO I = 0,3
              Q1(I) = Q(J,I)
           ENDDO
        endif
        j1 = j
*
*  NEXT-TO-LEADING ENERGETIC PHOTON IS EXTRACTED
*
      ENHARD = 0.D0
      J = 0
* 
      DO I = 1,40
         IF (I.NE.JL) THEN
*
            ENPHOT = Q(I,0)
            IF (ENHARD.GE.ENPHOT) THEN
               ENHARD = ENHARD
            ELSE
               J = I
               ENHARD = ENPHOT
            ENDIF
*     
         ENDIF
      ENDDO 
*   
      ENHARD = ENHARD
*
      if (j.gt.0) then
         DO I = 0,3
            Q2(I) = Q(J,I)
         ENDDO
      endif
      j2=j

      return
      end
***************************************************
      subroutine getmasstrack(ng,sqs,p1,p2,m2,m2trk)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3)
!     should be used only if the number of photons is >= 2!!
      
      m2trk = m2
      if (ng.lt.2) return

      p1m = sqrt(tridot(p1,p1))
      p2m = sqrt(tridot(p2,p2))

! by Stefan Mueller, from phokhara10
      squared = 0.5d0*((sqs-sqrt(tridot(p1+p2,p1+p2)))**2-p1m**2-p2m**2)
      m2trk = ( squared**2 - p1m**2*p2m**2 ) /
     .     ( 2.d0*squared + p1m**2 + p2m**2 )
      m2trk = abs(m2trk)
c      print*,m2trk,fforgetmasstrack(sqs,m2trk,p1,p2)
***   on gitlab
c      a = (sqs - sqrt(tridot(p1+p2,p1+p2)))**2
c      uu = 0.25d0*a - 0.5d0*(p1m**2+p2m**2)+0.25d0/a*(p1m**2-p2m**2)**2
c      print*,(uu-m2trk)/m2trk
      
      return            
c$$$      a = -min(p1m*p1m,p2m*p2m)
c$$$      b = sqs*sqs
c$$$      
c$$$      print*,sqrt(tridot(p1,p1)),sqrt(tridot(p2,p2))
c$$$      
c$$$      tol = 1d-10      
c$$$      do while (abs(a-b).gt.tol)
c$$$         c = (a+b)*0.5d0
c$$$         fc = fforgetmasstrack(sqs,c,p1,p2)
c$$$         fa = fforgetmasstrack(sqs,a,p1,p2)
c$$$         print*,fa,fc
c$$$         if (fa/fc.gt.0.d0) then
c$$$            a = c
c$$$         else
c$$$            b = c
c$$$         endif
c$$$      enddo
c$$$      m2 = c
c$$$      return
      end
      function fforgetmasstrack(sqs,m2,p1,p2)
      implicit double precision (a-h,m,o-z)
      dimension p1(0:3),p2(0:3)
      p1m = sqrt(tridot(p1,p1))
      p2m = sqrt(tridot(p2,p2))
      E1 = sqrt(p1m*p1m + m2)
      E2 = sqrt(p2m*p2m + m2)
c      fforgetmasstrack = sqs*sqs + 2.d0*m2 +
c     .     2d0*E1*E2 - 2d0*tridot(p1,p2)
c     .    -2.d0*sqs*(E1+E2)
c      f1 = fforgetmasstrack      
      fforgetmasstrack = (sqs - E1 - E2)**2 - tridot(p1+p2,p1+p2)     
      return
      end
***************************************************      
