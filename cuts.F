#if CMD
      subroutine cuts(p1,p2,qph,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/massainv/amassainvmin,amassainvmax
      common/ifirstcuts/ifirst
      data ifirst /0/
      common/idarkon/idarkon
      common/ecms/ecms,ecmsnom
      dimension ptvector(0:3)
#include "strong2020common.F"

      if (ifirst.eq.0) then         
         print*,'CMD SCENARIO'

         thavgmin = 1.d0
         thavgmax = pi - 1.d0

         dfmax  = 0.15d0
         csimax = 0.25d0

         ppmmin = 0.45d0*ecms*0.5d0
         
         ifirst = 1
      endif
      
      icut = 1 ! event rejected

      pm = sqrt(tridot(p1,p1))
      pp = sqrt(tridot(p2,p2))

      if (pp.lt.ppmmin) return
      if (pm.lt.ppmmin) return

      thm = acos(p1(3)/pm)
      thp = acos(p2(3)/pp)
      
      z = abs(thp+thm-pi)
      if (abs(z).gt.csimax) return

      thavg = (thm - thp + pi)*0.5d0
      if (thavg.lt.thavgmin) return
      if (thavg.gt.thavgmax) return

c      phip = getphi_here(p1) !atan2(p1(2),p1(1))
c      phim = getphi_here(p2) !atan2(p2(2),p2(1))
      phip = atan2(p1(2),p1(1))
      phim = atan2(p2(2),p2(1))      
      dphi = abs(abs(phip-phim) - pi)
c      if (dphi.gt.0.d0) then
c      print*,abs(abs(atan2(p1(2),p1(1))-atan2(p2(2),p2(1))) - pi)/dphi
c      endif
      if (dphi.gt.dfmax) return      
*
************************************
*
      icut = 0 ! event accepted
      return
      end
#elif KLOEI
      subroutine cuts(p1,p2,qph,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/massainv/amassainvmin,amassainvmax
      common/ifirstcuts/ifirst
      data ifirst /0/
      common/idarkon/idarkon
      dimension ptvector(0:3)
#include "strong2020common.F"
      
      if (ifirst.eq.0) then         
         print*,'KLOEI (tagged) SCENARIO'

         thgtlow  =  50.d0/180.d0*pi ! used for tagged
         thgthigh = 130.d0/180.d0*pi

         thgulow  =  15.d0/180.d0*pi ! used for untagged 
         thguhigh = 165.d0/180.d0*pi

         mxx2min = 0.1d0           ! these for tagged
         mxx2max = 0.85d0
!         mxx2min = 0.35d0          ! these for untagged
!         mxx2max = 0.95d0

         pzmin    = 0.09d0
         pperpmin = 0.160d0

         egammamin = 0.020d0
         
         
         ifirst = 1
      endif
      
      icut = 1 ! event rejected

      pm = sqrt(tridot(p1,p1))
      pp = sqrt(tridot(p2,p2))
      
      c1 = p1(3)/pm
      c2 = p2(3)/pp
      thm = acos(c1)
      thp = acos(c2)
      if (thm.lt.thmin.or.thm.gt.thmax) return
      if (thp.lt.thmin.or.thp.gt.thmax) return

      pzm = abs(p1(3))
      pzp = abs(p2(3))
      pperpm = sqrt(p1(1)**2+p1(2)**2)
      pperpp = sqrt(p2(1)**2+p2(2)**2)

      if (pzm.gt.pzmin.or.pperpm.gt.pperpmin) then
         continue
      else
         return
      endif
      if (pzp.gt.pzmin.or.pperpp.gt.pperpmin) then
         continue
      else
         return
      endif

      mxx2 = dot(p1+p2,p1+p2)
      mxx  = sqrt(mxx2)
      if (mxx2.lt.mxx2min.or.mxx2.gt.mxx2max) return

      thereisoneinside = .false.
      qinside           = 0.d0
      if (qph(1,0).gt.0.d0) then
         j = 0
         do k = 1,40
            if (qph(k,0).gt.0.d0) then               
               thg = acos(qph(k,3)/qph(k,0))
               if (thg.gt.thgtlow.and.thg.lt.thgthigh
     .              .and.qph(k,0).gt.egammamin) then
                  j = j + 1
                  thereisoneinside = .true.
                  qinside(j,:) = qph(k,:)
               endif
            endif
         enddo
      endif
      if (.not.thereisoneinside) return
*
************************************
*
      icut = 0 ! event accepted
      return
      end
#elif KLOEII
      subroutine cuts(p1,p2,qph,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/massainv/amassainvmin,amassainvmax
      common/ifirstcuts/ifirst
      data ifirst /0/
      common/idarkon/idarkon
      dimension ptvector(0:3)
#include "strong2020common.F"
      
      if (ifirst.eq.0) then         
         print*,'KLOEII (untagged) SCENARIO'

         thgtlow  =  50.d0/180.d0*pi ! used for tagged
         thgthigh = 130.d0/180.d0*pi

         thgulow  =  15.d0/180.d0*pi ! used for untagged 
         thguhigh = 165.d0/180.d0*pi

!         mxx2min = 0.1d0           ! these for tagged
!         mxx2max = 0.85d0
         mxx2min = 0.35d0          ! these for untagged
         mxx2max = 0.95d0

         pzmin    = 0.09d0
         pperpmin = 0.160d0

         egammamin = 0.020d0
         ifirst = 1
      endif
      
      icut = 1 ! event rejected

      pm = sqrt(tridot(p1,p1))
      pp = sqrt(tridot(p2,p2))
      
      c1 = p1(3)/pm
      c2 = p2(3)/pp
      thm = acos(c1)
      thp = acos(c2)
      if (thm.lt.thmin.or.thm.gt.thmax) return
      if (thp.lt.thmin.or.thp.gt.thmax) return

      pzm = abs(p1(3))
      pzp = abs(p2(3))
      pperpm = sqrt(p1(1)**2+p1(2)**2)
      pperpp = sqrt(p2(1)**2+p2(2)**2)
      
      if ((pzm.gt.pzmin.or.pperpm.gt.pperpmin).and.
     .     (pzp.gt.pzmin.or.pperpp.gt.pperpmin)  ) then
         continue
      else
         return
      endif
c      if (pzp.gt.pzmin.or.pperpp.gt.pperpmin) then
c         continue
c      else
c         return
c      endif
      
      mxx2 = dot(p1+p2,p1+p2)
      mxx  = sqrt(mxx2)
      if (mxx2.lt.mxx2min.or.mxx2.gt.mxx2max) return

      qfake = -p1-p2
      thgtilda = acos(qfake(3)/sqrt(tridot(qfake,qfake)))
      if (thgtilda.ge.thgulow.and.thgtilda.le.thguhigh) return      
*
************************************
*
      icut = 0 ! event accepted
      return
      end
#elif BES3
      subroutine cuts(p1,p2,qph,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/massainv/amassainvmin,amassainvmax
      common/ifirstcuts/ifirst
      data ifirst /0/
      common/idarkon/idarkon
      dimension ptvector(0:3)
      logical condin
#include "strong2020common.F"
      
      if (ifirst.eq.0) then         
         print*,'BES3 SCENARIO'

         cthpmmin = -0.93d0
         cthpmmax =  0.93d0

         egammamin1 = 0.025d0
         egammamin2 = 0.050d0
         
         pperpmin = 0.300d0

         acthgmax1 = 0.80d0
         acthgmax2 = 0.92d0
         acthgmin2 = 0.86d0
                           
         ifirst = 1
      endif
      
      icut = 1 ! event rejected

      pm  = sqrt(tridot(p1,p1))
      pp  = sqrt(tridot(p2,p2))
      
      cm  = p1(3)/pm
      cp  = p2(3)/pp
      thm = acos(cm)
      thp = acos(cp)
      if (thm.lt.thmin.or.thm.gt.thmax) return
      if (thp.lt.thmin.or.thp.gt.thmax) return

      pperpm = sqrt(p1(1)**2+p1(2)**2)
      pperpp = sqrt(p2(1)**2+p2(2)**2)
      if (pperpm.lt.pperpmin) return
      if (pperpp.lt.pperpmin) return
      
      thereisoneinside = .false.
      qinside           = 0.d0
      if (qph(1,0).gt.0.d0) then
         j = 0
         do k = 1,40
            if (qph(k,0).gt.0.d0) then               
               eg   = qph(k,0)
               acth = abs(qph(k,3)/eg)
               condin = acth.lt.acthgmax1.and.eg.gt.egammamin1
               condin = condin.or.(
     .              (acth.gt.acthgmin2.and.acth.lt.acthgmax2).and.
     .              eg.gt.egammamin2
     .              )
*** new cut
               condin = condin.and.eg.gt.0.4d0
***
               if (condin) then
                  j = j + 1
                  thereisoneinside = .true.
                  qinside(j,:) = qph(k,:)
               endif
            endif
         enddo
      endif
      if (.not.thereisoneinside) return

      z    = abs(thp+thm-pi)
      phip = atan2(p1(2),p1(1))
      phim = atan2(p2(2),p2(1))      
      dphi = abs(abs(phip-phim) - pi)      
************************************
      icut = 0 ! event accepted
      return
      end
#elif B
      subroutine cuts(p1,p2,qph,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/massainv/amassainvmin,amassainvmax
      common/ifirstcuts/ifirst
      data ifirst /0/
      common/idarkon/idarkon
      dimension ptvector(0:3)
      logical condin
      double precision lambda
      external lambda
      common/printing/iprint
      character*2 fs
      common/finalstate/fs
#include "strong2020common.F"
      iprint = 0
      
      if (ifirst.eq.0) then         
         print*,'B SCENARIO'
         thmaxcutB = 2.75d0
         thmincutB = 0.65d0
         mxxgmin   = 8d0
         thghgfmax = 0.3d0
         egammamin = 3d0
         ifirst = 1
      endif
      
      icut = 1 ! event rejected

      pm = sqrt(tridot(p1,p1))
      pp = sqrt(tridot(p2,p2))

c$$$      mxx = sqrt(dot(p1+p2,p1+p2))
c$$$      ame  = 0.510998950d-3
c$$$      s   = dot(pin1+pin2,pin1+pin2)
c$$$      xxx = lambda(s,0.d0,(ame+ame)**2)/s
c$$$      xxx = 0.5d0*sqrt(xxx)
c$$$c      print*,xxx,s,ame
c$$$c      stop
c$$$c      mxx = xxx
c$$$      if (1d3*mxx.lt.16d0) then
c$$$         iprint = 1
c$$$         print*,' '
c$$$         print*,mxx*1d3
c$$$         print*,p1
c$$$         print*,p2
c$$$c         print*,sqrt(dot(p1+p2,p1+p2))/ame*0.5d0
c$$$         print*,qph(1,:)
c$$$      endif
         
      cm  = p1(3)/pm
      cp  = p2(3)/pp
      thm = acos(cm)
      thp = acos(cp)
      if (thm.lt.thmincutB.or.thm.gt.thmaxcutB) return
c      if (iprint.eq.1) print*,'after 1'
      if (pm.lt.1d0) return
c      if (iprint.eq.1) print*,'after 2'
      if (thp.lt.thmincutB.or.thp.gt.thmaxcutB) return
c      if (iprint.eq.1) print*,'after 3'
      if (pp.lt.1d0) return
c      if (iprint.eq.1) print*,'after 4'

      thereisoneinside = .false.
      qinside          = 0.d0
      if (qph(1,0).gt.0.d0) then
         j = 0
         do k = 1,40
            if (qph(k,0).gt.0.d0) then               
               eg  = qph(k,0)
               thg = acos(qph(k,3)/eg)
               condin = thg.gt.thgmin.and.thg.lt.thgmax.and.
     .              eg.gt.egammamin
               if (condin) then
                  j = j + 1
                  thereisoneinside = .true.
                  qinside(j,:) = qph(k,:)
               endif
            endif
         enddo
      endif
      if (.not.thereisoneinside) return
c      if (iprint.eq.1) print*,'after 5'

      mxxg = 0.d0
      call mostenergeticphoton(qinside,qhardest)
      mxxg = sqrt(dot(p1+p2+qhardest,p1+p2+qhardest))
      if (mxxg.lt.mxxgmin) return
c      if (iprint.eq.1) print*,'after 6'
            
      qfake  = -p1-p2
      qmf    = sqrt(tridot(qfake,qfake))
      ctest  = tridot(qfake,qhardest)/qmf/qhardest(0)
      thghgf = 0.d0
      if (abs(ctest).lt.1.d0) then
         thghgf = acos(ctest)
      endif

      if (thghgf.gt.thghgfmax) return
c      if (iprint.eq.1) print*,'after 7'

***   new cut for Bhabha
      if (fs.eq.'ee') then
         mxx = sqrt(dot(p1+p2,p1+p2))
         if (mxx.lt.0.3d0) return
      endif      
***      

      z    = abs(thp+thm-pi)
      phip = atan2(p1(2),p1(1))
      phim = atan2(p2(2),p2(1))      
      dphi = abs(abs(phip-phim) - pi)      

************************************
      icut = 0                  ! event accepted
c      if (iprint.eq.1) print*,icut
      return
      end      
#else
! default cuts
      subroutine cuts(p1,p2,qph,icut)
! p1(0...3): final state electron four-momentum
! p2(0...3): final state positron four-momentum
! qph(0...40,0...3): four-momenta of 40 emitted photons
! icut: icut = 1 event rejected, icut = 0 event accepted
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),ptmp(0:3)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      common/angularranges/thmine,thmaxe,thminp,thmaxp
      common/momentainitial/pin1(0:3),pin2(0:3)
      common/massainv/amassainvmin,amassainvmax
      common/ifirstcuts/ifirst
      data ifirst /0/
      common/idarkon/idarkon
      dimension ptvector(0:3)
      character*2 fs
      common/finalstate/fs
      common/ecms/ecms,ecmsnom

      if (ifirst.eq.0) then         
         
         ifirst = 1
      endif
      
      icut = 1 ! event rejected

      if(p1(0).lt.emin.or.p2(0).lt.emin) return

      z = acollinearityrad(p1,p2)
      if (z.gt.zmax) return

      c1 = p1(3)/sqrt(tridot(p1,p1))
      c2 = p2(3)/sqrt(tridot(p2,p2))
      th1 = acos(c1)
      th2 = acos(c2)

      thmine = thmin
      thmaxe = thmax
c      thminp = thmin
c      thmaxp = thmax      
      thminp = pi - thmax
      thmaxp = pi - thmin

      if (th1.lt.thmine.or.th1.gt.thmaxe) return
      if (th2.lt.thminp.or.th2.gt.thmaxp) return

c      pcut = ecms*0.5d0*0.45d0
c      p1m = sqrt(tridot(p1,p1))
c      if (p1m.lt.pcut) return
c      p2m = sqrt(tridot(p2,p2))
c      if (p2m.lt.pcut) return
************************************
      icut = 0 ! event accepted
      return
      end
#endif
**************************************************************************
      subroutine cutsgg(ng,p1,p2,qph,icut)
      implicit double precision (a-h,o-z)
      dimension p1(0:3),p2(0:3),qph(40,0:3),q(0:3),q1(0:3),q2(0:3)
      dimension icutvec((ng+2)*(ng+1)/2)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      common/expcuts/thmin,thmax,emin,zmax,egmin,thgmin,thgmax
      character*2 fs
      common/finalstate/fs
      integer ifirst
      common/uuucutsgg/ifirst
      data ifirst /0/
      
      double precision qphtot(40,0:3)
      double precision p1inside(0:3),p2inside(0:3)
      double precision p1rad(0:3),p2rad(0:3)
      integer inside(40)
      common/cutsggcmn/qphtot,p1inside,p2inside,p1rad,p2rad,
     .     nph,npairs,inside,n_inside,icutgg
c$$$      if (ng.eq.0) then
c$$$         call cuts(p1,p2,qph,icut)
c$$$         return
c$$$      endif
c$$$*     **************************
c$$$*     below only for ng > 0      
c$$$      icut   = 1 ! event rejected
c$$$      nph    = ng + 2
c$$$      npairs = (nph)*(nph-1)/2
c$$$      qphtot(1,:) = p1
c$$$      qphtot(2,:) = p2
c$$$      
c$$$      do k = 1,ng
c$$$         qphtot(k+2,:) = qph(k,:)
c$$$      enddo
c$$$      
c$$$      ip = 0
c$$$      icuttot = 0
c$$$      do ka = 1,nph-1
c$$$         do kb = ka+1,nph
c$$$            ip = ip + 1
c$$$            q1 = qphtot(ka,:)
c$$$            q2 = qphtot(kb,:)
c$$$            call cuts(q1,q2,qph,icuttmp)
c$$$            icuttot = icuttot + icuttmp
c$$$c            icutvec(ip) = icuttmp
c$$$
c$$$c            if (icuttmp.eq.0) then
c$$$c               icut = 0
c$$$c               return
c$$$c            endif
c$$$
c$$$        enddo
c$$$      enddo
c$$$c      if (ng.gt.0.and.icuttot.lt.npairs.and.icutvec(1).ne.0) 
c$$$c     .     print*,icutvec
c$$$  if (icuttot.lt.npairs) icut = 0

c      better to implement cuts on single photons, and then maybe on pairs

      if (ifirst.eq.0) then
         ifirst = 1
         inside = 0
         qphtot = 0.d0
      endif
      
      icut = 1 ! even rejected
      icutgg = 1
      
      qphtot = 0.d0
      
      nph    = ng + 2
      npairs = (nph)*(nph-1)/2
      qphtot(1,:) = p1
      qphtot(2,:) = p2
      
      if (ng.ge.1) then
         do k = 1,ng
            qphtot(k+2,:) = qph(k,:)
         enddo
      endif

      n_inside = 0
      inside(1:nph) = 0
      do k = 1,nph
         e  = qphtot(k,0)
         th = acos(qphtot(k,3)/e)
         if (e.ge.emin.and.th.ge.thmin.and.th.le.thmax) then
            inside(k) = 1
            n_inside = n_inside + 1
         endif
      enddo
      
      if (n_inside.ge.2) icut = 0
      icutgg = icut
      
      if (icut.gt.0) return
***   if it's not cut, I fill p[12]inside and p[12]rad
***   p[12]inside are the two most energetic photons inside the cuts
***   p[12]rad are the two most energetic photons of the rest
      jj = 0
      k = 0
      i1 = 0
      do while(jj.lt.2)
         k = k + 1
         if (inside(k).eq.1) then
            jj = jj + 1
            if (jj.eq.1) then
               k1in = k
               p1inside = qphtot(k,:)
            endif
            if (jj.eq.2) then
               k2in = k
               p2inside = qphtot(k,:)
            endif
         endif       
      enddo

      p1rad = 0.d0
      p2rad = 0.d0
      jj = 0
      k  = 0
      do while(jj.lt.2.and.k.le.nph)
         k = k + 1
         if (k.ne.k1in.and.k.ne.k2in) then
            jj = jj + 1
            if (jj.eq.1) p1rad = qphtot(k,:)
            if (jj.eq.2) p2rad = qphtot(k,:)
         endif
      enddo

      z = acollinearityrad(p1inside,p2inside)
      if (z.ge.zmax) then
         icut   = 1
         icutgg = icut
         return
      endif
      
c      if (nph.eq.3) then         
c         print*,' ',nph,k1in,k2in
c      print*,p1inside
c      print*,p2inside
c      print*,p1rad
c      print*,p2rad
c      endif
      return
      end
************************************************************************
