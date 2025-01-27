* COMPLEX VACUUM POLARIZATION: ALPHA(0)->ALPHA(0)*VPOLC(Q2)
      FUNCTION VPOLC(Q2) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 QIN4,ST24,DER4,ERRDER4,DEG4,ERRDEG4
      DIMENSION AMASSES(9)
      common/parameters/ame,ammu,convfac,alpha,pi     
      common/ialpharunning/iarun
      common/teubner/iteubn
      integer inskvp
      common/nskvp/inskvp
      common/resonances/ires,ionDHAD
      common/vpolcomc/ifirst
      common/jpsiparameters/amjpsi,gtjpsi,gejpsi,effcjpsi
      common/jpsiparameters2/amjpsi2s,gtjpsi2s,gejpsi2s,effcjpsi2s
      common/jpsiparameters3/amjpsi3s,gtjpsi3s,gejpsi3s,effcjpsi3s
      common/jpsiparameters4/amjpsi4s,gtjpsi4s,gejpsi4s,effcjpsi4s
      common/jpsiparameters5/amjpsi5s,gtjpsi5s,gejpsi5s,effcjpsi5s

      data ifirst /0/

      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt
      
      common/rescaledhaderr/scal_err_dhad

      double complex vpolcnsk,vpolchlmnt
      external vpolcnsk,vpolchlmnt

      double complex vpolc,deltahad,im,cder

      common/dalphacmn/dalphaQED,der,errder,toppart
      integer ioffhvp,inlovp
      common/offhadrVP_nlovp/ioffhvp,inlovp
      
      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
*     AM=HADRON MASSES, GT=TOTAL WIDTH, GE=PARTIAL E+E- OR Mu+Mu- WIDTH.
* 1=RHO, 2=OMEGA, 3=PHI, 4=J/PSI, 5=PSI(2S), 6=PSI(3770), 7=PSI(4040),
* 8=PSI(4160), 9=PSI(4415), 10=Y, 11=Y(2S), 12=Y(3S), 13=Y(4S), 14=Y(10860),
* 15=Y(11020).
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege
** errors on the previous parameters, from PDG
      
****  TO BE UPDATED (e.g. from PDG)!!
      DATA AM/ ! GeV
     >     0.769D0,    ! 1.  rho
     >     0.78265D0,  ! 2.  omega
     >     1.019455D0, ! 3.  phi
     >     3.096916D0, ! 4.  J/Psi
     >     3.686109D0, ! 5.  Psi(2S)
     >     3.77315D0,  ! 6.  Psi(3770)
     >     4.039D0,    ! 7.  Psi(4040)
     >     4.191D0,    ! 8.  Psi(4160)
     >     4.421D0,    ! 9.  Psi(4415)
     >     9.4603D0,   ! 10. Y
     >     10.02326D0, ! 11. Y(2S)
     >     10.35520,   ! 12. Y(3S)
     >     10.5794D0,  ! 13. Y(4S)
     >     10.876D0,   ! 14. Y(10860)
     >     11.019D0/   ! 15. Y(11020)
      DATA GT/ ! GeV
     >     0.149D0,   ! 1.  rho
     >     8.49D-3,   ! 2.  omega
     >     4.266D-3,  ! 3.  phi
     >     92.9D-6,   ! 4.  J/Psi
     >     299.D-6,   ! 5.  Psi(2S)
     >     27.2D-3,   ! 6.  Psi(3770)
     >     80.D-3,    ! 7.  Psi(4040)
     >     70.D-3,    ! 8.  Psi(4160)
     >     62.D-3,    ! 9.  Psi(4415)
     >     54.02D-6,  ! 10. Y
     >     31.98D-6,  ! 11. Y(2S)
     >     20.32D-6,  ! 12. Y(3S)
     >     20.5D-3,   ! 13. Y(4S)
     >     55.D-3,    ! 14. Y(10860)
     >     79.D-3/    ! 15. Y(11020)
      DATA GE/ ! GeV
     >     7.04D-6,   ! 1.  rho
     >     0.60D-6,   ! 2.  omega
     >     1.27D-6,   ! 3.  phi
     >     5.55D-6,   ! 4.  J/Psi
     >     2.36D-6,   ! 5.  Psi(2S)
     >     0.262D-6,  ! 6.  Psi(3770)
     >     0.86D-6,   ! 7.  Psi(4040)
     >     0.48D-6,   ! 8.  Psi(4160)
     >     0.58D-6,   ! 9.  Psi(4415)
     >     1.34D-6,   ! 10. Y
     >     0.612D-6,  ! 11. Y(2S)
     >     0.443D-6,  ! 12. Y(3S)
     >     0.272D-6,  ! 13. Y(4S)
     >     0.31D-6,   ! 14. Y(10860)
     >     0.130D-6/  ! 15. Y(11020)
*** ERRORS ON THE PREVIOUS PARAMETERS, from PDG
      DATA EAM/                  ! GeV
     >     0.0009D0,   ! 1.  rho
     >     0.00012D0,  ! 2.  omega
     >     0.000019D0, ! 3.  phi
     >     0.000011D0, ! 4.  J/Psi
     >     0.000013D0, ! 5.  Psi(2S)
     >     0.00033D0,  ! 6.  Psi(3770)
     >     0.001D0,    ! 7.  Psi(4040)
     >     0.005D0,    ! 8.  Psi(4160)
     >     0.004D0,    ! 9.  Psi(4415)
     >     0.00026D0,  ! 10. Y
     >     0.00031D0,  ! 11. Y(2S)
     >     0.0005D0,   ! 12. Y(3S)
     >     0.0012D0,   ! 13. Y(4S)
     >     0.011D0,    ! 14. Y(10860)
     >     0.008D0/    ! 15. Y(11020)
      DATA EGT/ ! GeV
     >     0.001D0,  ! 1.  rho
     >     0.08D-3,  ! 2.  omega
     >     0.031D-3, ! 3.  phi
     >     2.8D-6,   ! 4.  J/Psi
     >     8.0D-6,   ! 5.  Psi(2S)
     >     1.0D-3,   ! 6.  Psi(3770)
     >    10.D-3,    ! 7.  Psi(4040)
     >    10.D-3,    ! 8.  Psi(4160)
     >    20.D-3,    ! 9.  Psi(4415)
     >     1.25D-6,  ! 10. Y
     >     2.63D-6,  ! 11. Y(2S)
     >     1.85D-6,  ! 12. Y(3S)
     >     2.5D-3,   ! 13. Y(4S)
     >    28.D-3,    ! 14. Y(10860)
     >    16.D-3/    ! 15. Y(11020)
      DATA EGE/ ! GeV
     >     0.06D-6,  ! 1.  rho
     >     0.02D-6,  ! 2.  omega
     >     0.04D-6,  ! 3.  phi
     >     0.16D-6,  ! 4.  J/Psi
     >     0.04D-6,  ! 5.  Psi(2S)
     >     0.018D-6, ! 6.  Psi(3770)
     >     0.07D-6,  ! 7.  Psi(4040)
     >     0.22D-6,  ! 8.  Psi(4160)
     >     0.07D-6,  ! 9.  Psi(4415)
     >     0.018D-6, ! 10. Y
     >     0.011D-6, ! 11. Y(2S)
     >     0.008D-6, ! 12. Y(3S)
     >     0.029D-6, ! 13. Y(4S)
     >     0.07D-6,  ! 14. Y(10860)
     >     0.03D-6/  ! 15. Y(11020)

*     
* E,MU,TAU,TOP MASSES: 
*
      DATA AMASSES /0.510999D-3,.1056584D0,1.77686d0,172.5D0,
     . 0.066d0,! u
     . 0.066d0,! d
     . 1.2d0,  ! c
     . 0.150d0,! s 
     . 4.6d0  /! b
*
      double complex dalep,dahadr
      common/forfixedorder/dalep,dahadr

      if (ifirst.eq.0) then

c         print*,'resetting e mu masses in loops as the ones in input!'
         amasses(1) = ame
         amasses(2) = ammu

         dalep = (0.d0,0.d0)
         dhadr = (0.d0,0.d0)
         
! filling the common for importance sampling...
         ifirst = 1

         amjpsi   = am(4)
         gtjpsi   = gt(4)
         gejpsi   = ge(4)
         effcjpsi = 3.d0*gejpsi/amjpsi/alpha

         amjpsi2s   = am(5)
         gtjpsi2s   = gt(5)
         gejpsi2s   = ge(5)
         effcjpsi2s = 3.d0*gejpsi2s/amjpsi2s/alpha

         amjpsi3s   = am(10)
         gtjpsi3s   = gt(10)
         gejpsi3s   = ge(10)
         effcjpsi3s = 3.d0*gejpsi3s/amjpsi3s/alpha
         amjpsi4s   = am(11)
         gtjpsi4s   = gt(11)
         gejpsi4s   = ge(11)
         effcjpsi4s = 3.d0*gejpsi4s/amjpsi4s/alpha
         amjpsi5s   = am(12)
         gtjpsi5s   = gt(12)
         gejpsi5s   = ge(12)
         effcjpsi5s = 3.d0*gejpsi5s/amjpsi5s/alpha

         ires    = 0
         ionDHAD = 1
         if (ires.eq.1) ionDHAD = 0

      endif

      if (iarun.eq.0) then
         vpolc = 1.d0
         return
      endif

      if (inskvp.eq.1) then
         vpolc = VPOLCnsk(Q2)
         dalep = 1.d0-1.d0/vpolc
         if (ioffhvp.eq.1) then
            print*,'cannot switch off HVP in NSK. Stopping!'
            stop
         endif
         return
      endif

      if (iteubn.eq.1) then
         vpolc = VPOLChlmnt(Q2)
         dalep = 1.d0-1.d0/vpolc
         return
      endif
      
      SOMMA=0.D0
      DO I=1,4
         xxx  = SUMMA(AMASSES(I),Q2,I)
cc         print*,i,q2,xxx*alpha/pi
         SOMMA=SOMMA+xxx
      ENDDO

      dalphaQED = ALPHA/PI*SOMMA
      toppart = alpha/pi * xxx

      AMODQ2=ABS(Q2)      
      ST2 = 0.2322d0
      qin = q2/abs(q2) * sqrt(abs(q2))

      der = 0.d0
      deltahad = 0.d0      
ccc      goto 333
      
      if (ionDHAD.eq.1.or.qin.lt.0.d0) then
         if (iteubn.eq.0) then
c            CALL DHADR5n12(QIN,ST2,DER,ERRDER,DEG,ERRDEG)
c            der12 = der
c            errder12 = errder
            
c            call dhadr5n17(qin,st2,der,errdersta,errdersys,deg,
c     .           errdegsta,errdegsys) ! it's the 17 version
c            errder = sqrt(errdersta**2+errdersys**2)
c            err17 = errder
c            der17 = der
c            CALL DHADR5n16(QIN,ST2,DER,ERRDER,DEG,ERRDEG) ! it's the 16 version


            der = 0.d0
            errder = 0.d0
            errderst = 0.d0
            deg = 0.d0
            errdeg = 0.d0
            errdegst = 0.d0
            if (ioffhvp.eq.0) then
               CALL DHADR5x(QIN,ST2,DER,ERRDER,errderst,DEG,ERRDEG,
     .              errdegst)   ! it's the 23 version
            endif
c            errsum = errdersta+errdersys
c            print*,' '
c            print*,qin,der17/der,errsum/errder,err17/errder
c            print*,qin,der12/der,errder12/errder
c            print*,'     ',der,der12,der17
            
            if (qin.lt.0.d0) then
               der = der  + iswvperr(3) * scal_err_dhad * errder
c               der = 0.d0
            else
               der = der  + iswvperr(1) * scal_err_dhad * errder
            endif
         elseif (iteubn.eq.1) then
            das = 0.d0
            ddas = 0.d0
            dat = 0.d0
            ddat = 0.d0
            if (ioffhvp.eq.0) then            
               call d_alpha(dabs(qin),das,ddas,dat,ddat)
            endif
            if (qin.lt.0.d0) then
               der = alpha*(das + iswvperr(3) * scal_err_dhad * ddas)
            else
               der = alpha*(dat + iswvperr(1) * scal_err_dhad * ddat)
            endif
         endif
      endif
      
      if (ires.eq.1.and.qin.gt.0.d0) then ! ires = 0 by default
         im = (0.d0,1.d0)
         do j=1,nres
            AMC = AM(J)
            GEC = GE(J)
            GTC = GT(J)
            V   = 3.D0*GEC/AMC / ALPHA
            deltahad = deltahad + 
     .          V * Q2 / (Q2-AMC*AMC+im*GTC*AMC)
         enddo
      endif

c      if (qin.lt.0.d0) der = 0.d0

 333  continue
      
      cder = deltahad + der

c     vpolc = 1.d0/(1.d0-dalphaQED-cder)
c     vpolc = 1.d0 + cder
      
c     vpolc = 1.d0/(1.d0-dalphaQED) + cder ! this has the same peaks as NSK routine if using the resonances!

      vpolc = 1.d0/(1.d0 - dalphaQED - cder)

      dalep  = dalphaQED
      dahadr = cder

      RETURN
      END
******
******
* COMPLEX VACUUM POLARIZATION, T. Teubner & al.
      FUNCTION VPOLChlmnt(Q2)
      implicit double precision (A-H,L,O-Z)
      parameter (np  = 7186)
c      parameter (np=6883)
      parameter (nw  = 1582)

c      parameter (npt = np) ! or
c      parameter (npt = nw) ! or
      parameter (npt = nw+np)

      dimension s(0:np-1)
      dimension st(0:npt-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      double complex vpolcnsk,vpolchlmnt,deltahad,im,res
      double complex cP,ressp,restm,ctmp
      external vpolcnsk
      common/vpcnsk/s,pres,pre,pim,epres,epre,epim,cv,im,ifirstnsk

      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege

      dimension iresbh(8) ! the indexes of the resonances to be added by hand (Thomas said to include **only** 1-6)
      dimension prest(0:npt-1),pret(0:npt-1),pimt(0:npt-1)
      dimension eprest(0:npt-1),epret(0:npt-1),epimt(0:npt-1)
      common/vpcteub/st,prest,pret,pimt,eprest,epret,epimt,iresbh,
     .     ifirst,imidguess
      data ifirst/0/

      double precision wt(nw)
      common/wtabulation/wt

      double precision wtot(npt)
      integer indexes(npt)
      common/wtotcommon/wtot
      
      common/tabulationlimits/wmin,wmax,ic
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/ecms/ecms,ecmsnom

      integer iloc(1)

      dimension iswvperr(3),iswvperrteub(6)
      common/iswitchvperror/iswvperr,iswvperrteub

      common/rescaledhaderr/scal_err_dhad
      
*     using vp_hlmnt_v2_1.f, sent by T. Teubner, private communication
      
      if (ifirst.eq.0) then
***   filling the common by calling vpolcnsk once (I need the array s --> st[eubner])
         vpolchlmnt = vpolcnsk(q2)
****  filling the common wtabulation    
         wtot = 0.d0
         call fillwtabcommon
         wmin = wt(1)
         wmax = wt(nw)
         wtot(1:nw)     = wt**2
         wtot(nw+1:npt) = s
         call sort_vector_inverse(wtot,wtot,indexes,npt)
c         do k = 1,npt
c            print*,k,sqrt(wtot(k)),indexes(k)
c         enddo
c         stop
****  
c
c I solve for P, so that is the same quantity returned by NSK
c     I also use the same evaluation points as NSK, and then linear
c     interpolation between them
c     1/(1-P) = 1/(1-danr) + res ->
c     1-P = (1-danr)/(1+(1-danr)*res)
c       P = 1 - (1-danr)/(1+(1-danr)*res)
c         = (1+(1-danr)*res - (1-danr))/(1+(1-danr)*res)
c         = ( danr + res - danr*res) / (1 + res- danr*res)
******************
* filling st, from the nsk data or Teubner tabulation points
c         st(0:np-1) = s   ! or
c         do k = 1,nw
c            st(k-1) = wt(k)*wt(k)
c         enddo ! or
         st = wtot
******************
         prest  = 0.d0
         pret   = 0.d0
         pimt   = 0.d0
         eprest = 0.d0
         epret  = 0.d0
         epimt  = 0.d0
         iresbh(1) = 4
         iresbh(2) = 5
         iresbh(3) = 10
         iresbh(4) = 11
         iresbh(5) = 12
         iresbh(6) = 13
         iresbh(7) = 14
         iresbh(8) = 15
c         print*,'Evaluating Teubner et al. VP array '
c         open(10,file='cfrvp-0-v2-3',status='unknown')
         do k = 0,npt-1
            energy = sqrt(st(k))

c            hlmnt is tabulated up to wmax and from wmin
            if (energy.le.wmax.and.energy.ge.wmin) then
               ressp = 0.d0
               restm = 0.d0
! I do it below, every time the function is called, better for error handling...               
c$$$               do i=1,6
c$$$                  j = iresbh(i)
c$$$                  AMC = AM(J) !- EAM(J)
c$$$                  GEC = GE(J) !- EGE(J)
c$$$                  GTC = GT(J) !- EGT(J)
c$$$                  V   = 3.D0*GEC/AMC / ALPHA
c$$$                  restm = restm + 
c$$$     .                 V * st(k) / ( st(k)-AMC*AMC+im*GTC*AMC )
c$$$                  ressp = ressp - 
c$$$     .                 V * st(k) / (-st(k)-AMC*AMC+im*GTC*AMC )
c$$$               enddo
               nrflag = 1
               call vphlmntv2(energy,
     &              vprehadsp,vprehadtm,vpimhad,vprelepsp,
     &              vpreleptm,vpimlep,
     &              vpretopsp,vpretoptm,nrflag,evprehadsp,
     &              evprehadtm,evpimhad)
c               vprehadsp = vprehadsp - evprehadsp
c               vprehadtm = vprehadtm - evprehadtm
c               vpimhad   = vpimhad -  evpimhad
!     ! space like
               prest(k) = vprehadsp + vprelepsp + vpretopsp
               eprest(k) = evprehadsp
!     !     time like
               nrflag = 0
               call vphlmntv2(energy, 
     &        vprehadsp,vprehadtm,vpimhad,vprelepsp,vpreleptm,vpimlep,
     &        vpretopsp,vpretoptm,nrflag,evprehadsp,evprehadtm,evpimhad)

               pret(k)  = vprehadtm + vpreleptm + vpretoptm
               pimt(k)  = vpimhad   + vpimlep

               epret(k) = evprehadtm
               epimt(k) = evpimhad

c     [[ useless, but harmless               
c   this is 1/(1-P)
               ctmp = 1.d0/(1.d0 - pret(k)-im*pimt(k)) + restm ! in the time-like not resummed!
c  ctmp = 1/(1-P) -> 1-P = 1/ctmp  -> P = 1-1/ctmp
               pret(k) = DReal(1.d0-1.d0/ctmp)
               pimt(k) = DImag(1.d0-1.d0/ctmp)
c   ]]
c               write(10,*)energy,pret(k),pimt(k),pre(k),pim(k),
c     .              prest(k),pres(k)
c               write(10,*)energy,pret(k),pimt(k),prest(k)
            endif
         enddo
c         close(10)
c         print*,'done!'
c         print*,'stopping line ~420 sv.f'
c         stop
*** guessing where to start the binary search....         
         k = 0
         aq2 = ecmsnom**2
         do while(aq2.gt.st(k))
            k = k + 1         
         enddo
         k = k - 1
         imidguess = k ! this is correct for born s-channel and no energy spread
***************
         ic = 0
         ifirst=1
      endif
c      k   = 0
c      sk  = 0.d0
c      aq2 = abs(q2)
c      do while(aq2.gt.st(k))
c         k = k + 1         
c      enddo
c      k = k - 1
******** binary search ****
      aq2     = abs(q2)
      istart  = 0
      ifinish = npt-1
c      imid = (istart + ifinish)/2
      imid = imidguess
      
c     do while(ifinish.ge.istart)
      do while(.not.(st(imid+1).gt.aq2.and.st(imid).le.aq2))
         if (aq2.gt.st(imid)) then
            istart  = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
******************
      if (q2.gt.0.d0) then
         vpre = pret(k)+(pret(k+1)-pret(k))/(st(k+1)-st(k)) * (q2-st(k))
         vpim = pimt(k)+(pimt(k+1)-pimt(k))/(st(k+1)-st(k)) * (q2-st(k))

         evpre=epret(k)+(epret(k+1)-epret(k))/(st(k+1)-st(k))*(q2-st(k))
         evpim=epimt(k)+(epimt(k+1)-epimt(k))/(st(k+1)-st(k))*(q2-st(k))

         ctmp = 0.d0
         do i=1,6
            j    = iresbh(i)
            AMC  = AM(J) + EAM(J)*iswvperrteub(4)*scal_err_dhad
            GEC  = GE(J) + EGE(J)*iswvperrteub(5)*scal_err_dhad
            GTC  = GT(J) + EGT(J)*iswvperrteub(6)*scal_err_dhad
            V    = 3.D0*GEC/AMC / ALPHA
            ctmp = ctmp + 
     .           V * q2 / ( q2 - AMC*AMC+im*GTC*AMC )
         enddo

         vpre = vpre + iswvperrteub(1)*scal_err_dhad*evpre
         vpim = vpre + iswvperrteub(2)*scal_err_dhad*evpim
         
      else
         vpre = prest(k) +
     .        (prest(k+1)-prest(k))/(st(k+1)-st(k))*(aq2-st(k))
         evpre = eprest(k) +
     .        (eprest(k+1)-eprest(k))/(st(k+1)-st(k))*(aq2-st(k))

         vpim = 0.d0
         evpim = 0.d0

         ctmp = 0.d0

         vpre = vpre + iswvperrteub(3)*scal_err_dhad*evpre
      endif
*
! ctmp is 0 in the spacelike, not 0 in the timelike. It has not to be resummed
      VPOLChlmnt = 1.d0/(1.d0  - vpre - im*vpim) + ctmp
      return
      end
*
*     COMPLEX VACUUM POLARIZATION from NSBSK DATA: ALPHA(0)->ALPHA(0)*VPOL(Q2)
      FUNCTION VPOLCnsk(Q2) 
      implicit double precision (A-H,L,O-Z)
      parameter (np=7186)
c      parameter (np=6883)
      dimension s(0:np-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)
      double complex vpolcnsk,deltahad,im
      common/vpcnsk/s,pres,pre,pim,epres,epre,epim,cv,im,ifirst
      data ifirst/0/

      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt

      common/rescaledhaderr/scal_err_dhad
* from    http://cmd.inp.nsk.su/~ignatov/vpl/
      
      if (ifirst.eq.0) then
         im = (0.d0,1.d0)
c         open(45,file='vpol_novosibirsk_NEW.dat',status='old')
         open(45,file='vpol_novosibirsk.dat',status='old')         
         do k = 0,np-1
            read(45,*)s(k),pres(k),pre(k),pim(k),epres(k),epre(k),
     .           epim(k),cv(k)
         enddo
         close(45)
         ifirst=1
      endif

      k   = 0
      sk  = 0.d0
      aq2 = abs(q2)
      do while(aq2.gt.s(k))
         k = k + 1         
      enddo
      k = k - 1
      
      if (q2.gt.0.d0) then
         vpre = pre(k) + (pre(k+1)-pre(k))/(s(k+1)-s(k)) * (q2-s(k))
         vpim = pim(k) + (pim(k+1)-pim(k))/(s(k+1)-s(k)) * (q2-s(k))

         evpre = epre(k) + (epre(k+1)-epre(k))/(s(k+1)-s(k)) * (q2-s(k))
         evpim = epim(k) + (epim(k+1)-epim(k))/(s(k+1)-s(k)) * (q2-s(k))

         rho = cv(k) + (cv(k+1)-cv(k))/(s(k+1)-s(k)) * (q2-s(k))
         cov = rho * evpre*evpim
      else
         vpre = pres(k) + (pres(k+1)-pres(k))/(s(k+1)-s(k))*(aq2-s(k))
         vpim = 0.d0

         evpre = epres(k)+(epres(k+1)-epres(k))/(s(k+1)-s(k))*(aq2-s(k))
         evpim = 0.d0

         rho = 0.d0
         cov = 0.d0
      endif

c      print*,iswvperr
      
      if (abs(scal_err_dhad).gt.1d-4) then
c     vpre = vpre + scal_err_dhad * evpre
c     vpim = vpim + scal_err_dhad * evpim
*     
*     USING FORMALISM FROM GLEN COWAN, "STATISTICAL DATA ANALYSIS", p.22+, OXFORD SCIENCE PUBLICATIONS
         if (q2.lt.0.d0) then

            vpre = vpre + scal_err_dhad * evpre * iswvperr(3)
            
         elseif (iswvperr(1).ne.0.or.iswvperr(2).ne.0) then
            x1   = vpre
            x2   = vpim
            s12  = evpre*evpre
            s22  = evpim*evpim
            rho2 = rho * rho
            sq   = sqrt((s12+s22)**2-4.d0*(1.d0-rho2)*s12*s22)
            lp   = 0.5d0*(s12 + s22 + sq)
            lm   = lp - sq
            th   = 0.5d0 * atan(2.d0*rho*evpre*evpim/(s12-s22))
            cs   = cos(th)
            ss   = sin(th)
            y1   =  cs * x1 + ss * x2
            y2   = -ss * x1 + cs * x2
            
            dy1 = scal_err_dhad * iswvperr(1) * sqrt(lp)
            dy2 = scal_err_dhad * iswvperr(2) * sqrt(lm)
! changing y1 and y2...
            y1   = y1 + dy1
            y2   = y2 + dy2
! ... going back to x1 and x2
            x1   = cs * y1 - ss * y2
            x2   = ss * y1 + cs * y2
            
!     updating the variables
            vpre = x1
            vpim = x2
         endif
      endif
*********************************************************************************************
      VPOLCnsk = 1.d0/(1.d0  - vpre - im*vpim)
      return
      end
*==============================================================================================
*     LEPTONIC AND TOP CONTRIBUTION TO VACUUM POLARIZATION
** with correct expansion for |q2| <<< m2!!
      FUNCTION SUMMA(AM,Q2,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NC(9),QF2(9)
      common/ncqf2/nc,qf2
      data nc  /1.d0,1.d0,1.d0,3.d0,
     . 3.d0, !u
     . 3.d0, !d
     . 3.d0, !c
     . 3.d0, !s
     . 3.d0 /!b
      data qf2 /1.d0,1.d0,1.d0,0.44444444444444444444d0,
     . 0.44444444444444444444d0, !u
     . 0.11111111111111111111d0, !d
     . 0.44444444444444444444d0, !c
     . 0.11111111111111111111d0, !s
     . 0.11111111111111111111d0 /!b
* NC AND QF ARE COLOR FACTOR (1 FOR LEPTONS, 
* 3 FOR TOP) AND CHARGE
      AM2=AM**2

      if (Q2.eq.0.d0) then
         summa = 0.d0
         return
      endif

      z6lim = 1d-8
      
      IF (Q2.GT.0.D0.AND.Q2.LT.(4.D0*AM2)) THEN
         x = 4.D0*AM2/Q2
         z = sqrt(1.d0/x)
         z2 = z*z
         z4 = z2*z2
         z6 = z2*z4

         if (z6.gt.z6lim) then
            SQ = DSQRT(1d0 - z2)
            s1 = -5.D0/3.D0 - 1.d0/z2
            s2 = (1.d0 + z2- 2.D0*z4)/z/z2*1.d0/SQ*DATAN(z/SQ)
            s1ps2 = s1 + s2
         else
*     with wxmaxima
* s2 = taylor((1+z^2-2*z^4)/sqrt(1-z^2)/z^3*atan(z/sqrt(1-z^2)), z , 0 , 20);
            s1ps2 = -z2*
     .        (4.d0/5.d0+12.d0/35.d0*z2+64.d0/315.d0*z4+32.d0/231.d0*z6)
         endif
         SUMMA=1.d0/3.d0*NC(I)*QF2(I)*s1ps2
      ELSE
         x = -4.D0*AM2/Q2
         if (x.gt.0.d0) then    ! space-like
            z  = sqrt(1.d0/x)
            z2 = z*z
            z4 = z2*z2
            z6 = z2*z4
            
            if (z6.gt.z6lim) then
               SQ     = DSQRT(1.D0 + z2)
               ARGLOG = DABS( (z - sq)/(z+sq) )              
               s1 = -5.D0/3.D0 + 1.d0/z2
               s2 = 0.5d0*(1.d0 - z2 - 2.D0*z4)/SQ*DLOG(ARGLOG)/z2/z
               
               s1ps2 = s1 + s2
            else
*     with wxmaxima
* s2 = taylor(log( abs ((z - sqrt(1+z^2))/(z + sqrt(1+z^2)) ))/sqrt(1+z^2)*(1-z^2-2*z^4)/2/z^3,z,0,20); 
               s1ps2 = z2*
     .        (4.d0/5.d0-12.d0/35.d0*z2+64.d0/315.d0*z4-32.d0/231.d0*z6)
            endif
         else                   ! for x > 0 ie q2 > 4m^2
            x = -x
            z  = sqrt(1.d0/x)
            z2 = z*z
            z4 = z2*z2
            
            SQ     = DSQRT(z2-1.d0)
            ARGLOG = DABS((z - sq)/(z+sq) )            
            s1 = -5.D0/3.D0 - 1.d0/z2
            s2 =  0.5d0*(1.d0 + z2 - 2.D0*z4)/SQ*DLOG(ARGLOG)/z2/z
            s1ps2 = s1 + s2
         endif
         SUMMA = 1.d0/3.d0 * NC(I)*QF2(I)*s1ps2
      ENDIF
      RETURN
      END
*-----
* LEPTONIC AND TOP CONTRIBUTION TO VACUUM POLARIZATION
      FUNCTION SUMMA_OLD(AM,Q2,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NC(9),QF2(9)
      common/ncqf2/nc,qf2
      data nc  /1.d0,1.d0,1.d0,3.d0,
     . 3.d0, !u
     . 3.d0, !d
     . 3.d0, !c
     . 3.d0, !s
     . 3.d0 /!b
      data qf2 /1.d0,1.d0,1.d0,0.44444444444444444444d0,
     . 0.44444444444444444444d0, !u
     . 0.11111111111111111111d0, !d
     . 0.44444444444444444444d0, !c
     . 0.11111111111111111111d0, !s
     . 0.11111111111111111111d0 /!b

* NC AND QF ARE COLOR FACTOR (1 FOR LEPTONS, 
* 3 FOR TOP) AND CHARGE
*
      AM2=AM**2
      IF (Q2.GE.0.D0.AND.Q2.LT.(4.D0*AM2)) THEN
         SQ=SQRT(4.D0*AM2/Q2-1.D0)
         SUMMA_OLD=NC(I)*QF2(I)*(-5.D0/9.D0-(4.D0/3.D0)*(AM2/Q2)+
     >        (4.D0/3.D0*(AM2/Q2)**2+1.D0/3.D0*AM2/Q2-1.D0/6.D0)*
     >        4.D0/SQ*ATAN(1.D0/SQ))
        ELSE
           SQ=SQRT(1.D0-4.D0*AM2/Q2)
           ARGLOG=ABS((1.D0-SQ)/(1.D0+SQ))
*
           SUMMA_OLD=NC(I)*QF2(I)*(-5.D0/9.D0-(4.D0/3.D0)*(AM2/Q2)+
     >          (4.D0/3.D0*(AM2/Q2)**2+1.D0/3.D0*AM2/Q2-1.D0/6.D0)*
     >          2.D0/SQ*LOG(ARGLOG))
       ENDIF
*
       RETURN
       END
***********************************************************************
      subroutine svfactor(model,ng,ecms,p1,p2,eps,sv,dsv)
#ifdef COLLIER
      use collier
#endif      
      implicit double precision (a-h,l,o-z)
      dimension p1(0:3),p2(0:3),pin1(0:3),pin2(0:3),p1b(0:3),p2b(0:3)
      character*10 model
      character*6 ord
      double precision masses(4),charges(4)
      character*2 fs
      common/finalstate/fs
      common/qedORDER/ord
      common/parameters/ame,ammu,convfac,alpha,pi
      common/hparam/ampi
      common/momentainitial/pin1,pin2
      common/reducedtoborn/p1b,p2b,iref
      common/momentainitialred/pin1b(0:3),pin2b(0:3),pin1r(0:3)
     >     ,pin2r(0:3)
      common/forborncrosssection/phsp2b,flux2b,bornme,bornmeq2,bornmez,
     >     bornmenovp
      common/ialpharunning/iarun
      common/icountsv/icsvf,ifirst
      data ifirst,icsvf /0,0/
      integer irecola
      common/recolamode/irecola
      common/various/beta,betafs
      common/softphotonmass/phmass
      character*100  path
      common/cmnpath/path
      double precision lambda2,mudim,mu2dim,m2min
      common/duvmudim/deltauv,mudim,mu2dim,m2min,pi216m1,ec6
      common/gammamass/lambda2
      integer ioffhvp ,inlovp
      common/offhadrVP_nlovp/ioffhvp,inlovp   
   
      integer iffpi
      character*10 what_ffpi
      common/ffpioption/iffpi,what_ffpi 
    
      real*8 virtplussoftpions
      external virtplussoftpions  
#ifdef LOOPTOOLS
#include "looptools.h"
#endif
      icsvf = icsvf + 1

      dsv = 0.d0
      sv  = 1.d0

      if (ord.eq.'born'.or.ord.eq.'struct') then
         sv = 1.d0
         return
      endif
      
      if (ifirst.eq.0) then
         ifirst = 1

         phmass = 1d-8
         deltauv= 0.d0
         mu2dim= 1.d0
         lambda2= phmass*phmass
         m2min= 0.d0

         pi216m1 = pi*pi*16.d0
         pi216m1 = 1.d0/pi216m1

         ec6 = 4.d0*pi*alpha
         ec6 = ec6*ec6*ec6
         
#ifdef LOOPTOOLS         
         call ltini
         call setdelta(deltauv)
         call setmudim(mu2dim)
         call setlambda(lambda2)
         call setminmass(m2min)
#endif
#ifdef COLLIER
         call Init_cll(4,2,path(1:len_trim(path))//'cll_out')
         call SetMode_cll(1)    ! 1 is COLI, 2 DD, 3 both. 1 is faster
         call InitCacheSystem_cll(1,4)
         call SetDeltaUV_cll(deltauv)
         call SetMuUV2_cll(mu2dim)
         call SetmuIR2_cll(phmass*phmass)
#endif         
      endif

      if (ord.eq.'alpha') then
         if (ng.ge.1) then
            sv = 1.d0
            return
         endif         
         if (ng.eq.0) then
            s = ecms**2
            c = p1(3)/dsqrt(tridot(p1,p1))
            t = -s/2.d0*(1.d0 - c)            

            sigma0  = bornme*phsp2b/flux2b
            sigmaq2 = bornmeq2*phsp2b/flux2b
**
            if (model.ne.'ps') then
**
            if (fs.eq.'gg') then
               if (irecola.eq.0) then
c     sv = aintsvgg(s,c,eps)/sigma0
c i.e.                  
c                  sv1 = aintsvgg(s,c,eps)/sigma0

                  BORNBK=2.D0*PI*ALPHA*ALPHA/S*BETA*BETA/
     >                 (1.D0+BETA*C)/(1.D0-BETA*C)*(1.D0+C*C)
                  BORNBK = BORNBK*S*128d0*pi

                  svgg = aintsvgg(s,c,eps)/bornbk*bornme ! to restore full mass effects in the IR
c                  phmass = 1.d-10
c                ssf=soft_integral_gg(eps*sqrt(s)*0.5d0,pin1,pin2,p1,p2)
c                ssf = ssf /128.d0/s/pi
c                  print*,svgg-ssf*bornme
c                stop

                  sv = svgg/bornme/phsp2b*flux2b
               else
                  call recola_ampl(0,pin1,pin2,p1,p2,p2,brcl,svrcl)
                  sv=svrcl/brcl
c                  print*,svrcl
c                  stop                  
               endif
            elseif(fs.eq.'pp'.or.fs.eq.'pi') then
	       if(iffpi.le.2) then
	           sv = aintsv_pions(pin1,pin2,p1,p2,eps)/sigma0
               elseif(iffpi.eq.3) then
                   sv = 1d0+virtplussoftpions(pin1,pin2,p1,p2,eps)
               endif
               smtob = 0.d0
               if (inlovp.eq.1) then
                  call get_twiceVPtimesTreeOverBorn(pin1,pin2,p1,p2,
     .                 smtob)
               endif
               sv = sv + smtob
               dsv= sv-1.d0     !order alpha correction
            else
c               sv = aintsv(s,t,eps)/sigma0
               sv = aintsv_new(pin1,pin2,p1,p2,eps)/sigma0
               smtob = 0.d0
               if (inlovp.eq.1) then
                  call get_twiceVPtimesTreeOverBorn(pin1,pin2,p1,p2,
     .                 smtob)
               endif
                  sv =  sv + smtob
            endif
************* TEST 2 LOOP**********************
c$$$            aieps = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
c$$$            aiepsnir = -1.5d0+2.d0*eps-0.5d0*eps**2
c$$$
c$$$            dp1p2 = dlog(2.d0*dot(pin1,pin2)/ame**2)
c$$$            dp1p3 = dlog(2.d0*dot(pin1,p1)/ame**2)
c$$$            dp1p4 = dlog(2.d0*dot(pin1,p2)/ame**2)
c$$$            dp2p3 = dlog(2.d0*dot(pin2,p1)/ame**2)
c$$$            dp2p4 = dlog(2.d0*dot(pin2,p2)/ame**2)
c$$$            dp3p4 = dlog(2.d0*dot(p1,p2)/ame**2)
c$$$            
c$$$            lcollnew = dp1p2 + dp1p3 - dp1p4 - dp2p3 + dp2p4 + dp3p4
c$$$            lcollnew = 2.d0 * lcollnew - 4.d0
c$$$            lcollnew = lcollnew/4.d0
c$$$            arg  = 2.d0*alpha/pi * lcollnew * aieps
c$$$            argnir  = 2.d0*alpha/pi * lcollnew * aiepsnir
c$$$            vips = - arg
c$$$            vipsnir = - argnir
c$$$
c$$$            call aintsv_work(s,t,eps,
c$$$     .           svs,svt,svst,bs,bt,bst,das,dat)
c$$$            dast = das+dat
c$$$
c$$$            svt = svt/sigma0
c$$$            svs = svs/sigma0
c$$$            svst = svst/sigma0
c$$$
c$$$            bs   = bs/sigma0
c$$$            bt   = bt/sigma0
c$$$            bst  = bst/sigma0
c$$$
c$$$            btot = bs+bt+bst
c$$$            
c$$$            Ds  = svs - vips*bs
c$$$            Dt  = svt - vips*bt
c$$$            Dst = svst - vips*bst
c$$$
c$$$            Dsnir  = Ds
c$$$            Dtnir  = Dt
c$$$            Dstnir = Dst
c$$$
c$$$            sv1l  = svt+svst+svs+bs*2.d0*das+bt*2.d0*dat+bst*dast
c$$$
c$$$            pippo  = svt+svst+svs
c$$$            ciccio = vips * (bs+bt+bst)
c$$$
c$$$            a1 = vips
c$$$            b1 = Ds+Dt+Dst
c$$$            c1 = 2.d0*das*bs+2.d0*dat*bt+(das+dat)*bst
c$$$            a1nir = vipsnir
c$$$            b1nir = Dsnir+Dtnir+Dstnir
c$$$            c1nir = 2.d0*das*bs+2.d0*dat*bt+(das+dat)*bst
c$$$
c$$$            a2 = vips**2/2.d0
c$$$            b2 = 2.d0*Ds*das+2.d0*Dt*dat+Dst*(das+dat)
c$$$            c2 = 3.d0*bs*das**2+3.d0*bt*dat**2+
c$$$     .           bst*(das**2+dat**2+das*dat)
c$$$
c$$$            a2nir = vipsnir**2/2.d0
c$$$            b2nir = 2.d0*Dsnir*das+2.d0*Dtnir*dat+Dstnir*(das+dat)
c$$$            c2nir = 3.d0*bs*das**2+3.d0*bt*dat**2+
c$$$     .           bst*(das**2+dat**2+das*dat)
c$$$
c$$$            a1c1 = a1*c1
c$$$            a1b1 = a1*b1
c$$$            b1c1 = b1*c1
c$$$
c$$$            a1nirc1nir = a1nir*c1nir
c$$$            a1nirb1nir = a1nir*b1nir
c$$$            b1nirc1nir = b1nir*c1nir
c$$$
c$$$            sv2l = a2+b2+c2+a1b1+a1c1+b1c1
c$$$
c$$$            sv2lnir = a2nir+b2nir+c2nir+a1nirb1nir+a1nirc1nir+b1nirc1nir
c$$$
c$$$            a2Lsv  = 2.d0*alpha/pi * (lcollnew + 0.d0)*
c$$$     .           (1.5d0-2.d0*eps+0.5d0*eps**2) * b1
c$$$
c$$$            sv = a2Lsv
c$$$
c$$$c            sv = svt+svs+svst
c$$$c            sv = a1
c$$$            
c$$$*** only photonic, to compare to Penin
c$$$            das = 0.d0
c$$$            dat = 0.d0
c$$$            dast = das+dat
c$$$            a1 = vips
c$$$            b1 = Ds+Dt+Dst
c$$$            c1 = 2.d0*das*bs+2.d0*dat*bt+(das+dat)*bst
c$$$
c$$$            a2 = vips**2/2.d0
c$$$            b2 = 2.d0*Ds*das+2.d0*Dt*dat+Dst*(das+dat)
c$$$            c2 = 3.d0*bs*das**2+3.d0*bt*dat**2+
c$$$     .           bst*(das**2+dat**2+das*dat)
c$$$            a1c1 = a1*c1
c$$$            a1b1 = a1*b1
c$$$            b1c1 = b1*c1
c$$$
c$$$c            print*,a2,b2,c2,a1b1,a1c1,b1c1
c$$$
c$$$            sv2lphot = a2+b2+c2+a1b1+a1c1+b1c1
c$$$
c$$$            sv2lphotnir = a2nir+a1nir*b1nir
c$$$
c$$$c            sv2lphotnira2L = a1nir*b1nir
c$$$c            sv = sv2lphotnira2L
c$$$
c$$$            sv = sv2lphot
c$$$            
c$$$c     print*,sv2lphotnir
c$$$** virtual pair
c$$$            ammu = ame
c$$$            L   = dlog(-t/ammu/ammu)
c$$$            L2  = L*L
c$$$            L3  = L2*L
c$$$            cf  = 383.d0/108.d0-11.d0/6.d0*pi*pi/6.d0
c$$$            z3  = 1.2020569032d0
c$$$            cfp = -z3/3.d0+3355.d0/1296.d0-19.d0/18.d0*pi**2/6.d0 
c$$$            
c$$$            fmu = -1.d0/36.d0*L3
c$$$     .           +19.d0/72.d0*L2
c$$$     .           +(-1.d0/36.d0*pi*pi-265.d0/216.d0)*L+
c$$$     .           cf + 19.d0/75.d0*pi*pi
c$$$            
c$$$            v = 4.d0*(alpha/pi)**2*fmu ! full virtual for t channel only
c$$$
c$$$            svpair = bt*v
c$$$
c$$$            oursvnf1 = sv2l-sv2lphot+svpair
c$$$            oursvnf1nir = sv2lnir-sv2lphotnir+svpair
c$$$
c$$$
c$$$            energia = ecms/2.d0
c$$$            de = 0.2d0 * ecms
c$$$            pl = -dlog(de/energia)
c$$$            cl =  dlog(abs(t)/ame/ame)
c$$$
c$$$            pallarealpairs = 2.d0*
c$$$     .           alpha**2/pi/pi * 1.d0/18.d0 *
c$$$     .           (        (cl - 2.d0*pl)**3
c$$$     .            -5.d0 * (cl - 2.d0*pl)**2
c$$$     .           +12.d0 * (cl - 2.d0*pl) *(14.d0/9.d0 - pi*pi/12.d0)
c$$$     .            )
c$$$
c$$$            pallarealpairs = bt * pallarealpairs
c$$$
c$$$C coppie spente !!!!!!!!
c$$$c            sv2l = sv2l + v + pair
c$$$*************
c$$$               sv = sv2lphot !- a1nir*b1nir*0.5d0
c$$$            esoft = eps*pin1(0)
c$$$            call twoloopvirtualsoft(pin1,pin2,p1,p2,esoft,
c$$$     .           unloop,dueloop,pair,penin)
c$$$
c$$$            boncianinf1 = dueloop - pair - penin
c$$$
c$$$            sv = sv2l-dueloop
c$$$            
c$$$
c$$$            sv = penin - sv2lphot !+ a1nir*b1nir*0.5d0
c$$$            
c$$$c            print*,penin,sv2lphot,penin/sv2lphot
c$$$
c$$$c            sv = oursvnf1nir
c$$$c            sv = boncianinf1 - oursvnf1
c$$$c            sv = pallarealpairs
c$$$cc            sv = dueloop - pair - (sv2l + svpair)
c$$$c            sv = a2Lsv
c$$$c            sv = a1+b1+c1 - unloop
************* End TEST **********************
         elseif (model.eq.'ps') then
               aieps     = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
               aiepsscal = -2.d0*dlog(eps)-2.d0+2.d0*eps
! see below
               dp1p2 = dlog(2.d0*dot(pin1,pin2)/ame**2)
               dp1p3 = dlog(2.d0*dot(pin1,p1)/ame**2)
               dp1p4 = dlog(2.d0*dot(pin1,p2)/ame**2)
               dp2p3 = dlog(2.d0*dot(pin2,p1)/ame**2)
               dp2p4 = dlog(2.d0*dot(pin2,p2)/ame**2)
               dp3p4 = dlog(2.d0*dot(p1,p2)/ame**2)

               if (fs.eq.'ee'.or.fs.eq.'mm'
     +              .or.fs.eq.'pp'.or.fs.eq.'pi') then
                  lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
                  lcollnew = 2.d0 * lcollnew - 4.d0
                  lcollnew = lcollnew/4.d0
                  arg = 2.d0*alpha/pi * lcollnew * aieps
**** NEW ***************
                  charges(1) = -1.d0
                  charges(2) = -1.d0
                  charges(3) = -1.d0
                  charges(4) = -1.d0
                  masses(1)  = ame
                  masses(2)  = ame
                  masses(3)  = ame
                  masses(4)  = ame
                  if (fs.eq.'mm') then
                     masses(3)  = ammu
                     masses(4)  = ammu
                  endif
                  if (fs.eq.'pp') then
                     masses(3)  = ampi
                     masses(4)  = ampi
                  endif
                  if (fs.eq.'pi') then
                     charges(3) = 0.d0
                     charges(4) = 0.d0
                     masses(3)  = ampi
                     masses(4)  = ampi
                  endif
                  lcoll1=eikonalintegral(pin1,pin2,p1,p2,charges,masses)
                  arg2  = 2.d0*alpha/pi * lcoll1 * aieps 
                  if (fs.eq.'pp') then
                     arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
                  endif
                  arg = arg2
**************************
               elseif (fs.eq.'gg') then
                  lcollnew = dp1p2
                  lcollnew = lcollnew - 1.d0
                  arg = alpha/pi * lcollnew * aieps
               endif
               sv = 1.d0 - arg
            endif
            return
         endif
      endif

      if (ord.eq.'exp') then
         aieps     = -2.d0*dlog(eps)-1.5d0+2.d0*eps-0.5d0*eps**2
         aiepsscal = -2.d0*dlog(eps)-2.d0+2.d0*eps
! with reduced to born momenta...
         aopi = alpha/pi
         s = 4.d0*pin1b(0)*pin1b(0)

! for gg not-reduced center of mass..
ccccc         if (fs.eq.'gg') s = ecms**2

         p1mod = sqrt(tridot(p1b,p1b))
         c = p1b(3)/p1mod
         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)

         ame2mu = 1.d0/ame/ame
         sdif = bornme*phsp2b/flux2b

** VP STUFF
*** changed for no VP in 1-loop
         sdif = bornmenovp*phsp2b/flux2b
***         
         if (fs.eq.'ee'.or.fs.eq.'mm'
     .        .or.fs.eq.'pp'.or.fs.eq.'pi') then
            scale = s*t/u
            lcoll = dlog(scale*ame2mu) - 1.d0
            arg = 2.d0 * aopi * lcoll * aieps

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mm') then
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            if (fs.eq.'pp') then
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            if (fs.eq.'pi') then
               charges(3) = 0.d0
               charges(4) = 0.d0
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            lcoll1 = eikonalintegral(pin1b,pin2b,p1b,p2b,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps 
            if (fs.eq.'pp') then
               arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
            endif
            arg = arg2
**************************
** VP STUFF
ccc   changed for no VP in virtual
            iarunbck = iarun ! added
            iarun = 0           ! added ! ON OFF
            if (model.ne.'ps') then
               if (fs.eq.'pp'.or.fs.eq.'pi') then
                 if(iffpi.le.2) then 
                 svpureoal=aintsv_pions(pin1b,pin2b,p1b,p2b,eps)/sdif
     .                 -1.d0
                 elseif(iffpi.eq.3) then
                 svpureoal=virtplussoftpions(pin1b,pin2b,p1b,p2b,eps)   
                 endif
               else
                  svpureoal = aintsv_new(pin1b,pin2b,p1b,p2b,eps)/sdif
     .                 -1.d0
               endif
            else 
            svpureoal = -arg
            endif
            iarun = iarunbck    ! added

         elseif (fs.eq.'gg') then
            scale = s
            lcoll = dlog(scale*ame2mu) - 1.d0
            arg = 1.d0 * aopi * lcoll * aieps
            svpureoal = aintsvgg(s,c,eps)/sdif - 1.d0
         endif
         deltasv = svpureoal + arg !!! NB

*** NEW SFF, formula 2.23 PhD. Thesis. in the limit dot(pi,pj) >> m**2
         dp1p2 = dlog(2.d0*dot(pin1,pin2)*ame2mu)
         dp1p3 = dlog(2.d0*dot(pin1,p1)*ame2mu)
         dp1p4 = dlog(2.d0*dot(pin1,p2)*ame2mu)
         dp2p3 = dlog(2.d0*dot(pin2,p1)*ame2mu)
         dp2p4 = dlog(2.d0*dot(pin2,p2)*ame2mu)
         dp3p4 = dlog(2.d0*dot(p1,p2)*ame2mu)

cc         dp1p2 = dlog(2.d0*dot(pin1b,pin2b)*ame2mu)
cc         dp1p3 = dlog(2.d0*dot(pin1b,p1b)*ame2mu)
cc         dp1p4 = dlog(2.d0*dot(pin1b,p2b)*ame2mu)
cc         dp2p3 = dlog(2.d0*dot(pin2b,p1b)*ame2mu)
cc         dp2p4 = dlog(2.d0*dot(pin2b,p2b)*ame2mu)
cc         dp3p4 = dlog(2.d0*dot(p1b,p2b)*ame2mu)

         argbck = arg
         
         if (fs.eq.'ee'.or.fs.eq.'mm'
     .        .or.fs.eq.'pp'.or.fs.eq.'pi') then
            lcollnew = dp1p2+dp1p3-dp1p4-dp2p3+dp2p4+dp3p4
            lcollnew = 2.d0 * lcollnew - 4.d0
            lcollnew = lcollnew/4.d0
            arg = 2.d0 * aopi * lcollnew * aieps

**** NEW ***************
            charges(1) = -1.d0
            charges(2) = -1.d0
            charges(3) = -1.d0
            charges(4) = -1.d0
            masses(1)  = ame
            masses(2)  = ame
            masses(3)  = ame
            masses(4)  = ame
            if (fs.eq.'mm') then
               masses(3)  = ammu
               masses(4)  = ammu
            endif
            if (fs.eq.'pp') then
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            if (fs.eq.'pi') then
               charges(3) = 0.d0
               charges(4) = 0.d0
               masses(3)  = ampi
               masses(4)  = ampi
            endif
            lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,masses)
            arg2   = 2.d0*aopi * lcoll1 * aieps 
            if (fs.eq.'pp') then
               arg2 = alpha/pi*lcoll1*(aieps+aiepsscal)
            endif
            arg = arg2
**************************
         elseif (fs.eq.'gg') then
            lcollnew = dp1p2
            lcollnew = lcollnew - 1.d0
            arg = 1.d0 * aopi * lcollnew * aieps
         endif

!         deltasv = svpureoal + arg      !!! NB
         ffs = dexp(-arg)
         if (model.eq.'ps'.or.model.eq.'eikonal') then 
            deltasv = 0.d0
         endif
         sv = ffs * (1.d0 + deltasv)
         dsv = deltasv         
         return
      endif

      if (ord.eq.'alpha2') then
         if (fs.ne.'ee') then
            print*,'alpha2 only for bhabaha! stopping...'
            stop
         endif
         
         aieps = -2.d0*dlog(eps) -1.5d0+2.d0*eps-0.5d0*eps**2
         s = (pin1b(0) + pin2b(0))**2
         p1mod = sqrt(tridot(p1b,p1b))
         c = p1b(3)/p1mod

         t = -s/2.d0*(1.d0 - c)
         u = -s/2.d0*(1.d0 + c)
c         scale = s*t/u
c         lcoll = dlog(scale/ame**2) - 1.d0
c         arg = 2.d0*alpha/pi * lcoll * aieps

         sdif = bornme*phsp2b/flux2b

         svpureoal = aintsv(s,t,eps)/sdif - 1.d0
c        svpureoal = aintsv_new(pin1b,pin2b,p1b,p2b,eps)/sdif - 1.d0
cc         xxx = aintsv(s,t,eps)
c     deltasv = svpureoal + arg      !!! NB
*** NEW SFF, formula 2.23 PhD. Thesis. in the limit dot(pi,pj) >> m**2
c         dp1p2 = dlog(2.d0*dot(pin1,pin2)/ame**2)
c         dp1p3 = dlog(2.d0*dot(pin1,p1)/ame**2)
c         dp1p4 = dlog(2.d0*dot(pin1,p2)/ame**2)
c         dp2p3 = dlog(2.d0*dot(pin2,p1)/ame**2)
c         dp2p4 = dlog(2.d0*dot(pin2,p2)/ame**2)
c         dp3p4 = dlog(2.d0*dot(p1,p2)/ame**2)
c         lcollnew = dp1p2 + dp1p3 - dp1p4 - dp2p3 + dp2p4 + dp3p4
c         lcollnew = 2.d0 * lcollnew - 4.d0
c         lcollnew = lcollnew/4.d0
c         arg = 2.d0*alpha/pi * lcollnew * aieps
**** NEW ***************
         charges(1) = -1.d0
         charges(2) = -1.d0
         charges(3) = -1.d0
         charges(4) = -1.d0
         masses(1)  = ame
         masses(2)  = ame
         masses(3)  = ame
         masses(4)  = ame
         if (fs.eq.'mm') then
            masses(3)  = ammu
            masses(4)  = ammu
         endif
         lcoll1 = eikonalintegral(pin1b,pin2b,p1b,p2b,charges,masses)
         arg2   = 2.d0*alpha/pi * lcoll1 * aieps
         arg = arg2         
         deltasv = svpureoal + arg !!! NB
!     I used arg with reduced momenta to define deltasv, now I recalculate it
! with true external momenta...
c         lcoll1 = eikonalintegral(pin1,pin2,p1,p2,charges,masses)
c         arg   = 2.d0*alpha/pi * lcoll1 * aieps
**************************
c         ffs = dexp(-arg)
         if (model.eq.'ps') deltasv = 0.d0

         if (ng.le.2) then
            sv = 1.d0
         endif
         if (ng.le.1) then
            sv = sv - arg + deltasv
         endif
         if (ng.eq.0) then
            sv = sv + arg**2*0.5d0 - arg*deltasv
         endif
         return
      endif
      return
      end
******
* COMPLEX VACUUM POLARIZATION, T. Teubner & al.
      FUNCTION VPOLChlmnt_BACKUP(Q2)
      implicit double precision (A-H,L,O-Z)
      parameter (np  = 7186)
c      parameter (np=6883)
      parameter (nw  = 1582)

c      parameter (npt = np) ! or
c      parameter (npt = nw) ! or
      parameter (npt = nw+np)

      dimension s(0:np-1)
      dimension st(0:npt-1)
      dimension pres(0:np-1),pre(0:np-1),pim(0:np-1)
      dimension epres(0:np-1),epre(0:np-1),epim(0:np-1),cv(0:np-1)

      dimension iswvperr(3),iswvperrt(6)
      common/iswitchvperror/iswvperr,iswvperrt

      double complex vpolcnsk,vpolchlmnt,deltahad,im,res
      double complex cP,ressp,restm,ctmp
      external vpolcnsk
      common/vpcnsk/s,pres,pre,pim,epres,epre,epim,cv,im,ifirstnsk

      common/rescaledhaderr/scal_err_dhad

      parameter (nres=15)
      double precision AM(nres),GT(nres),GE(nres)
      common/resonancesparams/am,gt,ge
      double precision EAM(nres),EGT(nres),EGE(nres)
      common/resonancesparamserrors/eam,egt,ege

      dimension iresbh(8) ! the indexes of the resonances to be added by hand (Thomas said to include **only** 1-6)
      dimension prest(0:npt-1),pret(0:npt-1),pimt(0:npt-1)
      dimension eprest(0:npt-1),epret(0:npt-1),epimt(0:npt-1)
      common/vpcteubBCK/st,prest,pret,pimt,eprest,epret,epimt,iresbh,
     .     ifirst,imidguess
      data ifirst/0/

      double precision wt(nw)
      common/wtabulation/wt

      double precision wtot(npt)
      integer indexes(npt)
      common/wtotcommon/wtot
      
      common/tabulationlimits/wmin,wmax,ic
      
      common/parameters/ame,ammu,convfac,alpha,pi
      common/ecms/ecms,ecmsnom

      integer iloc(1)
      
*     using vp_hlmnt_v2_1.f, sent by T. Teubner, private communication
      
      if (ifirst.eq.0) then
***   filling the common by calling vpolcnsk once (I need the array s --> st[eubner])
         vpolchlmnt = vpolcnsk(q2)
****  filling the common wtabulation        
         wtot = 0.d0
         call fillwtabcommon
         wmin = wt(1)
         wmax = wt(nw)
         wtot(1:nw)     = wt**2
         wtot(nw+1:npt) = s
         call sort_vector_inverse(wtot,wtot,indexes,npt)
c         do k = 1,npt
c            print*,k,sqrt(wtot(k)),indexes(k)
c         enddo
c         stop
****  
c
c I solve for P, so that is the same quantity returned by NSK
c     I also use the same evaluation points as NSK, and then linear
c     interpolation between them
c     1/(1-P) = 1/(1-danr) + res ->
c     1-P = (1-danr)/(1+(1-danr)*res)
c       P = 1 - (1-danr)/(1+(1-danr)*res)
c         = (1+(1-danr)*res - (1-danr))/(1+(1-danr)*res)
c         = ( danr + res - danr*res) / (1 + res- danr*res)

******************
* filling st, from the nsk data or Teubner tabulation points
c         st(0:np-1) = s   ! or
c         do k = 1,nw
c            st(k-1) = wt(k)*wt(k)
c         enddo ! or
         st = wtot
******************
         prest  = 0.d0
         pret   = 0.d0
         pimt   = 0.d0
         eprest = 0.d0
         epret  = 0.d0
         epimt  = 0.d0
         iresbh(1) = 4
         iresbh(2) = 5
         iresbh(3) = 10
         iresbh(4) = 11
         iresbh(5) = 12
         iresbh(6) = 13
         iresbh(7) = 14
         iresbh(8) = 15
         print*,'Evaluating Teubner et al. VP array '
         open(10,file='cfrvp-0-v2-3',status='unknown')

         do k = 0,npt-1
            energy = sqrt(st(k))

c            hlmnt is tabulated up to wmax and from wmin
            if (energy.le.wmax.and.energy.ge.wmin) then
               ressp = 0.d0
               restm = 0.d0
               do i=1,6
                  j = iresbh(i)
                  AMC = AM(J) !- EAM(J)
                  GEC = GE(J) !- EGE(J)
                  GTC = GT(J) !- EGT(J)
                  V   = 3.D0*GEC/AMC / ALPHA
                  restm = restm + 
     .                 V * st(k) / ( st(k)-AMC*AMC+im*GTC*AMC )
                  ressp = ressp - 
     .                 V * st(k) / (-st(k)-AMC*AMC+im*GTC*AMC )
               enddo

               print*,'HERE, CHECK WHATs done below...'
               
               nrflag = 1
c               print*,'attenzione'
               call vphlmntv2(energy,
     &              vprehadsp,vprehadtm,vpimhad,vprelepsp,
     &              vpreleptm,vpimlep,
     &              vpretopsp,vpretoptm,nrflag,evprehadsp,
     &              evprehadtm,evpimhad)

c               vprehadsp = vprehadsp - evprehadsp
c               vprehadtm = vprehadtm - evprehadtm
c               vpimhad   = vpimhad -  evpimhad
!     ! space like
               prest(k) = vprehadsp + vprelepsp + vpretopsp
               prest(k) = prest(k)  ! + ressp ! in the space-like always resummed!
c               prest(k) = vprehadsp + ressp
!     ! time like
               nrflag = 0
               call vphlmntv2(energy, 
     &        vprehadsp,vprehadtm,vpimhad,vprelepsp,vpreleptm,vpimlep,
     &        vpretopsp,vpretoptm,nrflag,evprehadsp,evprehadtm,evpimhad)

               pret(k)  = vprehadtm + vpreleptm + vpretoptm
               pimt(k)  = vpimhad   + vpimlep

c   this is 1/(1-P)
               ctmp = 1.d0/(1.d0 - pret(k)-im*pimt(k)) + restm ! in the time-like not resummed!
c   ctmp = 1/(1-P) -> 1-P = 1/ctmp  -> P = 1-1/ctmp
               pret(k) = DReal(1.d0-1.d0/ctmp)
               pimt(k) = DImag(1.d0-1.d0/ctmp)

c               write(10,*)energy,pret(k),pimt(k),pre(k),pim(k),
c     .              prest(k),pres(k)
               write(10,*)energy,pret(k),pimt(k),prest(k)
               
            endif
         enddo
         close(10)
         print*,'done!'
         print*,'stopping line ~420 sv.f'
         stop
*** guessing where to start the binary search....         
         k = 0
         aq2 = ecmsnom**2
         do while(aq2.gt.st(k))
            k = k + 1         
         enddo
         k = k - 1
         imidguess = k ! this is correct for born s-channel and no energy spread
***************
         ic = 0
         ifirst=1
      endif

c      k   = 0
c      sk  = 0.d0
c      aq2 = abs(q2)
c      do while(aq2.gt.st(k))
c         k = k + 1         
c      enddo
c      k = k - 1
******** binary search ****
      aq2     = abs(q2)
      istart  = 0
      ifinish = npt-1
c      imid = (istart + ifinish)/2
      imid = imidguess

c     do while(ifinish.ge.istart)
      do while(.not.(st(imid+1).gt.aq2.and.st(imid).le.aq2))
         if (aq2.gt.st(imid)) then
            istart  = imid + 1
         else
            ifinish = imid - 1
         endif
         imid = (istart + ifinish)/2
      enddo
      k = imid
******************
      
      if (q2.gt.0.d0) then
         vpre = pret(k)+(pret(k+1)-pret(k))/(st(k+1)-st(k)) * (q2-st(k))
         vpim = pimt(k)+(pimt(k+1)-pimt(k))/(st(k+1)-st(k)) * (q2-st(k))

c         evpre = epret(k)+ (epret(k+1)-epret(k))/(st(k+1)-st(k))*(q2-st(k))
c         evpim = epimt(k)+ (epimt(k+1)-epimt(k))/(st(k+1)-st(k))*(q2-st(k))
      else
         vpre = prest(k) +
     .        (prest(k+1)-prest(k))/(st(k+1)-st(k))*(aq2-st(k))
         vpim = 0.d0

c         evpre=eprest(k)
c     .        +(eprest(k+1)-eprest(k))/(st(k+1)-st(k))*(aq2-st(k))
c         evpim = 0.d0
      endif
*
      VPOLChlmnt = 1.d0/(1.d0  - vpre - im*vpim)
      return
      end
