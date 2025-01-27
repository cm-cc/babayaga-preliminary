      program combinefiles
      implicit none
      character*50 input1,input2,input3,output
      character*50 distrib,stringadd
      integer len1
      integer k,ifile,nfiles,ios,nfilesdim
      integer maxlines
      parameter (maxlines=400,nfiles=80,nfilesdim=80)
      character*(100) files(nfilesdim)
      character*(100) line(maxlines,nfilesdim)
      integer nlines(maxlines)
      real * 8 v1,v2,v3,y,err
      integer ilength
      external ilength
      character*1 slash
      integer i,j
      character*20 str
      external str

      slash= '/'

      
      do j=1,2
      if(j.eq.1) distrib= 'cmd_pp_thav_EXP.txt'
      if(j.eq.2) distrib= 'cmd_pp_mxx_EXP.txt'

         output=trim(distrib)
*
         open(unit=14,file=output,status='new')
*
         do ifile=1,(nfiles)
            files(ifile)= trim(str(ifile))//trim(slash)//trim(distrib)
            print*,'files(ifile)= ',files(ifile)
            open(unit=11,file=files(ifile),status='old')
            do k=1,maxlines+1
               read(unit=11,fmt='(a)',end=111) line(k,ifile)
               if(k.eq.maxlines+1) then
                write(*,*) ' too many lines in file, increase maxlines'
                call exit(-1)
               endif
               goto 12
 111           nlines(ifile)=k-1
               goto 11
 12            continue
            enddo
 11         continue
            print*,trim(str(ifile)),nlines(ifile)
         enddo
         do ifile=1,(nfiles)
            if(nlines(ifile).ne.nlines(1)) then
               write(*,*) ' error: file', files(ifile),
     1              ' does not match in length'
               call exit(-1)
            endif
         enddo
         
         y= 0.d0
         err= 0.d0
         do k=1,nlines(1)
            print*,'k= ',k
            read(unit=line(k,1),fmt=*,iostat=ios) v1,v2,v3
            if(ios.ne.0) then
               write(14,'(a)') line(k,1)(1:ilength(line(k,1)))
            else
               y=v2
               err=v3**2
               do ifile=2,(nfiles)
                  read(unit=line(k,ifile),fmt=*,iostat=ios) v1,v2,v3
                  y=y+v2
                  err=err+v3**2
               enddo
               err=sqrt(err)
               write(14,'(3(1x,e14.8))')
     +              v1,y/(float(nfiles)),err/(float(nfiles))
            endif
         enddo
c      close(11)
c      close(12)
c      close(13)
         close(14)
      enddo
      stop
      end

      function ilength(line)
      integer ilength
      character *(*) line
      ilength=len(line)
      do j=ilength,1,-1
         if(line(j:j).ne.' ') then
            ilength=j
            return
         endif
      enddo
      ilength=0
      end

      character(len=20) function str(k)
!      "Convert an integer to string."
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
      end function str
