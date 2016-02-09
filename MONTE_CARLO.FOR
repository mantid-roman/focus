! MODULE MONTE_CARLO.FOR
! Monte Carlo on free parameter

!--------------------------------------------------------------------------  
! Monte Carlo on all free and limited parameters
!--------------------------------------------------------------------------  
      subroutine random_sel
        implicit none

        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:XR.INC'
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:LIMITPAR.INC'

        integer is_varied

        integer i,j,k,l
        real*8 x1,x2,random

c         x1             x2
c      ---------      ---------
c      ! x !   !      ! x !   !  => change component x := x1+x2 e [0,1]
c      ---------      ---------
c
c      ---------      ---------
c      ! x !   !      !   ! x !  => do not change component
c      ---------      ---------
c
c      ---------      ---------
c      !   ! x !      !x  !   !  => do not change component
c      ---------      ---------
c
c      ---------      ---------
c      !   ! x !      !   ! x !  => change component x:= x1+x2-2 e [-1,0]
c      ---------      ---------
c      0  0.5  1      0  0.5  1
c
c
      do i=1,np         ! these are the number of free parameters
           j = inp(i)   ! here they are coded       
           if(pmax(j).gt.pmin(j)) then  ! vary only those which have limits set
                x1 = random()
                x2 = random()
                is_varied = 0
                if(x1.le.0.5d0.and.x2.le.0.5d0) then
                    is_varied = 1
                    p(j) = x1+x2       ! e [0,1]
                else if(x1.gt.0.5d0.and.x2.gt.0.5d0) then
                    is_varied = 1
                    p(j) = x1+x2-2.0d0 ! e [-1,0]
                end if
                if(is_varied.eq.1) then
                    p(j) = 0.5*(up_limit+low_limit)*p(j)+0.5*(up_limit-low_limit)
                end if
           end if
      end do
      return
      end 
!-------------------------------------------------------------------------------  
!create a random number between 0 and 1
!-------------------------------------------------------------------------------  
      real*8 function random()
         implicit none
         include 'CF_SOURCES:RNUMBER.INC'

         real*8 ran1
         random = ran1(r0_number)
      return
      end
c-------------------------------------------------------------------------------
c random number generator between 0 and 1
c-------------------------------------------------------------------------------
      double precision function RAN1(idum)
      integer idum, ia, im, iq, ir, ntab, ndiv
      real am, eps, rnmx
      parameter (ia=16807, im=2147483647, am=1.0/im, iq=127773, ir=2836,
     + ntab=32, ndiv=1+(im-1)/ntab, eps=1.2e-7, rnmx=1.0-eps)
C
C [Should be] identical to 2nd ed. routine
C Only the final result is double precision
C
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      if (idum .le. 0 .or. iy .eq. 0) then
        idum=max(-idum,1)
        do 11 j=ntab+8,1,-1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          if (idum .lt. 0) idum=idum+im
          if (j .le. ntab) iv(j)=idum
 11     continue
        iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum .lt. 0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      ran1=dble(min(am*iy,rnmx))
      return
      end
c---------------------------------------------------------------------------
c  store the best MAX_CHI2-th parameter
c---------------------------------------------------------------------------
      subroutine mc_store
        implicit none
c       --------------------------------------------------------------------        
	include 'CF_SOURCES:FRILLS.INC'
 	include 'CF_SOURCES:FOCUS.PAR'
 	include 'CF_SOURCES:MONTE.INC'
 	include 'CF_SOURCES:OVERALL.INC'
c       ---------------------------------------------------------------------        
	external M01DAF   ! NAG - sorts real numbers 
	external M01ZAF 
	external M01CAF
	integer ifail,ind(MAX_CHI2)
        real*8  w_par(MAX_CHI2,MAX_PAR), w_chi2(MAX_CHI2), w_ratio(MAX_CHI2)
        real*8  w_min(MAX_CHI2,MAX_PAR), w_max(MAX_CHI2,MAX_PAR)
c       --------------------------------------------------------
        integer i, j, k, l, n, is_stored, is_first
c       --------------------------------------------------------------------        
c       find out how many are already stored
c       --------------------------------------------------------------------        
        n=0
        do i=1,MAX_CHI2
           call ifnull(s_chisq(i))
           if(s_chisq(i).ne.0.0d0) n=n+1
        end do
c       --------------------------------------------------------------------        
c       order by ascending Chi**2 values
c       --------------------------------------------------------------------        
        if(n.eq.0) then
           s_chisq(1) = chisq
           s_ratio(1) = ratioenergy
           do i=1,MAX_PAR
              s_p   (1,i) = p(i)
              s_pmin(1,i) = pmin(i)
              s_pmax(1,i) = pmax(i)
           end do
        else
c         -------------------------------------
c         store old values
c         -------------------------------------
          do i=1,n
             w_chi2(i)  = s_chisq(i)
             w_ratio(i) = s_ratio(i)
             do j=1,MAX_PAR
                w_par(i,j) = s_p(i,j)
                w_min(i,j) = s_pmin(i,j)
                w_max(i,j) = s_pmax(i,j) 
             end do
          end do
c         ---------------------------------------
c         sort Chi2**2 values in ascending order
c         ---------------------------------------
          ifail=1
	  call M01DAF(w_chi2,1,n,'A',ind,ifail)
	  call M01ZAF(ind   ,1,n        ,ifail)
	  call M01CAF(w_chi2,1,n,'A'    ,ifail)      
c         ---------------------------------------
c         store 
c         ---------------------------------------
          k=0
          is_first = 0
          do i=1,n
             j=ind(i)
             if(chisq.lt.s_chisq(j).and. is_first.eq.0) then
                  k = k+1
                  if(k.le.MAX_CHI2) then
                     s_chisq(k) = chisq
                     s_ratio(k) = ratioenergy
                     do l=1,nptot
                        s_p   (k,l) = p(l)
                        s_pmin(k,l) = pmin(l)
                        s_pmax(k,l) = pmax(l)
                     end do
                  end if
                  is_first = 1
             end if
             k = k+1
             if(k.le.MAX_CHI2) then 
                s_chisq(k) = w_chi2(i)
                s_ratio(k) = w_ratio(j) ! because w_ratio is not sorted 
                do l=1,nptot
                   s_p   (k,l) = w_par(j,l)
                   s_pmin(k,l) = w_min(j,l)
                   s_pmax(k,l) = w_max(j,l)
                end do
             end if
          end do
          if(is_first.eq.0) then
               is_first = 1
               k = k+1
               if(k.le.MAX_CHI2) then
                  s_chisq(k) = chisq
                  s_ratio(k) = ratioenergy
                  do l=1,nptot
                     s_p   (k,l) = p(l)
                     s_pmin(k,l) = pmin(l)
                     s_pmax(k,l) = pmax(l)
                  end do
               end if
          end if

       end if

      return
      end
c---------------------------------------------------------------------------
c perform monte_carlo
c---------------------------------------------------------------------------
      subroutine monte_carlo (fun, line, l)
        implicit none

	include 'CF_SOURCES:FRILLS.INC'
 	include 'CF_SOURCES:FOCUS.INC'
 	include 'CF_SOURCES:MONTE.INC'
 	include 'CF_SOURCES:XR.INC'
 	include 'CF_SOURCES:CONTROL.INC'
 	include 'CF_SOURCES:OVERALL.INC'
 	include 'CF_SOURCES:DEGENERATION.INC'

        external fun   

        real*8 mc_chisq, mc_ratio, mc_p(MAX_PAR)
        character step*6, chi2*8, c, ic*3
        integer o_setbkq, o_setvkq, o_setxr
        integer n0, no_0_par
        integer sitime(7)
        integer eitime(7)
        external x05aaf  ! found in NAG library
c       ---------------------------------------------------------------------
	character word*12, line*80
	logical fitted
	integer*4 i, line_len, getlin, l, word_len, getwrd, j, ctoi, n
        real*4    ctor
        real*8    chi2max, chisq_old, po(MAX_PAR),pmino(MAX_PAR),pmaxo(MAX_PAR),ratio_old              
	character chint*8, cpmin*13, cpmax*13



        real*8  pminl, pmaxl
        integer flag, same, dim
        external same

             n0 = no_0_par()     

             n=0
             do i=1,np
                j = inp(i)
                if(pmax(j).gt.pmin(j)) n=n+1
             end do
                 
             if(np.lt.1) then
                  call remark ('ERROR: Can not carry out Monte-Carlo, because there are no FREE parameter.')
                  call remark ('       Free at least one parameter. Use the MULTI_FRILLS command C P [no.].')
             else if(n.eq.0) then
                  call remark ('ERROR: Can not carry out Monte-Carlo, because there are no free parameters')
                  call remark ('       with finite limits. Give limits to at least one free parameter.')
                  call remark ('       Remeber that the MULTI_FRILLS command LIMIT supports this task.')
             else
c               --------------------------------
c               look if there are no zero errors
c               --------------------------------
 	        do i = 1,ns
 	       	     do j = 1,nv(i)
 	       		  if (vsig(inv(j,i),i) .le. 0.0) then
 				call remark ('ERROR: Zero error in point '//chint(inv(j,i)))
 				return
 			  end if
 	             end do
 	        end do

c               --------------------------------
c               ask for input
c               --------------------------------
                call remark(' ')
  901           call prompt(' Enter the number of Monte-Carlo steps [100=default]:')
                line_len = getlin(line)
                if(line_len .gt. 0) then
                   l = 1
                   n = abs(ctoi(line(1:line_len),l))
                   if(n.eq.0) then 
                      call remark('ERROR: The number of Monte-Carlo steps must be greater than 0.')
                      goto 901
                   end if
                   if(n.gt.999999) then
                      call remark ('ERROR: The number of Monte-Carlo steps must be lower than 1e6.')
                      goto 901
                   end if
                else
                   n = 100
                end if

  902           call prompt(' Enter the maximum value of Chi**2     [100=default]:')
                line_len = getlin(line)
                if(line_len .gt. 0) then
                   l = 1
                   chi2max = ctor(line(1:line_len),l)
                   call ifnull(chi2max)
                   if(chi2max .eq. 0.0d0) then
                      call remark('ERROR: Please set maximum value of Chi**2 greater than 0.')
                      goto 902
                   end if
                else 
                   chi2max = 100.0d0
                end if

c               ----------------------------------------
c               relimit the R/W parameter if not so in the 
c               case of a fixed crystal field level
c               ----------------------------------------
                if( setoverall.eq.1 .and. (setxr.eq.1 .or. setwx.eq.1) ) then
                   flag=0
                   n0 = no_0_par()
                   do i=1,np
                      j = inp(i)
                      if( (setxr.eq.1.and.j.eq.indvold(1)).or.(setwx.eq.1.and.j.eq.n0+1) ) flag=1
                   end do
                   if(flag.eq.0) then 
                      np = np + 1
                      if(setxr.eq.1) inp(np) = indvold(1)
                      if(setwx.eq.1) inp(np) = n0+1
                   end if
                end if

                if(setoverall.eq.1) then
                   n0 = no_0_par()
                   do i=1,np
                      j = inp(i)
                      if( (setxr.eq.1.and.j.eq.indvold(1)).or.(setwx.eq.1.and.j.eq.n0+1) )then
                        p(j)    = 0
                        if(pmax(j).eq.pmin(j)) then
                             pmin(j) = -1.0d0
                             pmax(j) =  1.0d0
                        else
                             if(pmin(j).lt.0.0d0) pmin(j) = pmin(j)/abs(pmin(j))
                             if(pmax(j).gt.0.0d0) pmax(j) = pmax(j)/abs(pmax(j))

                             if(pmin(j).lt.0 .and. pmax(j).lt.0) pmax(j)=0.0d0
                             if(pmin(j).gt.0 .and. pmax(j).gt.0) pmin(j)=0.0d0
                         end if
                        if(setxr.eq.1) call remark (' FOCUS relimits the R-parameter due to a fixed CEF level.')
                        if(setwx.eq.1) call remark (' FOCUS relimits the W-parameter due to a fixed CEF level.')
                      end if
                   end do
                end if
c               --------------------------------------
c               prints out which components are varied
c               --------------------------------------
                call remark (' The following parameters will be varied by Monte-Carlo:')
  100           format(1x,a20,a3,a13,a3,a13,a3)
                n0     = no_0_par()
                flag = 0
                do i=1,np                      
                   l = inp(i)
                   if( pmax(l).gt. pmin(l)) then
                      write(cpmin, fmt='(g13.5)' ) pmin(l)
                      write(cpmax, fmt='(g13.5)' ) pmax(l)
!                     ---------------------------------------
!                     print out the new boundaries
!                     ---------------------------------------
                      if(setoverall.eq.1 .and. l.ge.n0-5) then
                         if(setxr.eq.1 .and. same(l,indvold(1)).eq.1) then
                             if(l.ne.indvold(1)) then
                                call remark (' '//nam(l)//' [ '//cpmin//' , '//cpmax//' ] ')  !x-par
                             else 
                                call remark (' '//nam(l)//' [ '//cpmin//' , '//cpmax//' ] * energy scaling factor')
                                flag=1
                             end if
                         else if(setwx.eq.1 .and. l.gt.n0) then
                             if(l.eq.n0+2) then 
                                call remark (' '//nam(l)//' [ '//cpmin//' , '//cpmax//' ] ')  !x-par
                             else 
                                call remark (' '//nam(l)//' [ '//cpmin//' , '//cpmax//' ] * energy scaling factor')
                                flag=1
                             end if
                         else
                             call remark (' '//nam(l)//' [ '//cpmin//','//cpmax//' ] * energy scaling factor')
                             flag=1
                         end if
                      else
                         call remark (' '//nam(l)//' [ '//cpmin//' , '//cpmax//' ] ')
                      end if
                   end if
                end do       
                if(setoverall.eq.1 .and. flag.ne.0) then 
                   call remark (' The  energy scaling factor  ensures that the fixed  CEF level  stays fixed') 
                   call remark (' during the Monte Carlo. To free the fixed CEF level use the command @ A O.') 
                end if
c               --------------------------------
c               Monte carlo start
c               --------------------------------
                do i=1,MAX_CHI2
                   s_chisq(i) = 0.0d0
                end do

                call x05aaf(sitime)
                call remark (' Monte-Carlo starts on '//ic(sitime(4))//':'//ic(sitime(5))//':'//ic(sitime(6)) )
  
                call prtin

c               --------------------------------
c               store old situation
c               --------------------------------
               	do i = 1,nptot  !Store old parameters 
	           po(i)    = p(i)
                   pmino(i) = pmin(i)
                   pmaxo(i) = pmax(i)
	        end do 
                ratio_old = ratioenergy              
                chisq_old = chisq

                mc_chisq = 1000000.0d0
                if_monte = 0
                write(*,*)
  910           format('+',x,x,a10,i6)
  911           format('+',x,x,a10,a6)
                do j=1,n
                  write(*,910) 'mc_step = ',j

!                 --------------------------------------------------------------------------
!                 1. transform the limit parameters pj e [pminj,pmaxj] to angles [0,pi] 
!                 2. random selection of the limit parameters
!                 3. calculate chi**2
!                 4. back transformation of parameters
!                 --------------------------------------------------------------------------
  799      	  call ptrans (np, p, inp, pmin, pmax) ! Transform all limited parameters into circular coordinates
                  call random_sel                          ! produce a random parameters
                  setnowarnings=1 ! give no warnings out
                  setmonte=1
                  call chisq_calc (fun)                    ! calculate initial Chi**2 
                  setmonte=0
                  setnowarnings=0 ! give warnings out
        	  call pback (np, nptot, p, inp, pmin, pmax, inr, rp) ! back transformation of ptrans
!                 --------------------------------------------------------------------------
!                 make sure that the total splitting is not zero
!                 --------------------------------------------------------------------------
                  dim = dimj
                  if(energy(dim).le.de) goto 799

                  if( chisq .le. chi2max ) then
                       write(*,911) '          ','      '
                       write(step, fmt = '(i6)'  ) j
                       write(chi2, fmt = '(f8.1)') chisq
                       call remark ( step//'. Monte-Carlo Step    Chi**2 = '//chi2)
                       write(*,*) '     '
c                      -------------------------------------------------------------
c                      store best parameter situation
c                      -------------------------------------------------------------
                       if( chisq .lt. mc_chisq ) then
                            if_monte = 1
                            mc_chisq = chisq
                            mc_ratio = ratioenergy
                            do i=1,nptot
                               mc_p(i) = p(i)
                            end do
                       end if
c                      -------------------------------------------------------------
c                      store parameter
c                      -------------------------------------------------------------
                       call mc_store
                  end if
                end do
 	        fitted = .TRUE.  

                call x05aaf(eitime)
                call remark (' Monte-Carlo starts on '//ic(sitime(4))//':'//ic(sitime(5))//':'//ic(sitime(6)) )
                call remark (' Monte-Carlo ends   on '//ic(eitime(4))//':'//ic(eitime(5))//':'//ic(eitime(6)) )
c               --------------------------------
c               Monte carlo end
c               --------------------------------
c               --------------------------------------------------------------------
c               restore the best MC-parameter or the old situation before Monte Carlo
c               --------------------------------------------------------------------
                if(if_monte .eq. 1) then
                    call remark (' ')
                    call prompt (' Do you wish to store the best Monte-Carlo parameters [Y]? ')
                    line_len = getlin(line)
                    c = 'y'
                    if(line_len .ne. 0) c = line(1:1)

                    if(c.eq.'N'.or.c.eq.'n') then
                       goto 903                   !Restore the old parameters
                    else
                       chisq       = mc_chisq     !Restore chisq
                       ratioenergy = mc_ratio     !Restore ratio important if setoverall=1
                       do i = 1,nptot             !Restore best Monte Carlo parameters 
 	                   p(i) = mc_p(i)
	               end do               
                    end if
                else
                    write(chi2, fmt = '(f8.1)') chi2max
                    call remark (' Sorry, could not find a parameter set with Chi**2 lower than or equal to '//chi2//'.')
  903               chisq       = chisq_old
                    ratioenergy = ratio_old
                    do i = 1,nptot  !Restore old parameters 
      	                p(i)    = po(i)
                        pmin(i) = pmino(i)
                        pmax(i) = pmaxo(i)
	            end do
                end if

             end if

      return
      end
c---------------------------------------------------------------------------
c display monte carlo CHI**2
c---------------------------------------------------------------------------
      subroutine dis_monte (fun, line, l)
        implicit none

	include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.PAR'
        include 'CF_SOURCES:MONTE.INC'

        external fun   

	character line*80
	integer*4 line_len, getlin,l, word_len, getwrd, ctoi
        real*4    ctor
c       ---------------------------------------------------------------------        
        real*8  chi2max
        integer i, j, k, n, is_displayed
        character*8 chint
c       ---------------------------------------------------------------------        
c       find out how many are already stored
c       --------------------------------------------------------------------        
        n=0
        do i=1,MAX_CHI2
           call ifnull(s_chisq(i))
           if(s_chisq(i).ne.0.0d0) n=n+1
        end do
c       --------------------------------------------------------------------        
c       ask for the maximun Chi**2 
c       --------------------------------------------------------------------        
 902    call prompt(' Enter the maximum value of Chi**2     [100=default]:')
        line_len = getlin(line)
        if(line_len .gt. 0) then
           l = 1
           chi2max = ctor(line(1:line_len),l)
           call ifnull(chi2max)
           if(chi2max .eq. 0.0d0) then
                call remark('ERROR: Please set maximum value of Chi**2 greater than 0.')
                goto 902
           end if
        else 
               chi2max = 100.0d0
        end if
c       --------------------------------------------------------------------        
c       display Chi**2 <= chi2max
c       --------------------------------------------------------------------        
        is_displayed = 0
        do i=1,n 
           if(s_chisq(i).le.chi2max) then
                is_displayed = 1
                write(*,100) i,'. stored Monte-Carlo Chi**2 is ', s_chisq(i)
           end if
        end do
  100   format(x,x,i6,a31,f8.1)
        if(is_displayed.eq.0) then
           write(chint, fmt = '(f8.1)') chi2max
           call remark ('Sorry, no Chi**2 found which are lower than or equal to '//chint)
        end if
      return 
      end 
c---------------------------------------------------------------------------
c performs a load on selected parameters
c---------------------------------------------------------------------------
      subroutine load_monte (fun, line, l)
        implicit none

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:FOCUS.PAR'
	include 'CF_SOURCES:MONTE.INC'
	include 'CF_SOURCES:CONTROL.INC'
	include 'CF_SOURCES:OVERALL.INC'

        external fun   

	character word*12, line*80, cint*10
	integer*4 i, line_len, getlin, l, word_len, getwrd, j, ctoi, n, nchi2
        integer*4 par_len, itoc, getari
	character*5 par_no
        real*4    ctor
        logical fitted
c       ---------------------------------------------------------------------        
c       find out how many are already stored
c       --------------------------------------------------------------------        
        n=0
        do i=1,MAX_CHI2
           call ifnull(s_chisq(i))
           if(s_chisq(i).ne.0.0d0) n=n+1
        end do
        if(n.eq.0) then
           call remark ('Sorry, there are no Chi**2 defined. Please try again @ M N.')
           return
        end if
c       --------------------------------------------------------------------        
c       ask for the Monte Carlo number 
c       --------------------------------------------------------------------        
        if(nset.eq.0) then
           if(if_plot.eq.0) then
              call remark ('Enter below the number of the Monte-Carlo-set you want to load.')
              call remark ('REMARK: The number of the Monte-Carlo-sets are displayed by @ m d.')
 10           call prompt ('Which set do you want to load [1=default]? ')
           else
              call remark ('Enter below the number of the Monte-Carlo-set you want to plot.')
              call remark ('REMARK: The number of the Monte-Carlo-sets are displayed by @ m d.')
 20           call prompt ('Which set do you want to plot [1=default]? ')
           end if
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l   = 1
              nchi2 = ctoi(line(1:line_len),l)
              if(nchi2.lt.1.or.nchi2.gt.n) then
                   write(cint,fmt='(i6)') n+1
                   call remark ('ERROR: The number of the Monte-Carlo-set must be greater than 1 ')
                   call remark ('       and lower than '//cint)
                   if(if_plot.eq.0) then
                      goto 10
                   else
                      goto 20
                   end if 
              end if
           else
               nchi2 = 1
           end if
           nset = 1
           n_mc_par = nchi2
        else
           if(if_plot.eq.0) n_mc_par = imc_array(1)
        end if

c       -------------------------
c       store specified parameter
c       -------------------------
        chisq       = s_chisq(n_mc_par)
        ratioenergy = s_ratio(n_mc_par)
        do i=1,MAX_PAR
           p(i)    = s_p   (n_mc_par,i)
           pmin(i) = s_pmin(n_mc_par,i)
           pmax(i) = s_pmax(n_mc_par,i)
        end do

      return 
      end 
c---------------------------------------------------------------------------
c plot monte carlo parameters
c---------------------------------------------------------------------------
      subroutine plot_monte (fun, line, l)
        implicit none

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:FOCUS.PAR'
	include 'CF_SOURCES:MONTE.INC'
	include 'CF_SOURCES:CONTROL.INC'
	include 'CF_SOURCES:OVERALL.INC'

        external fun   

	character word*12, line*80
	integer*4 i, line_len, getlin, l, word_len, getwrd, j, ctoi, n, ii
        real*4    ctor
c       ---------------------------------------------------------------------        
        logical fitted
        integer np_old, inp_old(MAX_PAR), ratio_old 
        real*8  p_old(MAX_PAR), p_min(MAX_PAR), p_max(MAX_PAR)
        character*8 cminx, cmaxx, cminy, cmaxy
        character*13 cpmin, cpmax
        character cint*6, c
        integer*4 nii
c-----------------------------------------
c store old situation
c-----------------------------------------
        np_old = np

        do i=1,np
           inp_old(i) = inp(i)
        end do

        do i=1,nptot
           p_old(i) = p(i)
           p_min(i) = pmin(i)
           p_max(i) = pmax(i)
        end do
        ratio_old = ratioenergy
c-----------------------------------------
c mark the plot mode
c-----------------------------------------
        if_plot = 1
c-----------------------------------------
c load parameters in  plot mode
c-----------------------------------------
        nii = nset
        if(nset.eq.0) nii=1
        
        do ii=1, nii
                n_mc_par = imc_array(ii) 
        	call load_monte (fun, line, l)
c-----------------------------------------
c ask for plot boundaries
c-----------------------------------------
        	write(cminx, fmt='(f8.1)' ) min_x
        	write(cmaxx, fmt='(f8.1)' ) max_x
        	write(cminy, fmt='(f8.1)' ) min_y
        	write(cmaxy, fmt='(f8.1)' ) max_y

        	if(if_first.eq.0 .or. nset.eq.1) then
                   call prompt(' Enter min energy       ['//cminx//']:')
        	   line_len = getlin(line)
        	   if(line_len .gt. 0) then
           		l = 1
           		min_x = ctor(line(1:line_len),l)
        	   end if

                   call prompt(' Enter max energy       ['//cmaxx//']:')
        	   line_len = getlin(line)
           	   if(line_len .gt. 0) then
           		l = 1
           		max_x = ctor(line(1:line_len),l)
        	   end if

                   call prompt(' Enter min intensity    ['//cminy//']:')
        	   line_len = getlin(line)
        	   if(line_len .gt. 0) then
           		l = 1
           		min_y = ctor(line(1:line_len),l)
        	   end if

           	   call prompt(' Enter max intensity    ['//cmaxy//']:')
        	   line_len = getlin(line)
        	   if(line_len .gt. 0) then
           		l = 1
           		max_y = ctor(line(1:line_len),l)
        	   end if
                end if

        	if_first = 1
c-----------------------------------------
c fix all free parameters 
c-----------------------------------------
        	np =0
        	do i=1,np_old
           	   inp(i) = 0
        	end do
c-----------------------------------------
c do a go 
c-----------------------------------------
                setnowarnings=1
        	call prtin
        	call pmfit (fun)
        	call prtout
                setnowarnings=0
        	fitted = .TRUE.
c--------------------------------------------
c print out what have been varied, pmin, pmax 
c--------------------------------------------
        	call remark    (' The following parameters have been varied by Monte-Carlo:')
 100    	format(1x,a20,a3,a13,a3,a13,a3)
       	 	do i=1,np_old
           		j = inp_old(i)
           		if( pmax(j).gt. pmin(j)) then
              			write(cpmin, fmt='(g13.5)' ) pmin(j)
             			write(cpmax, fmt='(g13.5)' ) pmax(j)
              			call remark (' '//nam(j)//' [ '//cpmin//' , '//cpmax//' ] ')
           		end if
        	end do       
c-----------------------------------------------
c call plot fit
c-----------------------------------------------
        	write(cminx, fmt='(f8.1)' ) min_x
        	write(cmaxx, fmt='(f8.1)' ) max_x
        	write(cminy, fmt='(f8.1)' ) min_y
        	write(cmaxy, fmt='(f8.1)' ) max_y
        	write(cint,  fmt='(i6)'   ) n_mc_par

        	line = 'p f'//cminx//cmaxx//cminy//cmaxy
        	l=2

        	call remark ('Plot Monte-Carlo-set '//cint)
        	call remark ('Please wait ...')
        	call plot (line, l, fitted, fun)
        	call remark ('Done.')


                if(ii.lt.nii) then
                     call prompt ('Continue plotting [Y]? ')
                     line_len = getlin(line)
                     c = 'y'
                     if( line_len .gt. 0 ) then
                         c = line(1:1)
                     end if
                     if(c.ne.'Y'.and.c.ne.'y') goto 21
                     
                end if
        end do
        call remark ('All numbers are plotted.')
        goto 22
c-----------------------------------------
c restore old situation
c-----------------------------------------
 21     call remark ('Plotting aborted.')    
 22     if_plot = 0

        np = np_old

        do i=1, np
           inp(i) = inp_old(i)
        end do

        do i=1,nptot
           p   (i) = p_old(i)
           pmin(i) = p_min(i)
           pmax(i) = p_max(i)
        end do

        ratioenergy = ratio_old
      return 
      end 

