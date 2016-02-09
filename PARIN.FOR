!	Subroutine PARIN - Reads in parameters and limits from standard output file produced by PAROUT or FITOUT

	subroutine parin (file_name)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:DEGENERATION.INC'
        include 'CF_SOURCES:RESOLUTION.INC'
        include 'CF_SOURCES:INSTRUMENT.INC'
        include 'CF_SOURCES:SELF.INC'
        include 'CF_SOURCES:OVERALL.INC'
        include 'CF_SOURCES:MONTE.INC'


	character*(*) file_name
	character*132 line
	logical found
	integer*4 file_len, i, l, n, getarr, k
	real*4 array(4)
        integer nnre,ssymmetry
	integer n0,no_0_par,no_cf_par,ncftot

        character*64 text0,text1,text2,text3,text4,text5,text6,text7,c7*7
        character*64 text10,text11,text12,text13
        character*64 text20,text21,text22,text23,text24,text30
   	real*8      c_fmevkelvin  
   	external    c_fmevkelvin 

        integer*4 n_monte, n_par, is_ratio
 

           data text0  /' FOCUS does not calculate the resolution function.              '/
           data text1  /' FOCUS calculates the resolution function.                      '/
           data text2  /' FOCUS sets elastic vanadium line width to.....: 0000.0000 meV. '/
           data text3  /' FOCUS sets chopper contribution to............: 0000.0000 meV. '/
           data text4  /' FOCUS sets decay constant to..................: 0000.0000 meV. '/
           data text5  /' FOCUS sets moderator-chopper distance to......: 0000.0000 m.   '/
           data text6  /' FOCUS sets chopper-sample    distance to......: 0000.0000 m.   '/
           data text7  /' FOCUS sets sample-detector   distance to......: 0000.0000 m.   '/

           data text10 /' FOCUS does not correct for selfshielding and absorption.       '/
           data text11 /' FOCUS corrects for selfshielding and absorption.               '/
           data text12 /' FOCUS sets mu*t (scattering) to...............: 000000.00.     '/
           data text13 /' FOCUS sets mu*t (absorption) to...............: 000000.00.     '/


           data text20 /' FOCUS does not fix the overall CEF splitting.                  '/
           data text21 /' FOCUS fix the overall CEF splitting.                           '/
           data text22 /' FOCUS sets the fixed CEF level to.............: 00.            '/
           data text23 /' FOCUS sets the desired energy  to.............: 000000.00 meV  '/
           data text24 /' FOCUS sets the CEF energy scaling factor to...: 000000.00      '/

           data text30 /' FOCUS has Monte Carlo parameter loaded. See @ m d.             '/


!Check that input file exists	
	call str$trim (file_name, file_name, file_len)
	inquire (file = file_name(1:file_len), exist = FOUND)
	if (.not.found) then
		call upcase (file_name(1:file_len))
		call remark ('"'//file_name(1:file_len)//'" does not exist')
		return
	end if
!       ---------------------------------------------------------------------
!       Open file, read in which kind of parameter was stored
!       and swap to the kind of parameters 
!       ---------------------------------------------------------------------
	open (unit=43, file = file_name(1:file_len), status = 'OLD', err = 200)
  101   format(1x,3i3,i5,8i3)
  102   format(1x,i5,i5,2f16.8)
  103   format(1x,a20)
  104   format(1x,f16.8)
  105   format(1x,8f16.8)
  106   format(1x,2f16.8)
  107   format(1x,i2,2f16.8)
  108   format(1x,3i2,5i6)
  109   format(1x,4f16.8)
!       --------------------------------------------------------------------------------
!       read control parameter in
!       focus version 1.0 does not have setakq,setwx,nnre,ssymmetry
!       --------------------------------------------------------------------------------
        read(43,101) setbkq,setvkq,setxr,nindvold, setakq,setwx,nnre,ssymmetry,
     *               setres,setself,setoverall,if_monte
        if(setakq.eq.1) type = 'Akq'
        if(setbkq.eq.1) type = 'Bkq'
        if(setvkq.eq.1) type = 'Vkq'
        if(setxr. eq.1) type = 'XR'
        if(setwx. eq.1) type = 'Wx'
	
!       --------------------------------------------------------------------------------
!       reads in the resolution function values
!       --------------------------------------------------------------------------------
        if(setres.eq.1) then
           read(43,105) fwhm_van_el,fwhm_chop_el,decay,lmc,lcs,lsd,de,di
           call remark(text1)
           write(text2(50:58),fmt='(f9.4)') fwhm_van_el
           write(text3(50:58),fmt='(f9.4)') fwhm_chop_el
           write(text4(50:58),fmt='(f9.4)') decay          
           write(text5(50:58),fmt='(f9.4)') lmc
           write(text6(50:58),fmt='(f9.4)') lcs
           write(text7(50:58),fmt='(f9.4)') lsd
           call remark(text2)
           call remark(text3)
           call remark(text4)
           call remark(text5)
           call remark(text6)
           call remark(text7)
        else
           write(*,*)
           call remark(text0)
           fwhm_van_el  = 0.0d0
           fwhm_chop_el = 0.0d0
           decay = 0
           lmc = lmc_het
           lcs = lcs_het
           lsd = lsd_het
           de  = de_default*c_fmevkelvin()
        end if
        write(c7,fmt='(f7.4)') de
        call remark (' FOCUS sets the energy resolution to de =  '//c7//' Kelvin.')
!       --------------------------------------------------------------------------------
!       reads in the selfshielding values
!       --------------------------------------------------------------------------------
        if(setself.eq.1) then
           read(43,106) muts,muta
           write(*,*)
           call remark(text11)
           write(text12(50:58),fmt='(f9.2)') muts
           write(text13(50:58),fmt='(f9.2)') muta
           call remark(text12)
           call remark(text13)
        else
           write(*,*)
           call remark(text10)
           muts = 0.0d0
           muta = 0.0d0
        end if
!       --------------------------------------------------------------------------------
!       reads in the overall values
!       --------------------------------------------------------------------------------
        if(setoverall.eq.1) then
           read(43,107) levelnumber, levelenergy, ratioenergy
           write(*,*)
           call remark(text21)
           write(text22(50:51),fmt='(i2)'  ) levelnumber
           write(text23(50:58),fmt='(f9.2)') levelenergy
           write(text24(50:58),fmt='(f9.2)') ratioenergy
           call remark(text22)
           call remark(text23)
           call remark(text24)
        else
           write(*,*)
	   call remark(text20)
           levelnumber = 2
           levelenergy = 1.0d0
           ratioenergy = 1.0d0
        end if
!       --------------------------------------------------------------------------------
!       reads in the monte carlo values
!       --------------------------------------------------------------------------------
        if(if_monte.eq.1) then
           write(*,*)
           call remark (text30)
           do i=1,MAX_CHI2
              s_chisq(i) = 0.0d0
           end do
           read(43,108) if_monte, if_first, if_plot, n_mc_par, nset, n_monte, n_par, is_ratio
           read(43,109) min_x,max_x, min_y,max_y
           do i=1,n_monte
              read(43,104) s_chisq(i)
              if(is_ratio.ne.0) read(43,104) s_ratio(i)
              do k=1,n_par
                 read(43,104)  s_p(i,k)
                 read(43,104)  s_pmin(i,k)
                 read(43,104)  s_pmax(i,k)
              end do
           end do
        end if
!       --------------------------------------------------------------------------------
!       read which parameter are fixed and which limits they have
!       --------------------------------------------------------------------------------
        if(nindvold.ne.0) then
          do i=1,nindvold
            read(43,102) fixold(i), indvold(i), pminold(i), pmaxold(i)
          end do
        end if
 
!       --------------------------------------------------------------------------------
!       looks if their is a new focus version parameter file
!       if so sets the new rare earth and the symmetry value
!       --------------------------------------------------------------------------------
        if(nnre    .ge.1 .and. nnre    .le.13 .and. 
     *     symmetry.ge.0 .and. symmetry.le. 8 .and.
     *     setwx   .ge.0 .and. setwx   .le. 1 .and.
     *     setakq  .ge.0 .and. setakq  .le. 1       ) then
           nre      = nnre
	   symmetry = ssymmetry
           do i=2,6,2
              do k=0,i
                 read(43,104) sbkq(i,k)
	      end do
	   end do
	end if

!       --------------------------------------------------------------------------------
!       find out how many fit parameters there are in total
!       --------------------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)
        n0     = no_0_par()
        nptot  = n0 + ncftot 
	

!       --------------------------------------------------------------------------------
!       read the names of the new parameter from the parameter file
!       --------------------------------------------------------------------------------
        do i=1,nptot
          read(43,103) nam(i)
        end do
 

!       --------------------------------------------------------------------------------
!       search for first line of parameter values
!       --------------------------------------------------------------------------------
	line = ' '
	do while (line(1:18) .ne. '     Set Parameter')
		read (43, '(a18)', end = 100) line
	end do

!       --------------------------------------------------------------------------------
!       Loop over number of expected parameters
!       --------------------------------------------------------------------------------
	np = 0
	do i = 1,nptot
		read (43, '(a)', end = 100) line
		l = 30
		n = getarr (line, l, array, 4)		!Read line for list of values to be stored in ARRAY (N = no. of values)
		if (n .gt. 0) then
			p(i) = array(1)			!Parameter value
			inr(i) = 0
		end if
		if (n .eq. 4 .and. l .gt. 0) then	!N=4 - Ignore 2nd value which is parameter error in FITOUT
			pmin(i) = array(3)
			pmax(i) = array(4)
			np = np + 1
			inp(np) = i
		else if (n .eq. 1 .and. l .lt. 0) then	!N=1 - Non-numeric characters found on line i.e. "Fixed" or "Bound"
			l = 56
			n = getarr (line, l, array, 2)	!Read end of line for PMIN and PMAX
			if (n .eq. 2) then
				pmin(i) = array(1)
				pmax(i) = array(2)
			end if
			if (line(45:47) .eq. 'Bnd') read (line(52:54), *) inr(i) !Parameter is bound by fixed ratio to P(INR(i))
		end if
	end do

!       --------------------------------------------------------------------------------
!       Loop over parameters to determine ratios of bound parameters
!       --------------------------------------------------------------------------------
	do i = 1,nptot
		if (inr(i) .gt. 0) then
			if (p(inr(i)) .ne. 0.0) then
				rp(i) = p(i) / p(inr(i))
			else
				inr(i) = 0
			end if
		end if
	end do
	close (43)

	return
!End of file found before expected number of parameters read in
  100	call remark ('WARNING: '//file_name(1:file_len)//' contains too few parameters')
	close (43)
	return
!Error in OPEN command
  200	call remark ('ERROR: Unable to open '//file_name(1:file_len))
	return  
	end
