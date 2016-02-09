!	Subroutine PAROUT	-	Outputs parameter values to file

	subroutine PAROUT (line, l)

	include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:DEGENERATION.INC'
        include 'CF_SOURCES:RESOLUTION.INC'
        include 'CF_SOURCES:SELF.INC'
        include 'CF_SOURCES:OVERALL.INC'
        include 'CF_SOURCES:MONTE.INC'

	character*60 file_name
	character*(*) line
	integer*4 l, line_len, file_len, i, j, k
	real*8 sigma /0.0D0/

        integer*4 n_monte, n_par, is_ratio

	call str$trim (line, line, line_len)
	if (l .lt. line_len) then
		file_name = line(l:line_len)
		file_len = line_len - l + 1
	else
		file_name = 'SYS$SCRATCH:PAR.OUT'
		file_len = 19
	end if
	open (unit=42, file = file_name(1:file_len), status = 'NEW', err = 500)
	inquire (unit=42, name = file_name)
	call str$trim (file_name, file_name, file_len)
  101   format(1x,3i3,i5,8i3)
        write(42,101) setbkq,setvkq,setxr,nindvold, setakq,setwx,nre,symmetry,
     *                setres,setself,setoverall,if_monte

  105   format(1x,8f16.8)
        if(setres.eq.1) then
           write(42,105) fwhm_van_el,fwhm_chop_el,decay,lmc,lcs,lsd,de,di
        end if

  106   format(1x,2f16.8)
        if(setself.eq.1) then
           write(42,106) muts,muta
        end if

  107   format(1x,i2,2f16.8)
        if(setoverall.eq.1) then
           write(42,107) levelnumber,levelenergy,ratioenergy
	end if


  108   format(1x,3i2,5i6)
  109   format(1x,4f16.8)
        if(if_monte.eq.1) then
           n_monte = 0
           do i=1,MAX_CHI2
              if(s_chisq(i).ne.0.0d0) n_monte = n_monte + 1
           end do
           if(n_monte.gt.10) n_monte = 10
           n_par    = nptot
           is_ratio = 1
           write(42,108) if_monte, if_first, if_plot, n_mc_par, nset , n_monte, n_par, is_ratio
           write(42,109) min_x,max_x, min_y,max_y
           do i=1,n_monte
              write(42,104) s_chisq(i)
              if(is_ratio.ne.0) write(42,104) s_ratio(i)
              do k=1,n_par
                 write(42,104)  s_p(i,k)
                 write(42,104)  s_pmin(i,k)
                 write(42,104)  s_pmax(i,k)
              end do
           end do
        end if



  102   format(1x,i5,i5,2f16.8)
        do i=1,nindvold
          write(42,102) fixold(i), indvold(i), pminold(i), pmaxold(i)
        end do

  104   format(1x,f16.8)
        do i=2,6,2
           do k=0,i
              write(42,104) sbkq(i,k)
	   end do
	end do


  103   format(1x,a20)
        do i=1,nptot
          write(42,103) nam(i)
        end do

	write (42, 100) 
	j = 1
	do i = 1,nptot
		if (i .eq. inp(j) .and. j .le. np) then
			write (42, 200) i, ins(i), nam(i), p(i), sigma, pmin(i), pmax(i)
			j = j + 1
		else if (inr(i) .gt. 0) then
			write (42, 300) i, ins(i), nam(i), p(i), inr(i), pmin(i), pmax(i)
		else
			write (42, 400) i, ins(i), nam(i), p(i), pmin(i), pmax(i)
		end if
	end do	  
	close (42)
	call remark ('Parameters stored in '//file_name(1:file_len))
  100	format (5x, 'Set Parameter', 15x, 'Value', 7x, 'Sigma', 8x, 'Pmin', 9x, 'Pmax')
  200	format (1x, i3, 1x, i3, 1x, a20, 1x, 2g12.5, 2g13.5)
  300	format (1x, i3, 1x, i3, 1x, a20, 1x, g12.5, 2x, 'Bnd to ', i<int(log10(real(inr(i))))+1>, t55, 2g13.5)
  400	format (1x, i3, 1x, i3, 1x, a20, 1x, g12.5, 3x, 'Fixed', t55, 2g13.5)
	return  

  500	call remark ('ERROR: Unable to open '//file_name(1:file_len))
	return
	end
