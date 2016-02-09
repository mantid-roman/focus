!	Subroutine FITDIS	-	Display fitted parameters on screen
!       13/4/96    P.Fabi  modified for use of FOCUS

	subroutine FITDIS 

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:OVERALL.INC'

        real*8   pfi, pfsigi
        integer  n0, no_0_par, same

	integer*4 i, j, k
!       -------------------------------------------------------------
	integer*4 iiarray(MAX_PAR), iinf(MAX_PAR)
	integer*4  iarray(MAX_PAR)
!       -------------------------------------------------------------

	call newlin (1)
	call remark ('           '//message)
	call newlin (1)
	write (*, 100) chisq
	call newlin (1)
	call remark ('           Parameter set             P            Sigma')
!       -------------------------------------------------------------
        do i=1,nptot
           iarray(i) = i
        end do
        call array(iarray,iiarray,nptot)                ! defined in pardis.for 
        call a_inf(iinf)
!       -------------------------------------------------------------
	j = 1
        n0     = no_0_par()
	do k = 1,nptot
!       -------------------------------------------------------------
                 i = iiarray(k)
!       -------------------------------------------------------------
                pfi    = pf(i)
                pfsigi = pfsig(i)

                if(     setxr.eq.1 .and. i.ge.n0-5 ) then
                     if( i.eq.indvold(1) .or. same(i,indvold(1)).eq.0 ) then
                          pfi    = pfi   *ratioenergy
                          pfsigi = pfsigi*ratioenergy
                     end if
                else if(setwx.eq.1 .and. i.ge.n0-5 ) then
                     if( (i.ge.n0-5.and.i.le.n0) .or. i.eq.n0+1 ) then
                          pfi    = pfi   *ratioenergy
                          pfsigi = pfsigi*ratioenergy
                     end if
                else if( i.ge.n0-5 ) then
                    pfi    = pfi   *ratioenergy
                    pfsigi = pfsigi*ratioenergy
                end if

		if (i.eq. iinf(j) .and. j .le. nf) then
			write (*, 200) i, ins(i), nam(i), pfi, pfsigi
			j = j + 1
		else if (inrf(i) .gt. 0) then
			write (*, 300), i, ins(i), nam(i), pfi, inrf(i)
		else 
			write (*, 400) i, ins(i), nam(i), pfi
		end if
	end do
	call newlin (1)

  100	format (24x, 'Chi ** 2', g13.5)
  200	format (1x, i5, 1x, i3, 2x, a20, g13.5, 2x, g13.5)
  300	format (1x, i5, 1x, i3, 2x, a20, g13.5, 6x, 'Bound to parameter ', i<int(log10(real(inrf(i))))+1>)
  400	format (1x, i5, 1x, i3, 2x, a20, g13.5, 6x, 'Fixed')
  500	format (2x,  a5, ':', 5g13.5)
  600	format (2x, '+/-  :', 5g13.5)
  700	format (2x, 'Fit  :', 5g13.5)

        setpardis = 1 ! subroutine pardis/fitdis have been called

	return
	end
