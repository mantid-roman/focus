!	Subroutine PRTOUT	-	Outputs results of fit to output file
!
!  Modified:
!    19/9/94	T.G.Perring	FUN now called through the subroutine FUNCTION_CALL
!    13/4/96    P.Fabi          modified for use of FOCUS

	subroutine PRTOUT (fun)

	include 'CF_SOURCES:FRILLS.INC'
	external fun
	integer*4 i, j, k
	real*8 yf(MAX_DAT), dif, pcdif
!       -------------------------------------------------------------
	integer*4 iiarray(MAX_PAR)
	integer*4  iarray(MAX_PAR)
!       -------------------------------------------------------------

	call newpag (41)
	write (41, 100) message, chisq, ((i, title(i)), i = 1,ns) 
!       -------------------------------------------------------------
        do i=1,nptot
           iarray(i) = i
        end do
        call array(iarray,iiarray,nptot)
!       -------------------------------------------------------------
	j = 1
	do k = 1,nptot
!       -------------------------------------------------------------
                i = iiarray(k) 
!       -------------------------------------------------------------
		if (i .eq. inf(j) .and. j .le. nf) then
			write (41, 200) i, ins(i), nam(i), pf(i), pfsig(i), pfmin(i), pfmax(i)
			j = j + 1
		else if (inrf(i) .gt. 0) then
			write (41, 300) i, ins(i), nam(i), pf(i), inrf(i), pfmin(i), pfmax(i)
		else
			write (41, 400) i, ins(i), nam(i), pf(i), pfmin(i), pfmax(i)
		end if
	end do	  
	write (41,500) fit_time
	if (output .ge. 10) then
		do i = 1,ns
			call function_call (i, nv(i), inv(1,i), x(1,i), fun, yf, pf)
!			call fun (i, nv(i), inv(1,i), x(1,i), yf, pf)
			call newpag (41)
			write (41, 600) i, title(i), xlab(i), ylab(i)
			do j = 1,nv(i)
				k = inv(j,i)
				dif = yf(k) - v(k,i)
				if (v(k,i) .eq. 0.0) then
					pcdif = 0.0
				else
					pcdif = dif * 100. / v(k,i)
				end if
				write (41,700) k, x(k,i), v(k,i), vsig(k,i), yf(k), dif, pcdif
			end do
		end do
	end if

  100	format (//, 9x, a50, //, 20x, 'Chi ** 2', g13.5, //, <ns>(5x, 'Dataset ', i2, 5x, a80, /), /, 
	1	5x, 'Set Final Parameters', 8x, 'Value', 7x, 'Sigma', 8x, 'Pmin', 9x, 'Pmax', /,
	2	5x, '=== ================', 8x, '=====', 7x, '=====', 8x, '====', 9x, '====')
  200	format (1x, i3, 1x, i3, 1x, a20, 1x, 2g12.5, 2g13.5)
  300	format (1x, i3, 1x, i3, 1x, a20, 1x, g12.5, 2x, 'Bnd to ', i<int(log10(real(inr(i))))+1>, t55, 2g13.5)
  400	format (1x, i3, 1x, i3, 1x, a20, 1x, g12.5, 3x, 'Fixed', t55, 2g13.5)
  500	format (/, 9x, 'Fit time : ', a20, /)
  600	format (//, 5x, 'Dataset ', i2, 5x, a80, //, 
	1	'     N    ',a5,'        ',a5,'        +/-          Fit         Delta        %Delta', /)
  700	format (1x, i4, 5g13.5, f10.2) 

	return
	end
