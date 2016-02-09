!	Subroutine PARDIS	-	Displays parameters on screen

	subroutine PARDIS (line, l)

	include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.PAR'
	include 'CF_SOURCES:XR.INC'
        include 'CF_SOURCES:MONTE.INC'
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:OVERALL.INC'


	character line*(*), chint*8, cint*10
	logical*1 fixed
	integer*4 l, n, getari, iarray(MAX_PAR), i, j, k, n_on_call

        real*8   pj, pminj, pmaxj
        integer  n0, no_0_par, same

!       -------------------------------------------------------------
	integer*4 iiarray(MAX_PAR)
!       -------------------------------------------------------------

	n         = getari (line, l, iarray, MAX_PAR)
        n_on_call = n
	if (n .eq. 0) then
		n = 0
		do i = 1,nptot
			if (ins(i) .eq. 0 .or. ins(i) .eq. is) then
				n = n + 1
				iarray(n) = i
			end if
		end do
		call remark ('           Dataset '//chint(is))
		call remark ('           ==========')
	end if	
        if(if_monte.ne.0 .and. n_mc_par.ne.0) then
            write(cint, fmt='(i6)' ) n_mc_par
            call remark ('           Monte Carlo set '//cint)
        end if
	call remark ('           Parameter set             P             Pmin         Pmax')
!       -------------------------------------------------------------
        if(n_on_call.eq.0) call array(iarray,iiarray,n)
!       -------------------------------------------------------------
	do i = 1,n
!       -------------------------------------------------------------
		if(n_on_call.eq.0) then
                   j = iiarray(i)      
                else               
                   j = iarray(i)
                end if
!       -------------------------------------------------------------
		if (j .le. 0 .or. j .gt. nptot) then 
			call remark ('ERROR: Parameter no. outside range')
			return
		end if
		fixed = .TRUE.
		do k = 1,np
			if (j .eq. inp(k)) fixed = .FALSE.
		end do

                pj    = p(j)
                pminj = pmin(j)
                pmaxj = pmax(j)
                n0    = no_0_par()

                if(     setxr.eq.1 .and. j.ge.n0-5 ) then
                     if( j.eq.indvold(1) .or. same(j,indvold(1)).eq.0 ) then
                             pj    = pj   *ratioenergy
                             pminj = pminj*ratioenergy
                             pmaxj = pmaxj*ratioenergy
                     end if
                else if(setwx.eq.1 .and. j.ge.n0-5 ) then
                     if( (j.ge.n0-5.and.j.le.n0) .or. j.eq.n0+1 ) then
                          pj    = pj   *ratioenergy
                          pminj = pminj*ratioenergy
                          pmaxj = pmaxj*ratioenergy
                     end if
                else if( j.ge.n0-5 ) then
                    pj    = pj   *ratioenergy
                    pminj = pminj*ratioenergy
                    pmaxj = pmaxj*ratioenergy
                end if

		if (.not. fixed) then
			if (pmaxj .gt. pminj ) then
				write (*, 100) j, ins(j), nam(j), pj, pminj, pmaxj
			else
				write (*, 200) j, ins(j), nam(j), pj
			end if
		else if (inr(j) .gt. 0) then
			write (*, 300) j, ins(j), nam(j), pj, inr(j)
		else 
			write (*, 400) j, ins(j), nam(j), pj
		end if
	end do	  
			  
  100	format (1x, i5, 1x, i3, 2x, a20, g13.5,  2x,  2g13.5)
  200	format (1x, i5, 1x, i3, 2x, a20, g13.5)
  300	format (1x, i5, 1x, i3, 2x, a20, g13.5,  6x, 'Bound to parameter ',  i<int(log10(real(inr(j))))+1>)
  400	format (1x, i5, 1x, i3, 2x, a20, g13.5, 13x, 'Fixed')

        setpardis = 1 ! subroutine pardis/fitdis have been called

	return  
	end


        subroutine array(iarray,iiarray,n)
          implicit none
	  include 'CF_SOURCES:FRILLS.INC'
	  include 'CF_SOURCES:XR.INC'
	  include 'CF_SOURCES:CONTROL.INC'

          integer*4 iarray(MAX_PAR),iiarray(MAX_PAR),n
          integer*4 i,k,ii,j

          if(setxr.eq.0) then
              do i=1,n
                 iiarray(i) = iarray(i)
              end do
          else
              ii=0
              do i=1,n
                 j = iarray(i)
                 do k=1,nindvold
                    if(j.eq.indvold(k)) goto 10 
                 end do
                 ii = ii + 1
                 iiarray(ii) = j
  10          end do
              do i=1,nindvold
                 ii = ii + 1
                 iiarray(ii) = indvold(i)
              end do
          endif
	return  
	end

        subroutine a_inf(iinf)
          implicit none
	  include 'CF_SOURCES:FRILLS.INC'
	  include 'CF_SOURCES:XR.INC'
	  include 'CF_SOURCES:CONTROL.INC'

          integer*4 iinf(MAX_PAR),n
          integer*4 i,k,ii,j

          if(setxr.eq.0) then
              do i=1,nf
                 iinf(i) = inf(i)
              end do
          else
              ii=0
              do i=1,nf
                 do k=1,nindvold
                    if(indvold(k).eq.inf(i)) goto 10 
                 end do
                 ii = ii + 1
                 iinf(ii) = inf(i)
  10          end do

              do k=1,nindvold
                 do i=1,nf
                    if( indvold(k).ne.inf(i) ) goto 20
                        ii = ii + 1
                        iinf(ii) = inf(i)
  20             end do
              end do

          endif

	return  
	end

