!	Subroutine CLEPAR	-	Clears parameters that were formerly fixed or bound to another parameter

	subroutine clepar (line, l)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:CONTROL.INC'

        integer n0,no_0_par,nindv0 

	character*(*) line
	integer*4 n, getari, l, iarray(MAX_PAR), i,j,k, is_in 

	n = getari (line, l, iarray, MAX_PAR)
        nindv0 = nindvold
	if (n .eq. 0) then
		np = nptot
		do i = 1,np
			inp(i) = i
			inr(i) = 0
		end do
                if(setxr.eq.1) then
                   n0 = no_0_par()
                   do i=n0-5,nptot
                      is_in=0
                      do k=1,nindvold
                         l = indvold(k)
                         if(i.eq.l) is_in =1
                      end do
                      if(is_in.eq.0) then
                         nindvold = nindvold + 1
                         indvold(nindvold) = i
                      end if
                   end do
                end if
	else
                n0 = no_0_par()
		do i = 1,n
                     j =iarray(i)
	             inr(j) = 0
                     if(setxr.eq.1 .and. j.le.nptot .and. j.ge.n0-5) then
                        is_in=0
                        do k=1,nindvold
                           l = indvold(k)
                           if(j.eq.l) is_in=1
                        end do
                        if(is_in.eq.0) then
                           nindvold = nindvold + 1
                           indvold(nindvold) = j
                        end if
                     end if
		end do
		call selbac (nptot, np, n, inp, iarray) 
	end if
        if(setxr.eq.1) call xr_add(nindv0)

	return
	end
