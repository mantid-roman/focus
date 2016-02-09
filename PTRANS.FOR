!	Subroutine PTRANS	-	Transforms parameters to circular coordinates
!       P.Fabi asin changed to acos to be consistent with (x,R)-transformation

	subroutine PTRANS (np, p, inp, pmin, pmax)

	real*8 p(*), pmin(*), pmax(*)
	integer*4 np, inp(*), i, j, k, l

	include 'CF_SOURCES:FRILLS.PAR'

	do i = 1,np
		j = inp(i)
                if (pmax(j) .gt. pmin(j)) then
		    	  p(j) = acos ((2.0D0*p(j)-pmax(j)-pmin(j)) / (pmax(j)-pmin(j)))
                end if
   	end do

        return
	end
