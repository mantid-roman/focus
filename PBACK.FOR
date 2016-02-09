!	Subroutine PBACK	-	Reverses parameter transforms
!       P.Fabi sin changed to cos to be consistent with (x,R)-transformation

	subroutine PBACK (np, nptot, p, inp, pmin, pmax, inr, rp)

	real*8 p(*), pmin(*), pmax(*), rp(*)
	integer*4 np, nptot, inp(*), inr(*), i, j

	include 'CF_SOURCES:FRILLS.PAR'
        
	do i = 1,np
		j = inp(i)
		if (pmax(j) .gt. pmin(j)) then
			p(j) = 0.5 * (pmax(j) + pmin(j)) + 0.5 * (pmax(j) - pmin(j)) * cos (p(j))
		end if
  	end do

	do i = 1,nptot
		if (inr(i) .gt. 0) p(i) = p(inr(i)) * rp(i)
	end do

	return
	end
