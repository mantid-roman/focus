!	Subroutine SELBAC	-	Restores indices to selected list

	subroutine SELBAC (ntot, nsel, n, in, iarray) 

	include 'FRILLS_SOURCES:FRILLS.INC'
	logical*1 flag
	integer*4 ntot, nsel, n, in(*), iarray(*), i, j, k, im(MAX_DAT)

	i = 0
	do j = 1,ntot
		flag = .FALSE.
		do k = 1,nsel
			if (in(k) .eq. j) flag = .TRUE.
		end do
		do k = 1,n
			if (iarray(k) .eq. j) flag = .TRUE.
		end do
		if (flag) then
			i = i + 1
			im(i) = j
		end if
	end do
	nsel = i
	do i = 1,nsel
		in(i) = im(i)
	end do

	return
	end
