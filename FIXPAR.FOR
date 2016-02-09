!	Subroutine FIXPAR	-	Sets parameters to be fixed in fits

	subroutine fixpar (line, l)

	include 'FRILLS_SOURCES:FRILLS.INC'

	character*(*) line
	integer*4 n, getari, iarray(MAX_PAR), l, line_len, getlin

	n = getari (line, l, iarray, MAX_PAR)
	if (n .eq. 0) then
		call prompt ('Give parameters to be fixed :')
		line_len = getlin (line)
		if (line_len .gt. 0) then
			l = 1
			n = getari (line(1:line_len), l, iarray, MAX_PAR)
		end if
	end if
	if (n .gt. 0) call selout (nptot, np, n, inp, iarray)

	return
	end
