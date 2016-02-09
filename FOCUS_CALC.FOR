! Calls the user-defined function
	subroutine FOCUS_CALC (fun)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:CONTROL.INC'

	external  fun
	integer*4 j, nvtot
	real*8    c(MAX_DAT)

        setfit = 0
	call prtin
        j     = 1
        nvtot = nv(1)
	call function_call (j, nv(j), inv(1,j), x(1,j), fun, c, p)
	call prtout
	setfit = 1

	return
	end
