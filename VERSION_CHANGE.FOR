	program version_change 

	include 'FOCUS_REF:FRILLS.INC'
	include 'FOCUS_REF:FOCUS.INC'

	character*132 line
	logical found
	integer*4 file_len, i, l, n, getarr
	real*4 array(4)
        real*8 version_old

!Check that input file exists	
	call str$trim (focus_version, focus_version, file_len)
	inquire (file = focus_version(1:file_len), exist = FOUND)
	if (.not.found) then
		call upcase (focus_version(1:file_len))
		call remark ('"'//focus_version(1:file_len)//'" does not exist')
		stop
	end if
!       ---------------------------------------------------------------------
!       Open file, read in which kind of parameter was stored
!       and swap to the kind of parameters 
!       ---------------------------------------------------------------------

	open (unit=44, file = focus_version(1:file_len), status = 'OLD', err = 200)
        read(44,1) line
        read(44,2) version_old
	close (44)

  1     format(1x,a132)
  2     format(1x,f5.1)
	open (unit=44, file = focus_version(1:file_len), status = 'NEW', err = 200)
        write(44,1) line
        write(44,2) version_old+0.1
!       write(44,2) version_old
	close (44)

        goto 300
!Error in OPEN command

  200	call remark ('ERROR: Unable to open '//focus_version(1:file_len))
	  
  300	end
