	subroutine version_number(file_name) 

	include 'CF_SOURCES:FRILLS.INC'

	character*(*) file_name, c3*3
	character*132 line
	logical found
	integer*4 file_len, i, l, n, getarr
	real*4 array(4)
        real*8 ver

!Check that input file exists	
	call str$trim (file_name, file_name, file_len)
	inquire (file = file_name(1:file_len), exist = FOUND)
	if (.not.found) then
		call upcase (file_name(1:file_len))
		call remark ('"'//file_name(1:file_len)//'" does not exist')
		return
	end if
!       ---------------------------------------------------------------------
!       Open file, read in which kind of parameter was stored
!       and swap to the kind of parameters 
!       ---------------------------------------------------------------------
	open (unit=44, file = file_name(1:file_len), status = 'OLD', err = 200)

        read(44,*) 
        read(44,*) ver
	close (44)
        write(c3,'(f3.1)') ver
        write(*,*)
        call remark ( '            -------------------------------------------------- ')
        call remark ( '           |                                                  |')
        call remark ( '           |         *****   ***    ***   *   *   ***         |')     
        call remark ( '           |         *      *   *  *      *   *  *            |')   
        call remark ( '           |         ***    *   *  *      *   *   ***         |')   
        call remark ( '           |         *      *   *  *      *   *      *        |')      
        call remark ( '           |         *       ***    ***    ***    ***         |')
        call remark ( '           |      An Interactive Crystal Electric Field       |')
        call remark ( '           |         Parameter Fitting Package using          |')
        call remark ( '           |             Neutron Scattering Data              |')
        call remark ( '           |                  Dr Peter Fabi                   |')
        call remark ( '           |                Version : '//c3//'                     |')
        call remark ( '           | Home: Am Klostergarten 21                        |')
        call remark ( '           |       52499 Baesweiler                           |')
        call remark ( '           |       Germany                                    |')
        call remark ( '           | Tel : +49 2401 51844                             |')
        call remark ( '           |                                                  |')
        call remark ( '           | If you have detected a bug, have problems or     |')
        call remark ( '           | comments mail to PF01@isise.rl.ac.uk             |')
        call remark ( '            -------------------------------------------------- ')

	return

!Error in OPEN command

  200	call remark ('ERROR: Unable to open '//file_name(1:file_len))
	return  
	end
