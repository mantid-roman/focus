      subroutine MULTI_FRILLS (fun)

!	Based on: 	FRILLS: An Interactive Least Squares Fitting Package
!			R. Osborn, Rutherford Appleton Laboratory Report RAL 91-011 (1991)

!	Modified:	Adapted for fitting of multiple datasets in May 1993 by
!			R.Osborn, Materials Science Division, Argonne, National Laboratory, Argonne, IL 60439-4845, USA

!	Command List
!	------------
!	Help		- help provided on each command
!	Display		- display data, parameters, calculation or fit
!	Plot		- plot data, calculation, fit, or residuals
!	Output		- output parameters or fit to disk file
!	Input		- input parameters from disk file
!	Keep		- keep hardcopy for plotting on the laser printer
!	Use		- use selected set of data points in future commands
!	Modify		- modify data points
!	Remove		- remove data points
!	Set		- set parameters
!	Fix		- fix parameters
!	Bind		- bind parameters
!	Clear		- clear fixed parameters or modified data
!	Limit		- limit parameters
!	Alter		- alter fitting parameters, plot characteristics and terminal emulation
!	Title		- change fit title
!	Go		- fit functions to data
!	Jump		- spawns DCL sub-process
!	Exit		- exit from routine
!       @focus          - focus defined command


      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:CONTROL.INC'

      integer no_commands,no_options
      parameter(no_commands=21,no_options=5)

      external fun

      character word*12, line*80
      character command(no_commands)*8, option(no_options)*11
      character*8 chint
      logical fitted
      integer*4 i, line_len, getlin, l, word_len, getwrd, j

      data command /	'DISPLAY', 'PLOT', 'INPUT', 'OUTPUT', 'MODIFY', 'REMOVE', 'SET', 'FIX', 'BIND', 'CLEAR', 'LIMIT', 
     >			'ALTER', 'TITLE', 'GO', 'HELP', 'KEEP', 'JUMP', 'USE', 'EXIT', '@FOCUS', 'ZERO'/
      data option /'PARAMETERS', 'DATA', 'CALCULATION', 'FIT', 'MATRIX'/

!-----------------------------------------------------------------------------------------------------------------------------------
!	Initialization
  	if (ns .le. 0) then
 		call remark ('ERROR: No. of datasets <= 0')
 		return
 	end if
 	do i = 1,ns
 		if (nd(i) .le. 0) then
 			call remark ('ERROR: No. of data points in dataset '//chint(i)//' <= 0')
 			return
 		end if
 	end do
 	if (nptot .le. 0) then
 		call remark ('ERROR: No. of parameters <= 0')
 		return
 	end if
 	call setup 
 	fitted = .FALSE.
!-----------------------------------------------------------------------------------------------------------------------------------
!	Input Command
   1	call prompt ('#')
 	line_len = getlin (line)
 	if (line_len .eq. 0) go to 1
 	l = 1
 	word_len = getwrd (line(1:line_len), l, word)
 	if (word_len .gt. 0) then
 		call upcase (word(1:word_len))
 	else
 		go to 1
 	end if
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is DISPLAY
  	if (word(1:word_len) .eq. command(1)(1:word_len)) then
 		word_len = getwrd (line(1:line_len), l, word)
 		do while (word_len .eq. 0) 
 			call prompt ('Data, Parameters, Calculation, Fit, Matrix :')
 			line_len = getlin (line)
 			if (line_len .eq. 0) then
 				word_len = 0
 			else
 				l = 1
 				word_len = getwrd (line(1:line_len), l, word)
 			end if
 		end do
 		if (word_len .gt. 0) call upcase (word(1:word_len))
!	 Option is PARAMETERS
 		if (word(1:word_len) .eq. option(1)(1:word_len)) then
 			call pardis (line, l)
!	 Option is DATA
 		else if (word(1:word_len) .eq. option(2)(1:word_len)) then
 			call datdis (line, l)
!	 Option is CALCULATION
 		else if (word(1:word_len) .eq. option(3)(1:word_len)) then
 			call caldis (fun, line, l)
!	 Option is FIT
 		else if (word(1:word_len) .eq. option(4)(1:word_len)) then
 			if (fitted) then
 				call fitdis 
 			else
 				call remark ('ERROR: No fitted parameters available')
 			end if
!	 Option is MATRIX
 		else if (word(1:word_len) .eq. option(5)(1:word_len)) then
 			if (fitted) then
 				call print_corr (6)
 			else
 				call remark ('ERROR: Correlation matrix not yet determined')
 			end if
 		else
 			call remark ('ERROR: "'//word(1:word_len)//'" not recognized')
 		end if
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is PLOT
 	else if (word(1:word_len) .eq. command(2)(1:word_len)) then
 		call plot (line, l, fitted, fun)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is INPUT
 	else if (word(1:word_len) .eq. command(3)(1:word_len)) then
 		if (l .gt. line_len) then
 			call prompt ('Input filename :')
 			line_len = getlin (line)
 			l = 1
 		else
 			l = l + 1
 		end if
 		if (line_len .gt. 0 .and. l .le. line_len) call parin (line(l:line_len))
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is OUTPUT
 	else if (word(1:word_len) .eq. command(4)(1:word_len)) then
 		word_len = getwrd (line(1:line_len), l, word)
 		do while (word_len .eq. 0) 
 			call prompt ('Parameters or Fit :')
 			line_len = getlin (line)
 			l = 1
 			word_len = getwrd (line(1:line_len), l, word)
 		end do
 		if (word_len .gt. 0) call upcase (word(1:word_len))
!	 Option is PARAMETERS
 		if (word(1:word_len) .eq. option(1)(1:word_len)) then
 			call parout (line, l)
!	 Option is CALCULATION
 		else if (word(1:word_len) .eq. option(3)(1:word_len)) then
 			call calout (fun, line, l)
!	 Option is FIT
 		else if (word(1:word_len) .eq. option(4)(1:word_len)) then
 			if (fitted) then
 				call fitout (line, l)
 			else
 				call remark ('ERROR: No fitted parameters available')
 			end if
 		else
 			call remark ('ERROR: "'//word(1:word_len)//'" not recognized')
 		end if
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is MODIFY
 	else if (word(1:word_len) .eq. command(5)(1:word_len)) then
 		call modify (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is REMOVE
 	else if (word(1:word_len) .eq. command(6)(1:word_len)) then
 		call remove_data (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is SET
 	else if (word(1:word_len) .eq. command(7)(1:word_len)) then
 		call setpar (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is FIX
 	else if (word(1:word_len) .eq. command(8)(1:word_len)) then
 		call fixpar (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is BIND
 	else if (word(1:word_len) .eq. command(9)(1:word_len)) then
 		call bndpar (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is CLEAR
 	else if (word(1:word_len) .eq. command(10)(1:word_len)) then
 		word_len = getwrd (line(1:line_len), l, word)
 		if (word_len .gt. 0) call upcase (word(1:word_len))
 		if (word_len .eq. 0) then
 			call setup
!	 Option is PARAMETERS
 		else if (word(1:word_len) .eq. option(1)(1:word_len)) then
 			call clepar (line, l)
!	 Option is DATA
 		else if (word(1:word_len) .eq. option(2)(1:word_len)) then
 			call cledat (line, l)
 		else
 			call remark ('ERROR: "'//word(1:word_len)//'" not recognized')
 		end if
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is LIMIT
 	else if (word(1:word_len) .eq. command(11)(1:word_len)) then
 		call limpar (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is ALTER
 	else if (word(1:word_len) .eq. command(12)(1:word_len)) then
 		call alter (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is TITLE
 	else if (word(1:word_len) .eq. command(13)(1:word_len)) then
 		i = 80
 		do while (i .gt. 1 .and. title(is)(i:i) .eq. ' ')
 			i = i - 1
 		end do
 		call remark ('Old title	: '//title(is)(1:i))
 		call prompt ('Give new title :')
 		line_len = getlin (line)
 		if (line_len .gt. 0) title(is) = line(1:line_len)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is GO
  	else if (word(1:word_len) .eq. command(14)(1:word_len)) then
 		do i = 1,ns
 			do j = 1,nv(i)
 				if (vsig(inv(j,i),i) .le. 0.0) then
 					call remark ('ERROR: Zero error in point '//chint(inv(j,i)))
 					go to 1
 				end if
 			end do
 		end do
                setnowarnings=1
 		call prtin
  		call pmfit (fun)
  		call prtout
                setnowarnings=0
 		fitted = .TRUE.
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is HELP
 	else if (word(1:word_len) .eq. command(15)(1:word_len)) then
 		call help (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is KEEP
 	else if (word(1:word_len) .eq. command(16)(1:word_len)) then
 		call keep_plot (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is JUMP
 	else if (word(1:word_len) .eq. command(17)(1:word_len)) then
 		if (l .lt. line_len) then
 			call lib$spawn (line(l:line_len))
 		else
 			call remark ('# Spawning DCL sub-process - LOGOUT to return to FRILLS')
 			call lib$spawn ()
 		end if
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is USE
 	else if (word(1:word_len) .eq. command(18)(1:word_len)) then
 		call use_data (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is EXIT
 	else if (word(1:word_len) .eq. command(19)(1:word_len)) then
                setnowarnings=1
  		return                                                                     
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is @USER
 	else if (word(1:word_len) .eq. command(20)(1:word_len)) then
 		call user_command (fun, line, l)                                                                     
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is ZERO
 	else if (word(1:word_len) .eq. command(21)(1:word_len)) then
 		call unlimpar (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is unrecognized
 	else 
 		call remark ('ERROR: "'//word(1:word_len)//'" not recognized')
 	end if
      go to 1
      end
