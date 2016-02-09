!	Subroutine USER_HELP		-	Outputs help information to screen

	subroutine USER_HELP (line, l)

	character*1 letter
	character*8 word
	character*(*) line
	integer*4 l, len, getwrd
	
	equivalence (word, letter)

	len = getwrd (line, l, word)
	if (len .gt. 0) call upcase (word(1:len))

	if (len .eq. 0) then
	    	call remark ('Possible FOCUS commands are:')
		call remark (' @ Alter        - Alter Parameter, Resolution, Selfshielding, Overall splitting')
		call remark (' @ Display      - Display Angles, Cef symmetry, Degeneration, Excitations,')
                call remark ('                          IOn R3+, INtensities, Level scheme, MAtrix elements,')
                call remark ('                          MOment of R3+, Operator, Shape of excitations, ')
                call remark ('                          Wavefunctions')
		call remark (' @ Help         - Help Alter, Display, Help, Info, Montecarlo, Set')
		call remark (' @ Info         - Info Version, Scattering function, Lineshapes, Resolution,')
		call remark ('                  Intensities, Authors, Cef symmetries, Mu*t                ')
		call remark (' @ Montecarlo   - Monte-Carlo')
		call remark (' @ Set          - Set  Angles(Euler), Cef symmetry, Degeneration, Operator,')
		call remark ('                  Shape of line, Ion                     ')
                call remark (' Capital letters sign the minimal abbreviation.')   
		call remark (' Type @ HELP HELP for further informations.')       

	else if (letter .eq. 'A') then
		call remark (' @ Alter Parameter  [Bkq] - use the Bkq-parameter for fitting (default)')
		call remark ('                    [Akq] - use the Akq-parameter for fitting. They are defined ')
		call remark ('                            by: Akq = Bkq/Stevens_factor/r**k ')
		call remark ('                    [Vkq] - use the normalised Vkq-parameter for fitting')
		call remark ('                    [Xr]  - use n-dimensional sphere parameter')
		call remark ('                    [Wx]  - use the Lea-Leask-Wolf parameter x,W (only ')
		call remark ('                            possible if the crystal field is cubic)    ')
                call remark ('        Resolution        - change resolution parameters')
                call remark ('        Selfshielding     - change selfshielding parameters')
                call remark ('        Overall splitting - give a level to fix the overall CEF splitting')
                call remark (' Capital letters sign the minimal abbreviation.')   


	else if (letter .eq. 'D') then 
		call remark (' @ Display Angles...........- list the Euler-angles of the R3+ moment')
                call remark ('                              -> use @ s a to set the Euler-angles')
                call remark ('           Cef symmetry.....- list the CEF point symmetry the calculation is for')
                call remark ('           Degeneration.....- list if the degeneration flag is ON (degenerated')
                call remark ('                              energy levels are combined) or  OFF (degenerated ')
                call remark ('                              energy levels are NOT combined')
                call remark ('                              -> use @ s d to change the degeneration flag')
		call remark ('           Excitations......- list the calculated crystal field excitations')
		call remark ('           IOn..............- list the R3+ ion the calculation is for')
		call remark ('           INtensities......- list the CEF transition intensities between')
		call remark ('                              CEF energy levels ')
                call remark ('           Level scheme.....- list the CEF level scheme (CEF transition')
                call remark ('                              energies and transition matrix elements for')
                call remark ('                              a powdered sample)')
		call remark ('           MAtrix elements..- list the CEF transition matrix elements')
		call remark ('           MOment R3+.......- list the calculated magnetic moment of R3+')
                call remark ('           Operator.........- list if calculation uses stevens-operators (ON) ')
                call remark ('                              or 3j-symbols (OFF) ')
                call remark ('                              -> use @ s o to change the operator flag')
		call remark ('           Shape............- list the lineshape of a special excitation')
                call remark ('                              -> use @ s s to change the lineshape')
                call remark ('           Wavefunctions....- list the wave function of the R3+ ion')
                call remark (' Capital letters sign the minimal abbreviation.')   

        else if(letter .eq. 'H') then
        	call remark (' @ HELP is possible on the following COMMANDS:')
		call remark (' @ Help Alter        - Help provided on the ALTER        command')
		call remark (' @ Help Display      - Help provided on the DISPLAY      command')
		call remark (' @ Help Help         - Help provided on the HELP         command')
		call remark (' @ Help Info         - Help provided on the INFO         command')
		call remark (' @ Help Montecarlo   - Help provided on the MONTECARLO   command')
		call remark (' @ Help Set          - Help provided on the SET          command')
                call remark (' Capital letters sign the minimal abbreviation.')   

	else if (letter .eq. 'I') then 
		call remark (' @ Info Version.............- prints out which version of FOCUS you have')
		call remark ('        Scattering function - prints out which scattering functions FOCUS uses')
		call remark ('        Resolution..........- prints out resolution function definiton')  
		call remark ('        Intensities ........- prints out CEF intensity definition           ')
		call remark ('        Authors.............- prints out the origin of FOCUS libraries      ')
		call remark ('        Cef symmetries......- prints out the relationship between point     ')
		call remark ('                              symmetries and symmetry numbers               ')
		call remark ('        Mu*t................- prints out the definition of mu*t      ')
                call remark (' Capital letters sign the minimal abbreviation.')   

	else if (letter .eq. 'M') then 
		call remark (' @ Montecarlo  [New]  - Monte-Carlo on all parameters with finite boundaries.')
                call remark ('                                                                       ')
		call remark ('               [Display chi**2] - display all chi**2 values found by @ M N.  ')
                call remark ('                                                                           ')
                call remark ('               [Plot parameters] [n1 n2-n3 n4 ...] - plot parameter sets defined ') 
                call remark ('                                                     by @ M N over the neutron') 
                call remark ('                                                     data. n1,n2,... are the')
                call remark ('                                                     positions of the Chi**2')
                call remark ('                                                     values in the Chi**2 list.')
                call remark ('                                                     --> See @ M D.')
                call remark ('                                                                             ')
                call remark ('               [Load parameters] [n1] - load in a selected parameter set which') 
                call remark ('                                        has been defined by @ M N. n1 is the')
                call remark ('                                        position of a Chi**2 value in the Chi**2')
                call remark ('                                        list.')
                call remark ('                                        --> See @ M D.')
                call remark (' Capital letters sign the minimal abbreviation.')   

	else if (letter .eq. 'S') then
		call remark (' @ Set Angles             - set the Euler-angles alpha, beta and gamma')
		call remark ('       Cef symmetry [nr]  - set a CEF symmetry of symmetry number nr  ')
		call remark ('                            nr must be between 0 and 8.               ')
		call remark ('       Degeneration [OFf] - degenerated energy levels are not combined (default) ')
		call remark ('                    [ON]  - degenerated energy levels are combined ')
		call remark ('       Ion [R3+]          - set the R3+ ion. R=CE,ND,PR ... ')
		call remark ('       Operator     [OFf] - calculation uses 3j-symbols')                            
		call remark ('                    [ON]  - calculation uses steven-operators')       
                call remark ('       Shape [n1 n2-n3 n4 ...] - set a specified lineshape for a')
		call remark ('                                 special excitation. n1,n2,... are') 
                call remark ('                                 numbers of the excitations. See: @ D S')
                call remark (' Capital letters sign the minimal abbreviation.')   

        else
		call remark ('ERROR: @ h "'//letter//'" not recognized')
 
	end if	     

	return
	end
