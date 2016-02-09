!	Subroutine USER_COMMAND (fun, line, l)

!	Based on: 	MULTI_FRILLS: Adapted for fitting of multiple datasets in May 1993 by
!			R.Osborn, Materials Science Division, Argonne, National Laboratory, Argonne, IL 60439-4845, USA

!	Modified:       Adapted for crystal field fitting (FOCUS) by
!	                P. Fabi, ISIS Facility, Rutherford Appleton Laboratory, Chilton, Didcot, Oxon, OX11 0QX, UK 

!       

	subroutine user_command (fun, line, l)

	include 'CF_SOURCES:FRILLS.INC'
 	include 'CF_SOURCES:FOCUS.INC'
 	include 'CF_SOURCES:MONTE.INC'
 	include 'CF_SOURCES:XR.INC'
 	include 'CF_SOURCES:CONTROL.INC'

        external fun
                                                             
        integer out      ! messages a written to file out
        integer in       ! reading is from file in
        parameter(in=5,out=5)

	character word*12, line*80, cint*6, c
	character command(6)*20, option(12)*20
        character s_option(6)*20,ss_option(2)*20
	character a_option(4)*20,aa_option(5)*20
	character i_option(8)*20
	character m_option(4)*20
	character r_option(13)*20, rr_option(9)*20
      
	logical fitted
	integer*4 i, line_len, getlin, l, word_len, getwrd, j, ctoi, getari, nmax, nnre,ssymmetry

	data command   /'DISPLAY', 'HELP', 'SET', 'ALTER', 'INFO', 'MONTECARLO'/
	data  option   /'ION R3+', 'CEF SYMMETRY', 'EXCITATIONS', 'SHAPE OF LINES', 'MOMENT OF R3+','LEVEL SCHEME',
     *                 'MATRIX ELEMENTS', 
     *                 'WAVEFUNCTIONS', 'DEGENERATION', 'INTENSITIES', 'OPERATOR', 'ANGLES'/
        data s_option  /'SHAPE','DEGENERATION','OPERATOR','ANGLES','ION R3+','CEF SYMMETRY'/
        data ss_option /'OFF','ON'/
        data r_option  /'CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB'/
        data rr_option /'0','1','2','3','4','5','6','7','8'/
        data a_option  /'PARAMETER','RESOLUTION','SELFSHIELDING','OVERALL SPLITTING'/
        data aa_option /'BKQ','VKQ','XR','AKQ', 'WX'/
        data i_option  /'VERSION', 'SCATTERING FUNCTION', 'LINESHAPES', 'RESOLUTION', 'INTENSITIES',
     *                  'AUTHORS', 'CEF SYMMETRIES', 'MU*t'/
        data m_option  /'NEW', 'DISPLAY CHI**2', 'PlOT PARAMETERS', 'LOAD PARAMETERS'/
!-----------------------------------------------------------------------------------------------------------------------------------
!	Search for the user defined commands
        line_len = len(line)
 	word_len = getwrd (line(1:line_len), l, word)
 	if (word_len .gt. 0) then
 		call upcase (word(1:word_len))
 	else 
 		call user_help (line, l)
                return
 	end if
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is DISPLAY
	if (word(1:word_len) .eq. command(1)(1:word_len)) then
		word_len = getwrd (line(1:line_len), l, word)
		do while (word_len .eq. 0)
                        call remark (' Display what? Implemented are:') 
			call remark (' Angles, Cef symmetry , Degeneration, Excitations, ')
                        call remark (' INtensities  , IOn R3+, Level scheme, Operator, MAtrix elements, ')
                        call prompt (' MOment of R3+, Shape of lines, Wavefunctions          :')
			line_len = getlin (line)
			if (line_len .eq. 0) then
				word_len = 0
			else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			end if
		end do
		if (word_len .gt. 0) call upcase (word(1:word_len))       
!	Option is ION
		if (word(1:2) .eq. option(1)(1:2)) then
			call iondis 
!	Option is INTENSITY
		else if (word(1:2) .eq. option(10)(1:2)) then
			call focus_calc (fun)
              		call intensitydis
!	Option is CEF 
		else if (word(1:word_len) .eq. option(2)(1:word_len)) then
			call symdis 
!	Option is ANGLES 
		else if (word(1:word_len) .eq. option(12)(1:word_len)) then
			call eulerdis 
!	Option is DEGENERATION
		else if (word(1:word_len) .eq. option(9)(1:word_len)) then
			 if(setdegoff.eq.1) then
                            call remark (' The degeneration flag is OFF.')
                            call remark (' Degenerated energy levels are not combined!') 
                         else
                            call remark (' The degeneration flag is ON.')
                            call remark (' Degenerated energy levels are combined!') 
                         end if
                         call remark    (' Use @ s d  to change the degeneration flag.')
!	Option is OPERATOR
		else if (word(1:word_len) .eq. option(11)(1:word_len)) then
			 if(stdhamoff.eq.1) then
                            call remark (' The operator flag is OFF.')
                            call remark (' The hamiltonian is calculated by 3j-symbols!') 
                         else
                            call remark (' The operator flag is ON.')
                            call remark (' The hamiltonian is calculated by steven-operators!') 
                         end if
                         call remark    (' Use @ s o to change the operator flag.')
!	Option is EXCITATIONS
		else if (word(1:word_len) .eq. option(3)(1:word_len)) then
			call focus_calc (fun)
			call excdis 
!	Option is SHAPE
		else if (word(1:word_len) .eq. option(4)(1:word_len)) then
			call focus_calc (fun)
			call shapedis 
!	Option is MATRIXELEMENTS
!	Option is MOMENT
		else if (word(1:1) .eq. option(5)(1:1) .or. word(1:1) .eq. option(7)(1:1) ) then
                        if (word(1:2) .eq. option(5)(1:2)) then
			    call focus_calc (fun)
                            call momentdis
		        else if (word(1:2) .eq. option(7)(1:2)) then 
			    call focus_calc (fun)
                            call matrixelements
                        else
		            call remark ('ERROR: @ D "'//word(1:2)//'" not recognized')
                        end if
!	Option is LEVELSCHEME
		else if (word(1:word_len) .eq. option(6)(1:word_len)) then
			call focus_calc (fun)
			call levelscheme 
!	Option is WAVEFUNCTIONS
		else if (word(1:word_len) .eq. option(8)(1:word_len)) then
			call focus_calc (fun)
			call wavefunctions 
		else                                                                                       
			call remark ('ERROR: @ D "'//word(1:word_len)//'" not recognized')
		end if                                  
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is HELP
	else if (word(1:word_len) .eq. command(2)(1:word_len)) then
		call user_help (line, l)
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is SET
	else if (word(1:word_len) .eq. command(3)(1:word_len)) then
		word_len = getwrd (line(1:line_len), l, word)
		do while (word_len .eq. 0) 
                        call remark ('Set what? Implemented are:')
			call remark ('Angles(Euler), Cef Symmetry, Degeneration, Operator, ')
			call prompt ('Shape of line,Ion                                  :')
			line_len = getlin (line)
			if (line_len .eq. 0) then
				word_len = 0
			else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			end if
		end do
		if (word_len .gt. 0) call upcase (word(1:word_len))
!	Option is SHAPE
		if (word(1:word_len) .eq. s_option(1)(1:word_len)) then
			call focus_calc (fun)
			call setshape (line, l) 
!	Option is ANGLES
		else if (word(1:word_len) .eq. s_option(4)(1:word_len)) then
				call seteuler (line, l)                                                                                      
!	Option is DEGENERATION
		else if (word(1:word_len) .eq. s_option(2)(1:word_len)) then
 	        	word_len = getwrd (line(1:line_len), l, word)
	        	do while (word_len .eq. 0)
                           call remark ('Set Degeneration to what? Implemented are:') 
                           call prompt ('OFf, ON :')
		           line_len = getlin (line)
		           if (line_len .eq. 0) then
				word_len = 0
		           else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			   end if
	        	end do
	        	if (word_len .gt. 0) call upcase (word(1:word_len))       

                        if (word(1:2) .eq. ss_option(1)(1:2)) then
                            setdegoff=1
	        	else if (word(1:2) .eq. ss_option(2)(1:2)) then 
                            setdegoff=0
                        else
 		            call remark ('ERROR: @ S D "'//word(1:2)//'" not recognized')
                        end if            

!	Option is OPERATOR STEVENS
		else if (word(1:word_len) .eq. s_option(3)(1:word_len)) then

 	        	word_len = getwrd (line(1:line_len), l, word)
	        	do while (word_len .eq. 0)
                           call remark ('Set stevens-operator-flag to what? Implemented are:') 
                           call prompt ('OFf, ON :')
		           line_len = getlin (line)
		           if (line_len .eq. 0) then
				word_len = 0
		           else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			   end if
	        	end do
	        	if (word_len .gt. 0) call upcase (word(1:word_len))       

                        if (word(1:2) .eq. ss_option(1)(1:2)) then
                            stdhamoff=1
	        	else if (word(1:2) .eq. ss_option(2)(1:2)) then 
                            stdhamoff=0
                        else
 		            call remark ('ERROR: @ S O "'//word(1:2)//'" not recognized')
                        end if            

!	Option is ION R3+
		else if (word(1:word_len) .eq. s_option(5)(1:word_len)) then

 	        	word_len = getwrd (line(1:line_len), l, word)
	        	do while (word_len .eq. 0)
                           call remark ('Set which R3+ ion? Implemented are:') 
                           call prompt ('CE,PR,ND,PM,SM,EU,GD,TB,DY,HO,ER,TM,YB:')
		           line_len = getlin (line)
		           if (line_len .eq. 0) then
				word_len = 0
		           else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			   end if
	        	end do
	        	if (word_len .gt. 0) call upcase (word(1:word_len))       

                        nnre = 0
                        if (word(1:2) .eq. r_option(1)(1:2)) then
                             nnre = 1
	        	else if (word(1:2) .eq. r_option(2)(1:2)) then 
                             nnre = 2
	        	else if (word(1:2) .eq. r_option(3)(1:2)) then 
                             nnre = 3
	        	else if (word(1:2) .eq. r_option(4)(1:2)) then 
                             nnre = 4
	        	else if (word(1:2) .eq. r_option(5)(1:2)) then 
                             nnre = 5
	        	else if (word(1:2) .eq. r_option(6)(1:2)) then 
                             nnre = 6
	        	else if (word(1:2) .eq. r_option(7)(1:2)) then 
                             nnre = 7
	        	else if (word(1:2) .eq. r_option(8)(1:2)) then 
                             nnre = 8
	        	else if (word(1:2) .eq. r_option(9)(1:2)) then 
                             nnre = 9
	        	else if (word(1:2) .eq. r_option(10)(1:2)) then 
                             nnre = 10
	        	else if (word(1:2) .eq. r_option(11)(1:2)) then 
                             nnre = 11
	        	else if (word(1:2) .eq. r_option(12)(1:2)) then 
                             nnre = 12
	        	else if (word(1:2) .eq. r_option(13)(1:2)) then 
                             nnre = 13
                        else
 		            call remark ('ERROR: @ S I "'//word(1:2)//'" not recognized')
                        end if            

                        if(nnre.eq.6 .and. nnre.ne.nre) then
                           call remark('Sorry, you can not set the rare earth ion to Eu3+.') 
			end if
                        if( nnre.ge.1.and.nnre.le.13.and.nnre.ne.nre .and. nnre.ne.6 ) then

                            call remark(' WARNING: By changing to another rare earth ion some ')
                            call remark('          CEF-parameters may get lost.               ')  
 99                         call prompt(' Do you really want to change to another ion  [No=def]?')
                            line_len = getlin(line)
                            c = 'N'
                            if( line_len .gt. 0 ) then
                                c = line(1:1)
                            end if
                            if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 99
                            if(c.eq.'Y'.or.c.eq.'y') then
                                call r_init(symmetry,nnre) 
!                               nre   = nnre   ! already done by r_init
                                call remark(' FOCUS sets the rare earth ion to '//r_option(nre)(1:2)//'3+.')
                            end if
			end if

!	Option is CEF SYMMETRY
		else if (word(1:word_len) .eq. s_option(6)(1:word_len)) then

                        if( nre.eq.6 .or. nre.eq.7) then
			     if(nre.eq.7) call remark ('Sorry, you can not change the CEF symmetry for Gd3+.')
			     if(nre.eq.6) call remark ('Sorry, you can not change the CEF symmetry for Eu3+.')
                             symmetry = 8
                        else 
 	        	  word_len = getwrd (line(1:line_len), l, word)
                          if(word_len.eq.0) call p_symmetry_numbers(out)
	        	  do while (word_len .eq. 0)
                             call remark ( ' Which is the point symmetry you want to set? ')
                             call prompt ( ' Please type in the symmetry number, see table above (0-8)  :')
		             line_len = getlin (line)
		             if (line_len .eq. 0) then
			  	word_len = 0
		             else
			  	l = 1
			  	word_len = getwrd (line(1:line_len), l, word)
			     end if
	        	  end do
      	        	  if (word_len .gt. 0) call upcase (word(1:word_len))       
                        end if

                        ssymmetry= -1
                        if (     word(1:1) .eq. rr_option(1)(1:1)) then
                             ssymmetry = 0
	        	else if (word(1:1) .eq. rr_option(2)(1:1)) then 
                             ssymmetry = 1
	        	else if (word(1:1) .eq. rr_option(3)(1:1)) then 
                             ssymmetry = 2
	        	else if (word(1:1) .eq. rr_option(4)(1:1)) then 
                             ssymmetry = 3
	        	else if (word(1:1) .eq. rr_option(5)(1:1)) then 
                             ssymmetry = 4
	        	else if (word(1:1) .eq. rr_option(6)(1:1)) then 
                             ssymmetry = 5
	        	else if (word(1:1) .eq. rr_option(7)(1:1)) then 
                             ssymmetry = 6
	        	else if (word(1:1) .eq. rr_option(8)(1:1)) then 
                             ssymmetry = 7
	        	else if (word(1:1) .eq. rr_option(9)(1:1)) then 
                             ssymmetry = 8
                        end if            


                        if(ssymmetry.eq.-1 .and. nre.ne.6 .and. nre.ne.7) then
 		            call remark ('ERROR: @ S C "'//word(1:2)//'" not recognized')
                        end if            

                        if( ssymmetry.ge.0.and.ssymmetry.le.8.and.ssymmetry.ne.symmetry .and.
     *                      nre.ne.7 .and. nre.ne.6 ) then

                            call remark(' WARNING: By changing to another CEF symmetry some ')
                            call remark('          CEF-parameters may get lost.               ')  
 98                         call prompt(' Do you really want to change to another symmetry [No=def]?')
                            line_len = getlin(line)
                            c = 'N'
                            if( line_len .gt. 0 ) then
                                c = line(1:1)
                            end if
                            if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 98
                            if(c.eq.'Y'.or.c.eq.'y') then
                                call r_init(ssymmetry,nre) 
!                               symmetry   = ssymmetry   ! already done by r_init
                                call remark(' FOCUS sets the CEF symmetry to '//word(1:1)//'.')
                            end if
			end if

!       if the set command is not regognized
		else                                                                                       
			call remark ('ERROR: @ S "'//word(1:word_len)//'" not recognized')
		end if                                  

!       -------------------------------------------------------------------------------------
!	Command is ALTER
!       -------------------------------------------------------------------------------------
	else if (word(1:word_len) .eq. command(4)(1:word_len)) then
		word_len = getwrd (line(1:line_len), l, word)
		do while (word_len .eq. 0) 
                        call remark ('Alter what? Implemented are:')
			call prompt ('Parameter, Resolution, Selfshielding, Overall splitting:')
			line_len = getlin (line)
			if (line_len .eq. 0) then
				word_len = 0
			else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			end if
		end do
		if (word_len .gt. 0) call upcase (word(1:word_len))

!	Option is Parameter
		if (word(1:word_len) .eq. a_option(1)(1:word_len)) then
 	        	word_len = getwrd (line(1:line_len), l, word)
	        	do while (word_len .eq. 0)
                           call remark ('Alter Parameter to what? Implemented are:') 
                           call prompt ('Bkq, Akq, Vkq, Xr, Wx :')
		           line_len = getlin (line)
		           if (line_len .eq. 0) then
				word_len = 0
		           else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			   end if
	        	end do
	        	if (word_len .gt. 0) call upcase (word(1:word_len))       

                        if(nre.ne.6) then

                          if (word(1:1) .eq. aa_option(1)(1:1)) then
                              call a_to_bkq(setbkq,setvkq,setxr,setakq,setwx)  ! alter to Bkq-parameter   
    	               	  else if (word(1:1) .eq. aa_option(2)(1:1)) then 
                              call a_to_vkq(setbkq,setvkq,setxr,setakq,setwx)  ! alter to Vkq-parameter     
 	        	  else if (word(1:1) .eq. aa_option(3)(1:1)) then 
                              call a_to_xr(setbkq,setvkq,setxr,setakq,setwx)   ! alter to xR-parameter     
    	        	  else if (word(1:1) .eq. aa_option(4)(1:1)) then 
                              call a_to_akq(setbkq,setvkq,setxr,setakq,setwx)  ! alter to Akq-parameter     
    	        	  else if (word(1:1) .eq. aa_option(5)(1:1)) then 
                              call a_to_wx(setbkq,setvkq,setxr,setakq,setwx)  ! alter to Wx-parameter     
                          else
 		              call remark ('ERROR: @ A P "'//word(1:1)//'" not recognized')
                          end if            
                        else
			      call remark ('Sorry, no parameter transformation allowed for Eu3+.')
		        end if
		else if (word(1:word_len) .eq. a_option(2)(1:word_len)) then
                        call a_resolution(setres) 
		else if (word(1:word_len) .eq. a_option(3)(1:word_len)) then
                        call a_selfshielding(setself)
		else if (word(1:word_len) .eq. a_option(4)(1:word_len)) then
                        call a_overall_splitting
		else                                                                                       
			call remark (' ERROR: @ A "'//word(1:word_len)//'" not recognized')
		end if                                  
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is MONTECARLO
	else if (word(1:word_len) .eq. command(6)(1:word_len)) then
		word_len = getwrd (line(1:line_len), l, word)
		do while (word_len .eq. 0) 
                        call remark ('Which Monte-Carlo command do you want to use? Implemented are:')
			call prompt ('New, Display chi**2, Load parameters, Plot parameters :')
			line_len = getlin (line)
			if (line_len .eq. 0) then
				word_len = 0
			else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			end if
		end do
		if (word_len .gt. 0) call upcase (word(1:word_len))
!	Option is NEW 
		if (word(1:word_len) .eq. m_option(1)(1:word_len)) then
                        call monte_carlo (fun, line, l)
!	Option is DISPLAY CHI**2
		else if (word(1:word_len) .eq. m_option(2)(1:word_len)) then
                        if(if_monte.ne.0) then 
		            call dis_monte (fun, line, l)
                        else
                            call remark ('ERROR: Monte Carlo is not yet defined.')
                            call remark ('Define Monte Carlo by @ M N.')
                        end if                                                                                      
!	Option is PLOT PARAMETERS 
		else if (word(1:word_len) .eq. m_option(3)(1:word_len)) then
                        if(if_monte.ne.0) then 
	                    nset = getari (line, l, imc_array, MAX_CHI2)
	                    if (nset .lt. 0) then
		                   call remark ('ERROR: something wrong with the parameter number list.')
		                   return
	                    end if 
                            if(nset.ne.0) then ! test if all are in the chi**2 list
                               if_first = 0
                               nmax=0 ! find out how many are already stored
                               do j=1,MAX_CHI2
                                  call ifnull(s_chisq(j))
                                  if(s_chisq(j).ne.0.0d0) nmax=nmax+1
                               end do
                               do i=1, nset
                                  if(imc_array(i).lt.1 .or. imc_array(i).gt.nmax) then
                                     write(cint, fmt='(i6)') imc_array(i)
                                     call remark('ERROR: the number '//cint//' is not in the Chi**2 list.')
                                     return
                                  end if 
                               end do
                            end if
		            call plot_monte (fun, line, l)
                        else
                            call remark ('ERROR: Monte-Carlo is not yet defined.')
                            call remark ('Define Monte-Carlo by @ M N.')
                        end if                                                                                      
!	Option is LOAD PARAMETERS
		else if (word(1:word_len) .eq. m_option(4)(1:word_len)) then
                        if(if_monte.ne.0) then 
	                    nset = getari (line, l, imc_array, MAX_CHI2)
	                    if(nset .lt. 0) then
		                   call remark ('ERROR: something wrong with the parameter number list.')
		                   return
	                    end if                                                 
                            if(nset. gt. 1) then
                                   call remark ('ERROR: you can load only one number.')
                                   return
                            end if   
                            if(nset.ne.0) then ! test if all are in the chi**2 list
                               nmax=0 ! find out how many are already stored
                               do j=1,MAX_CHI2
                                  call ifnull(s_chisq(j))
                                  if(s_chisq(j).ne.0.0d0) nmax=nmax+1
                               end do
                               do i=1, nset
                                  if(imc_array(i).lt.1 .or. imc_array(i).gt.nmax) then
                                     write(cint, fmt='(i6)') imc_array(i)
                                     call remark('ERROR: the number '//cint//' is not in the Chi**2 list.')
                                     return
                                  end if 
                               end do
                            end if
		            call load_monte (fun, line, l)
                        else
                            call remark ('ERROR: Monte-Carlo is not yet defined.')
                            call remark ('Define Monte-Carlo by @ M N.')
                        end if                                                                                      

		else                                                                                       
			call remark ('ERROR: @ M "'//word(1:word_len)//'" not recognized')
		end if                                  
!-----------------------------------------------------------------------------------------------------------------------------------
!	Command is INFO
	else if (word(1:word_len) .eq. command(5)(1:word_len)) then
  		word_len = getwrd (line(1:line_len), l, word)
		do while (word_len .eq. 0) 
                        call remark ('Info what? Implemented are:')
			call remark ('Version, Scattering function, Lineshapes, Resolution,' )
			call prompt ('Intensities, Cef symmetries, Mu*t                   :' )
			line_len = getlin (line)
			if (line_len .eq. 0) then
				word_len = 0
			else
				l = 1
				word_len = getwrd (line(1:line_len), l, word)
			end if
		end do
		if (word_len .gt. 0) call upcase (word(1:word_len))
!	Option is VERSION
		if (word(1:word_len) .eq. i_option(1)(1:word_len)) then
                        call version_number(focus_version)

!	Option is SCATTERING FUNCTION
		else if (word(1:word_len) .eq. i_option(2)(1:word_len)) then
                        call i_scattering_function

!	Option is LINESHAPES            
		else if (word(1:word_len) .eq. i_option(3)(1:word_len)) then
                        call i_lineshapes           

!	Option is RESOLUTION            
		else if (word(1:word_len) .eq. i_option(4)(1:word_len)) then
                        call i_resolution           

!	Option is INTENSITIES           
		else if (word(1:word_len) .eq. i_option(5)(1:word_len)) then
                        call i_intensities           

!	Option is AUTHORS               
		else if (word(1:word_len) .eq. i_option(6)(1:word_len)) then
                        call i_authors               

!	Option is CEF SYMMETRIES               
		else if (word(1:word_len) .eq. i_option(7)(1:word_len)) then
                        call p_symmetry_numbers(out)               

!	Option is MU*T              
		else if (word(1:word_len) .eq. i_option(8)(1:word_len)) then
                        call i_mut           

		else                                                                                       
			call remark ('ERROR: @ I "'//word(1:word_len)//'" not recognized')
		end if                                  
!-----------------------------------------------------------------------------------------------------------------------------------
!!	Command is unrecognized
	else 
		call remark ('ERROR: @ "'//word(1:word_len)//'" not recognized')
	end if
	end
