!	Subroutine HELP		-	Outputs help information to screen

	subroutine HELP (line, l)

	character*1 letter
	character*8 word
	character*(*) line
	integer*4 l, len, getwrd
	
	equivalence (word, letter)

	len = getwrd (line, l, word)
	if (len .gt. 0) call upcase (word(1:len))
	if (len .eq. 0 .or. letter .eq. 'H') then
		call remark ('Possible commands are:')
		call remark (' Alter     - Alter fitting constants and output characteristics')	     
		call remark (' Bind      - Bind two parameters together in fits')
		call remark (' Clear     - Clear fixed Parameters or modified Data')
		call remark (' Display   - Display Data, Parameters, Calculation, Fit or Matrix')
		call remark (' Exit      - Exit from routine')
		call remark (' Fix       - Fix parameters')
		call remark (' Go        - Fit function to data')	     
		call remark (' Help      - Help provided on each command')
		call remark (' Input     - Input parameters from a file')
		call remark (' Jump      - Issue DCL command or spawn DCL sub-process')	     
		call remark (' Keep      - Store previous plot in a disk file')
		call remark (' Limit     - Limit parameters')
		call remark (' Modify    - Modify data points')
		call remark (' Output    - Output Parameters, Calculation or Fit to a disk file')
		call remark (' Plot      - Plot Data, Calculation, Fit, Residuals or Over')
		call remark (' Remove    - Remove data points')
		call remark (' Set       - Set parameters')
		call remark (' Title     - Set plot title')	     
		call remark (' Use       - Use specified dataset in subsequent commands')
		call remark (' Zero      - Unlimit parameters')
		call remark (' @focus    - FOCUS defined commands')
      		call remark ('Each command may be abbreviated by the first letter')
		call remark ('Type Help followed by the command for further information')       
	else if (letter .eq. 'D') then 
		call remark (' Display Parameters  [<n1,n2-n3,n4...>] - list specified parameters')
		call remark ('         Data        [<Xmin,Xmax>]      - list specified data points')
		call remark ('         Calculation [<Xmin,Xmax>]      - list calculated points')
		call remark ('         Fit                            - list fitted parameters')
		call remark ('         Matrix                         - list fitted correlation matrix')
	else if (letter .eq. 'I') then
		call remark (' Input <filename> - input parameters from a disk file')
		call remark ('                    in the format produced by the OUTPUT command')
	else if (letter .eq. 'O') then
		call remark (' Output Parameters  <filename> - output parameters to a disk file')
		call remark ('        Fit         <filename> - output fitted parameters to a disk file')
		call remark ('        Calculation <filename> - output data points and calculated values to a tab-delimited file')
		call remark (' Default filenames are SYS$SCRATCH:PAR.OUT, FIT.OUT and CAL.OUT respectively')
	else if (letter .eq. 'U') then
		call remark (' Use <dataset>   - Use specified dataset in subsequent plots etc')  
	else if (letter .eq. 'P') then
		call remark (' Plot [Over] Data [<Xmin,Xmax>,[<Ymin,Ymax>]] - plot data points')
		call remark ('             Calculation        "             - plot calculated points')
		call remark ('             Fit                "             - plot fitted points')
		call remark ('             Residuals          "             - plot residuals of fit')
		call remark ('             <n1,n2-n3,n4...>   "             - plot requested components')
		call remark (' Plot is in current box if command given with "Over" option')
	else if (letter .eq. 'K') then
		call remark (' Keep <filename> - Store previous plot in a disk file ')
		call remark ('                   (Default is SYS$SCRATCH:FRILLS.PS)')
	else if (letter .eq. 'M') then
		call remark (' Modify <n1,n2-n3,n4...> - modify data points')
	else if (letter .eq. 'R') then
		call remark (' Remove <X1,X2>          - remove data points between X limits')
	else if (letter .eq. 'S') then
		call remark (' Set <n1,n2-n3,n4...>    - set specified parameters')
	else if (letter .eq. 'F') then
		call remark (' Fix <n1,n2-n3,n4...>    - fix specified parameters')
	else if (letter .eq. 'B') then
		call remark (' Bind <n1,n2-n3,n4...>   - bind specified parameters')
	else if (letter .eq. 'C') then
		call remark (' Clear Data [X1,X2] - clear modified or removed data points between X limits')
		call remark ('       Parameter [n1,n2-n3,n4...] - clear fixed or bound parameters')
	else if (letter .eq. 'L') then
		call remark (' Limit <n1,n2-n3,n4...>  -  limit specified parameters')
	else if (letter .eq. 'Z') then
		call remark (' Zero  <n1,n2-n3,n4...>  -  unlimit specified parameters')
	else if (letter .eq. 'A') then
		call remark (' Alter  Fits            - Alter least squares fitting constants after prompts')	     
		call remark ('        Output <option> - Alter terminal output.  Type A O for list of options')     
		call remark ('        Markers <mk_type> <mk_size> <mk_coli)')
		call remark ('                        - Alter marker types, size and colour index')	     
		call remark ('        Lines   <ln_type> <ln_size> <ln_coli)')
		call remark ('                        - Alter line types, thickness and colour index')	     
		call remark ('        Smoothing       - Alter smoothing options after prompts')	     
		call remark ('        Workstation     - Alter GKS terminal type after prompts')	     
		call remark ('        Hardcopy        - Alter hardcopy output paper size after prompts')	     
	else if (letter .eq. 'T') then
		call remark (' Title - set plot title after prompt')
	else if (letter .eq. 'J') then
		call remark (' Jump [<DCL command>] - spawn DCL sub-process.')
		call remark ('      If no command is given, LOGOUT to return to FRILLS')
	else if (letter .eq. 'G') then
		call remark (' Go - fit parameters to data')
	else if (letter .eq. 'E') then
		call remark (' Exit - exit from routine')

	else if (letter .eq. '@') then
		call USER_HELP (line, l)

	end if	     

	return
	end
