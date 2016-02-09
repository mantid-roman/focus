      program focus
!  ---------------------------------------------------------------------    
! |                        module FOCUS.FOR                             |
! | a program for Fitting Of Crystal field parameters interactively     |
! | by Using neutron spectra                                            |
!  --------------------------------------------------------------------- 
! | see Frills: an interactive Least-Squares Fitting Package RAL-91-011 |
! | and multi_frills will be published soon                             |
!  ---------------------------------------------------------------------
! | this file is written by: Dr P O Fabi                                |
! |                          ISIS Facility                              |
! |                          Rutherford Appleton Laboratory             |
! |                          Chilton,Didcot,Oxfordshire                 |
! |                          OX11 0QX                                   | 
!  ---------------------------------------------------------------------
!
! include header 
     implicit none
      include 'CF_SOURCES:FRILLS.INC'
!-----------------------------------------------------------------------
! define common blocks
!-----------------------------------------------------------------------
      include 'CF_SOURCES:FOCUS.INC' 
      include 'CF_SOURCES:LIMIT.INC'
      include 'CF_SOURCES:LIMITPAR.INC'
      include 'CF_SOURCES:MONTE.INC' 
      include 'CF_SOURCES:XR.INC' 
      include 'CF_SOURCES:CONTROL.INC'
      include 'CF_SOURCES:DEGENERATION.INC'
      include 'CF_SOURCES:FORMFACTOR.INC'
      include 'CF_SOURCES:CANTED.INC'
      include 'CF_SOURCES:RNUMBER.INC'
      include 'CF_SOURCES:NO0PAR.INC'
      include 'CF_SOURCES:RESOLUTION.INC'
      include 'CF_SOURCES:SELF.INC'
      include 'CF_SOURCES:INSTRUMENT.INC'
      include 'CF_SOURCES:OVERALL.INC'

!-----------------------------------------------------------------------
! the name of the fitting subroutine
!-----------------------------------------------------------------------
      external CFFIT
!-----------------------------------------------------------------------
! Frills Parameter names
!-----------------------------------------------------------------------
!     common /fr_nams/ nptot,nam(50)
!     integer*4    nptot ! total no. of parameters
!     character*20 nam   ! names of parameter    
!-----------------------------------------------------------------------
! setting instrument parameter
! they are stored in INSTRUMENT.INC
!-----------------------------------------------------------------------
      data lmc_het   /10.000d0/, 
     >     lcs_het   / 1.820d0/,
     >     lsd_het   / 2.512d0/,
     >     lsd_het_h / 4.047d0/

      data lmc_mari  /10.000d0/,
     >     lcs_mari  / 1.739d0/,
     >     lsd_mari  / 4.042d0/ 
!-----------------------------------------------------------------------
! setting Frills parameter
! they are stored in FRILLS.INC
! it is possible to add a new parameter here
! -change in no0par.inc the number from 17 to 18
! -add the new name to the data nam below
! -make in the subroutine CFFIT the new parameter visible
!  by using it
!-----------------------------------------------------------------------
! setting frills parameter names
      data nam(1) /'Temperature (Kelvin)'/,ins(1) /1/,
     >     nam(2) /'Incoming Energy(meV)'/,ins(2) /1/,
     >     nam(3) /'Detector angle      '/,ins(3) /1/,
     >     nam(4) /'Bkgd Const       (b)'/,ins(4) /1/,
     >     nam(5) /'Bkgd slope   (b/meV)'/,ins(5) /1/,
     >     nam(6) /'Elastic Int      (b)'/,ins(6) /1/,
     >     nam(7) /'Elastic pos    (meV)'/,ins(7) /1/,
     >     nam(8) /'Elastic   FWHM (meV)'/,ins(8) /1/,
     >     nam(9) /'Quasielas FWHM (meV)'/,ins(9) /1/,
     >     nam(10)/'Inelastic FWHM (meV)'/,ins(10)/1/,
     >     nam(11)/'Intensity factor    '/,ins(11)/1/, 
     >     nam(12)/'B_ext(x)     (Tesla)'/,ins(12)/1/,   
     >     nam(13)/'B_ext(y)     (Tesla)'/,ins(13)/1/,   
     >     nam(14)/'B_ext(z)     (Tesla)'/,ins(14)/1/,   
     >     nam(15)/'B_mol(x)     (Tesla)'/,ins(15)/1/,   
     >     nam(16)/'B_mol(y)     (Tesla)'/,ins(16)/1/,   
     >     nam(17)/'B_mol(z)     (Tesla)'/,ins(17)/1/    

!     data cno_0_par/17/ ! 17 crystal field independent parameter
      data ns /1/        ! 1  spectrum to fit 
!--------------------------------------------------------------------------
! declarations needed for rare-earth and symmetry determination
!--------------------------------------------------------------------------
      integer no_0_par ! function gives _no_0_par back 

      integer out      ! messages a written to file out
      integer in       ! reading is from file in
      parameter(in=5,out=5)
      integer no_cf_par! function defines ncftot
      integer ncftot   ! total number of crystal field parameter
      integer n0       ! total number of non crystal field parameter
 
! the variables nre and symmetry are already declared in focus.inc
!     integer nre      ! number of the rare-earth ion R3+
!     integer symmetry ! number of symmetry which classifies the point
!                      ! symmetry of the R3+ ion
 
      external no_cf_par          ! to be found in cf_fabi.for
 
      external p_input            ! to be found in cf_fabi.for
      external p_rare_earth       ! to be found in cf_fabi.for
      external p_symmetry_numbers ! to be found in cf_fabi.for
      external p_hamiltonian      ! to be found in cf_fabi.for
      external p_energy_values    ! to be found in cf_fabi.for
      external p_wave_functions   ! to be found in cf_fabi.for 
      external p_matrix_elements  ! to be found in cf_fabi.for
      external p_occupation_factor! to be found in cf_fabi.for
      external p_intensities      ! to be found in cf_fabi.for
      external p_moments          ! to be found in cf_fabi.for
      external p_excitations      ! to be found in cf_fabi.for
 
      external transform             ! cf_fabi.for
      external t_orthonormalisation  ! cf_fabi.for
      external c_crystal_field       ! cf_fabi.for
      external c_excitations         ! cf_fabi.for
!--------------------------------------------------------------------------
! Definitions used to plot out crystal field data
!--------------------------------------------------------------------------
 	real*8 pv(MAX_PAR)
 	real*8 mx,my,mz  
 	real*8 c_x_moment,c_y_moment,c_z_moment   
 	real*8 c_fmevkelvin     
 	integer k,j,i
        character c,c5*5,c6*6,c7*7
        real*8 pi
 
 	external c_x_moment,c_y_moment,c_z_moment  ! cf_fabi.for
 	external c_fmevkelvin                      ! cf_fabi.for
        external initform                          ! cf_formfactor.for
!--------------------------------------------------------------------------
! settings for the ctoi and ctor routines
!--------------------------------------------------------------------------
      character line*80
      integer*4 l, line_len, getlin, ctoi
      real*4    ctor
!--------------------------------------------------------------------------
! for de estimation in case of a resolution calculation
!--------------------------------------------------------------------------
      real*8 sigma_chop, gamma_elas
!--------------------------------------------------------------------------
! set default values
!--------------------------------------------------------------------------
      ns       = 1 ! one data set to fit

      up_limit = 4.0d0*atan(1.0d0) ! upper limit for transformed limited parameters 
      low_limit= 0.0d0             ! lower limit for transformed limited parameters 

      de_default = 0.005d0                   ! in meV
      di_default = 0.01d0                    ! in barn
      de         = de_default*c_fmevkelvin() ! in Kelvin
      di         = di_default                ! in barn

      setbkq    = 1   ! use Bkq parameters
      setakq    = 0   ! do not use the  Akq-parameters
      setvkq    = 0   ! do not use the normalised Vkq-parameters
      setwx     = 0   ! do not use the Lea, Leask Wolff parameter W, x
      setxr     = 0   ! do not use the spherical xR-parameters
      setdegoff = 0   ! degenerated levels are combined
      setres    = 0   ! no resolution function is calculated
      setself   = 0   ! no selfshielding and absorption is calculated
      setfit    = 1   ! FOCUS calculates crystal field properties and fits the data
      setfirst  = 1   ! the next call of r_init is the first call
      setprint  = 1   ! prints a heading in subroutine p_excitations
      stdhamoff = 0   ! use Stevens operators for crystal field calculation
      iscalc    = 0   ! the subroutine CFFIT has not yet been called
      if_monte  = 0   ! Monte Carlo has not been called
      setmonte  = 0   ! chisq_calc() is not called by subroutine monte_carlo 
      if_plot   = 0   ! No Monte Carlo plot mode
      if_first  = 0   ! no first call of Monte Carlo plot 
      cf_out    = 42  ! channel of output	
      min_x     =   0.0 ! in plot mode 
      max_x     =  67.0 ! in plot mode 
      min_y     =   0.0 ! in plot mode 
      max_y     = 300.0 ! in plot mode 

      setnowarnings = 0   ! FOCUS is not in fitting mode
      setquick    = 0     ! No quick startup for FOCUS environment
      nre         = 7     ! choose Gd by default
      symmetry    = 8     ! choose cubic symmetry be default
      setoverall  = 0     ! do not fix the overall splitting
      setpardis   = 1     ! subroutine pardis/fitdis have been called
      ratioenergy = 1.0d0 ! the ratio between the desired and the calculated energy level is 1.0
      levelnumber = 2     ! the fixed energy level is by default the first excited state
!                         ! the ground state has level number 1
      levelenergy = 1.0d0 ! the desired energy for the first excited level is 1.0meV by default
!     -----------------------------------------------------------------------------
!     starting random number for random number generator
      call x05aaf(itime)
      r0_number     = -(itime(7)*10000+itime(6)*100 +itime(5))             
!     -----------------------------------------------------------------------------

      do i=1,MAX_PAR      ! define all parameter to
         p(i) = 0.0d0     ! be zero
      end do             

      do i=1,MAX_PAR      ! define all parameter to
         fixold(i) = i    ! be free parameter for
      end do              ! fitting
 
      do i=0,n_lineshape  ! set all lineshapes of the excitations 
  	lineshape(i) = 1  ! to their default value 1=(gauss)     
      end do
 
      call initform	  ! initialise the r3+ formfactors  
 
!-----------------------------------------------------------------------
! give a welcome message
!-----------------------------------------------------------------------
      call remark ( '            -------------------------------------------------- ')
      call remark ( '           |                    Welcome to                    |')
      call remark ( '           |         *****   ***    ***   *   *   ***         |')     
      call remark ( '           |         *      *   *  *      *   *  *            |')   
      call remark ( '           |         ***    *   *  *      *   *   ***         |')   
      call remark ( '           |         *      *   *  *      *   *      *        |')      
      call remark ( '           |         *       ***    ***    ***    ***         |')
      call remark ( '           |      An Interactive Crystal Electric Field       |')
      call remark ( '           |         Parameter Fitting Package using          |')
      call remark ( '           |             Neutron Scattering Data              |')
      call remark ( '           |                   Peter Fabi                     |')
      call remark ( '            -------------------------------------------------- ')


      call section (out, ' CHANGES TO FOCUS VERSION 1.0')
      call remark  ( ' There is a manual for FOCUS version 1.0 availuable. The manual has been ')
      call remark  ( ' published as a proper Rutherford Appleton report: RAL-TR-95-023. Ask the')
      call remark  ( ' RAL  library for a copy. Also the manual will be a good introduction to ')
      call remark  ( ' the FOCUS program, FOCUS has extended its possibilities since version 1.0.') 
 29   call prompt  ( ' Do you want to be informed about the changes to FOCUS version 1.0 [No]?')
      line_len = getlin(line)
      c = 'N'
      if( line_len .gt. 0 ) then
          c = line(1:1)
      end if
      if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 29
      if(c.eq.'Y'.or.c.eq.'y') call i_changes

      call section (out, ' QUICK STARTUP FOR FOCUS ENVIRONMENT')
      call remark  ( ' If you are not an experienced user of FOCUS hit the <return>-key only.')
      call remark  ( ' If you have stored a parameter file by  OUT P FILENAME  in a previous ')
      call remark  ( ' FOCUS  session, type  Y and hit the <return>-key afterwards. Now load ')
      call remark  ( ' your stored FOCUS session by typing IN FILENAME. ')
  30  call prompt  ( ' Do you want a quick startup to the FOCUS environment......[No=default]?')
      line_len = getlin(line)
      c = 'N'
      if( line_len .gt. 0 ) then
          c = line(1:1)
      end if
      if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 30
      if(c.eq.'Y'.or.c.eq.'y') setquick = 1

      
      if( setquick.ne.1 ) then
        call section (out, ' RESOLUTION FUNCTION ')
        call set_resolution(setres)
        call section (out, ' SELFSHIELDING AND ABSORPTION CORRECTION')
        call remark ( ' FOCUS has a build in calculation to correct for selfshielding and absorption')
        call remark ( ' effects in the sample. The correction  is for a plate normal to the incident')
        call remark ( ' neutron beam. ')
  40    call prompt ( ' Do you want a selfshielding and absorption calculation.....[No=default]?')
        line_len = getlin(line)
        c = 'N'
        if( line_len .gt. 0 ) then
           c = line(1:1)
        end if
        if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 40
        if(c.eq.'Y'.or.c.eq.'y') then
          setself= 1
          muts = 0.0d0
          call remark ( ' Press  @ I M  to inform you about the meaning of mu*t.')
  41      call prompt ( ' Enter mu*t (scattering).....[0=default]:')
          line_len = getlin(line)
          if(line_len .gt. 0) then
             l    = 1
             muts = ctor(line(1:line_len),l)
             if(muts.lt.0.0d0) then
                call remark ( ' ERROR: mu*t can not be negativ! ')
                goto 41
             end if
          end if
          muta = 0.0d0
  42      call prompt ( ' Enter mu*t (absorption).....[0=default]:')
          line_len = getlin(line)
          if(line_len .gt. 0) then
             l    = 1
             muta = ctor(line(1:line_len),l)
             if(muta.lt.0.0d0) then
                call remark ( ' ERROR: mu*t can not be negativ! ')
                goto 42
             end if
          end if
        end if
      end if


      if(setquick.eq.1) goto 11
      write(*,*)
      call p_rare_earth(out)
  10  call remark ( ' For which rare-earth element is the calculation for ? Please type in the ')
      call prompt ( ' rare earth number, see table above (1-13) :')        
      line_len = getlin(line)
      if(line_len .gt. 0) then
         l   = 1
         nre = ctoi(line(1:line_len),l)
         if(nre.lt.1.or.nre.gt.13) goto 10 
      else
         goto 10
      end if
      call remark ( ' ' )
 
  11  if(nre.eq.6 .or. nre.eq.7) then
         symmetry = 8
      else 
         call p_symmetry_numbers(out)
  20     call remark ( ' What is the point symmetry of the rare earth ion ? Please type in the ')
         call prompt ( ' symmetry number, see table above (0-8)  :')
         line_len = getlin(line)
         if(line_len .gt. 0) then
            l   = 1
            symmetry = ctoi(line(1:line_len),l)
            if(symmetry.lt.0.or.symmetry.gt.8) goto 20
         else
            goto 20
         end if
      end if
      write(*,*) 
 


      if(setquick.eq.1) goto 25
      call section (out, ' THE DIRECTION OF THE RARE EARTH MOMENT ')
  21  call prompt ( ' Shows the rare earth moment in the +c-direction [Y]?')  
      line_len = getlin(line)
      c = 'y'
      if( line_len .gt. 0 ) then
          c = line(1:1)
      end if
      if(c.ne.'Y'.and.c.ne.'y'.and.c.ne.'n'.and.c.ne.'N') goto 21
      if(c.eq.'N'.or.c.eq.'n') then

         call remark ( ' Please give in now the Euler-angles (in units of degrees).')
  22     call prompt ( ' Euler-angle ALPHA [0=default]:')
         line_len = getlin(line)
         if(line_len .gt. 0) then
            l   = 1
            alpha_euler = ctor(line(1:line_len),l)
            if(alpha_euler.lt.-360.0d0 .or. alpha_euler.gt.360.0d0) goto 22
         end if

  23     call prompt ( ' Euler-angle BETA  [0=default]:')
         line_len = getlin(line)
         if(line_len .gt. 0) then
            l   = 1
            beta_euler = ctor(line(1:line_len),l)
            if(beta_euler.lt.-360.0d0 .or. beta_euler.gt.360.0d0) goto 23
         end if

  24     call prompt ( ' Euler-angle GAMMA [0=default]:')
         line_len = getlin(line)
         if(line_len .gt. 0) then
            l   = 1
            gamma_euler = ctor(line(1:line_len),l)
            if(gamma_euler.lt.-360.0d0 .or. gamma_euler.gt.360.0d0) goto 24
         end if

         pi = 4.0d0*atan(1.0d0)  
         alpha_euler = alpha_euler * pi/180.0
         beta_euler  = beta_euler  * pi/180.0
         gamma_euler = gamma_euler * pi/180.0
 
         stdhamoff = 1  ! use 3j-symbols for crystal field calculation
      else
         alpha_euler = 0.0d0
         beta_euler  = 0.0d0
         gamma_euler = 0.0d0
      end if
 
     
 25   call r_init(symmetry,nre) ! sets the R3+
      setfirst = 0  ! the next r_init call is not the first call
      setprint = 0  ! no heading in subroutine p_excitations

	

      write(*,*) 
      call remark (' The CONTROL will now be given to MULTI_FRILLS. ')
      write(*,*)  
      call remark (' Type H   <return> for a list of the available MULTI_FRILLS COMMANDS.')
      call remark (' Type @   <return> for a list of the available FOCUS COMMANDS.') 
      call remark (' Type D P <return> to DISPLAY the fitting PARAMETERS.')
      call remark (' Type S P <return> to SET     the fitting PARAMETERS.')
      call remark (' Type GO  <return> to fit the data.')
      call remark (' Type EX  <return> to exit.')
      write(*,*)
      call remark (' Details of the fit           are stored in SYS$SCRATCH:FOCUS.LPT')  
      call remark (' Details of the crystal field are stored in SYS$SCRATCH:FOCUS.CEF')
      call remark (' Please wait ...')
 
      open(unit=41, file='SYS$SCRATCH:FOCUS.LPT', status='NEW')
      open(unit=cf_out, file='SYS$SCRATCH:FOCUS.CEF', status='NEW') 
 
      call multi_frills_genie_in     
      xcap(1) = 'Energy Transfer E in meV'
      ycap(1) = 'S(phi,E) in mb per sr.meV.R-atom'
      xlab(1) = 'E'
      ylab(1) = 'S(PE)'
 
      call multi_frills (CFFIT)
      call multi_frills_genie_out (CFFIT)
 
      close(unit=41)
 
!---------------------------------------------------------------------------------------
! Prepare output file FOCUS.CEF
!---------------------------------------------------------------------------------------
! transform the actual parameters to the Bkq-parameter
!---------------------------------------------------------------------------------------
      call a_to_bkq(setbkq,setvkq,setxr,setakq,setwx)     

!---------------------------------------------------------------------------------------
! degenerated levels are combined
!---------------------------------------------------------------------------------------
      setdegoff = 0  

!---------------------------------------------------------------------------------------
! set the temperature
!---------------------------------------------------------------------------------------
      temp    = p(1)

!---------------------------------------------------------------------------------------
! set the magnetic fields
!---------------------------------------------------------------------------------------
      n0 = no_0_par()
      bext(1) = p(n0-5)
      bext(2) = p(n0-4)
      bext(3) = p(n0-3)
      bmol(1) = p(n0-2)
      bmol(2) = p(n0-1)
      bmol(3) = p(n0  )
 
!---------------------------------------------------------------------------------------
! delete all Bkq parameter 
!---------------------------------------------------------------------------------------
        do k=0,6
 	   do j=0,6
     	       bkq(k,j) = 0.0d0
           end do
        end do                                           
 
!---------------------------------------------------------------------------------------
! find out how many bkq parameters there are
!---------------------------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)
        nptot  = n0 + ncftot 
  	do i=1,ncftot
           pv(i) = p(n0+i)
 	end do
        call transform(pv,bkq,symmetry,nre)
 
!---------------------------------------------------------------------------------------
! transform the Bkq from the meV scale to the Kelvin scale
!---------------------------------------------------------------------------------------
        do k=2,6,2
 	   do j=0,6
     	       bkq(k,j) = bkq(k,j)*c_fmevkelvin()  
           end do
        end do 
  
!---------------------------------------------------------------------------------------
! prints an overview over the symmetry numbers
!---------------------------------------------------------------------------------------
    	call p_symmetry_numbers(cf_out)

!---------------------------------------------------------------------------------------
! prints an overview over the hamiltonian
!---------------------------------------------------------------------------------------
   	call p_hamiltonian(cf_out,symmetry)

!---------------------------------------------------------------------------------------
! prints the input out 
!---------------------------------------------------------------------------------------
  	call p_input(cf_out,symmetry,nre,bkq,sbkq,bmol,bext,temp)

!---------------------------------------------------------------------------------------
! calculate the crystal field splitting
!---------------------------------------------------------------------------------------
   	call c_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,'Bkq')

!---------------------------------------------------------------------------------------
! prints out the transition energies in meV
!---------------------------------------------------------------------------------------
     	call p_energy_values(cf_out,energy,dimj,fmevkelvin)

!---------------------------------------------------------------------------------------
! prints the wave functions
!---------------------------------------------------------------------------------------
     	call p_wave_functions(cf_out,energy,wavefunction,dimj)

!---------------------------------------------------------------------------------------
! tests if the wavefunctions are really orthonormal
!---------------------------------------------------------------------------------------
c     	call t_orthonormalisation(cf_out,wavefunction,dimj)

!---------------------------------------------------------------------------------------
! prints out the transition matrix elements
!---------------------------------------------------------------------------------------
    	call p_matrix_elements(cf_out,dimj,energy,jx2mat,jy2mat,jz2mat,jt2mat)

!---------------------------------------------------------------------------------------
! prints 'occupation  factor'
!---------------------------------------------------------------------------------------
     	call p_occupation_factor(cf_out,occupation_factor,temp)

!---------------------------------------------------------------------------------------
! prints out the transition intensities for the given temperature temp
!---------------------------------------------------------------------------------------
     	call p_intensities(cf_out,dimj,energy,intensity,temp)

!---------------------------------------------------------------------------------------
! calculates the R3+ moment in myb
!---------------------------------------------------------------------------------------
     	mx = c_x_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
    	my = c_y_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
    	mz = c_z_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)

!---------------------------------------------------------------------------------------
! prints out the R3+ moment 
!---------------------------------------------------------------------------------------
    	call p_moments(cf_out,mx,my,mz,gj,dimj) 

!---------------------------------------------------------------------------------------
! calculates the excitations for a given temperature
!---------------------------------------------------------------------------------------
        call c_excitations(cf_out,nre,symmetry,bkq,sbkq,bmol,bext,temp,
     *                    e_excitation,i_excitation,n_excitation,'Bkq',de,di)   

!---------------------------------------------------------------------------------------
! prints out the calculated excitations for a given temperature temp
!---------------------------------------------------------------------------------------
  	call p_excitations(cf_out,de,di)
                                      
        close(cf_out)

      write(*,*)
      call remark (' Details of the fit           are stored in SYS$SCRATCH:FOCUS.LPT')  
      call remark (' Details of the crystal field are stored in SYS$SCRATCH:FOCUS.CEF')
      call remark (' FOCUS wish you a nice day.' )
      end
!-----------------------------------------------------
! delete all bkq-parameter which do not belong
! to the symmetry
!-----------------------------------------------------
      subroutine del_bkq(bkq,symmetry,nre)
      	integer nre,symmetry
	complex*16 bkq(0:6,0:6),null

 	integer j2,k,q
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6
 	external i_re_constants      

 	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
 	j2   = dimj-1                                                         
        null = cmplx(0.0d0,0.0d0)

        if(j2.lt.6 .or. nre.eq.7) then
	   do q=0,6
              bkq(6,q) = null
	   end do
	end if

        if(j2.lt.4 .or. nre.eq.7) then
	   do q=0,4
              bkq(4,q) = null
	   end do
	end if

        if(j2.lt.2 .or. nre.eq.7) then
	   do q=0,2
              bkq(2,q) = null
	   end do
	end if

        if(j2.lt.2 .or. nre.eq.7) return
!       ----------------------------------------------
!       from here j2.ge.2 .and. nre.ne.7
!       ----------------------------------------------
        if(symmetry.eq.8) bkq(2,0)=null
        if(symmetry.ge.1) bkq(2,1)=null
        if(symmetry.ge.3) bkq(2,2)=null

	if(j2.lt.4) return
!       ----------------------------------------------
!       from here j2.ge.4
!       ----------------------------------------------
        if(symmetry.ge.1                    ) bkq(4,1)=null
        if(symmetry.ge.3                    ) bkq(4,2)=null
        if(symmetry.ge.1 .and. symmetry.le.4) bkq(4,3)=null
        if(symmetry.ge.7                    ) bkq(4,3)=null
        if(symmetry.ge.5 .and. symmetry.le.7) bkq(4,4)=null

	if(j2.lt.6) return
!       ----------------------------------------------
!       from here j2.ge.6
!       ----------------------------------------------
        if(symmetry.ge.1                    ) bkq(6,1)=null
        if(symmetry.ge.3                    ) bkq(6,2)=null
        if(symmetry.ge.1 .and. symmetry.le.4) bkq(6,3)=null
        if(symmetry.ge.7                    ) bkq(6,3)=null
        if(symmetry.ge.5 .and. symmetry.le.7) bkq(6,4)=null
        if(symmetry.ge.1                    ) bkq(6,5)=null
        if(symmetry.ge.3 .and. symmetry.le.4) bkq(6,6)=null
        if(symmetry.ge.8                    ) bkq(6,6)=null

      return 
      end 
!-----------------------------------------------------
      subroutine del_lim(symmetry,nre)
      	integer nre,symmetry

 	integer  j2,k,q
 	real*8   gj,dimj,alphaj,betaj,gammaj,r2,r4,r6
        real*8   null
 	external i_re_constants      

 	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
 	j2   = dimj-1                                                         

        if(j2.lt.6 .or. nre.eq.7) then
	   do q=0,6
              call set0(6,q)
	   end do
	end if

        if(j2.lt.4 .or. nre.eq.7) then
	   do q=0,4
              call set0(4,q)
	   end do
	end if

        if(j2.lt.2 .or. nre.eq.7) then
	   do q=0,2
              call set0(2,q)
	   end do
	end if

        if(j2.lt.2 .or. nre.eq.7) return
!       ----------------------------------------------
!       from here j2.ge.2 .and. nre.ne.7
!       ----------------------------------------------
        if(symmetry.eq.8) call set0(2,0)
        if(symmetry.ge.1) call set0(2,1)
        if(symmetry.ge.3) call set0(2,2)

	if(j2.lt.4) return
!       ----------------------------------------------
!       from here j2.ge.4
!       ----------------------------------------------
        if(symmetry.ge.1                    ) call set0(4,1)
        if(symmetry.ge.3                    ) call set0(4,2)
        if(symmetry.ge.1 .and. symmetry.le.4) call set0(4,3)
        if(symmetry.ge.7                    ) call set0(4,3)
        if(symmetry.ge.5 .and. symmetry.le.7) call set0(4,4)

	if(j2.lt.6) return
!       ----------------------------------------------
!       from here j2.ge.6
!       ----------------------------------------------
        if(symmetry.ge.1                    ) call set0(6,1)
        if(symmetry.ge.3                    ) call set0(6,2)
        if(symmetry.ge.1 .and. symmetry.le.4) call set0(6,3)
        if(symmetry.ge.7                    ) call set0(6,3)
        if(symmetry.ge.5 .and. symmetry.le.7) call set0(6,4)
        if(symmetry.ge.1                    ) call set0(6,5)
        if(symmetry.ge.3 .and. symmetry.le.4) call set0(6,6)
        if(symmetry.ge.8                    ) call set0(6,6)

      return 
      end 
!-----------------------------------------------------
! define which part of the stevenoperator to be chosen
!-----------------------------------------------------
      subroutine s_sbkq(in,out,nre,symmetry,sbkq)
 	implicit none
 	integer in,out,nre,symmetry,no,k,q
 	real*8 sbkq(0:6,0:6)
 	integer j2
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6
 	external i_re_constants      

        integer*4    l, line_len, getlin, ctoi
        character*80 line
 
! choose real part of the steven's operators as default
 	do k=2,6,2
           do q=0,6
              sbkq(k,q) = 1
 	   end do
  	end do
  
        if(nre.eq.7 .or. nre.eq.6) return

 	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
 	j2 = dimj-1                                                         
 
      if(symmetry.eq.0) then
        if(j2.ge.2) then
  1        call prompt (' B(2,1)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0 
           end if
           if(no.ne.0.and.no.ne.1) goto 1                                            
           if(no.eq.1) sbkq(2,1)=-1
        endif
      else if(symmetry.eq.1) then
        if(j2.ge.2) then
  2	   call prompt (' B(2,2)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 2   
 	   if(no.eq.1) sbkq(2,2)=-1                                             
        endif
      else if(symmetry.eq.2) then
        if(j2.ge.2) then
  3        call prompt (' B(2,2)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 3   
 	   if(no.eq.1) sbkq(2,2)=-1 
        endif
 	if(j2.ge.4) then                                    
  4	   call prompt (' B(4,2)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 4   
 	   if(no.eq.1) sbkq(4,2)=-1 
 
  5	   call prompt (' B(4,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0 
           end if
           if(no.ne.0.and.no.ne.1) goto 5   
 	   if(no.eq.1) sbkq(4,4)=-1                                                    
        endif
 	if(j2.ge.6) then                                    
  6        call prompt (' B(6,2)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 6  
 	   if(no.eq.1) sbkq(6,2)=-1 
 
  7	   call prompt (' B(6,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0 
           end if
           if(no.ne.0.and.no.ne.1) goto 7   
 	   if(no.eq.1) sbkq(6,4)=-1                                                    
 
  8	   call prompt (' B(6,6)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 8   
 	   if(no.eq.1) sbkq(6,6)=-1 
        endif                                                                                                         
      else if(symmetry.eq.3) then
        if(j2.ge.4) then
  9	   call prompt (' B(4,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 9   
 	   if(no.eq.1) sbkq(4,4)=-1    
        endif                              
      else if(symmetry.eq.4) then
        if(j2.ge.4) then
 10	   call prompt (' B(4,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 10  
 	   if(no.eq.1) sbkq(4,4)=-1                                   
        endif
 	if(j2.ge.6) then
 11	   call prompt (' B(6,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 11  
 	   if(no.eq.1) sbkq(6,4)=-1                                    
        endif                                                                            
      else if(symmetry.eq.5) then
        if(j2.ge.4) then
 12	   call prompt (' B(4,3)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 12  
 	   if(no.eq.1) sbkq(4,3)=-1                                   
        endif
      else if(symmetry.eq.6) then
        if(j2.ge.4) then
 13	   call prompt (' B(4,3)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 13  
 	   if(no.eq.1) sbkq(4,3)=-1                                   
        endif
 	if(j2.ge.6) then
 14	   call prompt (' B(6,3)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 14  
 	   if(no.eq.1) sbkq(6,3)=-1                                           
     
 15	   call prompt (' B(6,6)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 15  
 	   if(no.eq.1) sbkq(6,6)=-1                                           
        endif                                                                           
      else if(symmetry.eq.7) then
        if(j2.ge.6) then
 16	   call prompt (' B(6,6)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 16  
 	   if(no.eq.1) sbkq(6,6)=-1                                           
        endif                                                                   
      else if(symmetry.eq.8 .and. nre.ne.6 .and. nre.ne.7 ) then
        if(j2.ge.4) then
  17	   call prompt (' B(4,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 17 
 	   if(no.eq.1) sbkq(4,4)=-1                                   
        endif
 	if(j2.ge.6) then
  18	   call prompt (' B(6,4)  real (0) or imaginary (1) part [0=default] :')
           line_len = getlin(line)
           if(line_len .gt. 0) then
              l  = 1
              no = ctoi(line(1:line_len),l)
           else
              no = 0
           end if
           if(no.ne.0.and.no.ne.1) goto 18
 	   if(no.eq.1) sbkq(6,4)=-1                                    
        endif
      endif                                       
      call remark (' ')	
      return
      end
!---------------------------------------------------
! defines the parameter names nam for each symmetry 
!---------------------------------------------------
      subroutine d_pnames(nre,symmetry,sbkq,nam,n0)
	implicit none
	integer nre,symmetry,n0
	character nam(*)*(*)
	real*8 sbkq(0:6,0:6)
	integer j2
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6      
	external i_re_constants

	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
	j2 = dimj-1                                                         
                                                                              

        if(     symmetry.eq.0) then
	    if(j2.ge.2) then
               nam(n0+ 1)='Re B20       (meV)  '
               nam(n0+ 2)='Re B21       (meV)  '
               if(sbkq(2,1).eq.-1) nam(n0+ 2)='Im B21       (meV)  '   
               nam(n0+ 3)='Re B22       (meV)  '
               nam(n0+ 4)='Im B22       (meV)  '
	    endif
	    if(j2.ge.4) then   
               nam(n0+ 5)='Re B40       (meV)  '
               nam(n0+ 6)='Re B41       (meV)  ' 
               nam(n0+ 7)='Im B41       (meV)  '
               nam(n0+ 8)='Re B42       (meV)  '      
               nam(n0+ 9)='Im B42       (meV)  '
               nam(n0+10)='Re B43       (meV)  ' 
               nam(n0+11)='Im B43       (meV)  '
               nam(n0+12)='Re B44       (meV)  '   
               nam(n0+13)='Im B44       (meV)  '
	    endif
	    if(j2.ge.6) then
               nam(n0+14)='Re B60       (meV)  '                  
               nam(n0+15)='Re B61       (meV)  '
               nam(n0+16)='Im B61       (meV)  '      
               nam(n0+17)='Re B62       (meV)  '
               nam(n0+18)='Im B62       (meV)  ' 
               nam(n0+19)='Re B63       (meV)  '
               nam(n0+20)='Im B63       (meV)  '   
               nam(n0+21)='Re B64       (meV)  '
               nam(n0+22)='Im B64       (meV)  ' 
               nam(n0+23)='Re B65       (meV)  '
               nam(n0+24)='Im B65       (meV)  '      
               nam(n0+25)='Re B66       (meV)  '
               nam(n0+26)='Im B66       (meV)  '
	    endif      

        else if(symmetry.eq.1) then
            if(j2.ge.2) then
               nam(n0+ 1)='Re B20       (meV)  '
               nam(n0+ 2)='Re B22       (meV)  ' 
               if(sbkq(2,2).eq.-1) nam(n0+ 2)='Im B22       (meV)  '  
	    endif
	    if(j2.ge.4) then 
               nam(n0+ 3)='Re B40       (meV)  '
               nam(n0+ 4)='Re B42       (meV)  '   
               nam(n0+ 5)='Im B42       (meV)  '
               nam(n0+ 6)='Re B44       (meV)  ' 
               nam(n0+ 7)='Im B44       (meV)  '
	    endif
	    if(j2.ge.6) then
               nam(n0+ 8)='Re B60       (meV)  '      
               nam(n0+ 9)='Re B62       (meV)  '         
               nam(n0+10)='Im B62       (meV)  ' 
               nam(n0+11)='Re B64       (meV)  '
               nam(n0+12)='Im B64       (meV)  '   
               nam(n0+13)='Re B66       (meV)  '
               nam(n0+14)='Im B66       (meV)  '                  
            endif
        else if(symmetry.eq.2) then
            if(j2.ge.2) then
               nam(n0+ 1)='Re B20       (meV)  '
               nam(n0+ 2)='Re B22       (meV)  ' 
               if(sbkq(2,2).eq.-1) nam(n0+ 2)='Im B22       (meV)  ' 
	    endif
	    if(j2.ge.4) then
               nam(n0+ 3)='Re B40       (meV)  '
               nam(n0+ 4)='Re B42       (meV)  '   
               if(sbkq(4,2).eq.-1) nam(n0+ 4)='Im B42       (meV)  ' 
               nam(n0+ 5)='Re B44       (meV)  '
            if(sbkq(4,4).eq.-1) nam(n0+ 5)='Im B44       (meV)  ' 
	    endif
	    if(j2.ge.6) then
               nam(n0+ 6)='Re B60       (meV)  '          
               nam(n0+ 7)='Re B62       (meV)  '
               if(sbkq(6,2).eq.-1) nam(n0+ 7)='Im B62       (meV)  ' 
               nam(n0+ 8)='Re B64       (meV)  '       
               if(sbkq(6,4).eq.-1) nam(n0+ 8)='Im B64       (meV)  ' 
               nam(n0+ 9)='Re B66       (meV)  '         
               if(sbkq(6,6).eq.-1) nam(n0+ 9)='Im B66       (meV)  ' 
            endif
        else if(symmetry.eq.3) then
            if(j2.ge.2) then
               nam(n0+ 1)='Re B20       (meV)  '
	    endif
	    if(j2.ge.4) then
               nam(n0+ 2)='Re B40       (meV)  ' 
               nam(n0+ 3)='Re B44       (meV)  '
               if(sbkq(4,4).eq.-1) nam(n0+ 3)='Im B44       (meV)  ' 
            endif
	    if(j2.ge.6) then
               nam(n0+ 4)='Re B60       (meV)  '   
               nam(n0+ 5)='Re B64       (meV)  '          
               nam(n0+ 6)='Im B64       (meV)  '          
            endif
        else if(symmetry.eq.4) then
            if(j2.ge.2) then 
               nam(n0+ 1)='Re B20       (meV)  '
	    endif
	    if(j2.ge.4) then
               nam(n0+ 2)='Re B40       (meV)  ' 
               nam(n0+ 3)='Re B44       (meV)  '
               if(sbkq(4,4).eq.-1) nam(n0+ 3)='Im B44       (meV)  '
	    endif 
	    if(j2.ge.6) then
               nam(n0+ 4)='Re B60       (meV)  '   
               nam(n0+ 5)='Re B64       (meV)  '          
               if(sbkq(6,4).eq.-1) nam(n0+ 5)='Im B64       (meV)  ' 
            endif
        else if(symmetry.eq.5) then
            if(j2.ge.2) then
               nam(n0+ 1)='Re B20       (meV)  '
	    endif
	    if(j2.ge.4) then
               nam(n0+ 2)='Re B40       (meV)  ' 
               nam(n0+ 3)='Re B43       (meV)  '
               if(sbkq(4,3).eq.-1) nam(n0+ 3)='Im B43       (meV)  ' 
	    endif
	    if(j2.ge.6) then
               nam(n0+ 4)='Re B60       (meV)  '   
               nam(n0+ 5)='Re B63       (meV)  '
               nam(n0+ 6)='Im B63       (meV)  '          
               nam(n0+ 7)='Re B66       (meV)  '
               nam(n0+ 8)='Im B66       (meV)  '       
            endif
        else if(symmetry.eq.6) then
            if(j2.ge.2) then 
               nam(n0+ 1)='Re B20       (meV)  '
	    endif
	    if(j2.ge.4) then
               nam(n0+ 2)='Re B40       (meV)  ' 
               nam(n0+ 3)='Re B43       (meV)  '
               if(sbkq(4,3).eq.-1) nam(n0+ 3)='Im B43       (meV)  '
	    endif 
            if(j2.ge.6) then
               nam(n0+ 4)='Re B60       (meV)  '   
               nam(n0+ 5)='Re B63       (meV)  '
               if(sbkq(6,3).eq.-1) nam(n0+ 5)='Im B63       (meV)  ' 
               nam(n0+ 6)='Re B66       (meV)  '          
               if(sbkq(6,6).eq.-1) nam(n0+ 6)='Im B66       (meV)  ' 
            endif
        else if(symmetry.eq.7) then
            if(j2.ge.2) then
               nam(n0+ 1)='Re B20       (meV)  '
	    endif
	    if(j2.ge.4) then
               nam(n0+ 2)='Re B40       (meV)  '
	    endif
	    if(j2.ge.6) then 
               nam(n0+ 3)='Re B60       (meV)  '
               nam(n0+ 4)='Re B66       (meV)  '  
               if(sbkq(6,6).eq.-1) nam(n0+ 4)='Im B66       (meV)  ' 
            endif
        else if(symmetry.eq.8 .and. nre.ne.6 .and. nre.ne.7 ) then
            if(j2.ge.4) then
               nam(n0+ 1)='Re B40    (meV)  R44'
               if(sbkq(4,4).eq.-1) nam(n0+ 1)='Re B40    (meV)  I44'            
	    endif
	    if(j2.ge.6) then
               nam(n0+ 2)='Re B60    (meV)  R64'
               if(sbkq(6,4).eq.-1) nam(n0+ 2)='Re B60    (meV)  I64' 
            endif
        endif                               
                                      
      return
      end
c-----------------------------------------------------------------------
c     make default values public
c-----------------------------------------------------------------------
      subroutine p_public(value,pos)
         implicit none
         real*8 value
         integer pos
      include 'CF_SOURCES:FRILLS.INC'
         p(pos) = value
      return
      end
c-----------------------------------------------------------------------
c     FOCUS defined fitting subroutine
c-----------------------------------------------------------------------
      subroutine CFFIT 

        include 'CF_SOURCES:FUNCTION.INC'
        include 'CF_SOURCES:FOCUS.INC' 
        include 'CF_SOURCES:XR.INC' 
        include 'CF_SOURCES:CONTROL.INC' 
        include 'CF_SOURCES:DEGENERATION.INC'
        include 'CF_SOURCES:SELF.INC' 
	include 'CF_SOURCES:RESOLUTION.INC'


 	external c_fmevkelvin
 	external c_excitations
 	external transform
 	external no_cf_par
 	external exp_
        external conv_shape

 	integer out
 	integer*4 i,j,k 
 	integer no_cf_par,ncftot,n0
 	real*8 background,slope,intcorr,width
 	real*8 el_width,el_intensity
        real*8 qe_width
 	real*8 c_fmevkelvin
 	real*8 shape,formfactor,formf,qnorm
 	real*8 epos,iint,pi
 	real*8 det_angle,phi,ei,ki2,kfdki
c       ---------------------------------------
        real*8  e_old(17*17),e_n,e_o
        integer l_old(0:17*17),n_old
        real*8  e_neg,e_pos,elasticpos
	integer ind(17*17),positiv_start  
c       ---------------------------------------
        real*8 xx(MAX_PAR),pv(MAX_PAR)
        real*8 xxx,www
c       ---------------------------------------
        real*8 r_chop, r_mod
        real*8 conv_shape
c       ---------------------------------------
        real*8 inelas
  	real*8 sigma_van_el, sigma_chop_el
  	real*8 sigma_mod,    sigma_chop
 	real*8 gamma_van_el
        real*8 gamma,sigma
c       ---------------------------------------
 100    format(1x,3x,a11,a20,a4,f16.5) 
  	p(1) = abs(p(1))
 	if(p(1).eq.0.0d0) then 
           p(1) = 4.2d0
           call p_public(p(1),1)
           write(*,100) 'FOCUS sets ','Temperature (Kelvin)',' TO ',p(1)
        end if
 	temp = p(1)  
 
        p(2) = abs(p(2))
 	if(p(2).eq.0.0d0) then 
           p(2)=67.0d0
           call p_public(p(2),2)
           write(*,100) 'FOCUS sets ','Incoming Energy(meV)',' TO ',p(2)
        end if
 	ei = p(2)
 
 	p(3)      = abs(p(3))
 	pi        = 4.0d0*atan(1.0d0)     
        phi       = p(3)
 	det_angle = phi*pi/180.0d0 
   
 
 	background = p(4)
 	slope      = p(5)   
 
 	p(6) = abs(p(6))
 	el_intensity = p(6)

        elasticpos   = p(7)

 	p(8) = abs(p(8))
 	if(p(8).eq.0.0d0 .and. setres.eq.0) then 
           p(8)=ei * 0.02
           call p_public(p(8),8)
           write(*,100) 'FOCUS sets ','Elastic   FWHM (meV)',' TO ',p(8)
        else if(p(8).ne.0.0d0 .and. setres.eq.1) then 
           p(8)=0.0d0
           call p_public(p(8),8)
           write(*,100) 'FOCUS sets ','Elastic   FWHM (meV)',' TO ',p(8)
        end if


 	el_width    = p(8)   
 

 	p(9) = abs(p(9))
 	if(p(9).eq.0.0d0 .and. setres.eq.0) then 
           p(9)=ei * 0.02
           call p_public(p(9),9)
           write(*,100) 'FOCUS sets ','Quasielas FWHM (meV)',' TO ',p(9)
        end if
 	qe_width    = p(9)   
 
 	p(10) = abs(p(10))
 	if(p(10).eq.0.0d0 .and. setres.eq.0) then 
           p(10)=p(9)
           call p_public(p(10),10)
           write(*,100) 'FOCUS sets ','Inelastic FWHM (meV)',' TO ',p(10)
         end if
 	width = p(10)
 
 	p(11) = abs(p(11))
 	if(p(11).eq.0.0d0) then 
           p(11) = 1.0d0
           call p_public(p(11),11)
           write(*,100) 'FOCUS sets ','Intensity factor    ',' TO ',p(11)
        end if
 	intcorr = p(11)
c       ----------------------------------------------------------        
!       if the fitting parameter are the spherical xR
!       transform them back to vkq-parameter 
c       ----------------------------------------------------------  
        if(setxr.eq.1) then
           do k=1,nindvold
              xx(k) = p(indvold(k))
           end do
           call t_x_v(xx,pv,nindvold)
           do k=1,nindvold
              p(indvold(k)) = pv(k)
           end do
           type = 'Vkq'
        end if
c       ----------------------------------------------------------        
!       store the magnetic fields
c       ----------------------------------------------------------        
        n0 = no_0_par()
 	bext(1) = p(n0-5)
 	bext(2) = p(n0-4)
 	bext(3) = p(n0-3)                      
 	bmol(1) = p(n0-2)
 	bmol(2) = p(n0-1)
 	bmol(3) = p(n0)
c       ----------------------------------------------------------        
!       delete all Bkq parameter 
c       ----------------------------------------------------------        
        do k=2,6,2
 	   do j=0,6
     	       bkq(k,j) = 0.0d0
           end do
        end do                                           
c       ----------------------------------------------------------        
!       find out how many bkq parameters there are
c       ----------------------------------------------------------        
        ncftot = no_cf_par(symmetry,nre)
c       ----------------------------------------------------------        
c       store crystal field parameter in an one dimensional array
c       ----------------------------------------------------------        
  	do i=1,ncftot
           pv(i) = p(n0+i)
 	end do
c       ----------------------------------------------------------        
c       transform the wx-parameter to the bkq-parameter
c       ----------------------------------------------------------        
        if(setwx.eq.1) then
           www = pv(1)
           xxx = pv(2)
           call t_xw(nre,xxx,www,bkq)
           pv(1) = bkq(4,0)
           pv(2) = bkq(6,0)
           type = 'Bkq'
        end if

        call transform(pv,bkq,symmetry,nre)
c       ----------------------------------------------------------        
c       transform the akq-parameter to the bkq-parameter
c       ----------------------------------------------------------        
        if(setakq.eq.1) then
           call c_akq_bkq(nre,bkq,bkq)
           type = 'Bkq'
        endif
c       ----------------------------------------------------------        
!       transform the Bkq from the meV scale to the Kelvin scale
c       ----------------------------------------------------------        
         do k=2,6,2
 	   do j=0,6
     	       bkq(k,j) = bkq(k,j)*c_fmevkelvin()  
           end do
        end do                                       
c       ----------------------------------------------------------
!       save the old excitation situation (for the lineshape only)
c       ----------------------------------------------------------        
        if(iscalc.ne.0 .and. setfit.ne.0) then
          n_old    = n_excitation
          l_old(0) = lineshape(0)
          do k=1,n_old
             e_old(k) = e_excitation(k)
             l_old(k) = lineshape(k)
          end do
        end if
c       ----------------------------------------------------------        
!       calculate the excitations
c       ----------------------------------------------------------        
  	out = 42
 	call c_excitations(out,nre,symmetry,bkq,sbkq,bmol,bext,temp,
     *       e_excitation,i_excitation,n_excitation,type,de,di)
c       ----------------------------------------------------------        
!       if the fitting parameter are originally the spherical xR
!       parameter then  transform v-parameter back to xr-parameter 
c       ----------------------------------------------------------  
        if(setxr.eq.1) then
           do k=1,nindvold
              pv(k) = p(indvold(k))
           end do
           call t_v_x(pv,xx,nindvold)
           do k=1,nindvold
              p(indvold(k)) = xx(k)
           end do
           type = 'XR'
        end if
c       ----------------------------------------------------------        
!       set the lineshapes to the old situation
c       ----------------------------------------------------------        
!       set first all excitations to gaussian line
        if(iscalc.ne.0 .and. setfit.ne.0) then
          do k=0,17*17
             lineshape(k) = 1
          end do
!         set the old excitations to their already defined shape
          lineshape(0) = l_old(0)
          do k=1,n_excitation
             e_n = e_excitation(k)
             do i=1,n_old
                e_o = e_old(i)
                if(abs(e_n-e_o).le. de) then
                   lineshape(k) = l_old(i)
                   goto 98
                end if
             end do
  98      end do

!         set the line shape of the negative excitations to the value of the positive ones
          call energy_sort(e_excitation,ind,n_excitation,positiv_start)
          if(positiv_start.gt.1) then
             do i=1,positiv_start-1
                e_neg = e_excitation(ind(i))
                do j=positiv_start,n_excitation
                   e_pos = e_excitation(ind(j))
                   if(abs(e_pos+e_neg) .le. de ) then
                      lineshape(ind(i)) = lineshape(ind(j))
                      goto 99
                   end if
                end do
  99         end do    
          end if
        end if
        if(setfit.ne.0) then
c       ----------------------------------------------------------        
        if(setres.ne.0) then
           sigma_van_el  = fwhm_van_el /sqrt(log(256.0))
           sigma_chop_el = fwhm_chop_el/sqrt(log(256.0))
           sigma_mod     = sqrt(sigma_van_el**2 - sigma_chop_el**2)
           sigma_chop    = sigma_chop_el
           gamma_van_el  = decay/2.0d0   ! decay is in FWHM
        end if
c       ----------------------------------------------------------        
 	nc = n_excitation
 	ki2= 0.4825924583d0*ei    
        do j=1,nv  
!          -------------------------------------------------------
!          BACKGROUND
!          -------------------------------------------------------
 	   yc(j,1)  = 1000*(background + slope*x(j))/4/pi
!          -------------------------------------------------------
!          ELASTIC LINE
!          -------------------------------------------------------
           epos = elasticpos
 	   iint = 1000*el_intensity*intcorr/4/pi  ! intensity in mb per sr
           yc(j,2) = 0.0d0
           if(setres.eq.0 .and. x(j).le.ei) then
 	      yc(j,2) = shape(lineshape(0),x(j),epos,el_width,iint,temp)
           else if(setres.eq.1 .and. x(j).le.ei) then
                 r_chop = (1.0d0 + (lmc+lcs)/lsd * (1-x(j)/ei)**1.5)/(1.0d0 + (lmc+lcs)/lsd)
                 r_mod  = (1.0d0 +      lcs /lsd * (1-x(j)/ei)**1.5)/(1.0d0 +      lcs /lsd)
!                --------------------------------------------------------------------------       
!                sets gamma and sigma as a function of the lineshape
!                --------------------------------------------------------------------------       
                 call setgs(gamma,sigma,lineshape(0),el_width)       

                 gamma_exp = gamma_van_el * r_mod 
                 sigma_exp = sqrt((sigma_mod *r_mod)**2 + (sigma_chop*r_chop)**2 + sigma**2)
                 yc(j,2)   = conv_shape(x(j)-epos,gamma,gamma_exp,sigma_exp,lineshape(0),iint,temp)
           end if
           if(setself.eq.1) yc(j,2) = yc(j,2)*self(ei,x(j),det_angle,muts,muta,0.0)
!          -------------------------------------------------------
!          EXCITATIONS
!          -------------------------------------------------------
 	   yc(j,4) = 0.0d0
 	   do k=1,nc
              epos = e_excitation(k)/c_fmevkelvin()
              if(x(j).le.ei ) then
                 kfdki= sqrt(1-x(j)/ei)       ! kfdki := k_f/k_i
                 qnorm= sqrt(ki2*(1+kfdki**2-2*kfdki*cos(det_angle)))  
  	         formf= formfactor(qnorm,nre)
 	         iint = 1000*i_excitation(k)*intcorr*formf**2 /4/pi ! intensity in mb per sr
  	         if(iscalc.eq.0.and.epos.eq.0.0d0) lineshape(k) = 3 ! Lorentz*Bose  
                 if(epos.eq.0.0d0) then
!                    -------------------------------------------------
!                    QUASI ELASTIC LINE
!                    -------------------------------------------------
                     yc(j,3) = 0.0d0
                     if(setres.eq.0) then
                        yc(j,3)=shape(lineshape(k),x(j),epos,qe_width,iint,temp) 
                     else
    	             	r_chop = (1.0d0 + (lmc+lcs)/lsd * (1-x(j)/ei)**1.5)/(1.0d0 + (lmc+lcs)/lsd)
        	        r_mod  = (1.0d0 +      lcs /lsd * (1-x(j)/ei)**1.5)/(1.0d0 +      lcs /lsd)
!                       --------------------------------------------------------------------------       
!                       sets gamma and sigma as a function of the lineshape
!                       --------------------------------------------------------------------------       
                        call setgs(gamma,sigma,lineshape(k),qe_width)

                	gamma_exp = gamma_van_el * r_mod 
                	sigma_exp = sqrt((sigma_mod *r_mod)**2 + (sigma_chop*r_chop)**2 + sigma**2)
                        yc(j,3)   = conv_shape(x(j)-epos,gamma,gamma_exp,sigma_exp,lineshape(k),iint,temp)
                     end if
                     if(setself.eq.1) yc(j,3) = yc(j,3)*self(ei,x(j),det_angle,muts,muta,0.0)
                 else
!                    -------------------------------------------------
!                    INELASTIC LINES
!                    -------------------------------------------------
                     if(setres.eq.0) then
   	                inelas = shape(lineshape(k),x(j),epos,width,iint,temp) 
                     else
    	             	r_chop = (1.0d0 + (lmc+lcs)/lsd * (1-x(j)/ei)**1.5)/(1.0d0 + (lmc+lcs)/lsd)
        	        r_mod  = (1.0d0 +      lcs /lsd * (1-x(j)/ei)**1.5)/(1.0d0 +      lcs /lsd)
!                       --------------------------------------------------------------------------       
!                       sets gamma and sigma as a function of the lineshape
!                       --------------------------------------------------------------------------       
                        call setgs(gamma,sigma,lineshape(k),width) 

                	gamma_exp = gamma_van_el * r_mod 
                	sigma_exp = sqrt((sigma_mod *r_mod)**2 + (sigma_chop*r_chop)**2 + sigma**2)
                        inelas    = conv_shape(x(j)-epos,gamma,gamma_exp,sigma_exp,lineshape(k),iint,temp)
                     end if
                     if(setself.eq.1) inelas = inelas*self(ei,x(j),det_angle,muts,muta,0.0)
                     yc(j,4) = yc(j,4) + inelas
                 end if
              end if
 	   end do	
!          -------------------------------------------------------
!           FIT    =    BG   +  EL  +   QE  +   INEL
!          -------------------------------------------------------
           ycal(j) = yc(j,1)+yc(j,2)+yc(j,3)+yc(j,4) 
        end do                   
 	nc = 4 !  1=BG  2=el  3=QE  4=in 
        end if !  if(setfit.ne.0) then

  	iscalc = 1
 
      return
      end
c     -----------------------------------------------------------------------
c     set the lineshape of a special excitation
c     -----------------------------------------------------------------------
      real*8 function shape(no_lineshape,x,x0,fwhm,intensity,temp) 
 	implicit none
 	integer no_lineshape
 	real*8 x,x0,fwhm,intensity,temp
 	real*8 gauss,lorentz,lorentz_bose
        integer no_min
        integer no_max

        data no_min /1/
        data no_max /3/

  	if(fwhm.eq.0.0d0) then
           shape = 0.0d0
 	   return
 	end if   


  	if(no_lineshape.lt.no_min .or. no_lineshape.gt.no_max) no_lineshape=1
 
 	if(no_lineshape.eq. 1) then
 		shape = gauss(x,x0,fwhm,intensity)
 	else if(no_lineshape.eq. 2) then
                shape = lorentz(x,x0,fwhm,intensity)
        else if(no_lineshape.eq. 3) then
                shape = lorentz_bose(x,fwhm,intensity,temp) 
 	end if
      return
      end
!     --------------------------------------------------------------------------       
!     sets gamma and sigma as a function of the lineshape
!     --------------------------------------------------------------------------       
      subroutine setgs(gamma,sigma,no_lineshape,fwhm) 
           real*8  gamma,sigma,fwhm
	   integer no_lineshape
   	   if(no_lineshape.eq.1) then      !Gauss
              gamma = 0.0d0
   	      sigma = fwhm / sqrt(log(256.0d0))
           else if(no_lineshape.eq.2) then !Lorentz             
              gamma = fwhm/2.0d0
   	      sigma = 0.0d0 
 	   else if(no_lineshape.eq.3) then !Lorentz*Bosefactor
              gamma = fwhm/2.0d0
   	      sigma = 0.0d0 
           end if
      return
      end
!     --------------------------------------------------------------------------       
!     define a string which says which lineshape it is
!     -------------------------------------------------------------------------- 
      character*18 function text_line(no_lineshape)
 	implicit none               
 	integer no_lineshape               
 	text_line='Gauss             '  
 	if(no_lineshape.eq.1) text_line='Gauss             '  
 	if(no_lineshape.eq.2) text_line='Lorentz           '  
 	if(no_lineshape.eq.3) text_line='Lorentz*Bosefactor'
      return
      end
c     -----------------------------------------------------------------------
      real*8 function exp_control(x)
 	implicit none
        real*8 x,exp_max
 
 	exp_max = 20.0d0
 
 	if(abs(x).gt.exp_max) x = sign(exp_max,x)
 	exp_control = exp(x)
      return
      end
c     -----------------------------------------------------------------------
      real*8 function gauss(x,x0,width_FWHM,intensity)
 	implicit none
 	real*8 x,x0,width_FWHM,intensity    
 	real*8 pi,sigma,exp_control
 
 	pi    = 4.0d0*atan(1.0d0)
 	sigma = width_FWHM/sqrt(log(256.0d0))	
        gauss = intensity*exp_control(-0.5d0*((x-x0)/sigma)**2)/(sqrt(2*pi)*sigma)
 
      return
      end
c     -------------------------------------------------------------------------- 
      real*8 function lorentz_bose(x,fwhm,intensity,temp) 
	implicit none
        real*8 x,x0,fwhm,intensity,temp
        real*8 pi,bose
 	
 	external c_fmevkelvin
 
        pi = 4.0d0*atan(1.0d0)
 	lorentz_bose = bose(x,temp)*intensity/pi*fwhm/2/((fwhm/2)**2+x**2)
      return
      end                                                                                           
c     --------------------------------------------------------------------------
c     define the bose factor
c     --------------------------------------------------------------------------    
      real*8 function bose(x,temp)
	implicit none
	real*8 x,temp
	real*8 c_fmevkelvin,exp_control,expp,xx
 
 	xx   =  x*c_fmevkelvin()/temp
 	expp =  exp_control(-xx)
 
 	if( expp.eq.1.0d0 ) then
 	    bose = 1.0d0
 	else
 	    bose = xx/(1-expp)
 	end if
 
      return
      end	
c     --------------------------------------------------------------------------
      real*8 function lorentz(x,x0,fwhm,intensity)
	implicit none
        real*8 x,x0,fwhm,intensity,temp
        real*8 pi,exp_control,c_fmevkelvin
 	
 	external c_fmevkelvin
 
        pi = 4.0d0*atan(1.0d0)
 	lorentz = intensity/pi*fwhm/2/((fwhm/2)**2+(x-x0)**2)
      return
      end                                                                                           
