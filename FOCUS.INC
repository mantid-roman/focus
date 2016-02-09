c-------------------------------------------------------------------- 
c include file for the FOCUS.FOR module
c
c--------------------------------------------------------------------
c
c written by :  Dr. Peter O. Fabi
c               Rutherford Appleton Laboratory
c               ISIS Facility
c               Chilton
c               Didcot
c               Oxon, OX11, 0QX
c               tel.: 0235-44-5428
c--------------------------------------------------------------------
c       INPUT VALUES FOR C_CRYSTAL_FIELD
c       INPUT VALUES FOR C_EXCITATIONS

        include 'FOCUS.PAR'

	integer nre             ! 1=Ce3+,2=Pr3+, ... ,13=Yb3+  
	integer symmetry        ! chooses the point symmetry of the R3+ ion                           
	complex*16 bkq          ! crystal field parameters in kelvin
        real*8     bmol         ! molecular magnetic field in tesla
        real*8     bext         ! external  magnetic field in tesla
        real*8     sbkq         ! choose the real part (=1) or the
c                               ! imaginary part(=-1) of the crystal
c                               ! field parameter bkq, if possible.  
	real*8      temp        ! temperature of the measured sample 
        real*8  temp_calc       ! temperature used in calculation
	character*3 type        ! specify if Bkq is on input (type = 'bkq' or 'BKQ' or 'Bkq')
	                        !         or Vkq is on input (type = 'vkq' or 'VKQ' or 'Vkq')    
	common /FOCUS/ bkq(0:6,0:6),sbkq(0:6,0:6),bmol(3),bext(3),
     *                 temp,type,nre,symmetry,cf_out,temp_calc
c
	integer cf_out          ! not for C_CRYSTAL_FIELD
c       ------------------------------------------------------------- 
c       OUTPUT VALUES FOR C_CRYSTAL_FIELD
	real*8 gj,myb,dimj,alphaj,betaj,gammaj,r2,r4,r6
	real*8 fmevkelvin,energy,eshift
	real*8 jx2mat,jy2mat,jz2mat,jt2mat
	real*8 intensity,occupation_factor
	complex*16 wavefunction
	common /CFOUTPUT/ gj,myb,dimj,alphaj,betaj,gammaj,r2,r4,r6,
     *                    fmevkelvin,eshift,energy(17),wavefunction(17,17),
     *                    jx2mat(17,17),jy2mat(17,17),jz2mat(17,17),
     *                    jt2mat(17,17),intensity(17,17),occupation_factor
c-----------------------------------------------------------------------
c       OUTPUT VALUES FOR C_EXCITATIONS
c       the excitations are ordered by decreasing intensity
	real*8  e_excitation        ! energy position [kelvin] of the excitation
	real*8  i_excitation        ! intensity [barn] of the excitation
        integer n_excitation        ! number of excitations found with
                                    ! an intensity greater equal 0.01 barn
	integer lineshape           ! lineshapes of the excitations + elastic line
	integer n_lineshape         ! no of lineshapes
	integer iscalc              ! =1 or 0 : if excitations are calculated or not
	parameter(n_lineshape=17*17)
	
        common /EXCITATIONS/ e_excitation(17*17),i_excitation(17*17),
     *                       n_excitation,lineshape(0:n_lineshape),
     *                       iscalc
c--------------------------------------------------------------------        
        character*19 focus_version
        data focus_version/'FOCUS_REF:FOCUS.VER'/
c--------------------------------------------------------------------        
c end of file
c--------------------------------------------------------------------
