cc-------------------------------------------------------------------- 
cc 
cc the CF_FABI.FOR module calculates the crystal field splitting of 
cc a rare earth ion R3+ for all 32 point symmetries.
cc
cc--------------------------------------------------------------------
cc
cc written by :  Dr. Peter O. Fabi
cc               Rutherford Appleton Laboratory
cc               ISIS Facility
cc               Chilton
cc               Didcot
cc               Oxon, OX11, 0QX
cc               tel.: 0235-44-5428
cc--------------------------------------------------------------------
cc
cc Input  :   - number of R3+
cc
cc            - number of point symmetry
cc         
cc            - complex B     parameter k=2,4,6 and |q|<= k  [in kelvin]
cc                       kq
cc            - molecular magnetic field [in tesla]
cc
cc            - external  magnetic field [in tesla]
cc
cc            - temperature of the sample [in kelvin]
cc
cc
cc Output:    - transition energies in kelvin
cc
cc            - wavefunctions
cc
cc            - transition matrix elements for a single crystal
cc
cc            - transition matrix elements for a powder
cc
cc            - transition intensities for a powder [in barn]
cc
cc            - the magnetic moment of R3+ [in myb]
cc
cc--------------------------------------------------------------------
      program cf
	implicit none
cc      --------------------------------
 	integer out
 	parameter(out=6)
cc      --------------------------------
        include 'CF_SOURCES:CF_FABI.INC'
cc      --------------------------------    
        common /DEGENERATION/ de,di
        real*8 de,di
c       --------------------------
        common /CONTROL/ setdegoff,stdhamoff
        integer setdegoff,stdhamoff
c       ---------------------------------------------------------------------
cc      for internal use
 	complex*16 akq(0:6,0:6)
        integer i,k
 	real*8 mx,my,mz  
 	real*8 c_x_moment,c_y_moment,c_z_moment   
 	real*8 c_fmevkelvin     
cc	real*8  e_excitation(17*17) ! energy position [kelvin] of the excit.
cc	real*8  i_excitation(17*17) ! intensity [barn] of the excitations
cc	integer n_excitation        ! number of excitations found          
cc      -----------------
cc      for testing of the subroutines c_xsort_v and c_v_xsort
 	integer n,ind(26)
 	real*8  v(26),x_sort(26)
cc      --------------
 	real*8 x,w   ! lea,leask and wolff parameter
cc      --------------
cc      temperature of the sample in Kelvin
 	temp = 0.5d0 
cc      -------------------
 	nre = 3
cc      -------------------
 	symmetry = 4 ! 
cc      -------------------
        de = 0.01d0  !Kelvin
        di = 0.01d0  !barn
cc      -------------------
c       [akq] = Kelvin/a0**k
c  	akq(2,0) =  295
c 	akq(2,2) = -454
c 	akq(4,0) = -12.3
c 	akq(4,2) =  0.
c 	akq(4,4) =  0.
c 	akq(6,0) = -1.84
c 	akq(6,2) =  9.8
c 	akq(6,4) = -15.9
c 	akq(6,6) =  0
c
c calculates the bkq from the given akq
c
c       call c_akq_bkq(nre,akq,bkq) 
c
 	sbkq(2,0) = 1 
	sbkq(2,1) = 1
  	sbkq(2,2) = 1

 	sbkq(4,0) = 1 
	sbkq(4,1) = 1
  	sbkq(4,2) = 1
	sbkq(4,3) = 1
 	sbkq(4,4) = 1 

  	sbkq(6,0) = 1 
	SBKQ(6,1) = 1
 	sbkq(6,2) = 1
	SBKQ(6,3) = 1 
 	sbkq(6,4) = 1
	SBKQ(6,5) = 1   
 	sbkq(6,6) = 1  
cc      --------------
cc      [bmol] = tesla
 	bmol(1) = 0.0d0 
 	bmol(2) = 0.0d0   
 	bmol(3) = 0.0d0  
cc      -----------------
cc      [bext] = tesla
 	bext(1) = 0.0d0
 	bext(2) = 0.0d0
 	bext(3) = 0.0d0
cc      -----------------
cc      parameter for NdCu2Si2
  	bkq(2,0) = -3.1d-2 * c_fmevkelvin()
 	bkq(4,0) =  1.1d-3 * c_fmevkelvin()
 	bkq(4,4) = 1.33d-3 * c_fmevkelvin()
 	bkq(6,0) =-3.17d-5 * c_fmevkelvin()
 	bkq(6,4) = 6.62d-4 * c_fmevkelvin()          
cc      -----------------
cc      parameter for ErCu2Si2
c 	bkq(2,0) = -1.4d-2 * c_fmevkelvin()
c 	bkq(4,0) = -1.2d-4 * c_fmevkelvin()
c 	bkq(4,4) = 3.10d-4 * c_fmevkelvin()
c 	bkq(6,0) = 5.90d-7 * c_fmevkelvin()
c 	bkq(6,4) = 2.90d-5 * c_fmevkelvin()          
cc      -----------------
cc      parameter for TmCo2
cc      x = 0.70
cc	w = 0.77 ! Kelvin
cc
cc	bext(1) = 80/sqrt(3.0d0)
cc	bext(2) = 80/sqrt(3.0d0)
cc	bext(3) = 80/sqrt(3.0d0) ! Tesla
cc
cc	bext(3) = 130
cc
cc	sbkq(4,4) = 1
cc	sbkq(6,4) = 1
cc
cctransform the x,w to bkq
cc
cc      call t_xW(nre,x,w,bkq)
cc
cc
c prints an overview over the symmetry numbers
cc
  	call p_symmetry_numbers(out)
c
c prints an overview over the hamiltonian
c
  	call p_hamiltonian(out,symmetry)
c
c prints the input out 
c
  	call p_input(out,symmetry,nre,bkq,sbkq,bmol,bext,temp)
c
c calculate the crystal field splitting
c
   	call c_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,'Bkq')
c
c prints out the transition energies in meV
c
    	call p_energy_values(out,energy,dimj,fmevkelvin)
c
c prints the wave functions
c
    	call p_wave_functions(out,energy,wavefunction,dimj)
c
c tests if the wavefunctions are really orthonormal
c
    	call t_orthonormalisation(out,wavefunction,dimj)
c
c prints out the transition matrix elements
c
    	call p_matrix_elements(out,dimj,energy,jx2mat,jy2mat,jz2mat,jt2mat)
c
c prints 'occupation  factor'
c
    	call p_occupation_factor(out,occupation_factor,temp)
c
c prints out the transition intensities for the given temperature temp
c
    	call p_intensities(out,dimj,energy,intensity,temp)
c
c calculates the R3+ moment in myb
c
    	mx = c_x_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
    	my = c_y_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
    	mz = c_z_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
c
c prints out the R3+ moment 
c
   	call p_moments(out,mx,my,mz,gj,dimj)
c
c calculates the excitations for a given temperature
c
        call c_excitations(out,nre,symmetry,bkq,sbkq,bmol,bext,temp,
     *                    e_excitation,i_excitation,n_excitation,'Bkq',de,di)   
c
c
c a test if the transformations are working
c     
c	n=3
c	v(1) = -1
c	v(2) = 1.0d-10
c	v(3) = 1
c
c	do i=1,n
c	x_sort(i)=0
c	end do
c	call c_xsort_v(v,x_sort,n,ind)
c	write(6,*) (v(i),i=1,n)
c	write(6,*) (x_sort(i),i=1,n)
c	write(6,*) (ind(i),i=1,n)
c	do i=1,n
c	v(i)=0
c	end do
c	call c_v_xsort(v,x_sort,n,ind)
c	write(6,*) (v(i),i=1,n)
c                                              
c
c	call p_operator_norm(out,1)
c
c-----------------------------------
c end of program cf
      end
c-----------------------------------
