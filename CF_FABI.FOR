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
c
c-----------------------------------------------
c returns the number of parameters which are
c the same for all symmetries
c-----------------------------------------------
      integer function no_0_par()
	implicit none
        include 'CF_SOURCES:no0par.inc'
        no_0_par = cno_0_par
      return
      end
c-----------------------------------------------
c returns the number of crystal field parameters
c-----------------------------------------------
      integer function no_cf_par(symmetry,nre)
	implicit none
	integer symmetry,nre
	integer p2,p4,p6
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6      

        if(symmetry.eq.0) then
	   p2 =  4
	   p4 =  9
	   p6 = 13
        else if(symmetry.eq.1) then
	   p2 =  2
	   p4 =  5
	   p6 =  7
        else if(symmetry.eq.2) then 
           p2 =  2
	   p4 =  3
	   p6 =  4
        else if(symmetry.eq.3) then 
	   p2 =  1
	   p4 =  2 
	   p6 =  3
        else if(symmetry.eq.4) then
	   p2 =  1
	   p4 =  2
	   p6 =  2
        else if(symmetry.eq.5) then 
	   p2 =  1
           p4 =  2
	   p6 =  5
        else if(symmetry.eq.6) then 
	   p2 =  1
	   p4 =  2
	   p6 =  3 
        else if(symmetry.eq.7) then
	   p2 =  1
	   p4 =  1
	   p6 =  2
        else if(symmetry.eq.8) then 
	   p2 =  0
	   p4 =  1
	   p6 =  1
	endif
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
	no_cf_par = 0
	if(dimj-1.ge.6.0d0) no_cf_par = no_cf_par + p6
	if(dimj-1.ge.4.0d0) no_cf_par = no_cf_par + p4
	if(dimj-1.ge.2.0d0) no_cf_par = no_cf_par + p2

	if(nre.eq.7) no_cf_par = 0
      return
      end                                        
c--------------------------------------------
c prints the section capitals in the file out
c--------------------------------------------
      subroutine section(out,s)
	implicit none
	integer out
	character s*(*)
        integer i

        write(out,*)
	write(out,*) ('-',i=1, 74   )  
        write(out,*) s
	write(out,*) ('-',i=1, 74   )  

      return
      end
c-----------------------------------------------------------
c transforms the x to v
c                -    -
c-----------------------------------------------------------
c v: output
c x: input
c n: input
      subroutine t_x_v(x,v,n)
	implicit none
	integer n,i,j
	real*8 x(*),v(*)

c       ----------------------------------------------------
c       set all v(i) to the overall splitting x(1) 
c       ----------------------------------------------------
	do i=1,n
	   v(i) = x(1)
	end do
	if(n.eq.1) return
c       ----------------------------------------------------
c       calculate the v(i) for n>1 
c       ----------------------------------------------------
	do j=n,2,-1
           do i=1,j-1
              v(i) = v(i) * sqrt( 1.0d0 - x(j)**2 )
           end do
           v(j) = v(j) * x(j)
        end do
      return
      end
c----------------------------------------------------
c transforms the v to x
c                -    -
c----------------------------------------------------
c v: input
c x: output
c n: input
      subroutine t_v_x(v,x,n)
	implicit none
	integer n,i,j
	real*8 x(*),v(*)
        real*8 i_norm,d
c       ---------------------------------------------
c       set x(i) to default values
c       ---------------------------------------------
	do i=1,n
	   x(i) = 0.0d0           
	end do
c       ---------------------------------------------
c       calculate the overall CEF splitting x(1)
c       ---------------------------------------------
	x(1) = sign(1.0d0,v(1)) * i_norm(v,n)
	call ifnull(x(1))
	if(n.eq.1.or.x(1).eq.0.0d0) return
c       ---------------------------------------------
c       calculate the the x(i) for i=2...n
c       ---------------------------------------------
	do i=n,2,-1
	   d = i_norm(v,i)
           call ifnull(d)
	   if(d.eq.0.0d0) then
              x(i+1) = sign(1.0d0,v(i+1)) * sign(1.0d0,v(1))           
              return
           end if
           x(i) = v(i)/d * sign(1.0d0,v(1))
	end do   
      return
      end
c----------------------------------
c calculates || v ||_i
c----------------------------------
      real*8 function i_norm(v,i)
	implicit none
	integer i,k
	real*8 v(*)

	i_norm = 0.0d0
	if(i.le.0) return
	do k=1,i
           i_norm = i_norm + v(k)*v(k)
	end do
	i_norm = sqrt(i_norm)
      return
      end
!-----------------------------------
! calculates the vkq from the v(i)
!-----------------------------------
      subroutine transform(v,vkq,symmetry,nre)
	implicit none
	integer symmetry,nre
	real*8 v(*)
	complex*16 vkq(0:6,0:6)
	integer k,q
	integer j2
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6      

	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
	j2 = dimj-1                                                         
                                                                              
! delete all Vkq
	do k=2,6,2
	   do q=0,k
	      vkq(k,q) = cmplx(0.0d0,0.0d0)
	   end do
	end do

	if(symmetry.eq.0) then
           if(j2.ge.2) then
		vkq(2,0) = v(1)
		vkq(2,1) = v(2)   
		vkq(2,2) = cmplx(v(3),v(4))
           endif
	   if(j2.ge.4) then
      		vkq(4,0) = v(5)
		vkq(4,1) = cmplx(v( 6),v( 7))     
		vkq(4,2) = cmplx(v( 8),v( 9))     
		vkq(4,3) = cmplx(v(10),v(11))     
		vkq(4,4) = cmplx(v(12),v(13))     
           endif
	   if(j2.ge.6) then                                   
		vkq(6,0) = v(14)
		vkq(6,1) = cmplx(v(15),v(16))     
		vkq(6,2) = cmplx(v(17),v(18))     
		vkq(6,3) = cmplx(v(19),v(20))     
		vkq(6,4) = cmplx(v(21),v(22))    
		vkq(6,5) = cmplx(v(23),v(24))     
		vkq(6,6) = cmplx(v(25),v(26))    
           endif                                    
	else if(symmetry.eq.1) then 
           if(j2.ge.2) then
		vkq(2,0) = v(1)
		vkq(2,2) = v(2)   
           endif
	   if(j2.ge.4) then
      		vkq(4,0) = v(3)
		vkq(4,2) = cmplx(v( 4),v( 5))     
		vkq(4,4) = cmplx(v( 6),v( 7))  
           endif
	   if(j2.ge.6) then
      		vkq(6,0) = v(8)
		vkq(6,2) = cmplx(v( 9),v(10))     
		vkq(6,4) = cmplx(v(11),v(12))  
		vkq(6,6) = cmplx(v(13),v(14))          
           endif                                                                                              
	else if(symmetry.eq.2) then 
           if(j2.ge.2) then
		vkq(2,0) = v(1)
		vkq(2,2) = v(2)   
           endif
	   if(j2.ge.4) then
      		vkq(4,0) = v(3)
		vkq(4,2) = v(4)     
		vkq(4,4) = v(5)  
           endif
	   if(j2.ge.6) then
      		vkq(6,0) = v(6)
		vkq(6,2) = v(7)     
		vkq(6,4) = v(8) 
		vkq(6,6) = v(9)          
           endif
	else if(symmetry.eq.3) then 
           if(j2.ge.2) then
		vkq(2,0) = v(1)   
           endif
	   if(j2.ge.4) then
      		vkq(4,0) = v(2)     
		vkq(4,4) = v(3)  
           endif
	   if(j2.ge.6) then
      		vkq(6,0) = v(4)   
		vkq(6,4) = cmplx(v(5),v(6))   
           endif
	else if(symmetry.eq.4) then          
           if(j2.ge.2) then
		vkq(2,0) = v(1)   
           endif
	   if(j2.ge.4) then
      		vkq(4,0) = v(2)     
		vkq(4,4) = v(3)  
           endif
	   if(j2.ge.6) then
      		vkq(6,0) = v(4)   
		vkq(6,4) = v(5)   
           endif
	else if(symmetry.eq.5) then 
           if(j2.ge.2) then
		vkq(2,0) = v(1)
           endif
	   if(j2.ge.4) then
		vkq(4,0) = v(2)
	        vkq(4,3) = v(3)
           endif
	   if(j2.ge.6) then
		vkq(6,0) = v(4)
		vkq(6,3) = cmplx(v(5),v(6))
		vkq(6,6) = cmplx(v(7),v(8)) 
           endif
	else if(symmetry.eq.6) then 
           if(j2.ge.2) then
		vkq(2,0) = v(1)
           endif
	   if(j2.ge.4) then
		vkq(4,0) = v(2)
	        vkq(4,3) = v(3)
           endif
	   if(j2.ge.6) then
		vkq(6,0) = v(4)
		vkq(6,3) = v(5)
		vkq(6,6) = v(6)             
           endif
	else if(symmetry.eq.7) then 
           if(j2.ge.2) then
		vkq(2,0) = v(1)
           endif
	   if(j2.ge.4) then
		vkq(4,0) = v(2)
           endif
	   if(j2.ge.6) then
		vkq(6,0) = v(3)
		vkq(6,6) = v(4)             
           endif
	else if(symmetry.eq.8) then 
	   if(j2.ge.4) then         
		vkq(4,0) =     v(1)
		vkq(4,4) =   5*v(1)
	   endif
	   if(j2.ge.6) then
	        vkq(6,0) =     v(2)
		vkq(6,4) = -21*v(2)
           endif
        endif
      return 
      end
!-----------------------------------
! calculates the vkq from the v(i)
! for the "limit"-vectors
!-----------------------------------
      subroutine lim_trans(symmetry,nre)
	implicit none
	integer symmetry,nre

	integer k,q
	integer j2
        complex*16 null
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6      

	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
	j2   = dimj-1                                                         
        null = cmplx(0.0d0,0.0d0)                                                                     

! delete all "Vkq"
	do k=2,6,2
	   do q=0,k
	      call set0(k,q,null)
	   end do
	end do

	if(symmetry.eq.0) then
           if(j2.ge.2) then
	        call rkq(2,0,  1)
	        call rkq(2,1,  2)   
	        call ckq(2,2,  3,4)
           endif
	   if(j2.ge.4) then
      		call rkq(4,0,  5)
		call ckq(4,1,  6, 7)     
		call ckq(4,2,  8, 9)     
		call ckq(4,3, 10,11)     
		call ckq(4,4, 12,13)     
           endif
	   if(j2.ge.6) then                                   
		call rkq(6,0, 14)
		call ckq(6,1, 15,16)     
		call ckq(6,2, 17,18)     
		call ckq(6,3, 19,20)     
		call ckq(6,4, 21,22)    
		call ckq(6,5, 23,24)     
		call ckq(6,6, 25,26)    
           endif                                    
	else if(symmetry.eq.1) then 
           if(j2.ge.2) then
		call rkq(2,0,  1)
		call rkq(2,2,  2)   
           endif
	   if(j2.ge.4) then
      		call rkq(4,0,  3)
		call ckq(4,2,  4,5)     
		call ckq(4,4,  6,7)  
           endif
	   if(j2.ge.6) then
      		call rkq(6,0,  8)
		call ckq(6,2,  9,10)     
		call ckq(6,4, 11,12)  
		call ckq(6,6, 13,14)          
           endif                                                                                             
	else if(symmetry.eq.2) then 
           if(j2.ge.2) then
		call rkq(2,0, 1)
		call rkq(2,2, 2)   
           endif
	   if(j2.ge.4) then
      		call rkq(4,0, 3)
		call rkq(4,2, 4)     
		call rkq(4,4, 5)  
           endif
	   if(j2.ge.6) then
      		call rkq(6,0, 6)
		call rkq(6,2, 7)     
		call rkq(6,4, 8) 
		call rkq(6,6, 9)          
           endif
	else if(symmetry.eq.3) then 
           if(j2.ge.2) then
		call rkq(2,0, 1)   
           endif
	   if(j2.ge.4) then
      		call rkq(4,0, 2)     
		call rkq(4,4, 3)  
           endif
	   if(j2.ge.6) then
      		call rkq(6,0, 4)   
		call ckq(6,4, 5,6)   
           endif
	else if(symmetry.eq.4) then          
           if(j2.ge.2) then
		call rkq(2,0, 1)   
           endif
	   if(j2.ge.4) then
      		call rkq(4,0, 2)     
		call rkq(4,4, 3)  
           endif
	   if(j2.ge.6) then
      		call rkq(6,0, 4)   
		call rkq(6,4, 5)   
           endif
	else if(symmetry.eq.5) then 
           if(j2.ge.2) then
		call rkq(2,0, 1)
           endif
	   if(j2.ge.4) then
		call rkq(4,0, 2)
	        call rkq(4,3, 3)
           endif
	   if(j2.ge.6) then
		call rkq(6,0, 4)
		call ckq(6,3, 5, 6)
		call ckq(6,6, 7, 8) 
           endif
	else if(symmetry.eq.6) then 
           if(j2.ge.2) then
		call rkq(2,0, 1)
           endif
	   if(j2.ge.4) then
		call rkq(4,0, 2)
	        call rkq(4,3, 3)
           endif
	   if(j2.ge.6) then
		call rkq(6,0, 4)
		call rkq(6,3, 5)
		call rkq(6,6, 6)             
           endif
	else if(symmetry.eq.7) then 
           if(j2.ge.2) then
		call rkq(2,0, 1)
           endif
	   if(j2.ge.4) then
		call rkq(4,0, 2)
           endif
	   if(j2.ge.6) then
		call rkq(6,0, 3)
		call rkq(6,6, 4)             
           endif
	else if(symmetry.eq.8) then 
	   if(j2.ge.4) then         
		call rkq( 4,0, 1)
		call rkqc(4,4, 1, 5)
	   endif
	   if(j2.ge.6) then
	        call rkq( 6,0, 2)
		call rkqc(6,4, 2, -21)
           endif
        endif
      return 
      end
!-----------------------------------
! calculates the v(i) from the vkq
!-----------------------------------
      subroutine itransform(v,vkq,symmetry,nre)
	implicit none
	integer symmetry,nre
	real*8 v(*)
	complex*16 vkq(0:6,0:6)
	integer k,q
	integer j2
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6      

	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
	j2 = dimj-1                                                         
                                                                              

	if(symmetry.eq.0) then
           if(j2.ge.2) then
!		vkq(2,0) = v(1)
!		vkq(2,1) = v(2)   
!		vkq(2,2) = cmplx(v(3),v(4))
                v(1) = real( vkq(2,0))
                v(2) = real( vkq(2,1))
                v(3) = real( vkq(2,2))
                v(4) = dimag(vkq(2,2))
           endif
	   if(j2.ge.4) then
!      		vkq(4,0) = v(5)
!		vkq(4,1) = cmplx(v( 6),v( 7))     
!		vkq(4,2) = cmplx(v( 8),v( 9))     
!		vkq(4,3) = cmplx(v(10),v(11))     
!		vkq(4,4) = cmplx(v(12),v(13))     
                v( 5)  = real( vkq(4,0))
                v( 6)  = real( vkq(4,1))
                v( 7)  = dimag(vkq(4,1))
                v( 8)  = real( vkq(4,2))
                v( 9)  = dimag(vkq(4,2))
                v(10)  = real( vkq(4,3))
                v(11)  = dimag(vkq(4,3))
                v(12)  = real( vkq(4,4))
                v(13)  = dimag(vkq(4,4))
           endif
	   if(j2.ge.6) then                                   
!		vkq(6,0) = v(14)
!		vkq(6,1) = cmplx(v(15),v(16))     
!		vkq(6,2) = cmplx(v(17),v(18))     
!		vkq(6,3) = cmplx(v(19),v(20))     
!		vkq(6,4) = cmplx(v(21),v(22))    
!		vkq(6,5) = cmplx(v(23),v(24))     
!		vkq(6,6) = cmplx(v(25),v(26))    
                v(14)  = real( vkq(6,0))
                v(15)  = real( vkq(6,1))
                v(16)  = dimag(vkq(6,1))
                v(17)  = real( vkq(6,2))
                v(18)  = dimag(vkq(6,2))
                v(19)  = real( vkq(6,3))
                v(20)  = dimag(vkq(6,3))
                v(21)  = real( vkq(6,4))
                v(22)  = dimag(vkq(6,4))
                v(23)  = real( vkq(6,5))
                v(24)  = dimag(vkq(6,5))
                v(25)  = real( vkq(6,6))
                v(26)  = dimag(vkq(6,6))
           endif                                    
	else if(symmetry.eq.1) then 
           if(j2.ge.2) then
!		vkq(2,0) = v(1)
!		vkq(2,2) = v(2)   
                v(1 )  = real(vkq(2,0))
                v(2 )  = real(vkq(2,2))
           endif
	   if(j2.ge.4) then
!     		vkq(4,0) = v(3)
!		vkq(4,2) = cmplx(v( 4),v( 5))     
!		vkq(4,4) = cmplx(v( 6),v( 7))  
                v( 3)  = real( vkq(4,0))
                v( 4)  = real( vkq(4,2))
                v( 5)  = dimag(vkq(4,2))
                v( 6)  = real( vkq(4,4))
                v( 7)  = dimag(vkq(4,4))
           endif
	   if(j2.ge.6) then
!     		vkq(6,0) = v(8)
!		vkq(6,2) = cmplx(v( 9),v(10))     
!		vkq(6,4) = cmplx(v(11),v(12))  
!		vkq(6,6) = cmplx(v(13),v(14))          
                v( 8)  = real( vkq(6,0))
                v( 9)  = real( vkq(6,2))
                v(10)  = dimag(vkq(6,2))
                v(11)  = real( vkq(6,4))
                v(12)  = dimag(vkq(6,4))
                v(13)  = real( vkq(6,6))
                v(14)  = dimag(vkq(6,6))
           endif                                                                                                 
	else if(symmetry.eq.2) then 
           if(j2.ge.2) then
!		vkq(2,0) = v(1)
!		vkq(2,2) = v(2)   
                v(1 )  = real(vkq(2,0))
                v(2 )  = real(vkq(2,2))
           endif
	   if(j2.ge.4) then
!      		vkq(4,0) = v(3)
!		vkq(4,2) = v(4)     
!		vkq(4,4) = v(5)  
                v( 3)  = real( vkq(4,0))
                v( 4)  = real( vkq(4,2))
                v( 5)  = real( vkq(4,4))
           endif
	   if(j2.ge.6) then
!      		vkq(6,0) = v(6)
!		vkq(6,2) = v(7)     
!		vkq(6,4) = v(8) 
!		vkq(6,6) = v(9)          
                v( 6)  = real( vkq(6,0))
                v( 7)  = real( vkq(6,2))
                v( 8)  = real( vkq(6,4))
                v( 9)  = real( vkq(6,6))
           endif
	else if(symmetry.eq.3) then 
           if(j2.ge.2) then
!		vkq(2,0) = v(1)   
                v(1 )  = real(vkq(2,0))
           endif
	   if(j2.ge.4) then
!      		vkq(4,0) = v(2)     
!		vkq(4,4) = v(3)  
                v( 2)  = real( vkq(4,0))
                v( 3)  = real( vkq(4,4))
           endif
	   if(j2.ge.6) then
!      		vkq(6,0) = v(4)   
!		vkq(6,4) = cmplx(v(5),v(6))   
                v( 4)  = real(  vkq(6,0))
                v( 5)  = real(  vkq(6,4))
                v( 6)  = dimag( vkq(6,4))
           endif
	else if(symmetry.eq.4) then          
           if(j2.ge.2) then
!		vkq(2,0) = v(1)   
                v(1 )  = real(vkq(2,0))
           endif
	   if(j2.ge.4) then
!      		vkq(4,0) = v(2)     
!		vkq(4,4) = v(3)  
                v( 2)  = real( vkq(4,0))
                v( 3)  = real( vkq(4,4))
           endif
	   if(j2.ge.6) then
!      		vkq(6,0) = v(4)   
!		vkq(6,4) = v(5)   
                v( 4)  = real(  vkq(6,0))
                v( 5)  = real(  vkq(6,4))
           endif
	else if(symmetry.eq.5) then 
           if(j2.ge.2) then
!		vkq(2,0) = v(1)
                v(1 )  = real(vkq(2,0))
           endif
	   if(j2.ge.4) then
!		vkq(4,0) = v(2)
!	        vkq(4,3) = v(3)
                v( 2)  = real( vkq(4,0))
                v( 3)  = real( vkq(4,3))
           endif
	   if(j2.ge.6) then
!		vkq(6,0) = v(4)
!		vkq(6,3) = cmplx(v(5),v(6))
!		vkq(6,6) = cmplx(v(7),v(8)) 
                v( 4)  = real( vkq(6,0))
                v( 5)  = real( vkq(6,3))
                v( 6)  = dimag(vkq(6,3))
                v( 7)  = real( vkq(6,6))
                v( 8)  = dimag(vkq(6,6))
           endif
	else if(symmetry.eq.6) then 
           if(j2.ge.2) then
!		vkq(2,0) = v(1)
                v(1 )  = real(vkq(2,0))
           endif
	   if(j2.ge.4) then
!		vkq(4,0) = v(2)
!	        vkq(4,3) = v(3)
                v( 2)  = real( vkq(4,0))
                v( 3)  = real( vkq(4,3))
           endif
	   if(j2.ge.6) then
!		vkq(6,0) = v(4)
!		vkq(6,3) = v(5)
!		vkq(6,6) = v(6)             
                v( 4)  = real( vkq(6,0))
                v( 5)  = real( vkq(6,3))
                v( 6)  = real( vkq(6,6))
           endif
	else if(symmetry.eq.7) then 
           if(j2.ge.2) then
!		vkq(2,0) = v(1)
                v(1 )  = real(vkq(2,0))
           endif
	   if(j2.ge.4) then
!		vkq(4,0) = v(2)
                v( 2)  = real( vkq(4,0))
           endif
	   if(j2.ge.6) then
!		vkq(6,0) = v(3)
!		vkq(6,6) = v(4)             
                v( 3)  = real( vkq(6,0))
                v( 4)  = real( vkq(6,6))
           endif
	else if(symmetry.eq.8) then 
	   if(j2.ge.4) then         
!		vkq(4,0) =     v(1)
!		vkq(4,4) =   5*v(1)
                v( 1)  = real( vkq(4,0))
	   endif
	   if(j2.ge.6) then
!	        vkq(6,0) =     v(2)
!		vkq(6,4) = -21*v(2)
                v( 2)  = real( vkq(6,0))
           endif
        endif
      return 
      end
!-----------------------------------
! calculates the v(i) from the vkq
!-----------------------------------
      subroutine ilim_trans(symmetry,nre)
	implicit none
	integer symmetry,nre

	integer k,q
	integer j2
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6      

	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
	j2 = dimj-1                                                         
                                                                              

	if(symmetry.eq.0) then
           if(j2.ge.2) then
!               v(1) = real( vkq(2,0))
!               v(2) = real( vkq(2,1))
!               v(3) = real( vkq(2,2))
!               v(4) = dimag(vkq(2,2))
                call vr(1, 2,0)
                call vr(2, 2,1)
                call vr(3, 2,2)
                call vi(4, 2,2)
               

           endif
	   if(j2.ge.4) then
!               v( 5)  = real( vkq(4,0))
!               v( 6)  = real( vkq(4,1))
!               v( 7)  = dimag(vkq(4,1))
!               v( 8)  = real( vkq(4,2))
!               v( 9)  = dimag(vkq(4,2))
!               v(10)  = real( vkq(4,3))
!               v(11)  = dimag(vkq(4,3))
!               v(12)  = real( vkq(4,4))
!               v(13)  = dimag(vkq(4,4))
                call vr( 5 , 4,0) 
                call vr( 6 , 4,1)
                call vi( 7 , 4,1) 
                call vr( 8 , 4,2) 
                call vi( 9 , 4,2) 
                call vr(10 , 4,3) 
                call vi(11 , 4,3) 
                call vr(12 , 4,4) 
                call vi(13 , 4,4) 

           endif
	   if(j2.ge.6) then                                   
!               v(14)  = real( vkq(6,0))
!               v(15)  = real( vkq(6,1))
!               v(16)  = dimag(vkq(6,1))
!               v(17)  = real( vkq(6,2))
!               v(18)  = dimag(vkq(6,2))
!               v(19)  = real( vkq(6,3))
!               v(20)  = dimag(vkq(6,3))
!               v(21)  = real( vkq(6,4))
!               v(22)  = dimag(vkq(6,4))
!               v(23)  = real( vkq(6,5))
!               v(24)  = dimag(vkq(6,5))
!               v(25)  = real( vkq(6,6))
!               v(26)  = dimag(vkq(6,6))
                call vr(14 , 6,0)
                call vr(15 , 6,1)
                call vi(16 , 6,1)
                call vr(17 , 6,2)
                call vi(18 , 6,2)
                call vr(19 , 6,3)
                call vi(20 , 6,3)
                call vr(21 , 6,4)
                call vi(22 , 6,4)
                call vr(23 , 6,5)
                call vi(24 , 6,5)
                call vr(25 , 6,6)
                call vi(26 , 6,6)

           endif                                    
	else if(symmetry.eq.1) then 
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
!               v(2 )  = real(vkq(2,2))
                call vr(1 , 2,0)
                call vr(2 , 2,2)
           endif
	   if(j2.ge.4) then
!               v( 3)  = real( vkq(4,0))
!               v( 4)  = real( vkq(4,2))
!               v( 5)  = dimag(vkq(4,2))
!               v( 6)  = real( vkq(4,4))
!               v( 7)  = dimag(vkq(4,4))
                call vr( 3 , 4,0) 
                call vr( 4 , 4,2)
                call vi( 5 , 4,2)
                call vr( 6 , 4,4)
                call vi( 7 , 4,4)
           endif
	   if(j2.ge.6) then
!               v( 8)  = real( vkq(6,0))
!               v( 9)  = real( vkq(6,2))
!               v(10)  = dimag(vkq(6,2))
!               v(11)  = real( vkq(6,4))
!               v(12)  = dimag(vkq(6,4))
!               v(13)  = real( vkq(6,6))
!               v(14)  = dimag(vkq(6,6))
                call vr( 8 , 6,0)
                call vr( 9 , 6,2)
                call vi(10 , 6,2)
                call vr(11 , 6,4)
                call vi(12 , 6,4)
                call vr(13 , 6,6)
                call vi(14 , 6,6)
           endif                                                                                                 
	else if(symmetry.eq.2) then 
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
!               v(2 )  = real(vkq(2,2))
                call vr(1 , 2,0)
                call vr(2 , 2,2)
           endif
	   if(j2.ge.4) then
!               v( 3)  = real( vkq(4,0))
!               v( 4)  = real( vkq(4,2))
!               v( 5)  = real( vkq(4,4))
                call vr( 3 , 4,0)
                call vr( 4 , 4,2)
                call vr( 5 , 4,4)
           endif
	   if(j2.ge.6) then
!               v( 6)  = real( vkq(6,0))
!               v( 7)  = real( vkq(6,2))
!               v( 8)  = real( vkq(6,4))
!               v( 9)  = real( vkq(6,6))
                call vr( 6 , 6,0)
                call vr( 7 , 6,2)
                call vr( 8 , 6,4)
                call vr( 9 , 6,6)
           endif
	else if(symmetry.eq.3) then 
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
                call vr(1 , 2,0)
           endif
	   if(j2.ge.4) then
!               v( 2)  = real( vkq(4,0))
!               v( 3)  = real( vkq(4,4))
                call vr( 2 , 4,0)
                call vr( 3 , 4,4)
           endif
	   if(j2.ge.6) then
!               v( 4)  = real(  vkq(6,0))
!               v( 5)  = real(  vkq(6,4))
!               v( 6)  = dimag( vkq(6,4))
                call vr( 4 , 6,0)
                call vr( 5 , 6,4)
                call vi( 6 , 6,4)
           endif
	else if(symmetry.eq.4) then          
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
                call vr(1 , 2,0)
           endif
	   if(j2.ge.4) then
!               v( 2)  = real( vkq(4,0))
!               v( 3)  = real( vkq(4,4))
                call vr( 2 , 4,0)
                call vr( 3 , 4,4)
           endif
	   if(j2.ge.6) then
!               v( 4)  = real(  vkq(6,0))
!               v( 5)  = real(  vkq(6,4))
		call vr(4, 6,0)
		call vr(5, 6,4)
           endif
	else if(symmetry.eq.5) then 
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
                call vr(1 , 2,0)
           endif
	   if(j2.ge.4) then
!               v( 2)  = real( vkq(4,0))
!               v( 3)  = real( vkq(4,3))
                call vr( 2 , 4,0)
                call vr( 3 , 4,3)
           endif
	   if(j2.ge.6) then
!               v( 4)  = real( vkq(6,0))
!               v( 5)  = real( vkq(6,3))
!               v( 6)  = dimag(vkq(6,3))
!               v( 7)  = real( vkq(6,6))
!               v( 8)  = dimag(vkq(6,6))
                call vr( 4 , 6,0)
                call vr( 5 , 6,3)
                call vi( 6 , 6,3)
                call vr( 7 , 6,6)
                call vi( 8 , 6,6)
           endif
	else if(symmetry.eq.6) then 
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
                call vr(1 , 2,0)
           endif
	   if(j2.ge.4) then
!               v( 2)  = real( vkq(4,0))
!               v( 3)  = real( vkq(4,3))
                call vr( 2 , 4,0)
                call vr( 3 , 4,3)
           endif
	   if(j2.ge.6) then
!               v( 4)  = real( vkq(6,0))
!               v( 5)  = real( vkq(6,3))
!               v( 6)  = real( vkq(6,6))
                call vr( 4 , 6,0)
                call vr( 5 , 6,3)
                call vr( 6 , 6,6)
           endif
	else if(symmetry.eq.7) then 
           if(j2.ge.2) then
!               v(1 )  = real(vkq(2,0))
                call vr(1 , 2,0)
           endif
	   if(j2.ge.4) then
!               v( 2)  = real( vkq(4,0))
                call vr( 2 , 4,0)
           endif
	   if(j2.ge.6) then
!               v( 3)  = real( vkq(6,0))
!               v( 4)  = real( vkq(6,6))
                call vr( 3 , 6,0)
                call vr( 4 , 6,6)
           endif
	else if(symmetry.eq.8) then 
	   if(j2.ge.4) then         
!		vkq(4,0) =     v(1)
!		vkq(4,4) =   5*v(1)
!               v( 1)  = real( vkq(4,0))
                call vr( 1 , 4,0)
	   endif
	   if(j2.ge.6) then
!	        vkq(6,0) =     v(2)
!		vkq(6,4) = -21*v(2)
!               v( 2)  = real( vkq(6,0))
                call vr( 2 , 6,0)
           endif
        endif
      return 
      end
c---------------------------------------------------------------
c prints the calculated excitations out for a given temperature
c---------------------------------------------------------------
      subroutine p_excitations(out,de,di)
	implicit none

 	include 'CF_SOURCES:FOCUS.INC'
 	include 'CF_SOURCES:CONTROL.INC'

	integer out
	real*8  de,di
	real*8  c_fmevkelvin
	integer i
        character*18 text_line
        

  1     format(1x,a29,f7.2,a7)
  2     format(1x,a37,f4.2,a39) 
  3     format(1x,a7,f4.2,a62)	     
   
        if(setprint.eq.1) call section (out,' Crystal Electric Field Excitations ')
	write(out,1)' The sample temperature is : ',temp,' Kelvin' 
	write(out,2)' Energy levels which are closer than ',de,' Kelvin are assumed to  be degenerated.'  
 	write(out,*)' Only  those  excitations are printed out whose intensities are greater or'
        write(out,3)' equal ',di,' barn. The excitations are ordered by decreasing intensities.'                 
        write(out,*)
	
  4	format(1x, i3,' : ',f6.2,' meV',3x,f6.2,' barn',3x,a18)       
  5	format(1x, i3,' : ',f16.2,' meV',3x,f6.2,' barn',3x,a18)       
	do i=1,n_excitation
           if(e_excitation(i)/c_fmevkelvin().lt.1.0d3) then
              write(out,4) i,e_excitation(i)/c_fmevkelvin(),i_excitation(i), text_line(lineshape(i))
           else
              write(out,5) i,e_excitation(i)/c_fmevkelvin(),i_excitation(i), text_line(lineshape(i))
           end if
	end do               
        write(out,*)               
                       
      return
      end
c-----------------------------------------------------
c calculates the crystal field excitations for a given
c crystal field calculation
c the input parameter are the Vkq or the Bkq
c----------------------------------------------------- 
      subroutine c_excitations(out,nre,symmetry,bkq,sbkq,bmol,bext,temp,
     *                         eex,iex,nex,type,de,di)                         
	implicit none
	integer out
	real*8  eex(17*17),iex(17*17)
	integer nex                     
	include 'CF_SOURCES:CF_FABI.INC'   
c	real*8  e_excitation(17*17) ! energy position [kelvin] of the excit.
c	real*8  i_excitation(17*17) ! intensity [barn] of the excitations
c	integer n_excitation        ! number of excitations found          
	integer degeneration(17*17)  ! stores the degeneration of a level
 	integer n_energies        ! no. of degenerated energy levels
	real*8  e_energies(17*17) ! energy values of the degenerated energy levels
	real*8  i_energies(17,17) ! intensities of the degenerated energy levels      
                                                                                                             
        real*8 de,di

c       ------------------------------------------------- 
c	energy levels which are closer than de are assumed to be degenerated
c	only those excitations are taken into account whose intensities are 
c	greater equal than di
c       ------------------------------------------------- 
        integer i,k,ii,kk
	integer no       ! is a function which is called
c       ------------------------------------------------- 
	external M01DAF   ! NAG - sorts real numbers 
	external M01ZAF 
	external M01CAF
	integer ifail,ind(17*17)
c       --------------------------------------------------------

	call c_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,type) 
                                                                                
c       ------------------------------------------------- 
c       find out how many degenerated energy levels exists
c       store the intensities of the degenarated levels
c       ------------------------------------------------- 
        call deg_on(energy,intensity,dimj,degeneration,e_energies,i_energies,n_energies,de,di)   
c       ------------------------------------------------- 
c       calculate all possible excitations
c       ------------------------------------------------- 
	nex=0
	do i=1,n_energies
	   do k=1,n_energies
	      nex= nex + 1
              eex(nex) = e_energies(k)-e_energies(i)
	      iex(nex) = i_energies(i,k)             !I(i->k)
	      call ifnull(eex(nex))
 	   end do
	end do
c       ------------------------------------------------- 
c       print the result out
c       ------------------------------------------------- 
c	do i=1,nex
c	   write(6,*) i,eex(i),iex(i)
c	end do               
c	write(6,*)                                             
c       ------------------------------------------------- 
c       sort energies in ascending order
c       ------------------------------------------------- 
        ifail=1
	call M01DAF(eex,1,nex,'A',ind,ifail)
	call M01ZAF(ind,1,nex,ifail)
	call M01CAF(eex,1,nex,'A',ifail)      
c       ------------------------------------------------- 
c       print the result out
c       ------------------------------------------------- 
c	do i=1,nex
c	   write(6,*) i,eex(i),iex(ind(i))
c	end do                             
c       ------------------------------------------------- 
c       finds out how many non degenerated EXCITATIONS exists
c       ------------------------------------------------- 
	e_excitation(1)   = eex(1)
	degeneration(1) = 1
	k=1
	do i=2,nex
	   if( abs(e_excitation(k)-eex(i)).ge.de ) then
               k=k+1
	       e_excitation(k) = eex(i)
	       degeneration(k) = 1
	   else
	       degeneration(k)=degeneration(k)+1
	   endif
	end do
        n_excitation = k                   
c       ------------------------------------------------- 
c       store the intensities of the degenarated EXCITATIONS
c       ------------------------------------------------- 
	do i=1,17*17
           i_excitation(i) = 0.0d0
	end do
	do i=1,nex
	   ii =  no(i,degeneration,n_excitation)  
	   i_excitation(ii)=i_excitation(ii)+iex(ind(i))
	end do                                               
c       ------------------------------------------------- 
c       print the result out
c       ------------------------------------------------- 
c	do i=1,n_excitation
c	   write(6,*) i,e_excitation(i),i_excitation(i)
c	end do                             
c       ------------------------------------------------- 
c       sort intensities in descending order
c       ------------------------------------------------- 
        ifail=1
	call M01DAF(i_excitation,1,nex,'D',ind,ifail)
	call M01ZAF(ind,1,nex,ifail)
	call M01CAF(i_excitation,1,nex,'D',ifail)
c       ------------------------------------------------- 
c       print the result out
c       ------------------------------------------------- 
c	do i=1,n_excitation
c	   write(6,*) i,e_excitation(ind(i)),i_excitation(i)
c	end do                             
c       ------------------------------------------------- 
c       store the output values
c       ------------------------------------------------- 
	k=0	
        do i=1,n_excitation
	   if(i_excitation(i).ge.di .or. dimj.eq.1) then
           	k=k+1
	        iex(k)=i_excitation(i)
	        eex(k)=e_excitation(ind(i))
	   endif
	end do
        nex = k
c       ------------------------------------------------- 
c       print the result out
c       ------------------------------------------------- 
c	do i=1,nex
c	   write(6,*) i,eex(i),iex(i)
c	end do 
      return
      end
c-------------------------------------------------------------------------
c 
c-------------------------------------------------------------------------
      subroutine deg_on(energy,mat,dimj,degeneration,e_energies,i_energies,n_energies,de,di)                         
	implicit none
        real*8  energy(17)           ! already defined in CF_FABI.INC
        real*8  mat(17,17)
        real*8  dimj
	integer degeneration(17*17)  ! stores the degeneration of a level
 	integer n_energies           ! no. of degenerated energy levels
	real*8  e_energies(17*17)    ! energy values of the degenerated energy levels
	real*8  i_energies(17,17)    ! intensities of the degenerated energy levels      
        real*8 de,di

c	energy levels which are closer than de are assumed to be degenerated
c	only those excitations are taken into account whose intensities are 
c	greater equal than di

        integer i,k,ii,kk
	integer no       ! is a function which is called

                                                                                
c find out how many degenerated energy levels exists

	e_energies(1)   = 0.0d0
	degeneration(1) = 1
	k=1
	do i=2,dimj
	   if( abs(e_energies(k)-energy(i)).ge.de ) then
               k=k+1
	       e_energies(k)   = energy(i)
	       degeneration(k) = 1
	   else
	       degeneration(k)=degeneration(k)+1
	   endif
	end do
        n_energies = k                   

c store the intensities of the degenarated levels

	call clear(i_energies)
	do i=1,dimj
	   ii =  no(i,degeneration,n_energies)  
	   do k=1,dimj
	      kk = no(k,degeneration,n_energies) 
	      i_energies(ii,kk)=i_energies(ii,kk)+mat(i,k)	
	   end do
	end do                                                          
                                                            
      return
      end
c-------------------------
c transforms the indices
c-------------------------
      integer function no(i,d,n)  
	implicit none
	integer i,d(17*17),n,up,lo
	up = 0
	do no=1,n
	   lo = up+1
	   up = lo+d(no)-1
	   if(i.ge.lo.and.i.le.up) return
	end do
      return
      end

      integer function ro(i,d,n)  
	implicit none
	integer i,d(17*17),n,up,lo,noo
	up = 0
        ro = 0
	do noo=1,n
	   lo = up+1
	   up = lo+d(noo)-1
	   if(i.ge.lo.and.i.le.up) then 
              ro = i-lo+1
              return
           end if
	end do
      return
      end

c ---------------------------------------------------------
c  calculates the x and W from the given Bkq in kelvin
c ---------------------------------------------------------
c we have H=b40o40 + b44o44 + b60o60 + b64o64
c in the cubic case: b44=5*b40 and b64=-21*b60
c this gives       : H=b40(o40+5o44) + b60(o60-21o64)
c if we now define x and W values by:
c     B40 =       x*W/f4
c     B60 = (1-|x|)*W/f6
c
c then the representation after Lea,Leask and Wolf follows.
c
      subroutine c_xW(nre,x,w,bkq) 
	implicit none
        complex*16 bkq(0:6,0:6)
	real*8 x,w,f4,f6,c_f4,c_f6,b4,b6
	integer nre
	b4=bkq(4,0)
	b6=bkq(6,0)
	f4=c_f4(nre)
	f6=c_f6(nre)  
	call ifnull(b4)
	call ifnull(b6)
        x = 0.0d0 
	w = 0.0d0
	if(b4.eq.0.0d0.and.b6.eq.0.0d0) return
	x = 1.0d0                                                     
	w = b4*f4/x
	if(b4.ne.0.0d0.and.b6.eq.0.0d0) return
	x = 0.0d0
	W = B6*f6
	if(b4.eq.0.0d0.and.b6.ne.0.0d0) return
        if(b4*f4.eq.b6*f6) then
		x=0.5
	else
		if(b4*f4/b6/f6.lt.0.0d0) then
			x=b4*f4/(b6*f6-b4*f4)
                else
			x=b4*f4/(b6*f6+b4*f4)
                endif
	endif
	w = b4*f4/x               
      return
      end
c ---------------------------------------------------------
c  transformes the given x and W in Kelvin to Bkq in kelvin
c ---------------------------------------------------------
      subroutine t_xW(nre,x,w,bkq) 
	implicit none
        complex*16 bkq(0:6,0:6)
	real*8 x,w,f4,f6,c_f4,c_f6
	integer nre,k,q             
	real*8  gj,dimj,alphaj,betaj,gammaj,r2,r4,r6   
c       ------------------------------------------------- 
c       initialize some rare earth constants
c       ------------------------------------------------- 
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  

c       ------------------------------------------------- 
c       clear all Bkq's
c       ------------------------------------------------- 
 	do k=0,6
	   do q=0,k
	      bkq(k,q)=0.0d0
           end do
        end do

c       ------------------------------------------------- 
c       make sure that W is the only one parameter in a
c       cerium like case
c       ------------------------------------------------- 
        if( dimj-1.ge.4 .and. dimj-1.lt.6 ) x = 1.0d0

c       ------------------------------------------------- 
c       set the b40, b44, b60, b64
c       ------------------------------------------------- 
	f4 = c_f4(nre)
	bkq(4,0) = x*w/f4
	bkq(4,4) = 5*bkq(4,0)
	f6 = c_f6(nre)
	bkq(6,0) = (1-abs(x))*w/f6
	bkq(6,4) = -21*bkq(6,0)

      return
      end
c----------------------------------------------------------------------------
c just for an information: the f4,f6 values given below
c may be also defined by :
c
c  F(K) := F  (J) := Greatest_common_divisor({ <jM'|STEV  |JM> } )
c           k0                                                  M',M=-J...J
c
c  this definition differs from the usual one only for 
c  f4(Nd3+)=84 instead of 60, f4(Ho3+)=420 instead of 60
c  and f6(Nd3+)=5040 instead of 2540. 
c
c
c---------------------------
c returns f4 for a given nre 
c---------------------------
      real*8 function c_f4(nre)
	implicit none
	integer nre
	real*8 f4(13)
	f4( 1) = 60
	f4( 2) = 60
	f4( 3) = 60 ! 84
	f4( 4) = 60
	f4( 5) = 60
	f4( 6) =  1
	f4( 7) = 60
	f4( 8) = 60
	f4( 9) = 60
	f4(10) = 60 ! 420
	f4(11) = 60
	f4(12) = 60
	f4(13) = 60  
	c_f4 = f4(nre)
      return
      end  
c---------------------------
c returns f6 for a given nre 
c---------------------------
      real*8 function c_f6(nre)
	implicit none
	integer nre
	real*8 f6(13)
	f6( 1) =     1
	f6( 2) =  1260
	f6( 3) =  2520 ! 5040
	f6( 4) =  1260
	f6( 5) =     1
	f6( 6) =     1
	f6( 7) =  1260
	f6( 8) =  7560
	f6( 9) = 13860
	f6(10) = 13860
	f6(11) = 13860
	f6(12) =  7560
	f6(13) =  1260
	c_f6 = f6(nre)
      return
      end  
c--------------------------------------
c calculates the bkq from the given akq
c--------------------------------------
      subroutine c_akq_bkq(nre,akq,bkq) 
	implicit none
        complex*16 akq(0:6,0:6),bkq(0:6,0:6)
	integer nre  
	real*8  gj,dimj,alphaj,betaj,gammaj,r2,r4,r6   
c       ------------------------------------------------- 
c       initialize some rare earth constants
c       ------------------------------------------------- 
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  
c
	bkq(2,0) = alphaj*r2*akq(2,0)
	bkq(2,1) = alphaj*r2*akq(2,1)
	bkq(2,2) = alphaj*r2*akq(2,2)
c
	bkq(4,0) =  betaj*r4*akq(4,0) 
	bkq(4,1) =  betaj*r4*akq(4,1) 
	bkq(4,2) =  betaj*r4*akq(4,2) 
	bkq(4,3) =  betaj*r4*akq(4,3) 
	bkq(4,4) =  betaj*r4*akq(4,4) 
c
	bkq(6,0) = gammaj*r6*akq(6,0)
	bkq(6,1) = gammaj*r6*akq(6,1)
	bkq(6,2) = gammaj*r6*akq(6,2)
	bkq(6,3) = gammaj*r6*akq(6,3)
	bkq(6,4) = gammaj*r6*akq(6,4)
	bkq(6,5) = gammaj*r6*akq(6,5)
	bkq(6,6) = gammaj*r6*akq(6,6)                            
      return
      end                                           
c--------------------------------------
c calculates the akq from the given bkq
c--------------------------------------
      subroutine c_bkq_akq(nre,akq,bkq) 
	implicit none
        complex*16 akq(0:6,0:6),bkq(0:6,0:6)
	integer nre  
	real*8  gj,dimj,alphaj,betaj,gammaj,r2,r4,r6   
        integer k,l
c       ------------------------------------------------- 
c       initialize some rare earth constants
c       ------------------------------------------------- 
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)  

        if(alphaj.ne.0.0d0) then
        	akq(2,0) = bkq(2,0)/alphaj/r2
		akq(2,1) = bkq(2,1)/alphaj/r2
		akq(2,2) = bkq(2,2)/alphaj/r2
        end if
        if(betaj.ne.0.0d0) then
		akq(4,0) = bkq(4,0)/betaj/r4 
		akq(4,1) = bkq(4,1)/betaj/r4 
		akq(4,2) = bkq(4,2)/betaj/r4 
       		akq(4,3) = bkq(4,3)/betaj/r4
       		akq(4,4) = bkq(4,4)/betaj/r4
	end if
        if(gammaj.ne.0.0d0) then
		akq(6,0) = bkq(6,0)/gammaj/r6
      		akq(6,1) = bkq(6,1)/gammaj/r6
       		akq(6,2) = bkq(6,2)/gammaj/r6
       		akq(6,3) = bkq(6,3)/gammaj/r6
       		akq(6,4) = bkq(6,4)/gammaj/r6
       		akq(6,5) = bkq(6,5)/gammaj/r6
       		akq(6,6) = bkq(6,6)/gammaj/r6
        end if
      return
      end                                           
c----------------------------------------------------------------
c determines the name of the R3+ ion
c----------------------------------------------------------------
      character*4 function r_name(nre)
	implicit none
        integer nre
	if(nre.eq. 1) r_name = 'Ce3+'   
	if(nre.eq. 2) r_name = 'Pr3+'   
	if(nre.eq. 3) r_name = 'Nd3+'   
	if(nre.eq. 4) r_name = 'Pm3+'   
	if(nre.eq. 5) r_name = 'Sm3+' 
	if(nre.eq. 6) r_name = 'Eu3+'   
	if(nre.eq. 7) r_name = 'Gd3+'   
	if(nre.eq. 8) r_name = 'Tb3+'   
	if(nre.eq. 9) r_name = 'Dy3+'   
	if(nre.eq.10) r_name = 'Ho3+'    
	if(nre.eq.11) r_name = 'Er3+'   
	if(nre.eq.12) r_name = 'Tm3+'   
	if(nre.eq.13) r_name = 'Yb3+'  
      return
      end
c----------------------------------------------------------------
c chooses the kind of operator type 
c----------------------------------------------------------------
      character*24 function otype(type)
	implicit none
	real*8 type
	otype = 'imaginary part chosen of'
        if(type.eq.1.0d0)  otype = 'real      part chosen of'
      return
      end
c----------------------------------------------------------------
c prints the input out
c----------------------------------------------------------------
      subroutine p_input(out,symmetry,nre,bkq,sbkq,bmol,bext,temp)
	
	integer out,symmetry,nre
	complex*16 bkq(0:6,0:6)
	real*8    sbkq(0:6,0:6)
	real*8    bmol(3),bext(3),temp
	character*4 r_name
	character*24 otype
	real*8 c_fmevkelvin,c

	c=c_fmevkelvin()

	call section(out, ' THE INPUT OF CF_FABI IS: ')
 100  	format(     1x,a9,1pe14.6,' + ',1pe14.6,a10)     
	write(out,'(1x,a24,i2   )')' The symmetry number is :',symmetry 
	write(out,'(1x,a24,1x,a4)')' The calculation is for :',r_name(nre)
	write(out,*)
	write(out,100) ' B20  =  ',bkq(2,0)/c,'*i meV' 
	write(out,100) ' B21  =  ',bkq(2,1)/c,'*i meV' 
	write(out,100) ' B22  =  ',bkq(2,2)/c,'*i meV' 
	write(out,100) ' B40  =  ',bkq(4,0)/c,'*i meV' 
	write(out,100) ' B41  =  ',bkq(4,1)/c,'*i meV' 
	write(out,100) ' B42  =  ',bkq(4,2)/c,'*i meV' 
	write(out,100) ' B43  =  ',bkq(4,3)/c,'*i meV' 
	write(out,100) ' B44  =  ',bkq(4,4)/c,'*i meV' 
	write(out,100) ' B60  =  ',bkq(6,0)/c,'*i meV' 
	write(out,100) ' B61  =  ',bkq(6,1)/c,'*i meV' 
	write(out,100) ' B62  =  ',bkq(6,2)/c,'*i meV' 
	write(out,100) ' B63  =  ',bkq(6,3)/c,'*i meV' 
	write(out,100) ' B64  =  ',bkq(6,4)/c,'*i meV'   
	write(out,100) ' B65  =  ',bkq(6,5)/c,'*i meV' 
	write(out,100) ' B66  =  ',bkq(6,6)/c,'*i meV'   
        write(out,*)
 101  	format(1x,a9,f3.0,' : ',a24,' O(',i1,',',i1,')')
	write(out,101) ' sB20 =  ',sbkq(2,0),otype(sbkq(2,0)),2,0
	write(out,101) ' sB21 =  ',sbkq(2,1),otype(sbkq(2,1)),2,1 
 	write(out,101) ' sB22 =  ',sbkq(2,2),otype(sbkq(2,2)),2,2 
	write(out,101) ' sB40 =  ',sbkq(4,0),otype(sbkq(4,0)),4,0 
 	write(out,101) ' sB41 =  ',sbkq(4,1),otype(sbkq(4,1)),4,1 
  	write(out,101) ' sB42 =  ',sbkq(4,2),otype(sbkq(4,2)),4,2 
	write(out,101) ' sB43 =  ',sbkq(4,3),otype(sbkq(4,3)),4,3 
	write(out,101) ' sB44 =  ',sbkq(4,4),otype(sbkq(4,4)),4,4 
	write(out,101) ' sB60 =  ',sbkq(6,0),otype(sbkq(6,0)),6,0 
     	write(out,101) ' sB61 =  ',sbkq(6,1),otype(sbkq(6,1)),6,1 
	write(out,101) ' sB62 =  ',sbkq(6,2),otype(sbkq(6,2)),6,2 
	write(out,101) ' sB63 =  ',sbkq(6,3),otype(sbkq(6,3)),6,3 
	write(out,101) ' sB64 =  ',sbkq(6,4),otype(sbkq(6,4)),6,4  
	write(out,101) ' sB65 =  ',sbkq(6,5),otype(sbkq(6,5)),6,5 
	write(out,101) ' sB66 =  ',sbkq(6,6),otype(sbkq(6,6)),6,6  
        write(out,*)
	write(out,'(1x,a9,f9.3,a6)') ' Bmol_x :',bmol(1),' Tesla' 
	write(out,'(1x,a9,f9.3,a6)') ' Bmol_y :',bmol(2),' Tesla'     
	write(out,'(1x,a9,f9.3,a6)') ' Bmol_z :',bmol(3),' Tesla'     
        write(out,*)
	write(out,'(1x,a9,f9.3,a6)') ' Bext_x :',bext(1),' Tesla' 
	write(out,'(1x,a9,f9.3,a6)') ' Bext_y :',bext(2),' Tesla'     
	write(out,'(1x,a9,f9.3,a6)') ' Bext_z :',bext(3),' Tesla' 
	write(out,*)                                                                 
	write(out,'(1x,a28,f6.1,a7)') ' The sample temperature is: ',temp   ,' Kelvin'         
        write(out,*)
	call section(out, ' THE OUTPUT OF CF_FABI IS: ')
      return
      end
c-----------------------------------------------
c prints an overview of the implemented R3+ ions
c-----------------------------------------------
      subroutine p_rare_earth(out)
	implicit none
	integer out
	call section(out,' IMPLEMENTED RARE EARTH IONS ')
        write(out,*) '  ------------------------------------------------------------------------'
        write(out,*) ' |R3+|1=Ce|2=Pr|3=Nd|4=Pm|5=Sm|6=Eu|7=Gd|8=Tb|9=Dy|10=Ho|11=Er|12=Tm|13=Yb|'
        write(out,*) '  ------------------------------------------------------------------------'
        write(out,*) ' Table of implemented rare earth ions R3+ and their relationship to the '
        write(out,*) ' rare earth number.'                                                           
        write(out,*) 
      return
      end
c---------------------------------------------
c prints an overview over the symmetry numbers
c---------------------------------------------
      subroutine p_symmetry_numbers(out)
	implicit none
	integer out

	call section(out,' POINT SYMMETRIES OF CRYSTAL FIELD AND SYMMETRY NUMBERS ')
	write(out,*)'  --------------------------------------------------- ' 
	write(out,*)' | point symmetry of crystal field | symmetry number |'
	write(out,*)' |---------------------------------------------------|'
	write(out,*)' | Ci C1 ---------------------------------- 0 -------|'
	write(out,*)' | C2 Cs C2h ------------------------------ 1 -------|'
	write(out,*)' | C2v D2 D2h ----------------------------- 2 -------|'
	write(out,*)' | C4 S4 C4h ------------------------------ 3 -------|'
	write(out,*)' | D4 C4v D2d D4h ------------------------- 4 -------|'
	write(out,*)' | C3 S6 ---------------------------------- 5 -------|'
	write(out,*)' | D3 C3v D3d ----------------------------- 6 -------|'
	write(out,*)' | C6 C3h C6h D6 C6v D3h D6h -------------- 7 -------|'   
	write(out,*)' | T Td Th O Oh --------------------------- 8 -------|'
	write(out,*)'  --------------------------------------------------- '   
	write(out,*)' How symmetry numbers are connected to the point symmetry'
        write(out,*)' of the crystal field of the R3+ ion.'                                  
	write(out,*)
      return
      end
c---------------------------------------------
c prints an overview over the hamiltonian
c---------------------------------------------
      subroutine p_hamiltonian(out,symmetry)
	implicit none
	integer out,symmetry
	call section(out,' THE DEFINITION OF THE CRYSTAL FIELD HAMILTONIAN ')
	write(out,*)'    ----             ----     ----              '
	write(out,*)'    \       ^      1 \        \      *    ^     '
	write(out,*)' H = >  B   O(J) + -  >        >    B     O (J) '
	write(out,*)'    /    k0  k0    2 /        /      kq    kq   '
	write(out,*)'    ----             ----     ----              '
	write(out,*)'   k=2,4,6          k=2,4,6  0<|q|<=k           '
	write(out,*)                                                 
	write(out,*)                  '  --------------------------------------------------- '
	write(out,*)                  ' |S\B|20|21|22||40|41|42|43|44||60|61|62|63|64|65|66||'
	write(out,*)                  ' |---|--|--|--||--|--|--|--|--||--|--|--|--|--|--|--||'
	if(symmetry.eq.0) write(out,*)' | 0 |r |r |c ||r |c |c |c |c ||r |c |c |c |c |c |c ||' 
        if(symmetry.eq.1) write(out,*)' | 1 |r |  |r ||r |  |c |  |c ||r |  |c |  |c |  |c ||'
        if(symmetry.eq.2) write(out,*)' | 2 |r |  |r ||r |  |r |  |r ||r |  |r |  |r |  |r ||'
        if(symmetry.eq.3) write(out,*)' | 3 |r |  |  ||r |  |  |  |r ||r |  |  |  |c |  |  ||'
        if(symmetry.eq.4) write(out,*)' | 4 |r |  |  ||r |  |  |  |r ||r |  |  |  |r |  |  ||'
        if(symmetry.eq.5) write(out,*)' | 5 |r |  |  ||r |  |  |r |  ||r |  |  |c |  |  |c ||'
        if(symmetry.eq.6) write(out,*)' | 6 |r |  |  ||r |  |  |r |  ||r |  |  |r |  |  |r ||'
        if(symmetry.eq.7) write(out,*)' | 7 |r |  |  ||r |  |  |  |  ||r |  |  |  |  |  |r ||'
        if(symmetry.eq.8) write(out,*)' | 8 |  |  |  ||r |  |  |  |r ||r |  |  |  |r |  |  ||'

	write(out,*)'  --------------------------------------------------- '
	write(out,*)' B = Bkq-parameters and S = symmetry number'
	write(out,*)' Used System : crystallographic system     '
        write(out,*)

	write(out,*)'A table which shows you which crystal field parameter is complex (c) or real (r)'
  88    format(a34,i2,a1)
	write(out,88)'for the given symmetry number of ',symmetry,'.' 

	if(symmetry.eq.8) then 
		write(out,*)' The  B44  parameter is connected with the  B40  parameter and the B64'      
       		write(out,*)' is connected with the B60 one: B44 = 5*B40 and B64 = -21*B60.'
	end if
	write(out,*)
      return
      end 
c--------------------------
c prints out the R3+ moment
c--------------------------
      subroutine p_moments(out,mx,my,mz,gj,dimj)
	implicit none
	integer out
	real*8 mx,my,mz,gj,dimj,j,ms,mr
	real*8 m,phi,theta
 1      format(2x,'mx   = ',f6.2,' myb') 
 2      format(2x,'my   = ',f6.2,' myb') 
 3      format(2x,'mz   = ',f6.2,' myb') 
 21     format(2x,'mx/mr= ',f6.2,'%') 
 22     format(2x,'my/mr= ',f6.2,'%') 
 23     format(2x,'mz/mr= ',f6.2,'%') 
	
        call section (out, ' THE RARE EARTH MOMENT')
	write(out,1) mx
	write(out,2) my
	write(out,3) mz
	write(out,*) ' The R3+ moment in units of the full measurable R3+ moment mr = gj*myB*j'
	j = 0.5d0*(dimj-1)
	mr = gj*j
	if(mr.eq.0.0d0) mr=1.0d0
	write(out,21) mx/mr*100 
	write(out,22) my/mr*100  
	write(out,23) mz/mr*100
	write(out,*)  
        write(out,*) ' if we choose:  mx = m*sin(theta)*cos(phi)'
	write(out,*) '                my = m*sin(theta)*sin(phi)'
	write(out,*) '                mz = m*cos(theta)'
	write(out,*) ' we find:'
	call c_sphere(mx,my,mz,m,theta,phi)
 30     format(2x,'m = ',f6.2,' myb',5x,'theta = ',f7.2,' ',5x,'phi = ',f7.2,' ')
 40     format(2x,'    ',6x  ,'    ',5x,'        ',7x  ,'o',5x,'      ',7x, 'o')
	write(out,40)
  	write(out,30) m,theta,phi
	write(out,*) ' The x-,y- and z-axes are defined by the Euler-angles alpha,beta and gamma.'
        write(out,*) ' Type @ d a  to show what values they have.' 
      return
      end
c--------------------------------------
c calculates the spherical coordinates
c--------------------------------------
      subroutine c_sphere(x,y,z,r,theta,phi)
	implicit none
	real*8 x,y,z,r,theta,phi,pi,ccos,ssin,angle
	pi = 4.0d0*atan(1.0d0)
	r = sqrt(x*x+y*y+z*z)
	call ifnull(r)
	phi   = 0.0
	theta = 0.0
	if(r.eq.0.0d0) return
	theta = acos(z/r)*180/pi
	phi   = sin(theta*pi/180)
	call ifnull(phi)
	if(phi.eq.0.0d0) return
	ccos = x/r/sin(theta*pi/180)
	ssin = y/r/sin(theta*pi/180)
	phi  = angle(ccos,ssin)*180/pi
      return
      end
c--------------
c determine phi 
c--------------
      real*8 function angle(cphi,sphi)
	implicit none
	real*8 cphi,sphi,pi
	pi = 4.0d0*atan(1.0d0)
	call ifnull(sphi)
	angle = (1-sign(1.0d0,cphi))*pi/2
	if(sphi.eq.0.0d0) return
	angle=atan(cphi/sphi)
        if(cphi.ge.0.0d0.and.sphi.ge.0.0d0) return
	angle=angle+pi/2
	if(cphi.lt.0.0d0.and.sphi.ge.0.0d0) return
	angle=angle+pi/2
	if(cphi.lt.0.0d0.and.sphi.lt.0.0d0) return
	angle=angle+pi/2
      return
      end
c----------------------------------------------------------------------
c calculates mx = - gj* sum( exp(-E(i)/T)/T * <i|jx|i>, i=1 , 2j+1)
c in units of myB
c----------------------------------------------------------------------
      real*8 function c_x_moment(gj,ew,ev,dimj,z,t)
	implicit none
	real*8 gj,ew(17),dimj,z,t
	complex*16 ev(17,17),matjx,c_moment
	real*8 exp_
	integer i,dim
	dim = dimj
	c_moment = cmplx(0.0d0,0.0d0)  
	do i=1,dimj
           c_moment = c_moment + exp_(-ew(i)/t)/z * matjx(ev,i,i,dim)
	end do
	c_x_moment = -gj*c_moment
      return
      end     
c----------------------------------------------------------------------
c calculates my = - gj* sum( exp(-E(i)/T)/T * <i|jy|i>, i=1 , 2j+1)
c in units of myb
c----------------------------------------------------------------------
      real*8 function c_y_moment(gj,ew,ev,dimj,z,t)
	implicit none
	real*8 gj,ew(17),dimj,z,t
	complex*16 ev(17,17),matjy,c_moment
	real*8 exp_
	integer i,dim
	dim = dimj
	c_moment = cmplx(0.0d0,0.0d0)  
	do i=1,dimj
           c_moment = c_moment + exp_(-ew(i)/t)/z * matjy(ev,i,i,dim)
	end do
	c_y_moment = -gj*c_moment
      return
      end     
c----------------------------------------------------------------------
c calculates mz = - gj* sum( exp(-E(i)/T)/T * <i|jz|i>, i=1 , 2j+1)
c in units of myB
c----------------------------------------------------------------------
      real*8 function c_z_moment(gj,ew,ev,dimj,z,t)
	implicit none
	real*8 gj,ew(17),dimj,z,t
	complex*16 ev(17,17),matjz,c_moment
	real*8 exp_
	integer i,dim
	dim = dimj
	c_moment = cmplx(0.0d0,0.0d0)  
	do i=1,dimj
           c_moment = c_moment + exp_(-ew(i)/t)/z * matjz(ev,i,i,dim)
	end do
	c_z_moment = -gj*c_moment
      return
      end     
c-----------------------
c prints 'occupation  factor'
c-----------------------
      subroutine p_occupation_factor(out,occupation_factor,temp)
	implicit none
	integer out
	real*8 occupation_factor,temp
  1     format(1x,'Z(',f5.1,'K ) = ',f6.2)
	call section(out,' THE CALCULATED OCCUPATION FACTOR ')
	write(out,1) temp,occupation_factor 
	write(out,*)
      return
      end
c----------------------------------
c prints out the transition energies
c---------------------------------- 
      subroutine p_energy_values(out,energy,dimj,fmevkelvin)  
	implicit none
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:DEGENERATION.INC'

	integer out,i
 	real*8 energy(17),dimj,fmevkelvin

	integer degeneration(17*17)  ! stores the degeneration of a level
 	integer n_energies        ! no. of degenerated energy levels
	real*8  e_energies(17*17) ! energy values of the degenerated energy levels
	real*8  i_energies(17,17) ! intensities of the degenerated energy levels      

	call section(out,' THE CALCULATED ENERGY LEVELS ')

        if(setdegoff.eq.1) then
   1      format(x,x,i2,a2,f8.2,a4)
   2     format(x,x,i2,a2,f16.2,a4)
          do i=1,dimj
             if( energy(i)/fmevkelvin .lt.1.0d5) then
                write(out,1) i,' :',energy(i)/fmevkelvin,' meV'
             else
                write(out,2) i,' :',energy(i)/fmevkelvin,' meV'
             end if
          end do
        else
          call deg_on(energy,i_energies,dimj,degeneration,e_energies,i_energies,n_energies,de,di)   
   3      format(x,x,i2,a2,f8.2,a20,i3,a5)
   4      format(x,x,i2,a2,f16.2,a20,i3,a5)
          do i=1,n_energies
             if( e_energies(i)/fmevkelvin .lt.1.0d5) then
                write(out,3) i,' :',e_energies(i)/fmevkelvin, ' meV  Degeneration :',degeneration(i),'-fold'
             else
                write(out,4) i,' :',e_energies(i)/fmevkelvin, ' meV  Degeneration :',degeneration(i),'-fold'
             end if
          end do
        end if
      return
      end
c-------------------------
c prints the wavefunctions
c-------------------------
      subroutine p_wave_functions(out,energy,ev,dimj)
	implicit none
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:DEGENERATION.INC'

	integer out,s,t,u
        real*8 energy(17)
	complex*16 ev(17,17),i
	real*8 dimj,sj,tj,j,re,im

	integer degeneration(17*17)  ! stores the degeneration of a level
 	integer n_energies        ! no. of degenerated energy levels
	real*8  e_energies(17*17) ! energy values of the degenerated energy levels
	real*8  i_energies(17,17) ! intensities of the degenerated energy levels      

        integer no,ro

c use the columns of ev
	call section(out,' THE CALCULATED WAVEFUNCTIONS ')
  1     format(2x,'!',i2,'> = '   ,(1x,f7.3,' + ',f7.3,'*i'),1x,'!',f4.1,'>')
  2     format(2x,1x ,2x,2x,'+',1x,(1x,f7.3,' + ',f7.3,'*i'),1x,'!',f4.1,'>')
  3     format(2x,'!',i2,i3,'> = '   ,(1x,f7.3,' + ',f7.3,'*i'),1x,'!',f4.1,'>')
  4     format(2x,1x ,2x,2x,3x,'+',1x,(1x,f7.3,' + ',f7.3,'*i'),1x,'!',f4.1,'>')

        if(setdegoff.ne.1) then
           call deg_on(energy,i_energies,dimj,degeneration,e_energies,i_energies,n_energies,de,di)  
        end if

        i = cmplx(0.0d0 , 1.0d0)
	j = 0.5d0*(dimj-1)
	do s=1,dimj
           u = 0
 10        u  = u + 1
           if(u.gt.dimj) goto 999
	   tj = dble(u)-j-1
           re = 0.5d0*( ev(u,s) + conjg(ev(u,s)) )
           im = 0.5d0*(-ev(u,s) + conjg(ev(u,s)) )*i
           re = int(1000.0d0*re+0.5d0)/1000.0d0
           im = int(1000.0d0*im+0.5d0)/1000.0d0
           ev(u,s) = cmplx(re,im)
           call cifnull(ev(u,s))
	   if(ev(u,s).eq.cmplx(0.0d0,0.0d0)) goto 10 

           if(setdegoff.eq.1) then 
                write(out,1) s,ev(u,s),tj
           else               
                write(out,3) no(s,degeneration,n_energies),
     *                       ro(s,degeneration,n_energies),ev(u,s),tj
           end if

	   do t=u+1,dimj
	      tj = dble(t)-j-1
              call cifnull(ev(t,s))
	      if(ev(t,s).ne.cmplx(0.0d0,0.0d0)) then 
                if(setdegoff.eq.1) then
                     write(out,2)  ev(t,s),tj   
                else   
                     write(out,4)  ev(t,s),tj   
                end if
              endif
           end do
	   write(out,*)
        end do	
 999   return
      end
c-------------------------------------------------------
c tests if the matrix of wavefunctions is unitary
c-------------------------------------------------------
      subroutine t_orthonormalisation(out,ev,dim)
	implicit none
	integer out,i,k,s
	complex*16 ev(17,17),tik,ci
	real*8 dim,re,im
	call section(out,' TESTING THE UNITARY OF THE MATRIX OF WAVEFUNCTIONS U. ')
	write(out,*)
	write(out,*)' Only nonvanishing matrix elements (U*U+)ik are printed'
	ci = cmplx(0.0d0,1.0d0)
	do i=1,dim
	   do k=1,dim
	      tik = cmplx(0.0d0,0.0d0)
	      do s=1,dim
 	         tik =tik+ev(i,s)*conjg(ev(k,s))
              end do
	      re = 0.5d0*( tik+conjg(tik))
	      im = 0.5d0*(-tik+conjg(tik))*ci
c	      if(abs(re).lt.1.0d-6) re = 0.0d0
c	      if(abs(im).lt.1.0d-6) im = 0.0d0
	      tik = cmplx(re,im)
	      if(tik.ne.cmplx(0.0d0,0.0d0)) write(out,1) i,k,tik
           end do
        end do
  1     format(1x,' (U*U+)(',i2,',',i2,') = ',f8.6,' + i * ',f8.6)
	write(out,*)
      return
      end
c------------------------------------------
c prints out the transition matrix elements
c------------------------------------------
      subroutine p_matrix_elements(out,dimj,energy,jx2mat,jy2mat,jz2mat,jt2mat)
	implicit none
	integer out,i,k
	real*8 dimj,energy(17),jx2mat(17,17),jy2mat(17,17),jz2mat(17,17),jt2mat(17,17)
	real*8 sum
	call section(out, ' F O R  S I N G L E  C R Y S T A L  M E A S U R E M E N T S ')
	write(out,*) ' The transition matrix elements |<i|jx|k>|**2.'
        call table(out,dimj,energy,jx2mat,sum)
	write(out,*) ' The transition matrix elements |<i|jy|k>|**2.'
        call table(out,dimj,energy,jy2mat,sum)
	write(out,*) ' The transition matrix elements |<i|jz|k>|**2.'
        call table(out,dimj,energy,jz2mat,sum)

        call p_mat_powder(out,dimj,energy,jt2mat)

      return 
      end        
c-------------------------------------------------------
c prints out the transition matrix elements for a powder
c-------------------------------------------------------
      subroutine p_mat_powder(out,dimj,energy,jt2mat)
	implicit none
	integer out,i,k
	real*8 dimj,energy(17),jt2mat(17,17)
	real*8 sum
	call section(out, ' F O R  P O W D E R  M E A S U R E M E N T S ')
	write(out,*) ' The transition matrix elements |<i|jt|k>|**2.'
        call table(out,dimj,energy,jt2mat,sum)
        write(out,*) ' The sum of each row  i  is :     2/3*j*(j+1)*n_i '
        write(out,*) ' n_i  is  the  degeneration of the energy level i '
        write(out,*) ' which  belongs to  the  row  i. See  above : THE '
        write(out,*) ' CALCULATED ENERGY LEVELS  for the convention  of '
        write(out,*) ' labeling of the energy levels.' 
	write(out,*)
      return 
      end        
c---------------------------------------------------------------
c prints a table for max 10 rows
c---------------------------------------------------------------
      subroutine print_n(out,dimj,mat,sum0)
        implicit none
        integer out
        real*8 dimj,mat(17,17),sum0
        real*8 sum
        integer i,k
   3    format(1x,i2,18(1x,f6.2))   
   4    format(x,x,x,18(2x,i3,2x))

!---------------------------------
!       plot them out 
!---------------------------------
        write(out,4) (i,i=1,dimj)
        sum0 = 0.0d0
	do i=1,dimj
	   sum=0.0d0
	   do k=1,dimj
              sum=sum+mat(i,k)
           end do  
	   write(out,3) i,(mat(i,k),k=1,dimj),sum
           sum0 = sum0 + sum
        end do
        write(out,*)
      return
      end
c---------------------------------------------------------------
c prints a table for more than 10 rows
c---------------------------------------------------------------
      subroutine print_w(out,dimj,mat,sum0)
        implicit none
        integer out
        real*8 dimj,mat(17,17),sum0,mmat(17)        
        real*8 sum
        integer i,k
   2    format(x,x,i3,18(x,i3))
  21    format(x,x,3x,18(x,f3.2))
   1    format(x,x,3x,18(x,i2,x))

!---------------------------------
!       plot them out 
!---------------------------------
        write(out,1) (i,i=1,dimj)
        sum0=0.0d0
        do i=1,dimj
           sum=0.0d0
           do k=1,dimj
              sum     = sum+mat(i,k)
              mmat(k) = int(100*mat(i,k)+0.5d0)/100.0d0
           end do 
           sum = int(100*sum+0.5d0)/100.0d0
           write(out,2)  i,(int(mmat(k)),k=1,dimj),int(sum)
           write(out,21)   (mmat(k)-int(mmat(k)),k=1,dimj),sum-int(sum)
           sum0 = sum0 + sum
        end do
        write(out,*)
      return
      end
c---------------------------------------------------------------
c prints a table independent of rows
c---------------------------------------------------------------
      subroutine table(out,dimj,energy,mat,sum0)   
        implicit none
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:DEGENERATION.INC'

        integer out
        real*8  dimj,energy(17),mat(17,17),sum0        
        real*8  n_levels

	integer degeneration(17*17)  ! stores the degeneration of a level
 	integer n_energies        ! no. of degenerated energy levels
	real*8  e_energies(17*17) ! energy values of the degenerated energy levels
	real*8  i_energies(17,17) ! intensities of the degenerated energy levels      


        if(setdegoff.eq.1) then
           n_levels = dimj
           if(int(n_levels).le.9) then
              call print_n(out,n_levels,mat,sum0)
           else        
              call print_w(out,n_levels,mat,sum0)  
           end if
        else
           call deg_on(energy,mat,dimj,degeneration,e_energies,i_energies,n_energies,de,di)  
           n_levels = n_energies
           if(int(n_levels).le.9) then
              call print_n(out,n_levels,i_energies,sum0)
           else        
              call print_w(out,n_levels,i_energies,sum0)
           end if
        end if

      return
      end
c---------------------------------------------------------------
c prints the transition intensities for a given temperature temp
c---------------------------------------------------------------
      subroutine p_intensities(out,dimj,energy,intensity,temp)  
	implicit none
	integer out,i,k
	real*8 dimj,energy(17),intensity(17,17),temp
	real*8 sum,sum0

   6    format(x,x,a14,f7.2,a8)
	call section(out, ' F O R  P O W D E R  M E A S U R E M E N T S ')
	write(out,*) ' The transition intensities I(i->k) in barn.'
        write(out,6) 'temperature : ',temp,' Kelvin.'

        call table(out,dimj,energy,intensity,sum0)   

   5    format(x,x,a18,f6.2,a6)
	write(out,*) ' With const := 4*pi*(1/2*r0*gj)**2 we may write :'
	write(out,*) ' I(i->k) = const*exp(-E(i)/T)/Z(T)*|<i|jt|k>|**2.' 
	write(out,*) ' For each row the sum of the columns is equal to:'
	write(out,*) ' 2/3*j(j+1)*const*exp(-E(i)/T)/Z(T).'
        write(out,*)
	write(out,*) ' Summing additionally over all rows, we get the'
	write(out,*) ' total scattered magnetic intensity to :' 
	write(out,5) '2/3*const*j(j+1): ',sum0, ' barn.'
	write(out,*)

      return 
      end        
c---------------------------
c set some natural constants
c---------------------------
      subroutine s_natural_constants(pi,kb,hh,hq,ee,cc,me,na)  
	implicit none	
	real*8 pi,kb,hh,hq,me,ee,cc,na 
	pi = 4.0d0*atan(1.0d0)
	kb = 1.38062d0     ! x 10**(-23) J/K,   Boltzmann constant k_B
	hh = 6.626075540d0 ! x 10**(-34) J*sec, Planks constant h
	hq = hh/2/pi
	ee = 1.6021773349d0  ! x 10**(-19) Coulomb, electric charge
	cc = 2.99792458d0  ! x 10**(8) m/sec, speed of light
	me = 9.109389754d0 ! x 10**(-31) kg, electron mass
	na = 6.022045d0    ! x 10**(23) particles/mol, Avogardo constant  
      return
      end
c-----------------------
c calculates fmevkelvin 
c-----------------------
      real*8 function c_fmevkelvin()
	implicit none
	real*8 pi,kb,hh,hq,me,ee,cc,na 
        call s_natural_constants(pi,kb,hh,hq,ee,cc,me,na)  
c       within the above choose, the factor which connects
c       1meV and 1 Kelvin is fixed and given by:
c       10*ee/kb =: fmevkelvin = 11.6047...
c       this means 1 meV is nearly 11.6 K
	c_fmevkelvin = 10*ee/kb                          
      return
      end    
c--------------------------
c calculates Bohr radius a0 
c--------------------------
      real*8 function c_a0()
	implicit none
	real*8 pi,kb,hh,hq,me,ee,cc,na 
        call s_natural_constants(pi,kb,hh,hq,ee,cc,me,na)  
c       a0 is the Bohr radius in Angstroems
	c_a0 = 100*hq**2/me/ee**2/cc**2  
      return
      end    
c-----------------------------
c calculates Bohr magneton myb
c-----------------------------
      real*8 function c_myb()
	implicit none
	real*8 pi,kb,hh,hq,me,ee,cc,na 
	real*8 c_fmevkelvin
        call s_natural_constants(pi,kb,hh,hq,ee,cc,me,na)    
c       magneton of Bohr in kelvin per tesla
	c_myb = hq/me/2*c_fmevkelvin()               
      return
      end
c---------------------------------------------
c calculates the neutron scattering radius r0
c---------------------------------------------
      real*8 function c_r0()
	implicit none
	real*8 pi,kb,hh,hq,me,ee,cc,na 
        call s_natural_constants(pi,kb,hh,hq,ee,cc,me,na)    
c       magnetic neutron scattering radius
	c_r0 = -1.91*ee**2/me ! * 10**(-12) cm    
      return
      end
c--------------------------------------------------------
c calculate the crystal field splitting for a special R3+
c Moment is canted in the crystal system
c--------------------------------------------------------
      subroutine c_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,type)  
        implicit none

        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.PAR'
        include 'CF_SOURCES:CF_FABI.INC'
        include 'CF_SOURCES:XR.INC'
        include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:OVERALL.INC'
        include 'CF_SOURCES:DEGENERATION.INC'

	integer degeneration(17*17)  ! stores the degeneration of a level
 	integer n_energies        ! no. of degenerated energy levels
	real*8  e_energies(17*17) ! energy values of the degenerated energy levels
	real*8  i_energies(17,17) ! intensities of the degenerated energy levels      
        real*8  c_fmevkelvin

        integer  n0, ncftot, no_0_par, no_cf_par, same, i,j

        character*71 texte1, texte2, texte3, texte4
        data texte1 /' ERROR: Can not fix the 00. energy level to an energy of 000000.00 meV '/
        data texte2 /'        because the biggest levelnumber allowed is 00.                 '/

        data texte3 /' WARNING: FOCUS changed the CEF- and magnetic field parameter. Please  '/
        data texte4 /'          display the new parameter with the command D P.              '/

        integer k,dim
	real*8  pi, c_r0, r0, c_occupation_factor
                                                         
        if( stdhamoff.eq.1) then
            call dc_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,type)    
        else
            call sc_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,type)   
        end if

        if( setoverall.eq.1 ) then
c           ------------------------------------------------- 
c           find out how many degenerated energy levels exists
c           store the intensities of the degenerated levels
c           ------------------------------------------------- 
	    call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)           
            call deg_on(energy,intensity,dimj,degeneration,e_energies,i_energies,n_energies,de,di)   
            dim = dimj
c           ------------------------------------------------- 
c           release new parameter if levelnumber is not found
c           for the case of no monte carlo calculation
c           ------------------------------------------------- 
            if( levelnumber.gt.n_energies .and. setmonte.eq.0) then  
                write(texte1(25:26),fmt='(i2)')   levelnumber
                write(texte1(58:66),fmt='(f9.2)') levelenergy
                write(texte2(52:53),fmt='(i2)')   n_energies
                call remark( texte1)
                call remark( texte2)
!               ---------------------
!               release new parameter
!               ---------------------
                n0     = no_0_par()
                ncftot = no_cf_par(symmetry,nre)                                  
                do j=n0-5,n0+ncftot
                   if( setxr.eq.1 ) then
                        if( j.eq.indvold(1) .or. same(j,indvold(1)).eq.0 ) then
                             p(j)    = p(j)   *ratioenergy
                             pmin(j) = pmin(j)*ratioenergy
                             pmax(j) = pmax(j)*ratioenergy
                        end if
                   else if( setwx.eq.1 ) then
                        if( (j.ge.n0-5.and.j.le.n0) .or. j.eq.n0+1 ) then
                             p(j)    = p(j)   *ratioenergy
                             pmin(j) = pmin(j)*ratioenergy
                             pmax(j) = pmax(j)*ratioenergy
                        end if
                   else 
                       p(j)    = p(j)   *ratioenergy
                       pmin(j) = pmin(j)*ratioenergy
                       pmax(j) = pmax(j)*ratioenergy
                   end if
                end do
                setoverall  = 0
                setpardis   = 1
                ratioenergy = 1.0d0
                return
            end if

c           ----------------------------------------------------------------
c           calculate the ratio of energies
c           ----------------------------------------------------------------
            if(setmonte.eq.0 .or.(levelnumber.le.n_energies.and.setmonte.eq.1)) then
               ratioenergy = levelenergy*c_fmevkelvin()/e_energies(levelnumber) ! ground state is e_energies(1)
            else 
               ratioenergy = 1.0d0
            endif

            if(ratioenergy.ne.1.0d0) then
!              ------------------------------------------
!              calculate the new energies
!              ------------------------------------------
               dim = dimj
               do k=1,dim
                  energy(k) = energy(k) * ratioenergy
               end do
!              ------------------------------------------
!              calculates the new occupation_factor
!              ------------------------------------------
  	       occupation_factor = c_occupation_factor(energy,dimj,temp)
!              ------------------------------------------
!              calculates the new transition intensities for a powdered sample
!              ------------------------------------------
               pi = 4.0d0*atan(1.0d0)
               r0 = c_r0()
	       call intcalc(pi,r0,gj,occupation_factor,jt2mat,energy,intensity,dim,temp)      
               if( setpardis.eq.0 .and. setnowarnings.ne.1) then
                   call remark( texte3)
                   call remark( texte4)
               end if
            end if

           
        end if

      return
      end
c--------------------------------------------------------
c calculate the crystal field splitting for a special R3+
c Moment is lying on +c direction 
c calculation with steven operators
c--------------------------------------------------------
      subroutine sc_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,type)  
	implicit none
        include 'CF_SOURCES:CF_FABI.INC'
        include 'CF_SOURCES:CANTED.INC'

c       INTERNAL ARRAYS
	integer dim,m,n,k,q,qs
	real*8 nj,mj,j
	real*8 facmol,facext
	real*8 delta,jp,jm,jop

	complex*16 h(17,17)
        real*8       full_okq
        complex*16   dkq_star(6,-6:6), b_kq
        complex*16  rdkq_star(6,-6:6)
        complex*16   bex(1,-1:1)
        complex*16  rbex(1,-1:1)
        real*8 epsilon,omega
        real*8 a,b,c
        real*8 aj
        complex*16 ddrot
	complex*16 rbextp,rbextm,rbextz
	complex*16 bmolp,bmolm,bmolz
	complex*16 bextp,bextm,bextz,i
	complex*16 vec(17) ! for internal usage
	real*8 c_occupation_factor     
c       ------------------
	real*8 pi,c_fmevkelvin,c_myb,c_r0,r0
	real*8 c_operator_norm  ! is a function which is called
	real*8 norm
	real*8 f_bmol,f_bext
	complex*16 okq(17,17)
        real*8 ssbkq(0:6,0:6)
c       --------------------------------------------------------   
c       initialize some rare earth constants
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)           
                                                             
c       magneton of Bohr in kelvin per tesla         
	myb    = c_myb()
                                         
	if(type.eq.'Vkq' .or. type.eq.'vkq' .or. type.eq.'VKQ') then
c          the Vkq are in kelvin and are the parameters which uses the normalized 
c          operators Okq/||Okq|| and j+/||j+|| and j-/||j-|| and jz/||jz||
c          calculates the Bkq,B+,B- and Bz from the given Vkq,V+,V- and Vz 
	   f_bmol = 2*(gj-1)*myb
	   f_bext = gj*myb
	   norm = c_operator_norm(dimj,0,0,sbkq,okq)    ! 00 = j+ note:||j+||=||j-||
	   bmol(1) = 2*bmol(1)/norm/f_bmol !
	   bext(1) = 2*bext(1)/norm/f_bext !  V- / ||j+|| = 1/2*f*B-
	   bmol(2) = 2*bmol(2)/norm/f_bmol !   
	   bext(2) = 2*bext(2)/norm/f_bext !  V+ / ||j-|| = 1/2*f*B+
	   norm = c_operator_norm(dimj,0,1,sbkq,okq)    ! 01 = jz
	   bmol(3) = bmol(3)/norm/f_bmol      !
	   bext(3) = bext(3)/norm/f_bext      ! Vz / ||jz|| = f*Bz
	   do k=2,6,2
	      do m=0,k
	        call cifnull(bkq(k,m))
	        if(bkq(k,m).ne.0.0d0) then
		    ssbkq(k,m) = 1 ! note ||re(Okq)|| = ||im(Okq)|| for q<>0 
                    if(symmetry.eq.8) then
                       norm = 1.0d0
                       if(k.eq.4) norm = c_operator_norm(dimj,4,0,ssbkq,okq)
                       if(k.eq.6) norm = c_operator_norm(dimj,6,0,ssbkq,okq)
                    else
		        norm = c_operator_norm(dimj,k,m,ssbkq,okq)
                    end if
                    bkq(k,m) = bkq(k,m)/norm 
	        endif
	      end do
	   end do
        endif
c       ------------------
c       gives nearly 11.604... K
	pi = 4.0d0*atan(1.0d0)
	fmevkelvin = c_fmevkelvin()                        
c       ------------------
c       magnetic neutron scattering radius
	r0 = c_r0()
c       ------------------------------------------------------------
c       transform the Bkq with
c       H = sum_k=0 Bk0 Ok0 + sum_k>0_q>0  ReBkq ReOkq + ImBkq ImOkq
c       to a representation with
c                     *
c       H = sum_kq dkq  Okq
c
c       one finds:   dk0 = Bk0  and for q<>0: dkq = Bkq/2
c       ------------------------------------------------------------
	i = cmplx(0.0d0,1.0d0)
        do k=2,6,2
           do q=0,k
              b_kq = bkq(k,q)
              if(sbkq(k,q).eq.-1) b_kq = i*b_kq
              dkq_star(k,q) = b_kq
              if(q.ne.0) dkq_star(k,q) = dkq_star(k,q)/2
              dkq_star(k,-q) = conjg( dkq_star(k,q) )
           end do
        end do
c       -------------------------------------------------------------------
c       the parameters are   conjg(D_kq)
c       -------------------------------------------------------------------
        do k=2,6,2
           do q=-k,k
              dkq_star(k,q) = conjg( dkq_star(k,q) )
           end do
        end do
c       -------------------------------------------------------------------
        a = alpha_euler
        b = beta_euler
        c = gamma_euler
        do k=2,6,2
           do q=-k,k
              rdkq_star(k,q) = cmplx(0.0d0,0.0d0)
              do qs=-k,k
                rdkq_star(k,q)=rdkq_star(k,q) 
     *          + dkq_star(k,qs)
     *            * epsilon(k,qs)/epsilon(k,q) * omega(k,qs)/omega(k,q)
     *            * ddrot(k,q,qs,a,b,c)
              end do
           end do   
        end do
c       -------------------------------------------------------------------
c       rotate the external magnetic field
c       -------------------------------------------------------------------
c       aj  = sqrt(3/j/(j+1)/(2*j+1))
c
c       the following hamiltonian will be rotated
c
c        cry        -----    cry    cry
c       h     = f * >     Bex     V       with  <j|| V1 ||j> = sqrt(3)
c        mag        -----    1q    1q
c                     q
c

	bex(1, 1) = - (bext(1)-i*bext(2))/sqrt(2.0d0) ! /aj
	bex(1,-1) =   (bext(1)+i*bext(2))/sqrt(2.0d0) ! /aj
	bex(1, 0) =    bext(3)                        ! /aj

c       calculates Bex(1,q) for a canted moment
c
c       moment    cry       cry   cry     +
c      h       = R   (abc) h     R   (abc)
c       mag                 mag
c
        do q=-1,1
           rbex(1,q) = cmplx(0.0d0,0.0d0)
           do qs=-1,1
             rbex(1,q) = rbex(1,q) + bex(1,qs)*ddrot(1,q,qs,a,b,c)
           end do
        end do   

        rbextp =  rbex(1,-1)*sqrt(2.0d0)  !*aj
        rbextm = -rbex(1, 1)*sqrt(2.0d0)  !*aj
        rbextz =  rbex(1, 0)              !*aj
c       -------------------------------------------------------------------
c       magneton of Bohr in kelvin per tesla         
	myb    = c_myb()
	facmol = 2*(gj-1)*myb
	facext = gj*myb
	i = cmplx(0.0d0,1.0d0)
	bmolp = bmol(1)+i*bmol(2)
	bmolm = bmol(1)-i*bmol(2)
	bmolz = bmol(3)
	bextp = rbextp
	bextm = rbextm
	bextz = rbextz        
c       -----------------------------------------
c       set the hamiltonian matrix h to complex 0
c       -----------------------------------------
	call mx_clear(h)
c       -------------------------------------------------------------------
c       calculate the crystal field hamiltonian
c       define only the lower triangle of h(m,n)
c       -------------------------------------------------------------------
	dim = dimj
	j   = 0.5d0*(dimj-1) ! total momentum J of R3+    
	do m=1,dim
	   mj = dble(m) - j - 1
	   do n=1,m
	        nj = dble(n) - j - 1        
                h(m,n) = cmplx(0.0d0,0.0d0)
                do k=2,6,2
                   do q=-k,k
                      h(m,n) = h(m,n) + rdkq_star(k,q)*full_okq(k,q,mj,nj,j) 
                   end do   
                end do
           end do
        end do
c       -------------------------------------------------------------------
c       define only the lower triangle of h(m,n)
c       -------------------------------------------------------------------
	dim = dimj
	j   = 0.5d0*(dimj-1) ! total momentum J of R3+    
	do 10 m=1,dim
	   mj = dble(m) - j - 1
	   do 20 n=1,m
	        nj = dble(n) - j - 1        
c
c add the molecular field
c f*J*B = f*( 1/2*(J+ * B-  +  J- * B+) + Jz*Bz )
c
                h(m,n) = h(m,n) 
     *                 + 0.5d0*facmol*bmolm*delta(mj,nj+1,j)*jp(nj,j)
     *                 + 0.5d0*facmol*bmolp*delta(mj,nj-1,j)*jm(nj,j)
     *                 +       facmol*bmolz*delta(mj,nj,j)*nj                   
c
c add an external magnetic field
c
     *                 + 0.5d0*facext*bextm*delta(mj,nj+1,j)*jp(nj,j)
     *                 + 0.5d0*facext*bextp*delta(mj,nj-1,j)*jm(nj,j)
     *                 +       facext*bextz*delta(mj,nj,j)*nj             
 20	   continue
 10	continue  
c
c diagonalisation of the hamiltonian h
c only the lower triangle is needed
c
	call diagonal(h,energy,wavefunction,dim)
c
c shift the lowest energy level to 0
c
	eshift = energy(1)
	do 50 n=1,dim
	   energy(n) = energy(n) - eshift
 50     continue
c
c write the energies out (in meV)
c
c	write(6,*) (energy(m)/11.6,m=1,dim)
c
c test if the wavefunction matrix wavefunction(dim,dim) is unitary
c
c       call t_orthonormalisation(6,wavefunction,dimj)
c
c calculates the transition matrixelements for a single crystal and
c a powdered sample
c
	call matcalc(wavefunction,dim,jx2mat,jy2mat,jz2mat,jt2mat)
c
c calculates the sum over all occupation_factor
c
	occupation_factor = c_occupation_factor(energy,dimj,temp)
c
c calculates the transition intensities for a powdered sample
c
	call intcalc(pi,r0,gj,occupation_factor,jt2mat,energy,intensity,dim,temp)      
c
c end of subroutine cfcalc
c
      return
      end
c--------------------------------------------------------
c matrix element <j|| C_k(J) ||j>
c--------------------------------------------------------
      real*8 function redmat(j,k)
        implicit none
        real*8 j,k
        real*8 fac  ! look to 3j6j_8.for 
        redmat = sqrt(fac(2*j+k+1)/fac(2*j-k))/2**(int(k))
      return
      end
c--------------------------------------------------------
c binom (n over k)
c--------------------------------------------------------
      real*8 function binom(n,k)
        implicit none
        integer n,k
        real*8 fac  ! look to 3j6j_8.for 
        binom = fac(dble(n))/fac(dble(k))/fac(dble(n-k))
      return
      end
c--------------------------------------------------------
c              (k)
c calculates  D    (a,b,c)          
c              ms m
c
c see Lindner A, 'Drehimpulse in der Quantenmechanik',
c ISBN 3-519-03061-6                  (j)
c Stuttgart: Teubner, 1984  page 86 (d    (beta)=...)
c for equation (1)                    m ms
c
c see Buckmaster phys. stat. sol. (a) 13 (1972) 9
c for equation (2)
c
c--------------------------------------------------------
      complex*16 function ddrot(j,m,ms,a,b,c)
        implicit none
        integer j,ms,m,n
        real*8 a,b,c
        real*8 binom
        real*8 fac  ! look to 3j6j_8.for 
        real*8 d
        real*8 delta

c       equation (1)
        d = delta(dble(ms),dble(m),dble(j))
        call ifnull(b)
        if(b.ne.0.0d0) then
        d = 0.0d0
          do n=max(0,-(m+ms)),min(j-m,j-ms)
             d = d + (-1)**(j-ms-n)*binom(j-ms,n)*binom(j+ms,j-m-n)
     *             * cos(0.5*b)**(2*n+m+ms) * sin(0.5*b)**(2*j-2*n-m-ms)
          end do
        d = d * sqrt(fac(dble(j+m))/fac(dble(j+ms))*      
     *               fac(dble(j-m))/fac(dble(j-ms)))
        end if
c       equation (2)
        ddrot = cmplx(cos(m*a),-sin(m*a)) * d * cmplx(cos(ms*c),-sin(ms*c))
      return
      end
c--------------------------------------------------------
c calculate the crystal field splitting for a special R3+
c Moment is canted in the crystal system
c calculation with 3j-symbols 
c--------------------------------------------------------
      subroutine dc_crystal_field(nre,symmetry,bkq,sbkq,bmol,bext,temp,type)   
	implicit none
        include 'CF_SOURCES:CF_FABI.INC'
        include 'CF_SOURCES:CANTED.INC'

c       INTERNAL ARRAYS
	integer dim,m,n,k,q,qs
	real*8 nj,mj,j,kk,qq
	real*8 facmol,facext
	real*8 delta,jp,jm
	complex*16 h(17,17)
        complex*16   dkq_star(6,-6:6) ! H=sum_kq Dkq* x Ckq(J)
        complex*16  rdkq_star(6,-6:6)
        complex*16   bex(1,-1:1)
        complex*16  rbex(1,-1:1)
        complex*16  b_kq 
        real*8 dreij,redmat
        real*8 epsilon,omega
        real*8 theta_k
        real*8 a,b,c
        real*8 aj                                        
        complex*16 ddrot
	complex*16 bmolp,bmolm,bmolz
	complex*16 bextp,bextm,bextz,i
	complex*16 rbextp,rbextm,rbextz
	complex*16 vec(17) ! for internal usage
	real*8 c_occupation_factor     
c       ------------------
	real*8 pi,c_myb,c_r0,r0,c_fmevkelvin
	real*8 c_operator_norm  ! is a function which is called
	real*8 norm
        real*8 f_bmol,f_bext
	complex*16 okq(17,17)
        real*8 ssbkq(0:6,0:6)

c       --------------------------------------------------------
c       calculate the Bkq from the Vkq
c       --------------------------------------------------------
c       initialize some rare earth constants
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)           
                                                             
c       magneton of Bohr in kelvin per tesla         
	myb    = c_myb()
                                         
	if(type.eq.'Vkq' .or. type.eq.'vkq' .or. type.eq.'VKQ') then
c          the Vkq are in kelvin and are the parameters which uses the normalized 
c          operators Okq/||Okq|| and j+/||j+|| and j-/||j-|| and jz/||jz||
c          calculates the Bkq,B+,B- and Bz from the given Vkq,V+,V- and Vz 
	   f_bmol = 2*(gj-1)*myb
	   f_bext = gj*myb
	   norm = c_operator_norm(dimj,0,0,sbkq,okq)    ! 00 = j+ note:||j+||=||j-||
	   bmol(1) = 2*bmol(1)/norm/f_bmol !
	   bext(1) = 2*bext(1)/norm/f_bext !  V- / ||j+|| = 1/2*f*B-
	   bmol(2) = 2*bmol(2)/norm/f_bmol !   
	   bext(2) = 2*bext(2)/norm/f_bext !  V+ / ||j-|| = 1/2*f*B+
	   norm = c_operator_norm(dimj,0,1,sbkq,okq)    ! 01 = jz
	   bmol(3) = bmol(3)/norm/f_bmol      !
	   bext(3) = bext(3)/norm/f_bext      ! Vz / ||jz|| = f*Bz
	   do k=2,6,2
	      do m=0,k
	        call cifnull(bkq(k,m))
	        if(bkq(k,m).ne.0.0d0) then
		    ssbkq(k,m) = 1 ! note ||re(Okq)|| = ||im(Okq)|| for q<>0 
                    if(symmetry.eq.8) then
                       norm = 1.0d0
                       if(k.eq.4) norm = c_operator_norm(dimj,4,0,ssbkq,okq)
                       if(k.eq.6) norm = c_operator_norm(dimj,6,0,ssbkq,okq)
                    else
		        norm = c_operator_norm(dimj,k,m,ssbkq,okq)
                    end if
                    bkq(k,m) = bkq(k,m)/norm 
	        endif
	      end do
	   end do
        endif
c       --------------------------------------------------------
c       calculate the Dkq from the bkq(0:6,0:6) 
c       --------------------------------------------------------
c            ----   * 
c       H =  >     D   *  C  (J)
c            ----   kq     kq
c             kq
c
c
c       in the following only the Dkq (and not the Dkq*) are 
c       calculated from the Bkq. Remember that
c       bkq in the following has the meaning of Bkq==Re(Bkq)
c       or Bkq==Im(Bkq). This is actually how Bkq operates
c       in the programm and this is why below (for q not
c       equal to zero only) Dkq must be divided by 2.
c
c
	i = cmplx(0.0d0,1.0d0)
        do k=2,6,2
           do q=0,k
              b_kq = bkq(k,q)
c             in the case of sbkq(k,q)=-1 the real part of Bkq(k,q)
c             is to be interpreted as a imaginary part
              if( sbkq(k,q) .eq. -1) b_kq = i * b_kq
              dkq_star(k,q) = b_kq/epsilon(k,q)/omega(k,q)
              if(q.ne.0) dkq_star(k,q) = dkq_star(k,q)/2
              dkq_star(k,-q) = (-1)**q * conjg( dkq_star(k,q) )
           end do
        end do
c       -------------------------------------------------------------------
c       the parameters are   conjg(D_kq)
c       -------------------------------------------------------------------
        do k=2,6,2
           do q=-k,k
              dkq_star(k,q) = conjg( dkq_star(k,q) )
           end do
        end do
c       -------------------------------------------------------------------
c        cry     ----     *  cry         cry
c       h     =  >     [ D  ]    [ C (J)]      : hamiltonian in fixed
c                ----     kq        kq         : crystal system
c                 kq
c
c
c        moment      cry         cry    cry     +
c       h        =  R   (abc) * h   *  R   (abc)
c
c
c        cry             cry   cry     +    ----        moment    cry       cry          cry
c       R   (abc)  C  (J)     R   (abc)  =  >     C  (J)      *      <J qs|R   (abc)|q J>
c                   kq                      ----   k qs
c                                            qs
c
c
c
c
c       calculates D_kq for a canted moment
c
c                             (k)
c       drot(k,q,qs,a,b,c) = D    (abc) 
c                             q qs
c
c
c               ----   *         (k)
c       RD   =  >     D      *  D    (abc)
c         kq    ----   k qs      q qs
c                qs
c       -------------------------------------------------------------------
        a = alpha_euler
        b = beta_euler
        c = gamma_euler
        do k=2,6,2
           do q=-k,k
              rdkq_star(k,q) = cmplx(0.0d0,0.0d0)
              do qs=-k,k
                rdkq_star(k,q)=rdkq_star(k,q) 
     *                        + dkq_star(k,qs)*ddrot(k,q,qs,a,b,c)
              end do
           end do   
        end do
c       -------------------------------------------------------------------
c       rotate the external magnetic field
c       -------------------------------------------------------------------
c       aj  = sqrt(3/j/(j+1)/(2*j+1))
c
c       the following hamiltonian will be rotated
c
c        cry        -----    cry    cry
c       h     = f * >     Bex     V       with  <j|| V1 ||j> = sqrt(3)
c        mag        -----    1q    1q
c                     q
c

	bex(1, 1) = - (bext(1)-i*bext(2))/sqrt(2.0d0) ! /aj
	bex(1,-1) =   (bext(1)+i*bext(2))/sqrt(2.0d0) ! /aj
	bex(1, 0) =    bext(3)                        ! /aj

c       calculates Bex(1,q) for a canted moment
c
c       moment    cry       cry   cry     +
c      h       = R   (abc) h     R   (abc)
c       mag                 mag
c
        do q=-1,1
           rbex(1,q) = cmplx(0.0d0,0.0d0)
           do qs=-1,1
             rbex(1,q) = rbex(1,q) + bex(1,qs)*ddrot(1,q,qs,a,b,c)
           end do
        end do   

        rbextp =  rbex(1,-1)*sqrt(2.0d0)  !*aj
        rbextm = -rbex(1, 1)*sqrt(2.0d0)  !*aj
        rbextz =  rbex(1, 0)              !*aj

c       -------------------------------------------------------------------
c       set the hamiltonian matrix h to complex 0
c       -------------------------------------------------------------------
	call mx_clear(h)
c       -------------------------------------------------------------------
c       calculate the crystal field hamiltonian
c       define only the lower triangle of h(m,n)
c       -------------------------------------------------------------------
	dim = dimj
	j   = 0.5d0*(dimj-1) ! total momentum J of R3+    
	do m=1,dim
	   mj = dble(m) - j - 1
	   do n=1,m
	        nj = dble(n) - j - 1        
                h(m,n) = cmplx(0.0d0,0.0d0)
                do k=2,6,2
                   do q=-k,k
                      kk = dble(k)
                      qq = dble(q)
                      qs = int(j-mj)
                      h(m,n) = h(m,n) + rdkq_star(k,q)*(-1)**qs*
     *                                  dreij(j,kk,j,-mj,qq,nj)*redmat(j,kk)
                   end do   
                end do
           end do
        end do
c       -------------------------------------------------------------------
c       gives nearly 11.604... K
	pi = 4.0d0*atan(1.0d0)
	fmevkelvin = c_fmevkelvin()                        
c       -------------------------------------------------------------------
c       magneton of Bohr in kelvin per tesla         
	myb    = c_myb()
	facmol = 2*(gj-1)*myb
	facext = gj*myb
	i = cmplx(0.0d0,1.0d0)
	bmolp = bmol(1)+i*bmol(2)
	bmolm = bmol(1)-i*bmol(2)
	bmolz = bmol(3)
	bextp = rbextp
	bextm = rbextm
	bextz = rbextz        
c       -------------------------------------------------------------------
c       define only the lower triangle of h(m,n)
c       -------------------------------------------------------------------
	do 10 m=1,dim
	   mj = dble(m) - j - 1
	   do 20 n=1,m
	        nj = dble(n) - j - 1        
c
c               add the molecular field
c               f*J*B = f*( 1/2*(J+ * B-  +  J- * B+) + Jz*Bz )
c
                h(m,n) = h(m,n) 
     *                 + 0.5d0*facmol*bmolm*delta(mj,nj+1,j)*jp(nj,j)
     *                 + 0.5d0*facmol*bmolp*delta(mj,nj-1,j)*jm(nj,j)
     *                 +       facmol*bmolz*delta(mj,nj,j)*nj                   
c
c         add an external magnetic field
c
     *                 + 0.5d0*facext*bextm*delta(mj,nj+1,j)*jp(nj,j)
     *                 + 0.5d0*facext*bextp*delta(mj,nj-1,j)*jm(nj,j)
     *                 +       facext*bextz*delta(mj,nj,j)*nj             
 20	   continue
 10	continue  
c       -------------------------------------------------------------------
c       diagonalisation of the hamiltonian h
c       only the lower triangle is needed
	call diagonal(h,energy,wavefunction,dim)
c       -------------------------------------------------------------------
c       shift the lowest energy level to 0
	eshift = energy(1)
        do n=1,dim
	   energy(n) = energy(n) - eshift
        end do  
c       -------------------------------------------------------------------
c       test if the wavefunction matrix wavefunction(dim,dim) is unitary
c       call t_orthonormalisation(6,wavefunction,dimj)
c       -------------------------------------------------------------------
c       calculates the transition matrixelements for a single crystal and
c       a powdered sample
	call matcalc(wavefunction,dim,jx2mat,jy2mat,jz2mat,jt2mat)
c       -------------------------------------------------------------------
c       calculates the sum over all occupation_factor
	occupation_factor = c_occupation_factor(energy,dimj,temp)
c       -------------------------------------------------------------------
c       calculates the transition intensities for a powdered sample
c       magnetic neutron scattering radius
	pi = 4.0d0*atan(1.0d0)
	r0 = c_r0()
	call intcalc(pi,r0,gj,occupation_factor,jt2mat,energy,intensity,dim,temp)      
c       -------------------------------------------------------------------
      return
      end
c-------------------------------------
c calculation of the occupation factor
c-------------------------------------
      real*8 function c_occupation_factor(energy,dimj,temp)
	implicit none
	real*8 energy(17),dimj,exp_,temp
	integer s
	c_occupation_factor = 0.0d0
	if(temp.eq.0.0d0) temp=1.0d0  
	do 10 s=1,dimj
	   c_occupation_factor=c_occupation_factor+exp_(-energy(s)/temp)
 10     continue
      return
      end
c---------------------------------------------------------------
c calculates all transition matrix elements for a single crystal
c and a polycrystalline sample (powder)
c---------------------------------------------------------------
      subroutine matcalc(ev,dim,jx2,jy2,jz2,jt2)
	implicit none
	integer dim,i,k
	real*8 jx2(17,17),jy2(17,17),jz2(17,17),jt2(17,17)
	real*8 matjx2,matjy2,matjz2
	complex*16 ev(17,17)
	do 10 i=1,dim
	   do 20 k=1,dim
	      jx2(i,k) = matjx2(ev,i,k,dim)
	      jy2(i,k) = matjy2(ev,i,k,dim)
	      jz2(i,k) = matjz2(ev,i,k,dim)
	      jt2(i,k) = 2.0d0/3*(jx2(i,k)+jy2(i,k)+jz2(i,k))
 20     continue
 10     continue
      return
      end           
c------------------------------------------
c calculates all transition intensities for
c a polycrystalline sample (powder)
c------------------------------------------
      subroutine intcalc(pi,r0,gj,z,jt2,e,int,dim,temp)      
	implicit none
	integer dim,i,k
	real*8 pi,r0,gj,z,jt2(17,17),e(17),int(17,17),temp
	real*8 exp_,coeff,const
	const = 4.0*pi*(0.5*r0*gj)**2
	if(temp.eq.0.0d0) temp=1.0d0  

	do 10 i=1,dim
	   coeff = exp_(-e(i)/temp)/z * const
	   do 20 k=1,dim
              int(i,k) = coeff*jt2(i,k)
 20        continue
 10     continue
      return
      end                   
c-------------------------
c calculates |<i|jT|k>|**2
c-------------------------
      real*8 function matjt2(ev,i,k,dim)
	implicit none
	integer i,k,dim
	complex*16 ev(17,17)
	real*8 matjx2,matjy2,matjz2
	matjt2 = 2.0d0/3*( matjx2(ev,i,k,dim)+matjy2(ev,i,k,dim)
     *                    +matjz2(ev,i,k,dim))
      return
      end
c-------------------------
c calculates |<i|jx|k>|**2
c-------------------------
      real*8 function matjx2(ev,i,k,dim)
	implicit none
	integer i,k,dim
	complex*16 ev(17,17),matjx
	matjx2 = conjg(matjx(ev,i,k,dim))*matjx(ev,i,k,dim)
      return
      end
c-------------------------
c calculates |<i|jy|k>|**2
c-------------------------
      real*8 function matjy2(ev,i,k,dim)
	implicit none
	integer i,k,dim
	complex*16 ev(17,17),matjy
	matjy2 = conjg(matjy(ev,i,k,dim))*matjy(ev,i,k,dim)
      return
      end
c-------------------------
c calculates |<i|jz|k>|**2
c-------------------------
      real*8 function matjz2(ev,i,k,dim)
	implicit none
	integer i,k,dim
	complex*16 ev(17,17),matjz
	matjz2 = conjg(matjz(ev,i,k,dim))*matjz(ev,i,k,dim)
      return
      end 
c--------------------
c calculates <i|jx|k>
c--------------------
      complex*16 function matjx(ev,i,k,dim)
	implicit none
	integer i,k,dim
	complex*16 ev(17,17),matjp,matjm
	matjx = 0.5d0*(matjm(ev,i,k,dim)+matjp(ev,i,k,dim))
      return
      end
c--------------------
c calculates <i|jy|k>
c--------------------
      complex*16 function matjy(ev,i,k,dim)
	implicit none
	integer i,k,dim
	complex*16 ev(17,17),matjp,matjm,ci
	ci = cmplx(0.0d0,1.0d0)
	matjy = 0.5d0*ci*(matjm(ev,i,k,dim)-matjp(ev,i,k,dim))
      return
      end
c-------------------------------
c  calculates <i|j+|k>
c-------------------------------
      complex*16 function matjp(ev,i,k,dim)
	implicit none
	integer i,k,dim,s
	complex*16 ev(17,17),v(17),null
	real*8 j,sj,jp
	null= cmplx(0.0d0,0.0d0)
	j   = 0.5d0*(dble(dim)-1)
	do 10 s=dim,2,-1
	   sj   = dble(s)-j-1
           v(s) = ev(s-1,k)*jp(sj-1,j)
 10     continue
	v(1) = null
	matjp= null
	do 20 s=1,dim
	   matjp = matjp + conjg(ev(s,i))*v(s)
 20     continue 
      return
      end
c--------------------
c calculates <i|j-|k>
c--------------------
      complex*16 function matjm(ev,i,k,dim)
	implicit none
	integer i,k,dim,s
	complex*16 ev(17,17),v(17),null
	real*8 sj,j,jm
	null = cmplx(0.0d0,0.0d0)
	j    = 0.5d0*(dble(dim)-1)
	do 10 s=1,dim-1
	   sj   = dble(s)-j-1
	   v(s) = ev(s+1,k)*jm(sj+1,j)
 10     continue
	v(dim) = null
	matjm  = null
	do 20 s=1,dim
	   matjm = matjm + conjg(ev(s,i))*v(s)
 20     continue
      return
      end
c--------------------
c calculates <i|jz|k>
c--------------------
      complex*16 function matjz(ev,i,k,dim)
	implicit none
	integer i,k,dim,s
	complex*16 ev(17,17),v(17),null
	real*8 sj,j
	null = cmplx(0.0d0,0.0d0)
	j    = 0.5d0*(dble(dim)-1)
	do 10 s=1,dim
	   sj   = dble(s)-j-1
	   v(s) = ev(s,k)*sj
 10     continue
	matjz = null
	do 20 s=1,dim
	   matjz = matjz + conjg(ev(s,i))*v(s)
 20     continue
      return
      end
c------------------------------------------------------------
c make sure that no underflow and overflow appears within exp 
c------------------------------------------------------------
      real*8 function exp_(z)
	implicit none
	real*8 z,exp,zmax
	zmax = 71.0d0
	exp_ = exp(zmax)
	if(z.lt.0.0d0) exp_ = 0.0d0
	if(abs(z).lt.zmax) exp_ = exp(z)
      return
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c diagonalisation of a complex hermitian matrix
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine diagonal(a,eps,w,n)
	implicit none
	integer n
	complex*16 a(17,17),w(17,17),vec(17)
	real*8 eps(17)
c       ----------------
c       for internal use only
	real*8 ar(17,17),ai(17,17)
	real*8 vr(17,17),vi(17,17)
	real*8 wk1(17),wk2(17),wk3(17)
	integer i,j,ifail
	complex*16 ci
c       ----------------
	external F02AXF
c       ----------------
	call clear(ar)
	call clear(ai)
	call clear(vr)
	call clear(vi)
c
c set ar and ai (only the lower triangle)
c
	ci = cmplx(0.0d0,1.0d0)
	do 10 i=1,n
	   do 20 j=1,i 
		call cifnull(a(i,j))
		ar(i,j) = 0.5d0*( a(i,j)+conjg(a(i,j)))
		ai(i,j) = 0.5d0*(-a(i,j)+conjg(a(i,j)))*ci
 20        continue
 10     continue
c
c calls nag routine F02AXF
c
	ifail = 1
	call F02AXF(ar,17,ai,17,n,eps,vr,17,vi,17,wk1,wk2,wk3,ifail)
	if(ifail.ne.0) then
	   write(6,*) 'An error is detected by the subroutine DIAGONAL'
	   write(6,*) 'called by the subroutine CFCALC'
	   if(ifail.eq.1) then
	     write(6,*) 'more than ',30*n,
     *       ' iterations are required to isolate all eigenvalues'  
	   endif
	   if(ifail.eq.2) then
	     write(6,*) 'the diagonal elements of the complex matrix'
	     write(6,*) 'are not all zero,i.e. the matrix is not hermitian' 
	   endif
	endif
c
c set the eigenvectors into the complex matrix w
c
	do 30 i=1,n
	   do 40 j=1,n
	      w(i,j) = cmplx(vr(i,j),vi(i,j))
	      call cifnull(w(i,j))	
 40        continue
	   call ifnull(eps(i))
 30     continue

c
c orthonormalisation of the eigenvectors
c to make sure that even in the case of a degeneration 
c of energy values, the eigenvectors are orthonormal
c
	call orthonormal(w,vec,17,n)
c
c test if some components of the wavefunctions are zero 
c
	do 50 i=1,n
	   do 60 j=1,n
	      call cifnull(w(i,j))	
 60        continue
 50     continue

      return
      end
c-------------------------------------------
c orthonormalisation of all eigenvectors
c which are in the rows of mx
c-------------------------------------------
      subroutine orthonormal(mx,vec,lmx,nmx)
	implicit none
	integer lmx,nmx
	complex*16 mx(lmx,lmx),vec(lmx)
	integer i
	call normalisize(mx,lmx,nmx,1)
	do 10 i=2,nmx
	   call scalar(mx,lmx,nmx,i,vec,lmx)
	   call normalisize(mx,lmx,nmx,i)
 10     continue
      return
      end
c----------------------------------
c vec(s) = <w(s),v(sp)> i=1 to sp-1
c----------------------------------
      subroutine scalar(mx,lmx,nmx,sp,vec,dimvec)
	implicit none
	integer lmx,nmx,sp,dimvec,i,j
	complex*16 mx(lmx,lmx),vec(dimvec)
	do 10 i=1,sp-1
	   vec(i) = cmplx(0.0d0,0.0d0)
	   do 20 j=1,nmx
		vec(i)=vec(i)+mx(j,i)*conjg(mx(j,sp))
 20        continue
 10     continue
	do 30 i=1,sp-1
	   do 40 j=1,nmx
		mx(j,sp) = mx(j,sp) - mx(j,i)*vec(i)
40         continue 
30     continue
      return
      end
c-----------------------------------------------------------
c normalisize the sp-th row, per definition never equal to 0
c-----------------------------------------------------------
      subroutine normalisize(mx,lmx,nmx,sp)
	implicit none
	integer lmx,nmx,sp,i
	complex*16 mx(lmx,lmx)
	real*8 norm,c
	norm = 0.0d0
	do 10 i=1,nmx
	  c = mx(i,sp) * conjg(mx(i,sp))
	  norm = norm + c
 10     continue
        norm = sqrt(norm)
	do 20 i=1,nmx
	  mx(i,sp) = mx(i,sp)/norm
 20     continue
      return
      end
c--------------------------------
c define the delta function
c--------------------------------
      real*8 function delta(mj,nj,j)
	implicit none
	real*8 mj,nj,j
	delta = 0.0d0
	if(mj.eq.nj.and.abs(mj).le.j.and.abs(nj).le.j) delta = 1.0d0	
      return
      end
c-------------------------------
c clear the complex matrix mx
c-------------------------------
      subroutine mx_clear(mx)
	implicit none
	complex*16 mx(17,17)
	complex*16 null
	integer i,j
	null = cmplx(0.0d0,0.0d0)
	do 10 i=1,17
	   do 20 j=1,17
		mx(i,j) = null
 20	   continue
 10     continue
      return
      end                
c------------------------
c clears a real matrix mx
c------------------------
      subroutine clear(mx)
	real*8 mx(17,17)
	integer i,j
	do 10 i=1,17
	   do 20 j=1,17
		mx(i,j) = 0.0d0
 20        continue
 10     continue
      return
      end
c------------------------------
c  f(1) = 1, f(-1) = 1/i = (-i)
c------------------------------
      complex*16 function f(s)
	implicit none
	real*8 s
	f = cmplx((s+1)/2,(s-1)/2)
      return
      end
c-------------------------------
c ff(k,q,nj,j) := fkq(nj,j)
c-------------------------------
      real*8 function ff(k,q,nj,j)
	implicit none 
 	integer k,q,qq
 	real*8 nj,j
 	real*8 f20,f21,f22,f40,f41,f42,f43,f44,f60,f61,f62,f63,f64,f65,f66
 
        qq = abs(q)
	if(k.eq.2.and.qq.eq.0) then
		ff = f20(nj,j)   
	else if(k.eq.2.and.qq.eq.1) then  
		ff = f21(nj,j)   
	else if(k.eq.2.and.qq.eq.2) then  
		ff = f22(nj,j)   
	else if(k.eq.4.and.qq.eq.0) then  
		ff = f40(nj,j)   
	else if(k.eq.4.and.qq.eq.1) then     
		ff = f41(nj,j)   
	else if(k.eq.4.and.qq.eq.2) then  
		ff = f42(nj,j)   
	else if(k.eq.4.and.qq.eq.3) then     
		ff = f43(nj,j)   
	else if(k.eq.4.and.qq.eq.4) then  
		ff = f44(nj,j)   
	else if(k.eq.6.and.qq.eq.0) then     
		ff = f60(nj,j)   
	else if(k.eq.6.and.qq.eq.1) then     
		ff = f61(nj,j)   
	else if(k.eq.6.and.qq.eq.2) then  
		ff = f62(nj,j)   
	else if(k.eq.6.and.qq.eq.3) then     
		ff = f63(nj,j)   
	else if(k.eq.6.and.qq.eq.4) then  
		ff = f64(nj,j)   
	else if(k.eq.6.and.qq.eq.5) then     
		ff = f65(nj,j)   
	else if(k.eq.6.and.qq.eq.6) then 
		ff = f66(nj,j)   
	endif
      return
      end
c-------------------------------------------------
c calculates the norm ||A|| of a given operator A
c-------------------------------------------------
c k=0 q=0 : j+
c k=0 q=1 : jz
c k=0 q=2 : j-
c k=0 q=4 : o4 = o40+ 5*o44
c k=0 q=6 : o6 = o60-21*o64
      real*8 function c_operator_norm(dimj,k,q,sbkq,okq)
	implicit none
	integer k,q
	real*8 dimj,sbkq(0:6,0:6)
	complex*16 okq(17,17)
	integer m,n

	call s_okq(dimj,k,q,sbkq,okq)

	c_operator_norm = 0.0d0
	do m=1,dimj
	   do n=1,dimj
	      c_operator_norm = c_operator_norm + okq(m,n)*conjg(okq(m,n))	
           end do
	end do
	c_operator_norm = c_operator_norm/dimj
	c_operator_norm = sqrt(c_operator_norm)
      return
      end
c-----------------------------
c prints the operator norm out
c ----------------------------
      subroutine p_operator_norm(out,nre)
	implicit none
	integer out,nre
	real*8 norm,c_operator_norm,rnorm,inorm
	real*8 sbkq(0:6,0:6)
	complex*16 okq(17,17)
 	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6   
	integer k,q                            

        call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)       
  1	format(1x,'||j+|| = ',f16.6) 
  2	format(1x,'||jz|| = ',f16.6)  
  3	format(1x,'||j-|| = ',f16.6)  
  4	format(1x,'||o40 +  5 * o44|| = ',f16.6)  
  5	format(1x,'||o60 - 21 * o64|| = ',f16.6) 
  6     format(1x,'||re(O',i1,',',i1,')|| = ',f16.6,10x,'||im(O',i1,',',i1,')|| = ',f16.6)
	call section(out,' THE OPERATORNORM ')
	write(out,*) 'rare earth number :',nre
	norm = c_operator_norm(dimj,0,0,sbkq,okq) 
	write(out,1) norm       
	norm = c_operator_norm(dimj,0,1,sbkq,okq) 
	write(out,2) norm                       
	norm = c_operator_norm(dimj,0,2,sbkq,okq) 
	write(out,3) norm                          
	norm = c_operator_norm(dimj,0,4,sbkq,okq) 
	write(out,4) norm                          
	norm = c_operator_norm(dimj,0,6,sbkq,okq)  
	write(out,5) norm 
	write(out,*)
	do k=2,6,2
           do q=0,k
		sbkq(k,q) = 1
             	rnorm = c_operator_norm(dimj,k,q,sbkq,okq)   
		sbkq(k,q) = -1
             	inorm = c_operator_norm(dimj,k,q,sbkq,okq) 
		write(out,6) k,q,rnorm,k,q,inorm
	   end do
	   write(out,*)
	end do                          
      return
      end
c--------------------------------------------------
c set the full stevens' operators Okq to their R3+ value
c--------------------------------------------------
      subroutine s_okq(dimj,k,q,sbkq,okq)
    	implicit none
	integer k,q
	real*8 dimj,sbkq(0:6,0:6)
	complex*16 okq(17,17),o40(17,17),o44(17,17),o60(17,17),o64(17,17)     
  	real*8  mj,nj,j
	real*8  ff
	real*8  delta,d0,d1,d2,d3,d4,d5,d6
	real*8  jp,jp2,jp3,jp4,jp5,jp6,jm,jm2,jm3,jm4,jm5,jm6
	integer m,n,dim
	complex*16 i,f
	i = cmplx(0.0d0,1.0d0) 
	dim = dimj
	j   = 0.5d0*(dimj-1)
	do 10 m=1,dim
	   mj = dble(m)-j-1
	   do 20 n=1,dim
		nj = dble(n)-j-1
                if(k.eq.0) then
		   if(q.eq.0) then        !q=0 chose j+ operator
			okq(m,n) = delta(mj,nj+1,j)*jp(nj,j)
		   else if(q.eq.1) then   !q=1 chose jz operator
			okq(m,n) = delta(mj,nj,j)*nj
		   else if(q.eq.2) then   !q=2 chose j- operator
			okq(m,n) = delta(mj,nj-1,j)*jm(nj,j)
		   else if(q.eq.4) then   !q=4 chose o4 = o40 + 5*o44
      		        d0 = delta(mj,nj,j)   
	       	        o40(m,n) = d0*ff(4,0,nj,j) 
		        d4=delta(mj,nj+4,j)*jp4(nj,j)+delta(mj,nj-4,j)*jm4(nj,j)   
		        o44(m,n) = 0.25d0*(ff(4,4,mj,j)+ff(4,4,nj,j))*d4         
                        okq(m,n) = o40(m,n) +5*o44(m,n)
		   else if(q.eq.6) then   !q=6 chose o6 = o60 -21*o64
		        d0 = delta(mj,nj,j)   
		        o60(m,n) = d0*ff(6,0,nj,j) 
		        d4=delta(mj,nj+4,j)*jp4(nj,j)+delta(mj,nj-4,j)*jm4(nj,j)   
      		        o64(m,n) = 0.25d0*(ff(6,4,mj,j)+ff(6,4,nj,j))*d4 
		        okq(m,n) = o60(m,n) -21*o64(m,n)         
		   endif 
	        else
		   if(q.eq.0) then
	      	     	d0 = delta(mj,nj,j)
			okq(m,n) = d0*ff(k,q,nj,j)*(1+sbkq(k,q))/2
		   else if(q.eq.1) then 
			d1=delta(mj,nj+1,j)*jp(nj,j)+sbkq(k,q)*delta(mj,nj-1,j)*jm(nj,j)   
			okq(m,n) = 0.25d0*(ff(k,q,mj,j)+ff(k,q,nj,j))*d1*f(sbkq(k,q))                       
                   else if(q.eq.2) then 
                	d2=delta(mj,nj+2,j)*jp2(nj,j)+sbkq(k,q)*delta(mj,nj-2,j)*jm2(nj,j)
			okq(m,n) = 0.25d0*(ff(k,q,mj,j)+ff(k,q,nj,j))*d2*f(sbkq(k,q))                          
                   else if(q.eq.3) then
                	d3=delta(mj,nj+3,j)*jp3(nj,j)+sbkq(k,q)*delta(mj,nj-3,j)*jm3(nj,j)   
			okq(m,n) = 0.25d0*(ff(k,q,mj,j)+ff(k,q,nj,j))*d3*f(sbkq(k,q))          
		   else if(q.eq.4) then
			d4=delta(mj,nj+4,j)*jp4(nj,j)+sbkq(k,q)*delta(mj,nj-4,j)*jm4(nj,j)   
			okq(m,n) = 0.25d0*(ff(k,q,mj,j)+ff(k,q,nj,j))*d4*f(sbkq(k,q))          
		   else if(q.eq.5) then
			d5=delta(mj,nj+5,j)*jp5(nj,j)+sbkq(k,q)*delta(mj,nj-5,j)*jm5(nj,j)   
			okq(m,n) = 0.25d0*(ff(k,q,mj,j)+ff(k,q,nj,j))*d5*f(sbkq(k,q))           
		   else if(q.eq.6) then
			d6=delta(mj,nj+6,j)*jp6(nj,j)+sbkq(k,q)*delta(mj,nj-6,j)*jm6(nj,j)    
			okq(m,n) = 0.25d0*(ff(k,q,mj,j)+ff(k,q,nj,j))*d6*f(sbkq(k,q))        
		   endif
	        endif
 20	   continue
 10     continue
      return
      end  
c--------------------------------------------------
c <jm| jp^|q| or jm^|q|  |nj>
c--------------------------------------------------
      real*8 function jop(q,mj,nj,j)
	implicit none
        integer q
        real*8 mj,nj,j

	real*8 jp,jp2,jp3,jp4,jp5,jp6,jm,jm2,jm3,jm4,jm5,jm6
        real*8 delta

        if(q.eq.0) then
           jop = delta(mj,nj,j)
        else if(q.eq.1) then
           jop = delta(mj,nj+1,j)*jp(nj,j)
        else if(q.eq.2) then
           jop = delta(mj,nj+2,j)*jp2(nj,j)
        else if(q.eq.3) then  
           jop = delta(mj,nj+3,j)*jp3(nj,j)
        else if(q.eq.4) then
           jop = delta(mj,nj+4,j)*jp4(nj,j)
        else if(q.eq.5) then
           jop = delta(mj,nj+5,j)*jp5(nj,j)
        else if(q.eq.6) then
           jop = delta(mj,nj+6,j)*jp6(nj,j)
        else if(q.eq.-1) then
           jop = delta(mj,nj-1,j)*jm(nj,j)
        else if(q.eq.-2) then
           jop = delta(mj,nj-2,j)*jm2(nj,j)
        else if(q.eq.-3) then  
           jop = delta(mj,nj-3,j)*jm3(nj,j)
        else if(q.eq.-4) then
           jop = delta(mj,nj-4,j)*jm4(nj,j)
        else if(q.eq.-5) then
           jop = delta(mj,nj-5,j)*jm5(nj,j)
        else if(q.eq.-6) then
           jop = delta(mj,nj-6,j)*jm6(nj,j)
        end if
      return
      end
c--------------------------------------------------
c calculate the full stevens operators
c--------------------------------------------------
      real*8 function full_okq(k,q,mj,nj,j)
	implicit none
        real*8 jop,ff
        real*8 j,mj,nj
        integer k,q
        full_okq  = 0.5d0*jop(q,mj,nj,j)*( ff(k,q,mj,j)+ff(k,q,nj,j) )
      return
      end
c-----------------------------
      real*8 function f20(nj,j)
	implicit none
	real*8 nj,j
	f20 = 3*nj**2 - j*(j+1)
      return
      end   
c------------------------------
      real*8 function f21(nj,j)
	implicit none
	real*8 nj,j
	f21 = nj
      return
      end
c------------------------------
      real*8 function f22(nj,j)
	implicit none
	real*8 nj,j
	f22 = 1.0d0
      return
      end                 
c------------------------------
      real*8 function f40(nj,j)
	implicit none   
	real*8 nj,j
	f40 = 35*nj**4 - 30*j*(j+1)*nj**2 + 25*nj**2 - 6*j*(j+1) 
     *      + 3*j**2*(j+1)**2
      return
      end          
c------------------------------
      real*8 function f41(nj,j)
	implicit none
	real*8 nj,j
	f41 = 7*nj**3 - 3*j*(j+1)*nj - nj
      return
      end
c------------------------------
      real*8 function f42(nj,j)
	implicit none   
	real*8 nj,j
	f42 = 7*nj**2 - j*(j+1) - 5
      return
      end      
c------------------------------
      real*8 function f43(nj,j)
	implicit none
	real*8 nj,j
	f43 = nj
      return
      end
c------------------------------
      real*8 function f44(nj,j)
	implicit none   
	real*8 nj,j
	f44 = 1.0d0
      return
      end                      
c------------------------------
      real*8 function f60(nj,j)
	implicit none   
	real*8 nj,j
	f60 = 231*nj**6 - 315*j*(j+1)*nj**4 + 735*nj**4          
     *      + 105*j**2*(j+1)**2*nj**2 - 525*j*(j+1)*nj**2
     *      + 294*nj**2 - 5*j**3*(j+1)**3 + 40*j**2*(j+1)**2
     *      - 60*j*(j+1)
      return
      end          
c------------------------------
      real*8 function f61(nj,j)                                                    
	implicit none
	real*8 nj,j
	f61 = 33*nj**5 - ( 30*j*(j+1) - 15 )*nj**3 
     *      - 10*j*(j+1)*nj + 5*j**2*(j+1)**2*nj
     *      +  12*nj
      return
      end
c------------------------------
      real*8 function f62(nj,j)
	implicit none   
	real*8 nj,j
	f62 = 33*nj**4 - 18*j*(j+1)*nj**2 - 123*nj**2
     *      + j**2*(j+1)**2 + 10*j*(j+1) + 102                           
      return
      end      
c------------------------------     
      real*8 function f63(nj,j)
	implicit none
	real*8 nj,j
	f63 = 11*nj**3 - 3*j*(j+1)*nj - 59*nj 
      return
      end
c------------------------------  
      real*8 function f64(nj,j)
	implicit none   
	real*8 nj,j
	f64 = 11*nj**2 - j*(j+1) - 38          
      return
      end                      
c------------------------------
      real*8 function f65(nj,j)
	implicit none
	real*8 nj,j
	f65 = nj
      return
      end
c------------------------------  
      real*8 function f66(nj,j)
	implicit none   
	real*8 nj,j
	f66 = 1.0d0
      return
      end               
c-------------------------------------
c testing if the complex number c is 0
c-------------------------------------
      subroutine cifnull(c)
	implicit none
	complex*16 c,i
	real*8     re,im
	i = cmplx(0.0d0,1.0d0)
	re= 0.5d0*( c+conjg(c) )
	im= 0.5d0*(-c+conjg(c) )*i
	call ifnull(re)
	call ifnull(im)
	c=cmplx(re,im)
      return
      end
c----------------------------------
c testing if the real number r is 0
c----------------------------------
      subroutine ifnull(r)
	implicit none
	real*8 r,macheps
	macheps = 1.0d-14
	if(abs(r).le.macheps) r=0.0d0
      return
      end
c--------------------
c  jp(n) = <n+1|j+|n> 
c--------------------
      real*8 function jp(nj,j)
	real*8 nj,j
	jp = 0.0d0
	if(abs(nj).le.j) jp = sqrt(j*(j+1)-nj*(nj+1))
      return
      end
c--------------------
c  jm(n) = <n-1|j-|n>  
c--------------------
      real*8 function jm(nj,j)
	real*8 nj,j
	jm = 0.0d0
	if(abs(nj).le.j) jm = sqrt(j*(j+1)-nj*(nj-1))
      return
      end
c-----------------------
c  jp2(n) = <n+2|j+^2|n> 
c-----------------------
       real*8 function jp2(nj,j)
	real*8 nj,j,jp
	jp2 = jp(nj+1,j)*jp(nj,j)
      return
      end
c-----------------------
c  jm2(n) = <n-2|j-^2|n> 
c-----------------------
       real*8 function jm2(nj,j)
	real*8 nj,j,jm
	jm2 = jm(nj-1,j)*jm(nj,j)                                                                  
      return
      end                          
c-----------------------
c  jp3(n) = <n+3|j+^3|n> 
c-----------------------
       real*8 function jp3(nj,j)
	real*8 nj,j,jp
	jp3 = jp(nj+2,j)*jp(nj+1,j)*jp(nj,j)
      return
      end
c-----------------------
c  jm3(n) = <n-3|j-^3|n> 
c-----------------------
       real*8 function jm3(nj,j)
	real*8 nj,j,jm
	jm3 = jm(nj-2,j)*jm(nj-1,j)*jm(nj,j)
      return
      end                     
c-----------------------
c  jp4(n) = <n+4|j+^4|n> 
c-----------------------
       real*8 function jp4(nj,j)
	real*8 nj,j,jp
	jp4 = jp(nj+3,j)*jp(nj+2,j)*jp(nj+1,j)*jp(nj,j)
      return
      end
c-----------------------
c  jm4(n) = <n-4|j-^4|n> 
c-----------------------
       real*8 function jm4(nj,j)
	real*8 nj,j,jm
	jm4 = jm(nj-3,j)*jm(nj-2,j)*jm(nj-1,j)*jm(nj,j)
      return
      end              
c-----------------------
c  jp5(n) = <n+5|j+^5|n> 
c-----------------------
       real*8 function jp5(nj,j)
	real*8 nj,j,jp
	jp5 = jp(nj+4,j)*jp(nj+3,j)*jp(nj+2,j)*jp(nj+1,j)*jp(nj,j)
      return
      end
c-----------------------
c  jm5(n) = <n-5|j-^5|n> 
c-----------------------
       real*8 function jm5(nj,j)
	real*8 nj,j,jm
	jm5 = jm(nj-4,j)*jm(nj-3,j)*jm(nj-2,j)*jm(nj-1,j)*jm(nj,j)
      return
      end                        
c-----------------------
c  jp6(n) = <n+6|j+^6|n> 
c-----------------------
       real*8 function jp6(nj,j)
	real*8 nj,j,jp
	jp6 = jp(nj+5,j)*jp(nj+4,j)*jp(nj+3,j)*jp(nj+2,j)*jp(nj+1,j)*jp(nj,j)
      return
      end
c-----------------------
c  jm6(n) = <n-6|j-^6|n> 
c-----------------------
       real*8 function jm6(nj,j)
	real*8 nj,j,jm
	jm6 = jm(nj-5,j)*jm(nj-4,j)*jm(nj-3,j)*jm(nj-2,j)*jm(nj-1,j)*jm(nj,j)
      return
      end      
c---------------------------------------
      real*8 function epsilon(k,q)
c---------------------------------------
      implicit none
      integer k,q
      real*8 eps(49)

c ... epsilon_kq
      data eps/            1.00000000000000,     
     1  0.707106781186547, 1.00000000000000, -0.707106781186547,     
     2  0.612372435695795, 1.22474487139159,  0.500000000000000,     
     2 -1.22474487139159,  0.612372435695795,
     3  0.559016994374947, 1.36930639376292,  0.433012701892219,
     3  0.500000000000000,-0.433012701892219, 1.36930639376292,
     3 -0.559016994374947,
     4  0.522912516583797, 1.47901994577490,  0.395284707521047,     
     4  0.559016994374947, 0.125000000000000,-0.559016994374947,
     4  0.395284707521047,-1.47901994577490,  0.522912516583797,     
     5  0.496078370824611, 1.56873754975139,  0.369754986443726,     
     5  0.603807364424560, 0.114108866146910, 0.125000000000000,     
     5 -0.114108866146910, 0.603807364424560,-0.369754986443726,     
     5  1.56873754975139, -0.496078370824611,     
     6  0.474958879799083, 1.64530582263602,  0.350780380010057,     
     6  0.640434422872475, 0.106739070478746, 0.135015431216830,     
     6  0.062500000000000,-0.135015431216830, 0.106739070478746,     
     6 -0.640434422872475, 0.350780380010057,-1.64530582263602,     
     6  0.474958879799083/

        epsilon = eps(k*(k+1) + q+1)

      return
      end
c---------------------------------------
      real*8 function omega(k,q)
c---------------------------------------
      implicit none
      integer k,q
      real*8 oma(49)
c
c ... omega_kq
      data oma/                     1.0,
     1                         1.0, 1.0, 1.0,
     2                     1.0,1.0, 1.0, 1.0,1.0,
     3                 1.0,1.0,1.0, 1.0, 1.0,1.0,1.0,
     4             1.0,1.0,1.0,1.0, 1.0, 1.0,1.0,1.0,1.0,
     5         1.0,1.0,1.0,3.0,3.0, 1.0, 3.0,3.0,1.0,1.0,1.0,
     6     1.0,1.0,1.0,1.0,3.0,3.0, 1.0, 3.0,3.0,1.0,1.0,1.0,1.0 /

        omega = oma( k*(k+1) + q+1)

        return
        end
c---------------------------------------
c initialize some rare earth constants
c---------------------------------------
      subroutine i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)        
	implicit none
	integer nre
	real*8 gj,dimj,alphaj,betaj,gammaj,r2,r4,r6                             
	real*8 ggj(13),ddimj(13),aalphaj(13),bbetaj(13),ggammaj(13)
	real*8 rr2(13),rr4(13),rr6(13)
c       ------------------ 
	ggj( 1) = 6.0d0/7.
	ggj( 2) = 4.0d0/5.
	ggj( 3) = 8.0d0/11.
	ggj( 4) = 3.0d0/5.
	ggj( 5) = 2.0d0/7.
	ggj( 6) = 0.0d0
	ggj( 7) = 2.0d0
	ggj( 8) = 3.0d0/2.
	ggj( 9) = 4.0d0/3.
	ggj(10) = 5.0d0/4.
	ggj(11) = 6.0d0/5.
	ggj(12) = 7.0d0/6.
	ggj(13) = 8.0d0/7.
	gj = ggj(nre)
c       ------------------
	ddimj( 1) =  6.0d0
	ddimj( 2) =  9.0d0
	ddimj( 3) = 10.0d0
	ddimj( 4) =  9.0d0
	ddimj( 5) =  6.0d0
	ddimj( 6) =  1.0d0
	ddimj( 7) =  8.0d0
	ddimj( 8) = 13.0d0
	ddimj( 9) = 16.0d0
	ddimj(10) = 17.0d0
	ddimj(11) = 16.0d0
	ddimj(12) = 13.0d0
	ddimj(13) =  8.0d0
	dimj = ddimj(nre)
c       -----------------------------------------------
c       for the values of alphaj,betaj and gammaj
c       see: Abragam and Bleaney, Electronic Paramagnetic
c            Resonance of Transition Ions, 1970
c            appendix B, table 20, page 874-875
c
	aalphaj( 1) = -1.0d0 * 2/5/7
	aalphaj( 2) = -1.0d0 * 2*2*13/3/3/5/5/11
	aalphaj( 3) = -1.0d0 * 7/3/3/11/11
	aalphaj( 4) =  1.0d0 * 2*7/3/5/11/11
	aalphaj( 5) =  1.0d0 * 13/3/3/5/7
	aalphaj( 6) =  1.0d0 * 0
	aalphaj( 7) =  1.0d0 * 0
	aalphaj( 8) = -1.0d0 * 1/3/3/11
	aalphaj( 9) = -1.0d0 * 2/3/3/5/7
	aalphaj(10) = -1.0d0 * 1/2/3/3/5/5
	aalphaj(11) =  1.0d0 * 2*2/3/3/5/5/7
	aalphaj(12) =  1.0d0 * 1/3/3/11
	aalphaj(13) =  1.0d0 * 2/3/3/7
	alphaj = aalphaj(nre)
c       -------------------------------
	bbetaj( 1)  =  1.0d0 * 2/3/3/5/7
	bbetaj( 2)  = -1.0d0 * 2*2/3/3/5/11/11
	bbetaj( 3)  = -1.0d0 * 2*2*2*17/3/3/3/11/11/11/13
	bbetaj( 4)  =  1.0d0 * 2*2*2*7*17/3/3/3/5/11/11/11/13
	bbetaj( 5)  =  1.0d0 * 2*13/3/3/3/5/7/11
	bbetaj( 6)  =  1.0d0 * 0
	bbetaj( 7)  =  1.0d0 * 0
	bbetaj( 8)  =  1.0d0 * 2/3/3/3/5/11/11
	bbetaj( 9)  = -1.0d0 * 2*2*2/3/3/3/5/7/11/13
	bbetaj(10)  = -1.0d0 * 1/2/3/5/7/11/13
	bbetaj(11)  =  1.0d0 * 2/3/3/5/7/11/13
	bbetaj(12)  =  1.0d0 * 2*2*2/3/3/ 3/ 3/5/11/11
	bbetaj(13)  = -1.0d0 * 2/3/5/7/11
	betaj = bbetaj(nre)
c       -----------------------
	ggammaj( 1) =  1.0d0 * 0
	ggammaj( 2) =  1.0d0 * 2*2*2*2*17/3/3/3/3/5/7/11/11/13
	ggammaj( 3) = -1.0d0 * 5*17*19/3/3/3/7/11/11/11/13/13
	ggammaj( 4) =  1.0d0 * 2*2*2*17*19/3/3/3/7/11/11/11/13/13
	ggammaj( 5) =  1.0d0 * 0
	ggammaj( 6) =  1.0d0 * 0
	ggammaj( 7) =  1.0d0 * 0
	ggammaj( 8) = -1.0d0 * 1/3/3/3/3/7/11/11/13
	ggammaj( 9) =  1.0d0 * 2*2/3/3/3/7/11/11/13/13
	ggammaj(10) = -1.0d0 * 5/3/3/3/7/11/11/13/13
	ggammaj(11) =  1.0d0 * 2*2*2/3/3/3/7/11/11/13/13
	ggammaj(12) = -1.0d0 * 5/3/3/3/3/7/11/11/13
	ggammaj(13) =  1.0d0 * 2*2/3/3/3/7/11/13
	gammaj = ggammaj(nre)
c       ---------------
c       for the values of <r2>,<r4>,<r6> look to:
c       Freeman and Desclaux, Journal of Magnetism
c       and Magnetic Materials 12 (1979) 11 ff
c       [<r2>] = a0**2       
c
	rr2( 1) = 1.309
	rr2( 2) = 1.1963
	rr2( 3) = 1.114
	rr2( 4) = 1.0353
	rr2( 5) = 0.9743
	rr2( 6) = 0.9175
	rr2( 7) = 0.8671
	rr2( 8) = 0.8220
	rr2( 9) = 0.7814
	rr2(10) = 0.7446
	rr2(11) = 0.7111
	rr2(12) = 0.6804
	rr2(13) = 0.6522
	r2 = rr2(nre)
c       ---------------
c       [<r4>] = a0**4       
	rr4( 1) = 3.964
	rr4( 2) = 3.3335
	rr4( 3) = 2.910
	rr4( 4) = 2.5390
	rr4( 5) = 2.260
	rr4( 6) = 2.020
	rr4( 7) = 1.820
	rr4( 8) = 1.651
	rr4( 9) = 1.505
	rr4(10) = 1.379
	rr4(11) = 1.270
	rr4(12) = 1.174
	rr4(13) = 1.089
	r4 = rr4(nre)
c       ---------------
c       [<r6>] = a0**6     
	rr6( 1) = 23.31
	rr6( 2) = 18.353
	rr6( 3) = 15.03 
	rr6( 4) = 12.546
	rr6( 5) = 10.55
	rr6( 6) = 9.039
	rr6( 7) = 7.831
	rr6( 8) = 6.582
	rr6( 9) = 6.048
	rr6(10) = 5.379
	rr6(11) = 4.816
	rr6(12) = 4.340
	rr6(13) = 3.932
	r6 = rr6(nre)
      return
      end
