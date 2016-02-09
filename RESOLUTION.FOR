!-----------------------------------------------------------------------------
!Returns value of convolution of Lorentzian with an exponential decay 
!-----------------------------------------------------------------------------
! origin: het$disk:[hetmgr.library.peaks.ref]
! author: T G Perring
!         adapted for focus by P O FABI
!-----------------------------------------------------------------------------
	real*8 function CONV_SHAPE (hw, gamma, gamma_exp, sigma_exp,
     *                              no_lineshape, intensity, temperature)

        implicit none
	integer*4 i, j, ifail
	real*8  hw, gamma, gamma_exp, sigma_exp, xmin, dx, x(21), y(21), err
        integer no_lineshape
        real*8   intensity,temperature
	real*8   SHAPE,EXP_DECAY
        external SHAPE,EXP_DECAY,d01gaf
	
	if (gamma .le. 0.0) then
		CONV_SHAPE = intensity*EXP_DECAY (hw, gamma_exp, sigma_exp)
	else if (gamma_exp .le. 0.0 .and. sigma_exp .le. 0.0) then
		CONV_SHAPE = SHAPE (no_lineshape, hw, 0.0d0, 2.0*gamma, intensity, temperature) 
	else
		xmin = min(-4*gamma_exp, -2*sigma_exp)
		dx = (2*sigma_exp - xmin) / 20.0D0
		do i = 1,21
			x(i) = xmin + i * dx
			y(i) = SHAPE (no_lineshape, hw-x(i), 0.0d0, 2.0*gamma, intensity, temperature) 
     *                       * EXP_DECAY (x(i), gamma_exp, sigma_exp)
		end do
		ifail = 0
		call d01gaf (x, y, 21, CONV_SHAPE, err, ifail)
	end if
	return
	end
!-----------------------------------------------------------------------------
!Returns value of decaying exponential convolved with a Gaussian
!-----------------------------------------------------------------------------
! origin: het$disk:[hetmgr.library.peaks.ref]
! author: T G Perring
!         adapted for focus by P O FABI
!-----------------------------------------------------------------------------
	real*8 function EXP_DECAY (hw, gamma, sigma)

	real*8 hw, gamma, fwhm, sigma, GAUSS, x, ERFC, pi
	parameter (pi = 3.1415926535897932)
        external GAUSS


	if (sigma .le. 0.0d0 .and. gamma .le. 0.0d0) then
		EXP_DECAY = 0.0d0
	else if (gamma .le. 0.0d0) then
                fwhm      = sigma*sqrt(log(256.0d0))
		EXP_DECAY = GAUSS(hw, 0.0d0, fwhm, 1.0d0)
	else if (sigma .le. 0.0d0) then
		x = hw - gamma
		if (x .ge. 0.0d0) then
			EXP_DECAY = 0.0d0
		else if (x .lt. 50.0D0*gamma) then
			EXP_DECAY = exp_control(x / gamma) / gamma
		else
			EXP_DECAY = 0.0d0
		end if
	else
		x = hw - gamma
		if (x .lt. 5.0*sigma) then
			EXP_DECAY = exp_control(0.5D0*sigma**2/gamma**2 + x/gamma) 
     *                            * ERFC((sigma/gamma + x/sigma) / sqrt(2D0)) / (2D0*gamma)
		else
			EXP_DECAY = 0.0d0
		end if
	end if

	return
	end 
!-----------------------------------------------------------------------------
! Complementary Error Function
! Algorithm from Abramovitz and Stegun p.299
! Accurate to 3E-07 for both +ve and -ve x.
!-----------------------------------------------------------------------------
! origin: het$disk:[hetmgr.library.peaks.ref]
! author: T G Perring
!-----------------------------------------------------------------------------
	real*8 function ERFC (x)

	real*8 a1, a2, a3, a4, a5, a6, x, z
	data a1, a2, a3 /0.0705230784, 0.0422820123, 0.0092705272/
	data a4, a5, a6 /0.0001520143, 0.0002765672, 0.0000430638/

	z = abs(x)
	ERFC = 1.0 + a1*z + a2*z**2 + a3*z**3 + a4*z**4 + a5*z**5 + a6*z**6

	if (ERFC .ge. 10.0) then
		ERFC = 0.0
	else
		ERFC = 1.0 / ERFC**16
	end if
	if (x .lt. 0.0) ERFC = 2.0 - ERFC

	return
	end
