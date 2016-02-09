	real*8 function self (Ei, hw, phi, mutS, mutA, gamma)


	real*8 Ei, hw, phi, mutS, mutA, gamma, muti, mutf

	if (mutA .eq. 0.0 .and. mutS .eq. 0.0) then
		self = 1.0
	else if (gamma .lt. 89.9 .or. gamma .gt. 90.1) then
		muti = (mutS + mutA * sqrt(25.299/Ei)) / cosd(gamma)
		if (phi-gamma .lt. 89.9) then
			mutf = (mutS + mutA * sqrt(25.299/(Ei - hw))) / cosd(phi-gamma)
			self = (mutf - muti) / (exp(-muti) - exp(-mutf))
		else if (phi-gamma .gt. 90.1) then
			mutf = (mutS + mutA * sqrt(25.299/(Ei - hw))) / cosd(phi-gamma)
			self = (mutf - muti) / (exp(mutf-muti) - 1.0) 
		else
			self = 0.0
		end if
	else
		self = 0.0
	end if

	return
	end
