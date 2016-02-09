!	Subroutine PMFIT	-	Performs least squares fitting algorithm
!       26.4.95  Peter Fabi  adapted for x-parameter
!                            changes use ptrans and pback instead of ptrans and pback

	subroutine PMFIT (fun)

	external fun
	character*1 answer
	integer*4 n_iter, i, no_p, j, k, l

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:FOCUS.PAR'
	include 'CF_SOURCES:MONTE.INC'
	include 'CF_SOURCES:LIMITPAR.INC'
 
	real*8 po(MAX_PAR), delta_p(MAX_PAR), div, chisq_old
        real*8 pi, pold

!Store old parameters
	do i = 1,nptot
		po(i) = p(i)
	end do
!Transform limited parameters with parameter numbers lower than no_0_par()-5 into circular coordinates 
!and calculate initial Chi**2
	call ptrans (np, p, inp, pmin, pmax)   
	call chisq_calc (fun)
	chisq_old = chisq
!Loop over allowed number of iterations
	status = 0
	n_iter = 0
	do while (n_iter .lt. max_iter)
		n_iter = n_iter + 1
		do i = 1,nptot
			pf(i) = p(i)	 
		end do
		call dermat (fun)	!Calculate matrices of derivatives
		call delmat
		call hess		!Calculate Hessian matrix
		call invert (no_p)	!Invert Hessian matrix
		if (status .gt. 0) then
			if (status .eq. 2) then 
				message = 'ERROR: parameter xx affects none of the points'
				write (message(18:19), '(i2)') no_p
			else if (status .eq. 3) then
				message = 'ERROR: two or more parameters linearly dependent'
			end if
			go to 20
		end if
		call chp (delta_p)	!Calculate optimum parameter shifts
!Set new parameter values
		do i = 1,np
			j = inp(i)
			p(j) = p(j) + delta_p(i)
                        if(pmax(j).gt.pmin(j)) then                 ! added by P.fabi
                             if(p(j).gt.up_limit ) p(j) = up_limit  !
                             if(p(j).lt.low_limit) p(j) = low_limit !
                        end if                                      !
		end do
		div = 2D0
!Calculate new Chi**2
10		call chisq_calc (fun)
		if (output .ne. 0 .and. output .ne. 10) write (*, 100) n_iter, chisq_old, chisq
		write (41, 100) n_iter, chisq_old, chisq
		call pback (np, nptot, p, inp, pmin, pmax, inr, rp)
		if (output .eq. 2 .or. output .eq. 12) write (*, 200) (p(inp(i)), i = 1,np)
		write (41, 200) (p(inp(i)), i = 1,np)
		call ptrans (np, p, inp, pmin, pmax)
		if (abs(chisq-chisq_old) .lt. acc) then 			!Chi**2 converged
			status = 0
			message = 'Job completed: Agreement factor converged'
			go to 20
		else if (chisq .gt. chisq_old .and. div .lt. max_div) then      !Chi**2 diverging 
			!Attempt repair by reducing parameter step
			do i = 1,np
				j = inp(i)
				p(j) = p(j) - delta_p(i) / div
                                if(pmax(j).gt.pmin(j)) then               ! added by P.fabi
                                   if(p(j).gt.up_limit ) p(j) = up_limit  !
                                   if(p(j).lt.low_limit) p(j) = low_limit !
                                end if                                    !
			end do
			div = 2.0 * div
			go to 10
		else if (chisq .gt. chisq_old) then				!Chi**2 diverged
			!Reset parameters to those with lowest Chi**2
			chisq = chisq_old
			do i = 1,np
				j = inp(i)
				p(j) = pf(j)
			end do
			status = 4
			message = 'Job terminated: Agreement factor diverged'
			go to 20
		end if
		chisq_old = chisq
	end do
										!Chi**2 not converged after max. no. of iterations
	status = 1
	message = 'Job terminated: Max. iterations completed'

20	nf = np
                                             ! added by P.fabi to make sure 
	call chisq_calc (fun)                ! that the best parameter set
	                                     ! is taken for the function call FUN
!Calculate parameter errors from variance-covariance matrix
	do i = 1,np
		inf(i) = inp(i)
		j = inp(i)
		if (status .ne. 2 .and. status .ne. 3) then
			pfsig(j) = sqrt (-H(i,i)*chisq)
			if (pmax(j) .gt. pmin(j)) then
				pfsig(j) = - 0.5 * (pmax(j) - pmin(j)) * sin(p(j)) * pfsig(j)
			end if
		else
			pfsig(j) = 0.0
		end if
	end do
	if (status .ne. 2 .and. status .ne. 3) call print_corr (41)
	call pback (np, nptot, p, inp, pmin, pmax, inr, rp)
!Store fitted parameters 
	do i = 1,nptot
		pf(i) = p(i)
		pfmin(i) = pmin(i)
		pfmax(i) = pmax(i)
		rf(i) = rp(i)
		inrf(i) = inr(i)
	end do
	call lib$date_time (fit_time)
!Display fitted parameters at terminal
	call fitdis
!Restore old parameters if requested
        if(if_plot.ne.1) then
	   if (answer('Do you wish to restore old parameters [N]?') .eq. 'Y') then
		do i = 1,nptot
			p(i) = po(i)
		end do
	        call ptrans (np, p, inp, pmin, pmax) ! added by P.Fabi to make sure  
	        call chisq_calc (fun)                ! that the best parameter set
 	        call pback (np, nptot, p, inp, pmin, pmax, inr, rp)
	        chisq_old = chisq                    ! is taken for the function call FUN
	   end if
        end if
100	format (/, 3x, 'Iteration ', i4, / , 3x, 'Old Chi**2 = ', d20.8, '          New Chi**2 = ', d20.8)
200	format (3x, 'Parameters', <(np/5+1)>(/, 1x, 5d16.6))
		     
	return
	end
