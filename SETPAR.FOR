!	Subroutine SETPAR	-	Sets parameter values

	subroutine SETPAR (line, l)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:FOCUS.INC'
        include 'CF_SOURCES:OVERALL.INC'

        integer same, ncftot, no_cf_par, n0, no_0_par

	character*2 bind_no
	character*5 par_no
	character*(*) line, c3*3	
	integer*4 l, n, getari, iarray(MAX_PAR), i, j, par_len, itoc, bind_len, itocbf, line_len, getlin, k, getarr
	real*4 parray(3), pold, dummy 

        ncftot = no_cf_par(symmetry,nre)                                  
        n0     = no_0_par()
	n = getari (line, l, iarray, MAX_PAR)
	if (n .eq. 0) then
		n = nptot
		do i = 1,nptot
			iarray(i) = i
		end do
	else if (n .lt. 0) then
		call remark ('Error in parameter number list')
		return
	end if                                                 
	call remark ('Give the parameter value (with min and max if required)')
	do i = 1,n
		j       = iarray(i)
		par_len = itoc (j, par_no)
                pold    = p(j)
		if (j .lt. 1 .or. j .gt. nptot) then
			call remark ('"'//par_no(1:par_len)//'" is an invalid parameter number')
		else if (inr(j) .gt. 0) then
			bind_len = itoc (inr(j), bind_no)
			call remark ('Parameter '//par_no(1:par_len)//' is bound to parameter '//bind_no(1:bind_len))
		else
   			par_len = itocbf (j, par_no)
			call prompt (par_no//'  '//nam(j)//':')
			line_len = getlin (line)
			if (line_len .gt. 0) then
				l = 1
				k = getarr (line(1:line_len), l, parray, 3)

				if(k.eq.2) then
		                   par_len = itoc (j, par_no)
                                   call remark('ERROR: setup of parameter '//par_no(1:par_len)//' incorrect.')  
                                   call remark('FOCUS sets the parameter to its old values.')
                                   goto 99
				end if

				if (k .gt. 0) p(j) = parray(1)
				if (k .gt. 1) pmin(j) = parray(2)
				if (k .gt. 2) pmax(j) = parray(3)

				if(k.gt.2 .and. pmin(j).gt.pmax(j)) then
                                    dummy   = pmax(j)
				    pmax(j) = pmin(j)
				    pmin(j) = dummy
                                end if

                                if( (setwx.eq.1.and.j.eq.n0+1).or.(setxr.eq.1.and.j.eq.indvold(1)) )
     *                               ratioenergy = 1.0d0


                                if(setwx.eq.1 .and. j.eq.n0+2) then 

                                    if( p(j).lt.-1.0d0 .or. p(j).gt.1.0d0 ) then
		                        par_len = itoc (j, par_no)
                                        call remark('The x-parameter has to be between -1 and 1.')
                                        call remark('FOCUS  sets the parameter '//par_no(1:par_len)//'to the old value.')
                                        p(j) = pold
                                    end if
                                    if(pmin(j).ne.-1.0d0) then
		                       par_len = itoc (j, par_no)
                                       pmin(j) = -1.0d0
                                       call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//' to -1.')
                                    end if
                                    if(pmax(j).ne.1.0d0) then
	  	                       par_len = itoc (j, par_no)
                                       pmax(j) =  1.0d0
                                       call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//' to -1.')
                                    end if

                                end if

                                if(setxr.eq.1 .and. same(j,indvold(1)).eq.1 .and. j.ne.indvold(1)) then

                                    if(p(j).lt.-1.0d0 .or. p(j).gt.1.0d0 ) then
                                       p(j) = pold
           	                       par_len = itoc (j, par_no)
                                       call remark('The x-parameter has to be between -1 and 1.')
                                       call remark('FOCUS  sets the parameter '//par_no(1:par_len)//'to the old value.')
                                    end if
                                    if(pmin(j).ne.-1.0d0) then
          	                       par_len = itoc (j, par_no)
                                       pmin(j) = -1.0d0
                                       call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//' to -1.')
                                    end if
                                    if(pmax(j).ne.1.0d0) then
          	                       par_len = itoc (j, par_no)
                                       pmax(j) =  1.0d0
                                       call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//' to  1.')
                                    end if
                                end if

                                if( setakq.eq.1 .or. setbkq.eq.1 .or. setvkq.eq.1 .or.
     *                              (setxr.eq.1  .and. j.eq.indvold(1))           .or.
     *                              (setwx.eq.1  .and. j.eq.n0+1      )                ) then
                                    if( pmin(j).ne.pmax(j) ) then 
    			                if(p(j).lt.pmin(j)) then
                                           p(j) = pmin(j)
	                                   par_len = itoc (j, par_no)
                                           call remark('The parameter you want to set is outside the limits.')
                                           call remark('FOCUS sets the parameter '//par_no(1:par_len)//' to the minimal value.')
					end if 
					if(p(j).gt.pmax(j) ) then
                                           p(j) = pmax(j)
          	                           par_len = itoc (j, par_no)
                                           call remark('The parameter you want to set is outside the limits.')
                                           call remark('FOCUS sets the parameter '//par_no(1:par_len)//' to the maximal value.')
					end if 
				    end if
				end if

			end if
		end if
 99 	end do
	do i = 1,nptot
		if (inr(i) .gt. 0) p(i) = p(inr(i)) * rp(i)
	end do

	return
	end
