!	Subroutine LIMPAR	-	Sets parameter limits

	subroutine limpar (line, l)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:CONTROL.INC'
        include 'CF_SOURCES:FOCUS.INC'

        integer same, n0, no_0_par, ncftot, no_cf_par

	character*(*) line
	character*5 par_no
	integer*4 l, n, getari, iarray(MAX_PAR), line_len, getlin, par_len, itoc, i, j, k, getarr
	real*4 parray(2)

	n = getari (line, l, iarray, MAX_PAR)
	if (n .eq. 0) then
		call prompt ('Give parameters to be limited :')
		line_len = getlin (line)
		if (line_len .gt. 0) then
			l = 1
			n = getari (line(1:line_len), l, iarray, MAX_PAR)
		end if
	end if
	if (n .lt. 0) then
		call remark ('Error in parameter number list')
		return
	end if
	call remark ('Give the min and max value of each parameter')
        n0     = no_0_par()
        ncftot = no_cf_par(symmetry,nre)
	do i = 1,n
		j = iarray(i)
		par_len = itoc (j, par_no)
		if (j .lt. 1 .or. j .gt. nptot) then
			call remark ('"'//par_no(1:par_len)//'" is an invalid parameter number')
		else
 			call prompt (nam(j)//':')
			line_len = getlin (line)
			if (line_len .gt. 0) then
				l = 1
				k = getarr (line(1:line_len), l, parray, 2)
				if (k .gt. 0) pmin(j) = parray(1)
				if (k .gt. 1) pmax(j) = parray(2)

                                if(setakq.eq.1 .or. setbkq.eq.1 .or. setvkq.eq.1) then
				   if( p(j).lt.pmin(j) ) then
                                       pmin(j) = p(j)
                                       call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//
     *                                             ' to the parameter value.')
          	      		   end if
				   if( p(j).gt.pmax(j) ) then
                                       pmax(j) = p(j)
		                       call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//
     *                                             ' to the parameter value.')
          	                   end if
				else if(setxr.eq.1) then
				   if( p(j).lt.pmin(j) ) then
                                       pmin(j) = p(j)
		                       call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//
     *                                             ' to the parameter value.')
          	      		   end if
				   if( p(j).gt.pmax(j) ) then
                                       pmax(j) = p(j)
		                       call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//
     *                                             ' to the parameter value.')
          	      		    end if
                                   if(same(j,indvold(1)).eq.1 .and. j.ne.indvold(1)) then
                                     if(pmin(j).ne.-1.0d0) then
                                        pmin(j) = -1.0d0
                                        call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//
     *                                              ' to -1.')
                                     end if
                                     if(pmax(j).ne.1.0d0) then
                                        pmax(j) =  1.0d0
                                        call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//
     *                                              ' to 1.')
                                     end if
                                   else if( j.eq.indvold(1) ) then
				       if(p(j).lt.pmin(j) ) then
                                          pmin(j) = p(j)
		                          call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//
     *                                                ' to the parameter value.')
          	      		       end if
				       if(p(j).gt.pmax(j) ) then
                                          pmax(j) = p(j)
		                          call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//
     *                                                ' to the parameter value.')
          	      		       end if
                                   end if
                                else if(setwx.eq.1) then
                                  if(j.eq.n0+1) then
				       if(p(j).lt.pmin(j) ) then
                                          pmin(j) = p(j)
		                          call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//
     *                                                ' to the parameter value.')
          	      		       end if
				       if(p(j).gt.pmax(j) ) then
                                          pmax(j) = p(j)
		                          call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//
     *                                                ' to the parameter value.')
          	      		       end if
				  end if
                                  if(j.eq.n0+2) then
                                    if(pmin(j).ne.-1.0d0) then
                                       pmin(j) = -1.0d0
                                       call remark('FOCUS sets lower limit of parameter '//par_no(1:par_len)//
     *                                             ' to -1.')
                                    end if
                                    if(pmax(j).ne.1.0d0) then
                                       pmax(j) =  1.0d0
                                       call remark('FOCUS sets upper limit of parameter '//par_no(1:par_len)//
     *                                             ' to 1.')
                                    end if
                                  end if
                                  
                                end if
			end if
		end if
	end do

	return
	end
