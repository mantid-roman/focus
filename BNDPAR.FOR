!	Subroutine BNDPAR	-	Binds parameters together by fixed ratios

	subroutine BNDPAR (line, l)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:CONTROL.INC'

	character*(*) line
	integer*4 l, n, getari, line_len, getlin, iarray(MAX_PAR), i, j, k, jj, getarr
	real*4 parray(MAX_PAR)

        integer same, n0, no_0_par, ncftot, no_cf_par, no_r
        logical is_0, is_cf
        character*3 ic

	n = getari (line, l, iarray, MAX_PAR)
	if (n .eq. 0) then
		call prompt ('Give parameters to be bound :')
		line_len = getlin (line)
		if (line_len .gt. 0) then
			l = 1
			n = getari (line(1:line_len), l, iarray, MAX_PAR)
		end if
	end if
	call remark ('Give the binding parameter no. and the ratio')
        n0     = no_0_par()
        ncftot = no_cf_par(symmetry,nre)
	do i = 1,n
		j = iarray(i)
		call prompt (nam(j)//': ')
		line_len = getlin (line)
		if (line_len .gt. 0) then
			l = 1			  
			k = getarr (line(1:line_len), l, parray, MAX_PAR)

                        is_0  = .FALSE.
                        is_cf = .FALSE.
                        jj =  nint(parray(1))
                        if(j.ge.1   .and.j.le.n0-6  .and. jj.ge.1   .and.jj.le.n0-6 ) is_0  = .TRUE.
                        if(j.gt.n0-6.and.j.le.nptot .and. jj.gt.n0-6.and.jj.le.nptot) is_cf = .TRUE.

                        if( same(j,nint(parray(1))).eq.1 .and. (is_0.or.is_cf).and.
     *                      (   
     *                        (
     *                           setxr.eq.1 .and. j.ne.indvold(1) .and. jj.ne.indvold(1)
     *                        )
     *                        .or.(
     *                              setwx.eq.1 .and. j.eq.n0+1 
     *                            )
     *                      )  
     *                      .or. setbkq.eq.1 .or. setvkq.eq.1 .or. setakq.eq.1
     *                    ) then
            	                if (k .gt. 0) inr(j) = nint(parray(1))
       	     		        if (k .gt. 1) then
			        	rp(j) = parray(2)
			         	p(j) = p(inr(j)) * rp(j)
		         	else if (abs(p(inr(j))) .ne. 0) then
		        		rp(j) = p(j) / p(inr(j))
		        	else
			       	        call remark ('Specify ratio - binding parameter is currently zero')
				        iarray(j) = 0
			        end if
                        else
                                call remark ('ERROR: Impossible to bind the parameters. They are of different kind.')
                                is_0  = .FALSE.
                                is_cf = .FALSE.
                                if(j.ge.1   .and.j.le.n0-6)  is_0  = .TRUE.
                                if(j.gt.n0-6.and.j.le.nptot) is_cf = .TRUE.
                                if(is_0 ) then 
                                   call remark ('You can bind only to a parameter number lower than '//ic(n0-5))
                                   goto 99     
                                else if(is_cf) then  
                                   if(setxr.eq.1 .and. same(j,indvold(1)).eq.1 .and. j.ne.indvold(1)) then
                                      call remark ('You can bind only x-parmeters to x-parameters.')
                                      goto 99
                                   end if
                                   if(setxr.eq.1 .and. same(j,indvold(1)).eq.1 .and. j.eq.indvold(1)) then
                                      call remark ('The R parameter can not be bound to anything else.')
                                      goto 99
                                   end if
                                   if(setwx.eq.1 .and. j.eq.n0+2) then
                                      call remark ('You can not bind the x-parmeter to anything else.')
                                      goto 99
                                   end if
                                   if(setwx.eq.1 .and. j.eq.n0+1) then
                                      call remark ('You can not bind the W-parameter to anything else.')
                                      goto 99
                                   end if
                                   if(setxr.eq.1) then
                                     no_r = indvold(1)  
                                   else if(setwx.eq.1) then
                                     no_r = n0+1        
                                   else
                                     no_r = n0 + ncftot +1
                                   end if
                                   call remark ('You can bind only to a parameter number greater than '//ic(n0-6))   
                                   call remark ('and lower than '//ic(no_r))   
 99                             end if
                        end if
		end if
	end do
	if (n .gt. 0) call selout (nptot, np, n, inp, iarray)

	return
	end
c       -----------------------------------------------------------------------------------
c       finds out if both parameter are of the same kind
c       -----------------------------------------------------------------------------------
        integer function same(p1,p2)
            implicit none
	include 'FRILLS_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:XR.INC' 
            integer p1,p2,i,j,flag1,flag2
            flag1=0
            flag2=0
            do i=1, nindvold
               if(p1.eq.indvold(i)) flag1=1
               if(p2.eq.indvold(i)) flag2=1
            end do
            same = 0
            if( flag1.eq.0 .and. flag2.eq.0 ) same = 1
            if( flag1.eq.1 .and. flag2.eq.1 ) same = 1
	return
	end
