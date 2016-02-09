!	Subroutine DERMAT	-	Calculates matrix of parameter derivatives
!
!  Modified:
!    19/9/94	T.G.Perring	FUN now called through the subroutine FUNCTION_CALL
!

	subroutine DERMAT (fun)

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:LIMITPAR.INC'

	external fun
	integer*4 i, j, k, l, m, n
	real*8 pold, pp, pm, cp(MAX_DAT), cm(MAX_DAT)	 

	do i = 1,np
		j = inp(i)
		pold = p(j)
		if (abs(p(j)) .gt. 1.0D-10) then
			pp = p(j) * (1.0 + dp)
			pm = p(j) * (1.0 - dp)
                        if(pmax(j).gt.pmin(j)) then           ! added by P.Fabi
                           if(pp.gt.up_limit)  pp = up_limit  !
                           if(pm.gt.up_limit)  pm = up_limit  !
                           if(pp.lt.low_limit) pp = low_limit !
                           if(pm.lt.low_limit) pm = low_limit !
                        end if                                !
		else
			pp =  0.0001
			pm = -0.0001
		end if
                if(pp.eq.pm) then
                   pm = -pm
		   pp = max(pp,pm)
	 	   pm = min(pp,pm)
		end if
		if(abs(pp)+abs(pm).eq.0.0d0) then
		   pp =  0.0001
		   pm = -0.0001
		end if
		k = 0
		do l = 1,ns
			p(j) = pp
			call pback (np, nptot, p, inp, pmin, pmax, inr, rp)
			call function_call (l, nv(l), inv(1,l), x(1,l), fun, cp, p)
!			call fun (l, nv(l), inv(1,l), x(1,l), cp, p)
			call ptrans (np, p, inp, pmin, pmax)
			p(j) = pm
			call pback (np, nptot, p, inp, pmin, pmax, inr, rp)
			call function_call (l, nv(l), inv(1,l), x(1,l), fun, cm, p)
!			call fun (l, nv(l), inv(1,l), x(1,l), cm, p)
			call ptrans (np, p, inp, pmin, pmax)
			do m = 1,nv(l)
				k = k + 1
				n = inv(m,l)
				dv(i,k) = (cp(n) - cm(n)) / ((pp - pm) * vsig(n,l))
			end do
		end do
		p(j) = pold
	end do

	return
	end
