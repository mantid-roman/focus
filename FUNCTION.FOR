!
! module FUNCTION.FOR
!                                                                        
!-----------------------------------------------------
! delete all the limits/restrictions which do not belong
! to the symmetry
!-----------------------------------------------------
      subroutine rkq(k,q,s)
	implicit none
        integer k,q,s
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:LIMIT.INC'	
           pmikq( k,q) = cf_pmin(s)
           pmakq( k,q) = cf_pmax(s)
           pmiokq(k,q) = cf_pminold(s)
           pmaokq(k,q) = cf_pmaxold(s)
           indokq(k,q) = cf_indvold(s)
           fixokq(k,q) = cf_fixold(s)
           knpkq( k,q) = cf_knp(s)
           knrkq( k,q) = cf_knr(s)
           kkq(  k,q)  = cf_k(s)   
      return
      end

      subroutine rkqc(k,q,s, c)
	implicit none
        integer k,q,s, c
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:LIMIT.INC'	
        if(cf_pmin(s)*c.le.cf_pmax(s)*c) then
           pmikq(k,q) = cf_pmin(s)*c
           pmakq(k,q) = cf_pmax(s)*c
        else
           pmikq(k,q) = cf_pmax(s)*c
           pmakq(k,q) = cf_pmin(s)*c
	end if
        if(cf_pminold(s)*c.le.cf_pmaxold(s)*c) then
           pmiokq(k,q) = cf_pminold(s)*c
           pmaokq(k,q) = cf_pmaxold(s)*c
        else
           pmiokq(k,q) = cf_pmaxold(s)*c
           pmaokq(k,q) = cf_pminold(s)*c
	end if
        indokq(k,q) = cmplx( 0.0d0, 0.0d0)
        fixokq(k,q) = cmplx(-1.0d0,-1.0d0)
        knpkq( k,q) = cmplx( 0.0d0, 0.0d0)
        knrkq( k,q) = cmplx( 0.0d0, 0.0d0)
        rpkq(  k,q) = cmplx( 0.0d0, 0.0d0)
        kkq(  k,q)  = cmplx( 0.0d0, 0.0d0)
      return
      end

      subroutine ckq(k,q,s,t)
	implicit none
        integer k,q,s,t
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:LIMIT.INC'	
        pmikq( k,q) = cmplx(cf_pmin(s)   ,cf_pmin(t))
        pmakq( k,q) = cmplx(cf_pmax(s)   ,cf_pmax(t))
        pmiokq(k,q) = cmplx(cf_pminold(s),cf_pminold(t))
        pmaokq(k,q) = cmplx(cf_pmaxold(s),cf_pmaxold(t))
        indokq(k,q) = cmplx(cf_indvold(s),cf_indvold(t))
        fixokq(k,q) = cmplx(cf_fixold(s) ,cf_fixold(t))
        knpkq( k,q) = cmplx(cf_knp(s)    ,cf_knp(t))
        knrkq( k,q) = cmplx(cf_knr(s)    ,cf_knr(t))
        rpkq( k,q)  = cmplx(cf_rp(s)     ,cf_rp(t))
        kkq(  k,q)  = cmplx(cf_k(s)      ,cf_k(t))
      return
      end

      subroutine vr(s, k,q)
	implicit none
	integer s, k,q
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:LIMIT.INC'	
         cf_pmin(s)    = real( pmikq( k,q) )
         cf_pmax(s)    = real( pmakq( k,q) )
         cf_pminold(s) = real( pmiokq(k,q) )
         cf_pmaxold(s) = real( pmaokq(k,q) )
         cf_indvold(s) = real( indokq(k,q) )
         cf_fixold(s)  = real( fixokq(k,q) )
         cf_knp(s)     = real( knpkq( k,q) )
         cf_knr(s)     = real( knrkq( k,q) )
         cf_rp(s)      = real( rpkq(  k,q) )
         cf_k(s)       = real( kkq(   k,q) )
      return
      end

      subroutine vi(s, k,q)
	implicit none
	integer s, k,q
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:LIMIT.INC'	
         cf_pmin(s)    = dimag( pmikq( k,q) )
         cf_pmax(s)    = dimag( pmakq( k,q) )
         cf_pminold(s) = dimag( pmiokq(k,q) )
         cf_pmaxold(s) = dimag( pmaokq(k,q) )
         cf_indvold(s) = dimag( indokq(k,q) )
         cf_fixold(s)  = dimag( fixokq(k,q) )
         cf_knp(s)     = dimag( knpkq( k,q) )
         cf_knr(s)     = dimag( knrkq( k,q) )
         cf_rp(s)      = dimag( rpkq(  k,q) )
         cf_k(s)       = dimag( kkq(   k,q) )
      return
      end
!-----------------------------------------------------
      subroutine set0(k,q)
	implicit none
        integer k,q
	complex*16  null
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:LIMIT.INC'	
        null = cmplx(0.0d0,0.0d0)

        pmikq( k,q) = null
        pmakq( k,q) = null
        pmiokq(k,q) = null
        pmaokq(k,q) = null
        indokq(k,q) = null
        fixokq(k,q) = cmplx(-1.0d0,-1.0d0)   ! -1 means free
        knpkq( k,q) = null
        knrkq( k,q) = null
        rpkq(  k,q) = null
        kkq(   k,q) = null
      return
      end 
!    ---------------------------------------------------------------------------
!    initialisation for subroutine lim_trans
!    ---------------------------------------------------------------------------
        subroutine lim_init( symmetry, nre)
	    integer symmetry, nre

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:LIMIT.INC'

            integer k, ncftot, n0, no_0_par, no_cf_par
	    integer ik

c               ---------------------------------------------------------------
c               store the old situation: pmin, pmax, inp and inr
c               prepare a storage for indvold, fixold, pminold,pmaxold,
c               ---------------------------------------------------------------
      		n0     = no_0_par()                        
	        ncftot = no_cf_par(symmetry,nre)                                  

                do k=1,ncftot  ! store the old limits and restriction situation
                   cf_pmin(   k) = pmin(n0+k)
                   cf_pmax(   k) = pmax(n0+k)
		   cf_knp(    k) = 0
                   cf_knr(    k) = 0 
                   cf_rp(     k) = 0 
                   cf_k(      k) = 0 
	           cf_indvold(k) = 0
	           cf_fixold( k) = -1     ! by default free
	    	   cf_pminold(k) = 0.0d0
		   cf_pmaxold(k) = 0.0d0
		end do

                do k=1,6
		   mf_knp(    k) = 0
                   mf_knr(    k) = 0
                   mf_rp(     k) = 0
                   mf_k(      k) = 0
      	           mf_indvold(k) = 0
	           mf_fixold( k) = -1     ! by default free
	    	   mf_pminold(k) = 0.0d0
		   mf_pmaxold(k) = 0.0d0
		end do
c               ---------------------------------------------------------------
c               separate in magnetic and crystal field parameters
c               ---------------------------------------------------------------
                do k=1, nindvold
		   ik = indvold(k) 
                   if(ik.gt.n0 .and. ik.le.n0+ncftot) then
	              cf_indvold(ik-n0) = k
		      cf_fixold( ik-n0) = fixold(ik)/ik
		      cf_pminold(ik-n0) = pminold(ik)
		      cf_pmaxold(ik-n0) = pmaxold(ik)
                   else if(ik.gt.n0-6 .and. ik.le.n0) then
	              mf_indvold(ik-(n0-6)) = k
		      mf_fixold( ik-(n0-6)) = fixold(ik)/ik
		      mf_pminold(ik-(n0-6)) = pminold(ik)
		      mf_pmaxold(ik-(n0-6)) = pmaxold(ik)
		   end if
                end do

                do k=1, np
                   if(inp(k).gt.n0 .and. inp(k).le.n0+ncftot) then
                      cf_knp(inp(k)-n0) = k
                   else if(inp(k).gt.n0-6 .and. inp(k).le.n0) then
                      mf_knp(inp(k)-(n0-6)) = k
                   end if
                end do

                do k=n0-5, n0+ncftot 
                   if(inr(k).gt.n0 .and. inr(k).le.n0+ncftot) then
                      cf_knr(inr(k)-n0) = k
                      cf_rp( inr(k)-n0) = rp(k)
                   else if(inr(k).gt.n0-6 .and. inr(k).le.n0) then
                      mf_knr(inr(k)-(n0-6)) = k
                      cf_rp( inr(k)-(n0-6)) = rp(k)
                   end if
                end do

                do k=n0-5, n0+ncftot 
                   if(k.gt.n0 .and. k.le.n0+ncftot) then
                      cf_k( k-n0)    = k
                   else if(k.gt.n0-6 .and. k.le.n0) then
                      mf_k(k-(n0-6)) = k
                   end if
                end do

	return
	end
!    ---------------------------------------------------------------------------
!    initialisation for subroutine ilim_trans
!    ---------------------------------------------------------------------------
        subroutine ilim_init( symmetry, nre)
	    integer symmetry, nre
            integer k,i, ncftot, n0, no_0_par, no_cf_par, nn, kk, nk
	    real*8  dummy

	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:LIMIT.INC'

		integer n_nindvold, kcf, kmf, kold
		integer n_indvold(MAX_PAR), n_fixold( MAX_PAR)
		integer n_inp(MAX_PAR)    
		real*8  n_pminold(MAX_PAR), n_pmaxold(MAX_PAR)

      		n0     = no_0_par()                        
	        ncftot = no_cf_par(symmetry,nre)                                  
                nptot  = n0+ncftot
c               ---------------------------------------------------------------
c               store the new situation: inp and inr
c               ---------------------------------------------------------------
                do k=1, np
                   inp(k) = 0
                   do kcf=1,ncftot
		      if(cf_knp(kcf).eq.k .and. n0+kcf.le.nptot) inp(k) = kcf+n0
                   end do
		   do kmf=1,6
		      if(mf_knp(kmf).eq.k) inp(k) = kmf+n0-6 
		   end do
                end do

                do k=n0-5, nptot 
                   kk=0
                   do i=1,6
                      if(mf_k(i).eq.k) kk = i+n0-6
                   end do
                   do i=1,ncftot
                      if(cf_k(i).eq.k) kk = i+n0
                   end do
                   if(kk.ne.0) then
                      inr(kk)=0
                      do kcf=1,ncftot
      		         if(cf_knr(kcf).eq.k .and. n0+kcf.le.nptot) then 
                            inr(kk) = kcf+n0
                            rp(kk)  = cf_rp(kcf)
                         end if
                      end do
                      do kmf=1,6
		         if(mf_knr(kmf).eq.k ) then
                            inr(kk) = kmf+n0-6 
                            rp(kk)  = mf_rp(kmf)
                         end if
                      end do
                   end if
                end do
c               ---------------------------------------------------------------
c               compress new inp array -> determine new np
c               np = number of free parameter
c               ---------------------------------------------------------------
                nk = np
                np = 0
                do k=1,nk
                   if(inp(k).ne.0) then
                      np=np+1
                      n_inp(np) = inp(k)
                   end if
                end do
                do k=1,np
                   inp(k) = n_inp(k)
                end do
c               ---------------------------------------------------------------
c               store the new situation: pmin, pmax
c               ---------------------------------------------------------------
                do k=1,ncftot
		     pmin(n0+k) = cf_pmin(k) 
                     pmax(n0+k) = cf_pmax(k)
                     if(pmax(n0+k).lt.pmin(n0+k)) then
                        dummy      =  pmin(n0+k)
			pmin(n0+k) =  pmax(n0+k)
			pmax(n0+k) =  dummy
		     end if
		end do

c               ---------------------------------------------------------------
c               delete all arrays
c               ---------------------------------------------------------------
                do k=1, MAX_PAR
                   n_indvold(k) = 0
                   n_fixold(k)  = -1     ! -1 means free
		   n_pminold(k) = 0.0d0
		   n_pmaxold(k) = 0.0d0
                end do
c               ---------------------------------------------------------------
c               store the new situation in the old order
c               ---------------------------------------------------------------
                n_nindvold = 0
                do kold=1, nindvold
                   do kcf=1, ncftot
                      if(kold.eq.cf_indvold(kcf)) then
                         n_nindvold = n_nindvold + 1
			 n_indvold(n_nindvold) = n0+kcf
			 n_fixold( n0+kcf    ) = cf_fixold( kcf) * (n0+kcf)
			 n_pminold(n0+kcf    ) = cf_pminold(kcf)
			 n_pmaxold(n0+kcf    ) = cf_pmaxold(kcf)
                         cf_indvold(kcf) = - cf_indvold(kcf)  ! mark that this value is taken
                         goto 99
		      end if
                   end do

                   do kmf=1, 6
                      if(kold.eq.mf_indvold(kmf)) then
                         n_nindvold = n_nindvold +1
			 n_indvold(n_nindvold) = n0-6 + kmf
			 n_fixold( n0-6+kmf) = mf_fixold( kmf) * (n0-6+kmf)
			 n_pminold(n0-6+kmf) = mf_pminold(kmf)
			 n_pmaxold(n0-6+kmf) = mf_pmaxold(kmf)
                         mf_indvold(kmf) = - mf_indvold(kmf)  ! mark that this value is taken
                         goto 99
		      end if
		   end do
  99		end do
c               ---------------------------------------------------------------
c               store next crystal fields
c               ---------------------------------------------------------------
                do kcf=1,ncftot
                   if(cf_indvold(kcf).gt.0) then
                         n_nindvold = n_nindvold +1
			 n_indvold(n_nindvold) = n0+kcf
			 n_fixold( n0+kcf) = cf_fixold( kcf) * (n0+kcf)
			 n_pminold(n0+kcf) = cf_pminold(kcf)
			 n_pmaxold(n0+kcf) = cf_pmaxold(kcf)
                   end if
                end do
c               ---------------------------------------------------------------
c               store next molecular fields
c               ---------------------------------------------------------------
                do kmf=1,6
                   if(mf_indvold(kmf).gt.0) then
                         n_nindvold = n_nindvold +1
			 n_indvold(n_nindvold) = n0-6+kmf
			 n_fixold( n0-6+kmf) = mf_fixold( kmf) * (n0-6+kmf)
			 n_pminold(n0-6+kmf) = mf_pminold(kmf)
			 n_pmaxold(n0-6+kmf) = mf_pmaxold(kmf)
		   end if
                end do
c               ---------------------------------------------------------------
c               now store the new situation
c               ---------------------------------------------------------------
                nindvold = n_nindvold
                do k=1,nindvold
                   indvold(k)  = n_indvold(k)
	           nn          = indvold(  k)
                   fixold( nn) = n_fixold( nn)
                   pminold(nn) = n_pminold(nn)
                   pmaxold(nn) = n_pmaxold(nn)
                   if(pmaxold(nn).lt.pminold(nn)) then
                      dummy       =  pminold(nn)
		      pminold(nn) =  pmaxold(nn)
	              pmaxold(nn) =  dummy
		   end if
                end do

	return
	end
!    ---------------------------------------------------------------------------
!    set a new R3+ ion
!    ---------------------------------------------------------------------------
        subroutine r_init(new_symmetry, new_nre)
		integer new_symmetry, new_nre
                integer out      ! messages a written to file out
                integer in       ! reading is from file in
                parameter(in=5,out=5)

	include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:CONTROL.INC'
	include 'CF_SOURCES:XR.INC'
	include 'CF_SOURCES:LIMIT.INC'

	       	integer n0, no_0_par, no_cf_par, ncftot, k, kk, nn, new_nindvold		
                integer isbkq,isvkq,isxr,isakq,iswx
                real*8 pv(MAX_PAR)

		integer old_symmetry , old_nre
!               ---------------------------------------------------------------------------
!               store the old situation
!               ---------------------------------------------------------------------------
		old_nre      = nre
                old_symmetry = symmetry
		isakq = setakq
		isbkq = setbkq
                isvkq = setvkq
		isxr  = setxr
		iswx  = setwx

!               ---------------------------------------------------------------------------
!               alter in the old scheme to Bkq-parameter first
!               ---------------------------------------------------------------------------
		nre      = old_nre        ! to be stored on common block in FOCUS.INC
		symmetry = old_symmetry   !
      		n0     = no_0_par()                        
	        ncftot = no_cf_par(old_symmetry,old_nre)                                  
                nptot  = n0 + ncftot
                call a_to_akq(setbkq,setvkq,setxr,setakq,setwx)     
                if(setfirst.eq.0) then
                   call remark (' FOCUS transformed the actual parameter set to the Akq-parameter.')
		end if

!               ---------------------------------------------------------------------------
!               set all old bkq-parameters to 0 which do not belong to the new symmetry
!               ---------------------------------------------------------------------------
                do k=1,ncftot
                   pv(k) = p(n0+k)
                end do
                call transform(pv,bkq,old_symmetry,old_nre)


                call lim_init ( old_symmetry, old_nre )
                call lim_trans( old_symmetry, old_nre )

                call del_bkq(bkq,new_symmetry,new_nre)  ! adaptes CEF parameters to new R3+ & symmetry
                call del_lim(new_symmetry, new_nre )    ! adaptes limits to new R3+ & symmetry

                if(setfirst.eq.0 .and. (old_symmetry.ne.new_symmetry .or. old_nre.ne.new_nre)) then
                   if(new_symmetry.ne.old_symmetry) call remark (' FOCUS adapted the CEF parameters to the new symmetry.')
                   if(new_nre.ne.old_nre) call remark (' FOCUS adapted the CEF parameters to the new rare earth ion.')
		end if
!               ---------------------------------------------------------------------------
!               transform now back to the new sheme
!               ---------------------------------------------------------------------------
		nre      = new_nre        ! to be stored on common block in FOCUS.INC
		symmetry = new_symmetry   !
                call itransform(pv,bkq,new_symmetry,new_nre)

	        ncftot = no_cf_par(new_symmetry,new_nre)                                  
                nptot  = n0 + ncftot
                do k=1,ncftot
                   p(n0+k) = pv(k)
                end do

		call ilim_trans( new_symmetry, new_nre )
                call ilim_init ( new_symmetry, new_nre )

!               ---------------------------------------------------------------------------
!               define which part of the operator (real or imaginary) to be chosen
!               ---------------------------------------------------------------------------
                if(new_nre.ne.6 .and. new_nre.ne.7 .and. setfirst.eq.1 ) then   
 	                call p_hamiltonian(out,new_symmetry)
                end if
                call s_sbkq(in,out,new_nre,new_symmetry,sbkq)

!               ---------------------------------------------------------------------------
!               define the parameter names for each symmetry
!               ---------------------------------------------------------------------------
                if(new_nre.ne.6 .and. new_nre.ne.7) then   
                   call d_pnames(new_nre,new_symmetry,sbkq,nam,n0)
                end if

!               ---------------------------------------------------------------------------
!               transform back to the original parameter scheme
!               ---------------------------------------------------------------------------
                if(        isakq.eq.1) then
                   call a_to_akq(setbkq,setvkq,setxr,setakq,setwx) ! needed to get nam(i) right in # D P    
                else if(   isvkq.eq.1) then
                   call a_to_vkq(setbkq,setvkq,setxr,setakq,setwx)     
		else if(   isxr.eq.1) then
                   call a_to_xr(setbkq,setvkq,setxr,setakq,setwx)     
		else if(   isbkq.eq.1) then
                   call a_to_bkq(setbkq,setvkq,setxr,setakq,setwx)     
		else if(   iswx.eq.1) then
                   call a_to_wx(setbkq,setvkq,setxr,setakq,setwx)     
                end if

                if(setfirst.eq.0 .and. (old_symmetry.ne.new_symmetry .or. old_nre.ne.new_nre)) then
                   call remark (' FOCUS transformed the Akq-parameter to the actual parameter scheme.')
      	  	end if

	return
	end
!    ---------------------------------------------------------------------------
!    set the overall splitting
!    ---------------------------------------------------------------------------
        subroutine a_overall_splitting
	   include 'CF_SOURCES:FRILLS.INC'
	   include 'CF_SOURCES:FOCUS.INC'
	   include 'CF_SOURCES:OVERALL.INC'
	   include 'CF_SOURCES:CONTROL.INC'
	   include 'CF_SOURCES:XR.INC'
        
           character line*80, c  
           character*64 text1,text2,text3, text4, text5
           integer*4 l, line_len, getlin, ctoi, lneu, idimj
           real*4    ctor, r

           integer  n0, ncftot, no_cf_par, no_0_par, same, i,j

           data text1 /' Should FOCUS fix a CEF energy level to a given value...... [0]?'/
           data text2 /' Enter new level number (e.g. 2=second,3=third,.. )........[00]?'/
           data text3 /' Enter new energy value in meV......................[000000.00]:'/

           data text4 /' ERROR: The level number must be between 2 and 00 !             '/
           data text5 /'        The ground state has level number 1.                    '/

c          ---------------------------------------------------------------------
c          initialize some rare earth constants
c          ---------------------------------------------------------------------
	   call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)           
           idimj = dimj

           write(*,*)
           call remark ( ' Press the <enter>-key only if you do not want to change the value.' )
           if(setoverall.eq.1) then
              text1(62:62)='Y'
           else
              text1(62:62)='N'
           end if

  1        call prompt (text1)
           line_len = getlin(line)
           c = 'N'
           if(setoverall.eq.1) c = 'Y' 
           if( line_len .gt. 0 ) then
               c = line(1:1)
           end if
           if(c.ne.'Y'.and.c.ne.'y'.and.c.ne.'n'.and.c.ne.'N') goto 1
           if(c.eq.'y'.or.c.eq.'Y') setoverall = 1
           if(c.eq.'n'.or.c.eq.'N') setoverall = 0

           if(setoverall.eq.1) then

              setpardis = 0

              if(levelnumber.eq.0) levelnumber=2
              if(levelenergy.eq.0) levelenergy=1
              write(text2(61:62),fmt='(i2)'  ) levelnumber
              write(text3(54:62),fmt='(f9.2)') levelenergy


  2           call prompt (text2)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l   = 1
                 lneu = ctoi(line(1:line_len),l)
                 if(lneu.lt.2.or.lneu.gt.idimj) then
                    write(text4(48:49),fmt='(i2)'  ) idimj
                    call remark(text4)
                    call remark(text5)
                    goto 2
                 end if
                 levelnumber = lneu
              end if
            

  3           call prompt (text3)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.eq.0.0) then
                    call remark ( ' ERROR: the excited energy can not be zero! ')
                    goto 3
                 end if
                 levelenergy = abs(r)
              end if

        else
!            ---------------------
!            release new parameter
!            ---------------------
             n0     = no_0_par()
             ncftot = no_cf_par(symmetry,nre)                                  
             do j=n0-5,n0+ncftot
                if(  setxr.eq.1 ) then
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

        end if

        return
        end
!    ---------------------------------------------------------------------------
        subroutine a_selfshielding(setself)

	   integer setself

	include 'CF_SOURCES:SELF.INC'

           character line*80, c  
           character*64 text1,text2,text3
           integer*4 l, line_len, getlin, ctoi
           real*4    ctor, r

           data text1 /' Calculate the selfshielding and absorption correction..... [0]?'/
           data text2 /' Enter new value for   mu*t (scattering)............[000000.00]:'/
           data text3 /' Enter new value for   mu*t (absorption)............[000000.00]:'/

           write(*,*)
           call remark ( ' Press the <enter>-key only if you do not want to change the value.' )
           call remark ( ' Press  @ I M  to inform you about the meaning of mu*t.')
           if(setself.eq.1) then
              text1(62:62)='Y'
           else
              text1(62:62)='N'
           end if

  1        call prompt (text1)
           line_len = getlin(line)
           c = 'N'
           if(setself.eq.1) c = 'Y' 
           if( line_len .gt. 0 ) then
               c = line(1:1)
           end if
           if(c.ne.'Y'.and.c.ne.'y'.and.c.ne.'n'.and.c.ne.'N') goto 1
           if(c.eq.'y'.or.c.eq.'Y') setself = 1
           if(c.eq.'n'.or.c.eq.'N') setself = 0

           if(setself.eq.1) then

              write(text2(54:62),fmt='(f9.2)') muts
              write(text3(54:62),fmt='(f9.2)') muta


  2           call prompt (text2)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.lt.0.0) then
                    call remark ( ' ERROR: mu*t can not be negativ! ')
                    goto 2
                 end if
                 muts = r
              end if

  3           call prompt (text3)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.lt.0.0) then
                    call remark ( ' ERROR: mu*t can not be negativ! ')
                    goto 3
                 end if
                 muta = r
              end if

           end if
	return
	end
!    ---------------------------------------------------------------------------
!    sophisticated setup of the resolution function values
!    ---------------------------------------------------------------------------
        subroutine set_resolution(setres) 
           implicit none
           integer setres

           integer out      ! messages a written to file out
           integer in       ! reading is from file in
           parameter(in=5,out=5)

        include 'CF_SOURCES:INSTRUMENT.INC'
        include 'CF_SOURCES:DEGENERATION.INC'
        include 'CF_SOURCES:RESOLUTION.INC'

           character line*80, c,c6*6,c7*7,c5*5  
           character*64 text1,text2,text3,text4,text5,text6,text7
           integer*4 l, line_len, getlin, ctoi
           real*4    ctor, r
   	   real*8      c_fmevkelvin  
   	   external    c_fmevkelvin  
!          --------------------------------------------------------------------------
!          for de estimation in case of a resolution calculation
!          --------------------------------------------------------------------------
           real*8 sigma_chop, gamma_elas
!          --------------------------------------------------------------------------

      call remark ( ' FOCUS  has a build in resolution function calculation for a Fermi chopper')
      call remark ( ' time-of-flight neutron scattering instrument. The instrument should be of')
      call remark ( ' direct geometry. The Fermi chopper defines the energy of the incident ')
      call remark ( ' neutrons.')
      call remark ( ' ADVICE: A resolution function  calculation  is time consuming.  Therefore')
      call remark ( '         use this opition only to improve your  final  fit. To do this set')
      call remark ( '         the width of the elastic line to 0, and reduce the line widths of')
      call remark ( '         the quasielastic an the inelatic lines.                          ') 
  30  call prompt ( ' Do you want a resolution function calculation..............[No=default]?')
      line_len = getlin(line)
      c = 'N'
      if( line_len .gt. 0 ) then
          c = line(1:1)
      end if
      if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 30
      setres = 0
      if(c.eq.'Y'.or.c.eq.'y') then
          setres = 1
 999      format(1x,a47,f6.3,a3)
          call remark (' ')
  50      call prompt ( ' Enter the chopper instrument used....H=HET M=MARI O=OTHER...[H=default]?')
          line_len = getlin(line)
          c = 'H'
          if( line_len .gt. 0 ) then
              c = line(1:1)
          end if
          if(c.ne.'H'.and.c.ne.'h'.and.c.ne.'M'.and.c.ne.'m'.and.c.ne.'O'.and.c.ne.'o') goto 50
          if(c.eq.'H'.or.c.eq.'h') then
             lmc = lmc_het
             lcs = lcs_het
             lsd = lsd_het
             call remark (   ' HET default lengths chosen' )
             write (out,999) ' FOCUS sets the moderator-chopper distance to ',lmc,' m.'
             write (out,999) ' FOCUS sets the chopper-sample    distance to ',lcs,' m.'
             write (out,999) ' FOCUS sets the sample-detector   distance to ',lsd,' m.'
             call remark (' ')
             write(c5,fmt='(f5.3)') lsd_het_h
  51         call prompt (   ' Do you want to change the sample-detector distance to '//c5//'m [N=default]?')
             line_len = getlin(line)
             c = 'N'
             if( line_len .gt. 0 ) then
                 c = line(1:1)
             end if
             if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 51
             if(c.eq.'Y'.or.c.eq.'y') then
                lsd = lsd_het_h
                write (out,999) ' FOCUS sets the sample-detector   distance to ',lsd,' m.'
             end if
          else if(c.eq.'M'.or.c.eq.'m') then
            lmc = lmc_mari
            lcs = lcs_mari
            lsd = lsd_mari
            call remark (   ' MARI default lengths chosen' )
            write (out,999) ' FOCUS sets the moderator-chopper distance to ',lmc,' m.'
            write (out,999) ' FOCUS sets the chopper-sample    distance to ',lcs,' m.'
            write (out,999) ' FOCUS sets the sample-detector   distance to ',lsd,' m.'
          else
            lmc = lmc_het
            lcs = lcs_het
            lsd = lsd_het

            call remark (' ')
            write(c6,fmt='(f6.3)') lmc_het
  52        call prompt ( ' Enter the moderator-chopper distance in meter.....['//c6//']:')
            line_len = getlin(line)
            if(line_len .gt. 0) then
               l    = 1
               lmc = ctor(line(1:line_len),l)
               if(lmc.le.0.0d0) then
                  call remark ( ' ERROR: The distance must be positiv! ')
                  goto 52
               end if
            end if
            write(c5,fmt='(f5.3)') lcs_het
  53        call prompt ( ' Enter the chopper-sample    distance in meter......['//c5//']:')
            line_len = getlin(line)
            if(line_len .gt. 0) then
               l    = 1
               lcs = ctor(line(1:line_len),l)
               if(lcs.le.0.0d0) then
                  call remark ( ' ERROR: The distance must be positiv! ')
                  goto 53
               end if
            end if
            write(c5,fmt='(f5.3)') lsd_het
  54        call prompt ( ' Enter the sample-detector   distance in meter......['//c5//']:')
            line_len = getlin(line)
            if(line_len .gt. 0) then
               l    = 1
               lsd = ctor(line(1:line_len),l)
               if(lsd.le.0.0d0) then
                  call remark ( ' ERROR: The distance must be positiv! ')
                  goto 54
               end if
            end if
          end if

          write(*,*)
          call remark ( ' With a genie program called @G:RESFIT the elastic line width (FWHM)')
          call remark ( ' and the decay constant of monochromatic Vanadium can be determined.')  
          call remark ( ' The  monochromatic Vanadium  has to be measured under same chopper') 
          call remark ( ' conditions (i.e. incoming energy, chopper frequency) as the sample') 
          call remark ( ' measurement.')
  55      call prompt ( ' Do you have these data available......................[Y=default]?')
          line_len = getlin(line)
          c = 'Y'
          if( line_len .gt. 0 ) then
              c = line(1:1)
          end if
          if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 55
          if(c.eq.'N'.or.c.eq.'n') then
               setres = 0
               goto 9999
          end if
          call remark (' ')
  56      call prompt ( ' Enter the elastic width (FWHM) of Vanadium in meV.....[0=default]?')
          line_len = getlin(line)
          fwhm_van_el = 0.0d0
          if(line_len .gt. 0) then
             l    = 1
             fwhm_van_el = ctor(line(1:line_len),l)
             if(fwhm_van_el.lt.0.0d0) then
                call remark ( ' ERROR: a line width can not be negativ! ')
                goto 56
             end if
          end if
          call remark (' ')
  57      call prompt ( ' Enter the decay constant in meV of Vanadium...........[0=default]?')
          line_len = getlin(line)
          decay = 0.0d0
          if(line_len .gt. 0) then
             l    = 1
             decay = ctor(line(1:line_len),l)
             if(decay.lt.0.0d0) then
                call remark ( ' ERROR: The decay constant can not be negativ! ')
                goto 57
             end if
          end if
          write(*,*)
          call remark ( ' With a program called CHOP the chopper contribution to the elastic ')
          call remark ( ' Vanadium line width can be calculated.                             ')  
  58      call prompt ( ' Do you have this value  available.....................[Y=default]?')
          line_len = getlin(line)
          c = 'Y'
          if( line_len .gt. 0 ) then
              c = line(1:1)
          end if
          if(c.ne.'N'.and.c.ne.'n'.and.c.ne.'Y'.and.c.ne.'y') goto 58
          if(c.eq.'N'.or.c.eq.'n') then
               setres = 0
               goto 9999
          end if
          call remark (' ')
  59      call prompt ( ' Enter the chopper contribution (FWHM) in meV..........[0=default]?')
          line_len = getlin(line)
          fwhm_chop_el = 0.0d0
          if(line_len .gt. 0) then
             l    = 1
             fwhm_chop_el = ctor(line(1:line_len),l)
             if(fwhm_chop_el.lt.0.0d0) then
                call remark ( ' ERROR: The chopper contribution to the Vanadium width ')
                call remark ( '        can not be negativ! ')
                goto 59
             end if
             if(fwhm_chop_el.gt.fwhm_van_el) then
                 write(c6,fmt='(f6.2)') fwhm_van_el
                 call remark ( ' ERROR: The chopper contribution to the Vanadium width ')
                 call remark ( '        can not be greater than '//c6//' meV.')
                goto 59
             end if
          end if

          sigma_chop = fwhm_chop_el / sqrt(log(256.0d0))
          gamma_elas = decay        /2.0d0
          de = max(sigma_chop/5.0d0, gamma_elas/5.0d0)
          if(de.le.0.0d0) de = de_default
          de = de * c_fmevkelvin()

          write(c7,fmt='(f7.4)') de
          call remark (' FOCUS sets the energy resolution to de =  '//c7//' Kelvin.')

 9999 end if

      if(setres.eq.0) then
          de = de_default * c_fmevkelvin()
          write(c7,fmt='(f7.4)') de
          call remark (' FOCUS sets the energy resolution to de =  '//c7//' Kelvin.')

      end if
        return
	end
!    ---------------------------------------------------------------------------
        subroutine a_resolution(setres) 
           
           integer setres

        include 'CF_SOURCES:INSTRUMENT.INC'
        include 'CF_SOURCES:DEGENERATION.INC'
        include 'CF_SOURCES:RESOLUTION.INC'

           character line*80, c,c6*6,c7*7  
           character*64 text1,text2,text3,text4,text5,text6,text7
           integer*4 l, line_len, getlin, ctoi
           real*4    ctor, r
   	   real*8      c_fmevkelvin  
   	   external    c_fmevkelvin  
   
!--------------------------------------------------------------------------
! for de estimation in case of a resolution calculation
!--------------------------------------------------------------------------
           real*8 sigma_chop, gamma_elas
!--------------------------------------------------------------------------

           data text1  /' Calculate the resolution function..........................[0]?'/
           data text2  /' Enter new elastic vanadium line width (FWHM) in meV..[0000.00]:'/
           data text3  /' Enter new chopper contribution (FWHM) in meV.........[0000.00]:'/
           data text4  /' Enter new decay constant (FWHM) in meV...............[0000.00]:'/
           data text5  /' Enter new moderator-chopper distance in meter.......[000.0000]:'/
           data text6  /' Enter new chopper-sample    distance in meter.......[000.0000]:'/
           data text7  /' Enter new sample-detector   distance in meter.......[000.0000]:'/


           write(*,*)
           call remark (' You can have a sophisticated help in setting the resolution ')
           call remark (' function values.')
 10        call prompt (' Do you want a sophisticated help..................[No=default]?')
           line_len = getlin(line)
           c = 'N'
           if( line_len .gt. 0 ) then
               c = line(1:1)
           end if
           if(c.ne.'Y'.and.c.ne.'y'.and.c.ne.'n'.and.c.ne.'N') goto 10
           if(c.eq.'y'.or.c.eq.'Y') then
              call set_resolution(setres)
              return
           end if

           write(*,*)
           call remark ( ' Press the <enter>-key only if you do not want to change the value.' )
           if(setres.eq.1) then
              text1(62:62)='Y'
           else
              text1(62:62)='N'
           end if

  1        call prompt (text1)
           line_len = getlin(line)
           c = 'N'
           if(setres.eq.1) c = 'Y' 
           if( line_len .gt. 0 ) then
               c = line(1:1)
           end if
           if(c.ne.'Y'.and.c.ne.'y'.and.c.ne.'n'.and.c.ne.'N') goto 1
           if(c.eq.'y'.or.c.eq.'Y') setres = 1
           if(c.eq.'n'.or.c.eq.'N') setres = 0

           if(setres.eq.1) then

              write(text2(56:62),fmt='(f7.2)') fwhm_van_el
              write(text3(56:62),fmt='(f7.2)') fwhm_chop_el
              write(text4(56:62),fmt='(f7.2)') decay          

              if(lmc.eq.0.0d0) lmc = lmc_het
              if(lcs.eq.0.0d0) lcs = lcs_het
              if(lsd.eq.0.0d0) lsd = lsd_het

              write(text5(55:62),fmt='(f8.4)') lmc
              write(text6(55:62),fmt='(f8.4)') lcs
              write(text7(55:62),fmt='(f8.4)') lsd

  2           call prompt (text2)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.lt.0.0) then
                    call remark ( ' ERROR: A line width can not be negativ! ')
                    goto 2
                 end if
                 fwhm_van_el = r
              end if

  3           call prompt (text3)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.lt.0.0) then
                    call remark ( ' ERROR: A line width can not be negativ! ')
                    goto 3
                 end if
                 fwhm_chop_el = r
                 if(fwhm_chop_el.gt.fwhm_van_el) then
                     write(c6,fmt='(f6.2)') fwhm_van_el
                     call remark ( ' ERROR: The chopper contribution to the Vanadium width ')
                     call remark ( '        can not be greater than '//c6//' meV.')
                     goto 3
                 end if

              end if

  4           call prompt (text4)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.lt.0.0) then
                    call remark ( ' ERROR: The decay constant can not be negativ! ')
                    goto 4
                 end if
                 decay = r
              end if

  5           call prompt (text5)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.le.0.0) then
                    call remark ( ' ERROR: The  distance must be positiv! ')
                    goto 5
                 end if
                 lmc = r
              end if

  6           call prompt (text6)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.le.0.0) then
                    call remark ( ' ERROR: The  distance must be positiv! ')
                    goto 6
                 end if
                 lcs = r
              end if

  7           call prompt (text7)
              line_len = getlin(line)
              if(line_len .gt. 0) then
                 l = 1
                 r = ctor(line(1:line_len),l)
                 if(r.le.0.0) then
                    call remark ( ' ERROR: The  distance must be positiv! ')
                    goto 7
                 end if
                 lsd = r
              end if

           end if

           if(setres.eq.1) then
          	sigma_chop = fwhm_chop_el / sqrt(log(256.0d0))
          	gamma_elas = decay        /2.0d0
          	de = max(sigma_chop/5.0d0, gamma_elas/5.0d0)
		if (de.le.0.0d0) de = de_default
          	de = de * c_fmevkelvin()
          	write(c7,fmt='(f7.4)') de
          	call remark (' FOCUS sets the energy resolution to de = '//c7//' Kelvin.')
	   else
		de = de_default*c_fmevkelvin()
          	write(c7,fmt='(f7.4)') de
          	call remark (' FOCUS sets the energy resolution to de = '//c7//' Kelvin.')
	   end if
        return
        end
!    ---------------------------------------------------------------------------
!    FOCUS defined commands
!    ---------------------------------------------------------------------------    
 	subroutine SETSHAPE (line, l)
 
 	include 'CF_SOURCES:FRILLS.INC' 
 	include 'CF_SOURCES:FOCUS.INC'
 	include 'CF_SOURCES:DEGENERATION.INC'
 
 	character*5 par_no
 	character*(*) line	
 	integer*4 l, n, getari, iarray(n_lineshape+1), i, j, par_len, line_len, getlin, k, ctoi
 	real*4 parray(3)
 	real*8  c_fmevkelvin
        real*8  e_neg,e_pos
 	character*18 text_line                       
	integer ind(17*17),positiv_start,jj   
        logical first_time

        call energy_sort(e_excitation,ind,n_excitation,positiv_start)
 
 	n = getari (line, l, iarray, n_lineshape+1)
 	if (n .eq. 0) then
                n = (n_excitation-positiv_start+1)+1 
 		do i = 1,n
 		   iarray(i) = i-1
 		end do          
 	else if (n .lt. 0) then
  		call remark ('ERROR in excitation number list')
 		return
 	end if	                                                 
 
   1    format(2x,a6,1x,i4,1x,a32)
   2	format(x,x, i3,' : ',f6.2,' meV',3x,f6.2,' barn',3x,a18,':',$)       
   4	format(x,x, i3,' : ',f16.2,' meV',3x,f6.2,' barn',3x,a18,':',$)       
   3    format(x,x,a21,'            ',a18,':',$)

        first_time = .TRUE.
 	do i = 1,n
 		j = iarray(i)
 		if (j .lt. 0 .or. j .gt. n_excitation-positiv_start+1) then
                        write(*,1) 'ERROR:',j,' is an invalid excitation number'
 		else
                   if(first_time) then
                      first_time = .FALSE. 
 	              call remark ('Give the excitation value. For its definition see: @ D S')
                      call remark ('INPUT BELOW is exspected as follows:')
                      call remark ('   <return> : old value is not changed')
                      call remark (' 1 <return> : to choose GAUSSIAN   shape')
                      call remark (' 2 <return> : to choose LORENTZIAN shape')
                      call remark (' 3 <return> : to choose LORENTZIAN*BOSEFACTOR shape')
                   end if
 		   if(j.eq.0) then    
 	              write(*,3) '  0 :    elastic line',text_line(lineshape(0))           
                   else
                      jj = positiv_start + j  - 1
                      if(e_excitation(ind(jj))/c_fmevkelvin().lt.1.0d3) then 
                         write(*,2) j,e_excitation(ind(jj))/c_fmevkelvin(),i_excitation(ind(jj)), 
     *                              text_line(lineshape(ind(jj)))  
                      else
                         write(*,4) j,e_excitation(ind(jj))/c_fmevkelvin(),i_excitation(ind(jj)),
     *                              text_line(lineshape(ind(jj)))  
                      end if
 		   end if
 	           line_len = getlin (line)
 		   if (line_len .gt. 0) then
 			l = 1
 			k = ctoi (line(1:line_len), l)
                        if(j.eq.0) then
 			  lineshape(0) = k
                        else
 			  lineshape(ind(jj)) = k
                        end if
 		   end if
 		end if
 	end do

! set the line shape of the negative excitations to the value of the positive ones

        if(positiv_start.gt.1) then
           do i=1,positiv_start-1
              e_neg = e_excitation(ind(i))
              do j=positiv_start,n_excitation
                 e_pos = e_excitation(ind(j))
                 if(abs(e_pos+e_neg) .le. de ) then
                    lineshape(ind(i)) = lineshape(ind(j))
                    goto 99
                 end if
              end do
  99       end do
        end if

	return
	end
c       --------------------------------------------------------------------------    
	subroutine SETEULER (line, l)
        include 'CF_SOURCES:CANTED.INC'

        character line*80
        integer*4 l, line_len, getlin, ctoi
        real*4    ctor
        real*8 pi
                                
        call remark ( ' The R3+ moment is canted compared to the crystal system (a,b,c)')
        call remark ( ' by the Euler-angles.')
        call remark ( ' Please give in now the Euler-angles (in units of degrees).')
        call remark ( ' <return>-key only: old value will not be changed.')

  22    call prompt ( ' Euler-angle ALPHA :')
        line_len = getlin(line)
        if(line_len .gt. 0) then
           l   = 1
           alpha_euler = ctor(line(1:line_len),l)
           if(alpha_euler.lt.-360.0d0 .or. alpha_euler.gt.360.0d0) goto 22
        end if

  23    call prompt ( ' Euler-angle BETA  :')
        line_len = getlin(line)
        if(line_len .gt. 0) then
           l   = 1
           beta_euler = ctor(line(1:line_len),l)
           if(beta_euler.lt.-360.0d0 .or. beta_euler.gt.360.0d0) goto 23
        end if

  24    call prompt ( ' Euler-angle GAMMA :')
        line_len = getlin(line)
        if(line_len .gt. 0) then
           l   = 1
           gamma_euler = ctor(line(1:line_len),l)
           if(gamma_euler.lt.-360.0d0 .or. gamma_euler.gt.360.0d0) goto 24
        end if

        pi = 4.0d0*atan(1.0d0)  
        alpha_euler = alpha_euler * pi/180.0
        beta_euler  = beta_euler  * pi/180.0
        gamma_euler = gamma_euler * pi/180.0

      	return
	end
!     --------------------------------------------------------------------------    
!     sort energies in ascending order
!     --------------------------------------------------------------------------    
      subroutine energy_sort(energy,ind,nex,np)
       	implicit none
        real*8  energy(17*17),eex(17*17)
        integer ind(17*17),nex
        integer ifail,i,np
c       ------------------------------------------------- 
	external M01DAF   ! NAG - sorts real numbers 
	external M01ZAF 
	external M01CAF
c       ------------------------------------------------- 
        do i=1,nex
           eex(i) = energy(i)
        end do
        ifail=1
	call M01DAF(eex,1,nex,'A',ind,ifail)
	call M01ZAF(ind,1,nex,ifail)
	call M01CAF(eex,1,nex,'A',ifail)

        do i=1,nex
          if(eex(i).ge.0.0d0) then
             np = i
             return
          end if
        end do
      return
      end
!     --------------------------------------------------------------------------    
      subroutine shapedis
       	implicit none
        include 'CF_SOURCES:FOCUS.INC' 
        include 'CF_SOURCES:DEGENERATION.INC' 

 	real*8  c_fmevkelvin
 	integer out,i
 	character*18 text_line
	integer ind(17*17),npositiv   

        call energy_sort(e_excitation,ind,n_excitation,npositiv)
        
                                
  1     format(1x,a29,f7.2,a7)
  2     format(1x,a37,f4.2,a38) 
  3     format(1x,a1,f4.2,a6)	     
 
 	out=5	
 	write(out,*)
 	write(out,1)' The sample temperature is : ',temp,' Kelvin' 
 	write(out,2)' Energy levels which are closer than ',de,' Kelvin are assumed to be degenerated.'  
 	write(out,*)' Only those excitations are printed out whose intensities are greater or equal'
        write(out,3)' ',di,' barn.' 
        write(out,*)' Excitations with negative energies are not printed out.' 
        write(out,*)' The line shape of a negative excitation is determined by the positive one.'     
        write(out,*)' The excitations are ordered by ascending energies.'
        write(out,*)
 	
 
  4	format(1x, i3,' : ',f6.2,' meV',3x,f6.2,' barn',3x,a18)       
  6	format(1x, i3,' : ',f16.2,' meV',3x,f6.2,' barn',3x,a18)       
  5     format(1x,a21,'            ',a18)   
 	              

        write(*,5) '  0 :    elastic line',text_line(lineshape(0))           

 	do i=npositiv,n_excitation
           if(e_excitation(ind(i))/c_fmevkelvin().lt.1.0d3) then
               write(out,4) i-npositiv+1,e_excitation(ind(i))/c_fmevkelvin(),i_excitation(ind(i)), 
     *                      text_line(lineshape(ind(i)))
           else
               write(out,6) i-npositiv+1,e_excitation(ind(i))/c_fmevkelvin(),i_excitation(ind(i)), 
     *                      text_line(lineshape(ind(i)))
           end if
 	end do               
        write(out,*)               
      return
      end                                                         
!     -------------------------------------------------------------------------- 
      subroutine iondis
       	implicit none
 	external r_name ! cf_fabi.for
 	character*4 r_name
      include 'CF_SOURCES:FOCUS.INC'
        integer out
        real*8 j,a0

 	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)     
        j = (dimj-1)/2.0
        a0= 0.52917706d0 
        out=5 
        write(out,*)
   	write(out,'(a1,1x,a40,a4)'     ) ' ', 'The calculation is for rare earth ion : ',r_name(nre)
   	write(out,'(a1,1x,a40,i2)'     ) ' ', 'Number of 4f electrons in the shell   : ',nre
   	write(out,'(a1,1x,a40,f4.1)'   ) ' ', 'Total angular momentum J              : ',j
   	write(out,'(a1,1x,a40,f9.6)'   ) ' ', 'Lande-factor g_J                      : ',gj
   	write(out,'(a1,1x,a40,f9.6)'   ) ' ', 'Moment      g_J*J              [mu_B] : ',gj*j
   	write(out,'(a1,1x,a40,f9.6)'   ) ' ', 'Eff. moment g_J*sqrt{J(J+1)}   [mu_B] : ',gj*sqrt(j*(j+1))
   	write(out,'(a1,1x,a40,1pe11.4)') ' ', 'Radial Integral <r^2>   [Angstroem^2] : ',r2*a0**2
   	write(out,'(a1,1x,a40,1pe11.4)') ' ', 'Radial Integral <r^4>   [Angstroem^4] : ',r4*a0**4
    	write(out,'(a1,1x,a40,1pe11.4)') ' ', 'Radial Integral <r^6>   [Angstroem^6] : ',r6*a0**6
   	write(out,'(a1,1x,a40,1pe11.4)') ' ', 'Stevens factor alpha_J                : ',alphaj
   	write(out,'(a1,1x,a40,1pe11.4)') ' ', 'Stevens factor beta_J                 : ',betaj
   	write(out,'(a1,1x,a40,1pe11.4)') ' ', 'Stevens factor gamma_J                : ',gammaj
        write(out,*)
      return
      end
!     --------------------------------------------------------------------------    
      subroutine momentdis
       	implicit none
      include 'CF_SOURCES:FOCUS.INC' 
 	integer out
 	real*8 mx,my,mz  
 	real*8 c_x_moment,c_y_moment,c_z_moment        
 
 	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)     
                                                 
    	mx = c_x_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
    	my = c_y_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
    	mz = c_z_moment(gj,energy,wavefunction,dimj,occupation_factor,temp)
 
 	out=5
   	call p_moments(out,mx,my,mz,gj,dimj)                                                     
 	
      return
      end                                                         
!     --------------------------------------------------------------------------    
      subroutine symdis
       	implicit none
      include 'CF_SOURCES:FOCUS.INC' 
 	integer out
        out = 5
        write(out,*)
 	write(out,*)' point symmetry of crystal field | symmetry number '
 	if(symmetry.eq.0     ) then
                               write(out,*)' Ci C1 ---------------------------------- 0 -------'
 	else if(symmetry.eq.1) then
                               write(out,*)' C2 Cs C2h ------------------------------ 1 -------'
 	else if(symmetry.eq.2) then
                               write(out,*)' C2v D2 D2h ----------------------------- 2 -------'
 	else if(symmetry.eq.3) then
                               write(out,*)' C4 S4 C4h ------------------------------ 3 -------'
 	else if(symmetry.eq.4) then
                               write(out,*)' D4 C4v D2d D4h ------------------------- 4 -------'
 	else if(symmetry.eq.5) then
                               write(out,*)' C3 S6 ---------------------------------- 5 -------'
  	else if(symmetry.eq.6) then
                               write(out,*)' D3 C3v D3d ----------------------------- 6 -------'
 	else if(symmetry.eq.7) then
                               write(out,*)' C6 C3h C6h D6 C6v D3h D6h -------------- 7 -------'   
 	else 
                               write(out,*)' T Td Th O Oh --------------------------- 8 -------'
        end if
        write(out,*)
      return
      end            
!     --------------------------------------------------------------------------       
      subroutine excdis
 	implicit none
        include 'CF_SOURCES:FOCUS.INC'     
        include 'CF_SOURCES:DEGENERATION.INC' 

 	integer out
 
! prints out the calculated excitations for a given temperature temp
 
 	out = 5
 	call p_excitations(out,de,di)
      return
      end                                                         
!     --------------------------------------------------------------------------  
      subroutine intensitydis
        implicit none
        include 'CF_SOURCES:FOCUS.INC'

        integer out
        out=5     
        call p_energy_values(out,energy,dimj,fmevkelvin)
        call p_intensities(out,dimj,energy,intensity,temp)
      return
      end 
!     --------------------------------------------------------------------------  
      subroutine levelscheme
        implicit none
      include 'CF_SOURCES:FOCUS.INC'
        integer out
        out=5     
        call p_energy_values(out,energy,dimj,fmevkelvin)
        call p_mat_powder(out,dimj,energy,jt2mat)     
      return
      end 
!     --------------------------------------------------------------------------       
      subroutine wavefunctions
        implicit none
      include 'CF_SOURCES:FOCUS.INC'
        integer out
        out=5
        call p_wave_functions(out,energy,wavefunction,dimj)     
      return
      end 
!     --------------------------------------------------------------------------       
      subroutine matrixelements
        implicit none
      include 'CF_SOURCES:FOCUS.INC'
        integer out
        out=5     
        call p_energy_values(out,energy,dimj,fmevkelvin)
        call p_matrix_elements(out,dimj,energy,jx2mat,jy2mat,jz2mat,jt2mat)     
      return
      end 
!     --------------------------------------------------------------------------       
      subroutine eulerdis
        implicit none
        include 'CF_SOURCES:CANTED.INC'

        integer out
        real*8 pi
        pi = 4.0d0*atan(1.0d0)
        out=5
  1     format(x,x,a17,f8.2,a8)
        write(out,*)
        write(out,*) ' The R3+ moment is canted compared to the crystal system (a,b,c)'
        write(out,*) ' by the three following Euler-angles: '
        write(out,1) ' alpha (Euler) = ',alpha_euler*180.0d0/pi,' degrees'
        write(out,1) ' beta  (Euler) = ',beta_euler *180.0d0/pi,' degrees'
        write(out,1) ' gamma (Euler) = ',gamma_euler*180.0d0/pi,' degrees'
        write(out,*)
      return
      end
!     --------------------------------------------------------------------------  
!     transform the akq-parameter to bkq-parameter
!     --------------------------------------------------------------------------  
      subroutine t_akq_bkq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC'  
        integer i,ncftot,no_cf_par,n0,no_0_par
        n0 = no_0_par()
        call akq_bkq
        ncftot = no_cf_par(symmetry,nre)                                  
        do i=n0+1,n0+ncftot
             nam(i )(4:4) = 'B'
             nam(i )(10:20)= '    (meV)  '
        end do
      return
      end 
!     --------------------------------------------------------------------------  
!     transform the bkq-parameter to akq-parameter
!     --------------------------------------------------------------------------  
      subroutine t_bkq_akq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC'  
        integer i,ncftot,no_cf_par,n0,no_0_par
        n0 = no_0_par()
        call bkq_akq
        ncftot = no_cf_par(symmetry,nre)                                  
        do i=n0+1,n0+ncftot
             nam(i )(4:4) = 'A'
             nam(i )(10:20)= '(meV/a0**k)'
        end do
       return
      end 
!     --------------------------------------------------------------------------  
!     transform the bkq-parameter to wx-parameter
!     --------------------------------------------------------------------------  
      subroutine t_bkq_wx(setwx)
        implicit none
        integer setwx
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC'  
        integer i,ncftot,no_cf_par,n0,no_0_par
        call bkq_wx(setwx)
	if(setwx.eq.0) return
        ncftot = no_cf_par(symmetry,nre)                                  
        n0     = no_0_par()
        if(ncftot.ge.1) then               
             nam(n0+1 )(9:13) = ': W  '
             nam(n0+1 )(14:20)= '(meV)  '    
        end if
        if(ncftot.eq.2) then               
             nam(n0+2 )(9:13) = ': x  '
             nam(n0+2 )(14:20)= '       '    
        end if
       return
      end 
!     --------------------------------------------------------------------------  
!     transform the  wx-parameter to bkq-parameter
!     --------------------------------------------------------------------------  
      subroutine t_wx_bkq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC'  
        integer i,ncftot,no_cf_par,n0,no_0_par
        n0 = no_0_par()
        call wx_bkq
        ncftot = no_cf_par(symmetry,nre)                                  
        if(ncftot.ge.1) then               
             nam(n0+1 )(9:13) = '     '
             nam(n0+1 )(14:20)= '(meV)  '    
        end if
        if(ncftot.eq.2) then               
             nam(n0+2 )(9:13) = '     '
             nam(n0+2 )(14:20)= '(meV)  '    
        end if
       return
      end 
!     --------------------------------------------------------------------------  
!     transform the bkq-parameter to vkq-parameter
!     --------------------------------------------------------------------------  
      subroutine t_bkq_vkq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC'  
        integer i,ncftot,no_cf_par,n0,no_0_par
        n0 = no_0_par()
        call bkq_vkq
        do i=n0,n0-5,-1
             nam(i )(1:1) = 'V'
             nam(i )(10:20)= '    (meV)  '
        end do
        ncftot = no_cf_par(symmetry,nre)                                  
        do i=n0+1,n0+ncftot
             nam(i )(4:4) = 'V'
        end do
      return
      end 
!     --------------------------------------------------------------------------  
!     transform the vkq-parameter to bkq-parameter
!     --------------------------------------------------------------------------  
      subroutine t_vkq_bkq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC' 
        integer i,ncftot,no_cf_par,n0,no_0_par
        n0 = no_0_par()
        call vkq_bkq
        do i=n0,n0-5,-1
             nam(i )(1:1) = 'B'
             nam(i )(10:20)= '    (Tesla)'
        end do
        ncftot = no_cf_par(symmetry,nre)                                  
        do i=n0+1,n0+ncftot
             nam(i )(4:4) = 'B'
        end do
      return
      end 
c--------------------------------------------------------
c alter from akq to bkq parameter type
c       transform the Akq to Bkq
c--------------------------------------------------------
      subroutine akq_bkq
	implicit none
c       ---------------------------------------------------------------------
        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:LIMIT.INC'
c       ---------------------------------------------------------------------
        real*8 pv(MAX_PAR)
        integer no_cf_par,ncftot,n0,no_0_par,k
c       ---------------------------------------------------------------------
        type = 'Bkq'
c       ---------------------------------------------------------------------

c       ---------------------------------------------------------------------
c       transform parameters
c       ---------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
        n0     = no_0_par()
        do k=1,ncftot
           pv(k) = p(n0+k)
        end do
        call transform(pv,bkq,symmetry,nre)
        call c_akq_bkq(nre,bkq,bkq)
        call itransform(pv,bkq,symmetry,nre)
        do k=1,ncftot
           p(n0+k) = pv(k)
        end do
c       ---------------------------------------------------------------------
c       transform limits
c       ---------------------------------------------------------------------
        call lim_init(symmetry,nre)
        call lim_trans(symmetry,nre)
	call c_akq_bkq(nre,pmikq,pmikq)
	call c_akq_bkq(nre,pmakq,pmakq)
	call c_akq_bkq(nre,pmiokq,pmiokq)
	call c_akq_bkq(nre,pmaokq,pmaokq)
        call ilim_trans(symmetry,nre)
        call ilim_init(symmetry,nre)

      return
      end
c--------------------------------------------------------
c alter from bkq to akq parameter type
c       transform the Bkq to Akq
c--------------------------------------------------------
      subroutine bkq_akq
	implicit none
c       ---------------------------------------------------------------------
        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:LIMIT.INC'
c       ---------------------------------------------------------------------
        real*8 pv(MAX_PAR)
        integer no_cf_par,ncftot,n0,no_0_par,k
c       ---------------------------------------------------------------------
        type = 'Akq'
c       ---------------------------------------------------------------------
c       ---------------------------------------------------------------------
c       transform parameters
c       ---------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
        n0     = no_0_par()
        do k=1,ncftot
           pv(k) = p(n0+k)
        end do
        call transform(pv,bkq,symmetry,nre)
        call c_bkq_akq(nre,bkq,bkq)
        call itransform(pv,bkq,symmetry,nre)
        do k=1,ncftot
           p(n0+k) = pv(k)
        end do
c       ---------------------------------------------------------------------
c       transform limits
c       ---------------------------------------------------------------------
        call lim_init(symmetry,nre)
        call lim_trans(symmetry,nre)
	call c_bkq_akq(nre,pmikq,pmikq)
	call c_bkq_akq(nre,pmakq,pmakq)
	call c_bkq_akq(nre,pmiokq,pmiokq)
	call c_bkq_akq(nre,pmaokq,pmaokq)
        call ilim_trans(symmetry,nre)
        call ilim_init(symmetry,nre)

      return
      end
c--------------------------------------------------------
c alter from bkq to vkq parameter type
c--------------------------------------------------------
      subroutine bkq_vkq
	implicit none
c       ---------------------------------------------------------------------
        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:LIMIT.INC'
c       ---------------------------------------------------------------------
        integer k,m
	real*8 c_myb
	real*8 c_operator_norm  ! is a function which is called
	real*8 norm
        real*8 f_bmol,f_bext, dummy
	complex*16 okq(17,17)
        real*8 ssbkq(0:6,0:6)
        real*8 pv(MAX_PAR)
        integer no_cf_par,ncftot,n0,no_0_par
c       ---------------------------------------------------------------------
        type = 'Vkq'
c       ---------------------------------------------------------------------
c       initialize some rare earth constants
c       ---------------------------------------------------------------------
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)           
c       magneton of Bohr in kelvin per tesla         
	myb    = c_myb()
c       ---------------------------------------------------------------------
c       transform parameters Bext,Bmol to Vext,Vmol
c       ---------------------------------------------------------------------
c       the Vkq are in kelvin and are the parameters which uses the normalized 
c       operators Okq/||Okq|| and j+/||j+|| and j-/||j-|| and jz/||jz||
c       calculates the Bkq,B+,B- and Bz from the given Vkq,V+,V- and Vz 
        n0 = no_0_par()
        bext(1) = p(n0-5)
        bext(2) = p(n0-4)
        bext(3) = p(n0-3)
        bmol(1) = p(n0-2)
        bmol(2) = p(n0-1)
        bmol(3) = p(n0)
	f_bmol = 2*(gj-1)*myb
	f_bext = gj*myb
	norm = c_operator_norm(dimj,0,0,sbkq,okq)    ! 00 = j+ note:||j+||=||j-||
	bmol(1) = 0.5d0*bmol(1)*norm*f_bmol !
	bext(1) = 0.5d0*bext(1)*norm*f_bext !  V- / ||j+|| = 1/2*f*B-
	bmol(2) = 0.5d0*bmol(2)*norm*f_bmol !   
	bext(2) = 0.5d0*bext(2)*norm*f_bext !  V+ / ||j-|| = 1/2*f*B+
	norm = c_operator_norm(dimj,0,1,sbkq,okq)    ! 01 = jz
	bmol(3) = bmol(3)*norm*f_bmol      !
	bext(3) = bext(3)*norm*f_bext      ! Vz / ||jz|| = f*Bz
        p(n0-5) = bext(1)
        p(n0-4) = bext(2)
        p(n0-3) = bext(3)
        p(n0-2) = bmol(1)
        p(n0-1) = bmol(2)
        p(n0  ) = bmol(3)
c       ---------------------------------------------------------------------
c       transform limits of Bext,Bmol to Vext,Vmol
c       ---------------------------------------------------------------------
c       the Vkq are in kelvin and are the parameters which uses the normalized 
c       operators Okq/||Okq|| and j+/||j+|| and j-/||j-|| and jz/||jz||
c       calculates the Bkq,B+,B- and Bz from the given Vkq,V+,V- and Vz 
        n0 = no_0_par()
        bext(1) = pmin(n0-5)
        bext(2) = pmin(n0-4)
        bext(3) = pmin(n0-3)
        bmol(1) = pmin(n0-2)
        bmol(2) = pmin(n0-1)
        bmol(3) = pmin(n0)
	f_bmol = 2*(gj-1)*myb
	f_bext = gj*myb
	norm = c_operator_norm(dimj,0,0,sbkq,okq)    ! 00 = j+ note:||j+||=||j-||
	bmol(1) = 0.5d0*bmol(1)*norm*f_bmol !
	bext(1) = 0.5d0*bext(1)*norm*f_bext !  V- / ||j+|| = 1/2*f*B-
	bmol(2) = 0.5d0*bmol(2)*norm*f_bmol !   
	bext(2) = 0.5d0*bext(2)*norm*f_bext !  V+ / ||j-|| = 1/2*f*B+
	norm = c_operator_norm(dimj,0,1,sbkq,okq)    ! 01 = jz
	bmol(3) = bmol(3)*norm*f_bmol      !
	bext(3) = bext(3)*norm*f_bext      ! Vz / ||jz|| = f*Bz
        pmin(n0-5) = bext(1)
        pmin(n0-4) = bext(2)
        pmin(n0-3) = bext(3)
        pmin(n0-2) = bmol(1)
        pmin(n0-1) = bmol(2)
        pmin(n0  ) = bmol(3)

        bext(1) = pmax(n0-5)
        bext(2) = pmax(n0-4)
        bext(3) = pmax(n0-3)
        bmol(1) = pmax(n0-2)
        bmol(2) = pmax(n0-1)
        bmol(3) = pmax(n0)
	f_bmol = 2*(gj-1)*myb
	f_bext = gj*myb
	norm = c_operator_norm(dimj,0,0,sbkq,okq)    ! 00 = j+ note:||j+||=||j-||
	bmol(1) = 0.5d0*bmol(1)*norm*f_bmol !
	bext(1) = 0.5d0*bext(1)*norm*f_bext !  V- / ||j+|| = 1/2*f*B-
	bmol(2) = 0.5d0*bmol(2)*norm*f_bmol !   
	bext(2) = 0.5d0*bext(2)*norm*f_bext !  V+ / ||j-|| = 1/2*f*B+
	norm = c_operator_norm(dimj,0,1,sbkq,okq)    ! 01 = jz
	bmol(3) = bmol(3)*norm*f_bmol      !
	bext(3) = bext(3)*norm*f_bext      ! Vz / ||jz|| = f*Bz
        pmax(n0-5) = bext(1)
        pmax(n0-4) = bext(2)
        pmax(n0-3) = bext(3)
        pmax(n0-2) = bmol(1)
        pmax(n0-1) = bmol(2)
        pmax(n0  ) = bmol(3)

        do k=n0-5,n0
           if(pmax(k).lt.pmin(k)) then
              dummy   = pmin(k)
	      pmin(k) = pmax(k)
	      pmax(k) = dummy
	   end if
	end do
c       ---------------------------------------------------------------------
c       transform the parameters Bkq to Vkq
c       ---------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
        do k=1,ncftot
           pv(k) = p(n0+k)
        end do
        call transform(pv,bkq,symmetry,nre)
	do k=2,6,2
	   do m=0,k
	     call cifnull(bkq(k,m))
             if(bkq(k,m).ne.0.0d0) then
	        ssbkq(k,m) = 1 ! note ||re(Okq)|| = ||im(Okq)|| for q<>0 
		norm = c_operator_norm(dimj,k,m,ssbkq,okq)
                bkq(k,m) = bkq(k,m)*norm 
	     endif
	   end do
        end do
        call itransform(pv,bkq,symmetry,nre)
        do k=1,ncftot
           p(n0+k) = pv(k)
        end do
c       ---------------------------------------------------------------------
c       transform the limits of the crystal field parameters
c       ---------------------------------------------------------------------
        call lim_init(symmetry,nre)
	call lim_trans(symmetry,nre)
	do k=2,6,2
	   do m=0,k
	        ssbkq(k,m) = 1 ! note ||re(Okq)|| = ||im(Okq)|| for q<>0 
		norm = c_operator_norm(dimj,k,m,ssbkq,okq)
                pmikq( k,m) = pmikq( k,m)*norm 
                pmakq( k,m) = pmakq( k,m)*norm 
                pmiokq(k,m) = pmiokq(k,m)*norm 
                pmaokq(k,m) = pmaokq(k,m)*norm 
	   end do
        end do
	call ilim_trans(symmetry,nre)
	call ilim_init(symmetry,nre)
      return
      end
c--------------------------------------------------------
c alter from Vkq to Bkq parameter type
c--------------------------------------------------------
      subroutine vkq_bkq
	implicit none
c       ---------------------------------------------------------------------
        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.INC'
	include 'CF_SOURCES:LIMIT.INC'
c       ---------------------------------------------------------------------
        integer k,m
	real*8 c_myb
	real*8 c_operator_norm  ! is a function which is called
	real*8 norm
        real*8 f_bmol,f_bext
	complex*16 okq(17,17)
        real*8 ssbkq(0:6,0:6)
        real*8 pv(MAX_PAR)
        integer no_cf_par,ncftot,n0,no_0_par
c       ---------------------------------------------------------------------
        type = 'Bkq'
c       ---------------------------------------------------------------------
c       initialize some rare earth constants
c       ---------------------------------------------------------------------
	call i_re_constants(nre,gj,dimj,alphaj,betaj,gammaj,r2,r4,r6)           
c       magneton of Bohr in kelvin per tesla         
	myb    = c_myb()
c       ---------------------------------------------------------------------
c       transform the parameters Vext,Vmol to Bext,Bmol
c       ---------------------------------------------------------------------
c       the Vkq are in kelvin and are the parameters which uses the normalized 
c       operators Okq/||Okq|| and j+/||j+|| and j-/||j-|| and jz/||jz||
c       calculates the Bkq,B+,B- and Bz from the given Vkq,V+,V- and Vz 
        n0 = no_0_par()
        bext(1) = p(n0-5)
        bext(2) = p(n0-4)
        bext(3) = p(n0-3)
        bmol(1) = p(n0-2)
        bmol(2) = p(n0-1)
        bmol(3) = p(n0  )
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
        p(n0-5) = bext(1)
        p(n0-4) = bext(2)
        p(n0-3) = bext(3)
        p(n0-2) = bmol(1)
        p(n0-1) = bmol(2)
        p(n0  ) = bmol(3)
c       ---------------------------------------------------------------------
c       transform the limits 
c       ---------------------------------------------------------------------
        bext(1) = pmin(n0-5)
        bext(2) = pmin(n0-4)
        bext(3) = pmin(n0-3)
        bmol(1) = pmin(n0-2)
        bmol(2) = pmin(n0-1)
        bmol(3) = pmin(n0  )
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
        pmin(n0-5) = bext(1)
        pmin(n0-4) = bext(2)
        pmin(n0-3) = bext(3)
        pmin(n0-2) = bmol(1)
        pmin(n0-1) = bmol(2)
        pmin(n0  ) = bmol(3)

        bext(1) = pmax(n0-5)
        bext(2) = pmax(n0-4)
        bext(3) = pmax(n0-3)
        bmol(1) = pmax(n0-2)
        bmol(2) = pmax(n0-1)
        bmol(3) = pmax(n0  )
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
        pmax(n0-5) = bext(1)
        pmax(n0-4) = bext(2)
        pmax(n0-3) = bext(3)
        pmax(n0-2) = bmol(1)
        pmax(n0-1) = bmol(2)
        pmax(n0  ) = bmol(3)
c       ---------------------------------------------------------------------
c       transform the parameters Vkq to Bkq
c       ---------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
        do k=1,ncftot
           pv(k) = p(n0+k)
        end do
        call transform(pv,bkq,symmetry,nre)
	do k=2,6,2
	   do m=0,k
	     call cifnull(bkq(k,m))
             if(bkq(k,m).ne.0.0d0) then
	        ssbkq(k,m) = 1 ! note ||re(Okq)|| = ||im(Okq)|| for q<>0 
		norm = c_operator_norm(dimj,k,m,ssbkq,okq)
                if(norm.eq.0) norm=1
                bkq(k,m) = bkq(k,m)/norm 
	     endif
	   end do
        end do
        call itransform(pv,bkq,symmetry,nre)
        do k=1,ncftot
           p(n0+k) = pv(k)
        end do
c       ---------------------------------------------------------------------
c       transform the limits
c       ---------------------------------------------------------------------
        call lim_init(symmetry,nre)
	call lim_trans(symmetry,nre)
	do k=2,6,2
	   do m=0,k
	        ssbkq(k,m) = 1 ! note ||re(Okq)|| = ||im(Okq)|| for q<>0 
		norm = c_operator_norm(dimj,k,m,ssbkq,okq)
                if(norm.eq.0) norm=1
                pmikq( k,m) = pmikq( k,m)/norm
                pmakq( k,m) = pmakq( k,m)/norm
                pmiokq(k,m) = pmiokq(k,m)/norm
                pmaokq(k,m) = pmaokq(k,m)/norm
	   end do
        end do
	call ilim_trans(symmetry,nre)
	call ilim_init(symmetry,nre)
      return
      end
!     --------------------------------------------------------------------------  
!     converts  0<= i <=999 to '1  ' ... '10 ' ... '100'
!     --------------------------------------------------------------------------  
      character*3 function ic(n)
        implicit none
        integer n,i1,i2,i3
        character*1 c
        i1 = dble(n       )/100.0d0
        i2 = dble(n-i1*100)/10.0d0
        i3 = n-i1*100-i2*10
        if(i1.eq.0) then
           if(i2.eq.0) then
              ic=c(i3)//' '//' '
           else
              ic=c(i2)//c(i3)//' '
           end if
        else
           ic=c(i1)//c(i2)//c(i3)
        end if
      return
      end
      character*1 function c(i)
         implicit none
         integer i
         if(i.eq.0) c = '0'
         if(i.eq.1) c = '1'
         if(i.eq.2) c = '2'
         if(i.eq.3) c = '3'
         if(i.eq.4) c = '4'
         if(i.eq.5) c = '5'
         if(i.eq.6) c = '6'
         if(i.eq.7) c = '7'
         if(i.eq.8) c = '8'
         if(i.eq.9) c = '9'
      return
      end 
!     --------------------------------------------------------------------------  
!     alter parameter to bkq-parameter
!     --------------------------------------------------------------------------  
      subroutine a_to_bkq(setbkq,setvkq,setxr,setakq,setwx)
        implicit none
        integer  setbkq,setvkq,setxr,setakq,setwx

        if(setbkq.eq.1) then  ! to get the parameter names right in the case of @ s io .. or @ s c ..
             call t_bkq_akq
             call t_akq_bkq
        else if(setvkq.eq.1) then
             call t_vkq_bkq
             setvkq = 0
        else if(setakq.eq.1) then
             call t_akq_bkq
             setakq = 0
        else if(setwx.eq.1) then
             call t_wx_bkq
             setwx = 0
        else if(setxr.eq.1) then
             call t_xr_vkq
             call t_vkq_bkq
             setxr  = 0
        end if
        setbkq = 1
      return
      end
!     --------------------------------------------------------------------------  
!     alter parameter to akq-parameter
!     --------------------------------------------------------------------------  
      subroutine a_to_akq(setbkq,setvkq,setxr,setakq,setwx)
        implicit none
        integer  setbkq,setvkq,setxr,setakq,setwx

        if(setakq.eq.1) then   ! to get the parameter names right in the case of @ s io .. or @ s c ..
             call t_akq_bkq
             call t_bkq_akq
        else if(setvkq.eq.1) then
             call t_vkq_bkq
             call t_bkq_akq
             setvkq = 0
        else if(setbkq.eq.1) then
             call t_bkq_akq
             setbkq = 0
        else if(setwx.eq.1) then
             call t_wx_bkq
             call t_bkq_akq
             setwx = 0
        else if(setxr.eq.1) then
             call t_xr_vkq
             call t_vkq_bkq
             call t_bkq_akq
             setxr  = 0
        end if
        setakq = 1
      return
      end
!     --------------------------------------------------------------------------  
!     alter parameter to vkq-parameter
!     --------------------------------------------------------------------------  
      subroutine a_to_vkq(setbkq,setvkq,setxr,setakq,setwx)    
        implicit none
        integer  setbkq,setvkq,setxr,setakq,setwx
        if(setvkq.eq.1) then ! to get the parameter names right in the case of @ s io .. or @ s c ..
             call t_vkq_bkq
             call t_bkq_vkq
        else if(setbkq.eq.1) then
             call t_bkq_vkq
             setbkq = 0
        else if(setakq.eq.1) then
             call t_akq_bkq
             call t_bkq_vkq
             setakq = 0
        else if(setwx.eq.1) then
             call t_wx_bkq
             call t_bkq_vkq
             setwx = 0
        else if(setxr.eq.1) then
             call t_xr_vkq
             setxr  = 0
        end if
        setvkq = 1
      return
      end
!     --------------------------------------------------------------------------  
!     alter parameter to Wx-parameter
!     --------------------------------------------------------------------------  
      subroutine a_to_wx(setbkq,setvkq,setxr,setakq,setwx)
        implicit none
        integer  setbkq,setvkq,setxr,setakq,setwx

        if(setwx.eq.1) then ! to get the parameter names right in the case of @ s io .. or @ s c ..
             call t_wx_bkq
	     call t_bkq_wx(setwx)
        else if(setbkq.eq.1 ) then
             call t_bkq_wx(setwx)
             if(setwx.eq.1) then 
                  setbkq=0
             end if
        else if(setakq.eq.1) then
             call t_akq_bkq
             call t_bkq_wx(setwx)
             if(setwx.eq.1) then 
                  setakq=0
             else
                  call t_bkq_akq
             end if
        else if(setxr.eq.1) then
             call t_xr_vkq
             call t_vkq_bkq
             call t_bkq_wx(setwx)
             if(setwx.eq.1) then 
                  setxr=0
             else
                  call t_bkq_vkq
                  call t_vkq_xr(setxr)
             end if

        else if(setvkq.eq.1) then
             call t_vkq_bkq
             call t_bkq_wx(setwx)
             if(setwx.eq.1) then 
                  setvkq=0
             else
                  call t_bkq_vkq
             end if
        end if
      return
      end
!     --------------------------------------------------------------------------  
!     alter parameter to xr-parameter
!     --------------------------------------------------------------------------  
      subroutine a_to_xr(setbkq,setvkq,setxr,setakq,setwx)
        implicit none
        integer  setbkq,setvkq,setxr,setakq,setwx

        if(setxr.eq.1) then ! to get the parameter names right in the case of @ s io .. or @ s c ..
             call t_xr_vkq
	     call t_vkq_xr(setxr)
        else if(setbkq.eq.1 ) then
             call t_bkq_vkq
             call t_vkq_xr(setxr)
             if(setxr.eq.1) then 
                  setbkq=0
             else
                  call t_vkq_bkq
             end if

        else if(setwx.eq.1) then
             call t_wx_bkq
             call t_bkq_vkq
             call t_vkq_xr(setxr)
             if(setxr.eq.1) then 
                  setwx=0
             else
                  call t_vkq_bkq
                  call t_bkq_wx(setwx)
             end if
             
        else if(setakq.eq.1) then
             call t_akq_bkq
             call t_bkq_vkq
             call t_vkq_xr(setxr)
             if(setxr.eq.1) then 
                  setakq=0
             else
                  call t_vkq_bkq
                  call t_bkq_akq
             end if

        else if(setvkq.eq.1) then
             call t_vkq_xr(setxr)
             if(setxr.eq.1) setvkq=0
        end if
      return
      end
!     --------------------------------------------------------------------------  
!     correct xr-parameter after <clear> <p> [n1-n2] [n1] ... command 
!     --------------------------------------------------------------------------  
      subroutine xr_add(nindv0)
        implicit none
        integer nindv0
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC' 
      include 'CF_SOURCES:XR.INC' 
        integer ncftot,no_cf_par
        integer indv(MAX_PAR),n_indv
        integer no_f,no_fixed
        integer i,j,k
        real*8  rs,r,normi
        character*3 ic
c       ------------------------------------------------------------------------
c       look if there are new x-values          
c       ------------------------------------------------------------------------
        if(nindv0.eq.nindvold) return
c       ------------------------------------------------------------------------
c       correct for new value of the radius R
c       ------------------------------------------------------------------------
        r  = p(indvold(1))
        rs = r**2
        do i=nindv0+1,nindvold
           j = indvold(i)
           rs = rs + p(j)**2
        end do
        rs = sqrt(rs)*sign(1.0d0,r)
        p(indvold(1)) = rs
c       ------------------------------------------------------------------------
c       calculate the new x-values from the v-values
c       ------------------------------------------------------------------------
        do i=nindv0+1,nindvold
           normi = r**2
           do k=nindv0+1,i
              j = indvold(k)
              normi = normi + p(j)**2
           end do 
           normi = sqrt(normi)
           j = indvold(i)
           call ifnull(normi)
           if(normi.ne.0.0d0) then
              p(j) = p(j)/normi*sign(1.0d0,rs)
           else
              p(j) = 0.0d0               ! after the appendix in the manual see:
                                         ! The n-dimensional spherical parameter,
                                         ! x_j can be set to an arbitary value
                                         ! if the norm vanishes
           end if
              nam(j)(9:13) = ':x'//ic(i)
              nam(j)(14:20)= '       '
              pminold(j) = pmin(j)
              pmaxold(j) = pmax(j)
              pmin(j) = -1.0d0
              pmax(j) =  1.0d0
              fixold(j) = j
        end do
      return
      end
!     --------------------------------------------------------------------------  
!     transform the Bkq-parameter to Wx-parameter
!     --------------------------------------------------------------------------  
      subroutine bkq_wx(setwx)
        implicit none
        integer setwx
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC' 
      include 'CF_SOURCES:XR.INC' 
        integer ncftot,no_cf_par,n0,no_0_par
        real*8 pv(MAX_PAR),xxx,www
        integer k
c       --------------------------------------------------------------------------
c       test if transformation to Wx-parameter is possible
c       --------------------------------------------------------------------------
        setwx = 0
        if(symmetry.ne.8) then
           call remark ('ERROR: The transformation to the Lea-Leask-Wolf parameter W,x is not')
           call remark ('       possible because the symmetry is not cubic. ')
           return
        end if
c       --------------------------------------------------------------------------
c       mark that the transformation to Wx-parameter was successful
c       --------------------------------------------------------------------------
        setwx=1
        type = 'WX'
c       --------------------------------------------------------------------------
c       transform to wx-parameter                                        
c       --------------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
        n0 = no_0_par()
        do k=1,ncftot
           pv(k) = p(n0+k)
        end do
        bkq(4,0) = pv(1)
        bkq(6,0) = pv(2)
        call c_xW(nre,xxx,www,bkq)
        pv(1) = www
        pv(2) = xxx
        call transform(pv,bkq,symmetry,nre) ! to get the right bkq!
        call itransform(pv,bkq,symmetry,nre)
        do k=1,ncftot
           p(n0+k) = pv(k)
           pminold(n0+k) = pmin(n0+k)
           pmaxold(n0+k) = pmax(n0+k)
	   pmin(n0+k)    = 0.0d0
	   pmax(n0+k)    = 0.0d0
        end do
        if(ncftot.eq.2) then
           pmin(n0+2) = -1.0d0
           pmax(n0+2) =  1.0d0
        end if
      return
      end  
!     --------------------------------------------------------------------------  
!     transform the Wx-parameter to Bkq-parameter
!     --------------------------------------------------------------------------  
      subroutine wx_bkq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC' 
      include 'CF_SOURCES:XR.INC' 
        integer ncftot,no_cf_par,n0,no_0_par
        real*8 pv(MAX_PAR),xxx,www
        integer k
c       --------------------------------------------------------------------------
c       mark that the transformation to Wx-parameter was successful
c       --------------------------------------------------------------------------
        type = 'Bkq'
c       --------------------------------------------------------------------------
c       transform to bkq-parameter                                        
c       --------------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
        n0 = no_0_par()
        do k=1,ncftot
           pv(k) = p(n0+k)
        end do
        www = pv(1)
        xxx = pv(2)
        call t_xW(nre,xxx,www,bkq)
        pv(1) = bkq(4,0)
        pv(2) = bkq(6,0)
        call transform(pv,bkq,symmetry,nre) ! to get the right bkq!
        call itransform(pv,bkq,symmetry,nre)
        do k=1,ncftot
           p(n0+k) = pv(k)
           pmin(n0+k) = pminold(n0+k)
           pmax(n0+k) = pmaxold(n0+k)
        end do
     
      return
      end  
!     --------------------------------------------------------------------------  
!     transform the vkq-parameter to xr-parameter
!     --------------------------------------------------------------------------  
      subroutine t_vkq_xr(setxr)
        implicit none
        integer setxr
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC' 
      include 'CF_SOURCES:XR.INC' 
        integer ncftot,no_cf_par
        integer no_f,no_fixed
        real*8 pv(MAX_PAR),xx(MAX_PAR)
        integer i,j
        character*3 ic
c       --------------------------------------------------------------------------
c       test if transformation to xR-parameter is possible
c       --------------------------------------------------------------------------
        no_f = no_fixed(fixold,indvold,nindvold)
        ncftot = no_cf_par(symmetry,nre)                                  
        setxr = 0
        if(no_f.eq.6+ncftot) then 
           call remark ('ERROR: There are all magnetic field and crystal field parameters fixed.')
           call remark ('       The transformation to (x_i,R)-parameter is not possible.')
           return
        end if
c       --------------------------------------------------------------------------
c       mark that the transformation to xR-parameter was successful
c       --------------------------------------------------------------------------
        setxr=1
        type = 'XR'
c       --------------------------------------------------------------------------
c       transform to xR-parameter                                        
c       --------------------------------------------------------------------------
        do i=1,nindvold
           pv(i) = p(indvold(i))
        end do
        call t_v_x(pv,xx,nindvold)
        do i=1,nindvold
           j = indvold(i)
           p(j) = xx(i)
           if(i.eq.1) then
              nam(j)( 9:13) = ': R  '
              nam(j)(14:20) = '(meV)  '
              pminold(j) = pmin(j)
              pmaxold(j) = pmax(j)
              pmin(j) = 0.0d0
              pmax(j) = 0.0d0
           else
              nam(j)(9:13) = ':x'//ic(i)
              nam(j)(14:20)= '       '
              pminold(j) = pmin(j)
              pmaxold(j) = pmax(j)
              pmin(j) = -1.0d0
              pmax(j) =  1.0d0
           end if
           fixold(j) = j
        end do
      return
      end
!     --------------------------------------------------------------------------  
!     transform the xr-parameter to vkq-parameter
!     --------------------------------------------------------------------------  
      subroutine t_xr_vkq
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC' 
      include 'CF_SOURCES:XR.INC' 
        integer i,j
        real*8 xx(MAX_PAR),pv(MAX_PAR)
        character*3 ic
c       --------------------------------------------------------------------------
c       transform to Vkq-parameter                                        
c       --------------------------------------------------------------------------
        type = 'Vkq'
        do i=1,nindvold
           xx(i) = p(indvold(i))
        end do
        call t_x_v(xx,pv,nindvold)
        do i=1,nindvold
           j = indvold(i)
           p(j) = pv(i)
           if(i.eq.1) then
              nam(j)(9:9) = ' '
              nam(j)(10:20)= '    (meV)  '
              pmin(j) = pminold(j)
              pmax(j) = pmaxold(j)
           else
              nam(j)(9:9) = ' '
              nam(j)(10:20)= '    (meV)  '
              pmin(j) = pminold(j)
              pmax(j) = pmaxold(j)
           end if
           fixold(j) = j
        end do
      return
      end 
!     --------------------------------------------------------------------------  
!     determine which Vkq-parameter and magnetic field components are fixed 
!     and also how many
!     output: fix(i) for i=1...no_cf_par(symmetry,nre)
!     1: Bext(x), 2: Bext(y), 3: Bext(z), 4: Bmol(x), 5: Bmol(y), 6: Bmol(z)
!     7 ... no_cf_par(symmetry,nre): Bkq-parameter
!     -------------------------------------------------------------------------- 
      integer function no_fixed(fix,indv,n_indv)
        implicit none
      include 'CF_SOURCES:FRILLS.INC'
      include 'CF_SOURCES:FOCUS.INC'      
        integer fix(*),indv(*),n,n_indv
        integer i,n0,no_0_par,ncftot,no_cf_par
c       ------------------------------------------------------------------------
c       determine what parameter there are in total
c       stored as fix(j) = +/- j, negativ if fixed
c       ------------------------------------------------------------------------
        call which_fixed(fix)
c       ------------------------------------------------------------------------
c       determine how many parameter are independent of symmetry
c       ------------------------------------------------------------------------
        n0 = no_0_par()
c       ------------------------------------------------------------------------
c       determine how many parameter depend on symmetry
c       ------------------------------------------------------------------------
        ncftot = no_cf_par(symmetry,nre)                                  
c       ------------------------------------------------------------------------
c       n0-5,..,n0 : Bext(x),Bext(y),Bext(z),  Bmol(x),Bmol(y),Bmol(z), 
c       ------------------------------------------------------------------------
        n=0
        do i=n0-5,n0+ncftot
           if(fix(i).gt.0) then
               n = n + 1
               indv(n) = i
           end if
        end do 
        no_fixed = ncftot+6 - n
        n_indv = n
      return 
      end
!     --------------------------------------------------------------------------  
!     determine which parameter are fixed 
!     -------------------------------------------------------------------------- 
      subroutine which_fixed(fix)
        implicit none
        include 'CF_SOURCES:FRILLS.INC'
        include 'CF_SOURCES:FOCUS.INC'      

        integer fix(*)
        integer*4 iarray(MAX_PAR),i,j,k,n
        integer n0,no_0_par,same
c       ------------------------------------------------------------------------
c       determine how many parameter are independent of symmetry
c       ------------------------------------------------------------------------
        n0 = no_0_par()
c       ------------------------------------------------------------------------
c       fix(j) = j  for parameter which are free  in the fit
c       fix(j) =-j  for parameter which are fixed during fitting
c       ------------------------------------------------------------------------
        n=0
        do i=1,nptot
           if(ins(i).eq.0 .or. ins(i).eq.is) then
              n = n + 1
              iarray(n) = i
           end if
        end do
        do i=1,n
           j = iarray(i)
           if(j.le.0 .or. j.gt.nptot) then
              call remark ('ERROR: Parameter no. outside range')
              return
           end if
           fix(j) = -j
           do k=1,np
              if(j.eq.inp(k)) fix(j) = j
!              if(inr(i).ge.n0-5 .and. i.ge.n0-5) fix(j) = j         
!              if( same(i,inr(i)).eq.0 ) fix(j) = -j
           end do
        end do
      return
      end
c---------------------------------------------
c end of file
c---------------------------------------------
