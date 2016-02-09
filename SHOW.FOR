!
! show.for
!
      subroutine show
	include 'CF_SOURCES:FRILLS.INC'
	include 'CF_SOURCES:XR.INC'
        integer  k

        write(*,'(1x,a11,i2)') 'nindvold = ',nindvold
        do k=1,nindvold
           write(*,'(1x,a8,i2,a4,i2)')    'indvold(',k,') = ',indvold(k)
           write(*,'(1x,a8,i2,a4,i2)')    ' fixold(',k,') = ',fixold( indvold(k))
           write(*,'(1x,a8,i2,a4,f10.4)') 'pminold(',k,') = ',pminold(indvold(k))
           write(*,'(1x,a8,i2,a4,f10.4)') 'pmaxold(',k,') = ',pmaxold(indvold(k))
        end do

        write(*,'(1x,a5,i3)') 'np = ',np
        do k=1,np
           write(*,'(1x,a4,i2,a4,i2)') 'inp(',k,') = ',inp(k)
           write(*,'(1x,a4,i2,a4,i2)') 'inr(',k,') = ',inr(k)
	end do
      return
      end
