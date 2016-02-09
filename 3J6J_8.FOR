c***********************************************************************
c                                                                      *
c                   Function dreij                                     *
c                                                                      *
c nach Lindner: Drehimpulse in der Quantenmechanik  Seite 38-39        *
c***********************************************************************
      real*8 function dreij(j1,j2,j3,m1,m2,m3)
      implicit none
      real*8 j1,j2,j3,m1,m2,m3
      real*8 m,s,n(3,3),n11,n12,n13,n21,n22,n23,n31,n32,n33
      real*8 fac,min,s0,sum
      integer pz(3),ps(3),i,j,i0,j0
 
      m      =  m1+m2+m3
      s      =  j1+j2+j3
      n(1,1) = -j1+j2+j3
      n(1,2) =  j1-j2+j3
      n(1,3) =  j1+j2-j3
      n(2,1) =  j1-m1
      n(2,2) =  j2-m2
      n(2,3) =  j3-m3
      n(3,1) =  j1+m1
      n(3,2) =  j2+m2
      n(3,3) =  j3+m3
 
      dreij  =  0.0d0
 
c ... auswahlregeln beachten
      if(m.ne.0.0d0) return
 
c ... testen ob alle n(i,j) nicht negativ und ganz sind
      do 10,i=1,3
         do 20,j=1,3
            if( n(i,j)             .lt. 0.0d0 ) return
            if( n(i,j)-int(n(i,j)) .ne. 0.0d0 ) return
 20      continue
 10   continue
 
c ... position (i0,j0) des kleinsten elementes n(i,j) bestimmen
      i0  = 1
      j0  = 1
      min = n(i0,j0)
      do 30,i=1,3
         do 40,j=1,3
            if( n(i,j).lt.min ) then
                i0  = i
                j0  = j
                min = n(i0,j0)
            endif
 40      continue
 30   continue
 
c ... permutatonsvektoren pz(i) und ps(i) bestimmen, so dass
c     n(i0,j0) links oben in der matrix steht und das regge-symbol
c     dabei den wert nicht aendert
      call permut(pz,ps,3,3,i0,j0)
 
      n11 = n( pz(1),ps(1) )
      n12 = n( pz(1),ps(2) )
      n13 = n( pz(1),ps(3) )
      n21 = n( pz(2),ps(1) )
      n22 = n( pz(2),ps(2) )
      n23 = n( pz(2),ps(3) )
      n31 = n( pz(3),ps(1) )
      n32 = n( pz(3),ps(2) )
      n33 = n( pz(3),ps(3) )
 
      dreij = (-1)**( int(n23+n32) )
      dreij = dreij / sqrt( fac(s+1) )
      dreij = dreij * sqrt( fac(n12)/fac(n11) )
      dreij = dreij * sqrt( fac(n13)/fac(n22) )
      dreij = dreij * sqrt( fac(n21)/fac(n33) )
      dreij = dreij * sqrt( fac(n31)/fac(n23)/fac(n32) )
      s0 = (fac(n23)/fac(n23-n11)) * (fac(n32)/fac(n32-n11))
 
      sum = s0
      do 50,i=1,int(n11)
        s0=-s0*(n11+1-i)*(n22+1-i)*(n33+1-i)/(i*(n23-n11+i)*(n32-n11+i))
        sum = sum + s0
 50   continue
      dreij = dreij * sum
 
      return
      end
c***********************************************************************
c                                                                      *
c                   Function sechsj                                    *
c                                                                      *
c nach Lindner: Drehimpulse in der Quantenmechanik  Seite 49-50/58-59  *
c***********************************************************************
      real*8 function sechsj(j1,j2,j3,js1,js2,js3)
      implicit none
      real*8 j1,j2,j3,js1,js2,js3
      real*8 n(3,4)
      real*8 n11,n12,n13,n14
      real*8 n21,n22,n23,n24
      real*8 n31,n32,n33,n34
      real*8 ds,dss,d1,d2,d3
      real*8 fac,min,s0,sum
      integer pz(3),ps(4),i,j,i0,j0
 
      n11    = -j1+j2+j3
      ds     =  j1+js1-j2-js2
      dss    =  j1+js1-j3-js3
      d1     = -j2+js2-j3+js3
      d2     =  j1-js1-j3+js3
      d3     =  j1-js1-j2+js2
 
      n(1,1) = n11
      n(1,2) = n11 + d1
      n(1,3) = n11 + d2
      n(1,4) = n11 + d3
      n(2,1) = n11      + ds
      n(2,2) = n11 + d1 + ds
      n(2,3) = n11 + d2 + ds
      n(2,4) = n11 + d3 + ds
      n(3,1) = n11      + dss
      n(3,2) = n11 + d1 + dss
      n(3,3) = n11 + d2 + dss
      n(3,4) = n11 + d3 + dss
 
      sechsj =  0.0d0
 
c ... auswahlregeln beachten
c ... testen ob alle n(i,j) nicht negativ und ganz sind
      do 10,i=1,3
         do 20,j=1,4
            if( n(i,j)             .lt. 0.0d0 ) return
            if( n(i,j)-int(n(i,j)) .ne. 0.0d0 ) return
 20      continue
 10   continue
 
c ... position (i0,j0) des kleinsten elementes n(i,j) bestimmen
      i0  = 1
      j0  = 1
      min = n(i0,j0)
      do 30,i=1,3
         do 40,j=1,4
            if( n(i,j).lt.min ) then
                i0  = i
                j0  = j
                min = n(i0,j0)
            endif
 40      continue
 30   continue
 
c ... permutatonsvektoren pz(i) und ps(i) bestimmen, so dass
c     n(i0,j0) links oben in der matrix steht und das regge-symbol
c     dabei den wert nicht aendert
      call permut(pz,ps,3,4,i0,j0)
 
      n11 = n( pz(1),ps(1) )
      n12 = n( pz(1),ps(2) )
      n13 = n( pz(1),ps(3) )
      n14 = n( pz(1),ps(4) )
      n21 = n( pz(2),ps(1) )
      n22 = n( pz(2),ps(2) )
      n23 = n( pz(2),ps(3) )
      n24 = n( pz(2),ps(4) )
      n31 = n( pz(3),ps(1) )
      n32 = n( pz(3),ps(2) )
      n33 = n( pz(3),ps(3) )
      n34 = n( pz(3),ps(4) )
 
      sechsj = (-1)**(int(n12+n23+n34))
      sechsj = sechsj * sqrt(fac(n22)/fac(n11))
      sechsj = sechsj * sqrt(fac(n23)/fac(n12))
      sechsj = sechsj * sqrt(fac(n24)/fac(n13))
      sechsj = sechsj * sqrt(fac(n32)/fac(n14))
      sechsj = sechsj * sqrt(fac(n33)/fac(n21))
      sechsj = sechsj * sqrt(fac(n34)/fac(n31))
      sechsj = sechsj * sqrt(fac(n12+n23+n34+1)/fac(n13+n24+n31+1))
      sechsj = sechsj / sqrt(fac(n14+n21+n32+1)*fac(n11+n22+n33+1))
 
      s0 =      (fac(n12)/fac(n12-n11))
      s0 = s0 * (fac(n13)/fac(n13-n11))
      s0 = s0 * (fac(n14)/fac(n14-n11))
 
      sum = s0
      do 50,i=1,int(n11)
        s0=-s0*(n11+1-i)*(n21+1-i)*(n31+1-i)*(n12+n23+n34+1+i)
     .        /(i*(n12-n11+i)*(n13-n11+i)*(n14-n11+i))
        sum = sum + s0
 50   continue
      sechsj = sechsj * sum
 
      return
      end
c***********************************************************************
c                                                                      *
c                   Function  fac                                      *
c                                                                      *
c***********************************************************************
      real*8 function fac(n)
      implicit none
      real*8 n,i
      fac = 0.0d0
      if(n.lt.0.0) return
      fac = 1.0d0
      if(n.eq.0.0d0) return
      do 10,i=1.0d0,n,1.0d0
         fac = fac*i
 10   continue
      return
      end
c***********************************************************************
c                                                                      *
c                   Subroutine permut                                  *
c                                                                      *
c***********************************************************************
      subroutine permut(pz,ps,nz,ns,i0,j0)
      implicit none
      integer nz,ns,i0,j0,i
      integer pz(nz),ps(ns)
 
      pz(1) = i0
      do 10,i=2,nz
         pz(i) = mod(i0+i-1,nz)
         if(pz(i).eq.0) pz(i) = nz
 10   continue
 
      ps(1) = j0
      if(ns.lt.2) return
      do 20,i=2,ns
         ps(i) = mod(j0+i-1,ns)
         if(ps(i).eq.0) ps(i) = ns
 20   continue
      return
      end
c***********************************************************************
c                                                                      *
c                   Subroutine                                         *
c                                                                      *
c***********************************************************************
c     subroutine
c     implicit none
c     return
c     end
c***********************************************************************
c                                                                      *
c                   End of File                                        *
c                                                                      *
c***********************************************************************
