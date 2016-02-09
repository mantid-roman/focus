c include file for user's function
      include 'FRILLS.PAR'

      integer*4 is, nv, nc
      real*8    x,  ycal, p, yc

      common /fr_func/ is, nv, x(MAX_DAT), ycal(MAX_DAT), p(MAX_PAR)
      common /fr_cmps/ nc, yc(MAX_DAT,MAX_CMP)
