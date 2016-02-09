c---------------------------------------------------------------------
c FILE DEGENERATION.INC
c---------------------------------------------------------------------
c     all energy levels of the level scheme which are closer than DE Kelvin are
c     assumed to be degenerated. Only those excitations are plotted which have
c     a transition intensity greater equal than DI barn
      common /DEGENERATION/ de,di,de_default,di_default
      real*8 de,di,de_default,di_default
c---------------------------------------------------------------------
