!---------------------------------------------------------------------------------------------
! INFO.FOR
!---------------------------------------------------------------------------------------------
! INFO about changes to FOCUS version 1.0
!---------------------------------------------------------------------------------------------
      subroutine i_changes
      implicit none

      integer out      ! messages a written to file out
      integer in       ! reading is from file in
      parameter(in=5,out=5)

      character line*80, c*1 
      integer*4  line_len, getlin, ctoi
      real*4    ctor

      write(*,*)
      call section(out, ' WHAT IS NEW COMPARED TO FOCUS VERSION 1.0 ')
      write(*,*)
      call remark(' The resolution function calculation for a chopper type instrument. FOCUS ')
      call remark(' supports especially the chopper instruments HET and  MARI. Type @FOCUS ')
      call remark(' ALTER RESOLUTION (or @ A R) to set the resolution function values.')
      write(*,*)
      call remark(' The selfshielding and absorption calculation to correct for selfshielding')
      call remark(' and absorption effects in the sample. Type @FOCUS ALTER SELFSHIELDING (or')
      call remark(' @ A S) to set the selfshielding and absorption values.')
      write(*,*)
      call remark(' The possibility to fix any CEF energy level to a predefined value. This ')
      call remark(' option works within the fit and Monte-Carlo runs. Type @FOCUS ALTER')
      call remark(' OVERALL SPLITTING (or @ A O) to set the CEF level to the energy you want.')
      write(*,*)
      call remark(' Switching to another rare earth ion with the CEF parameters will be adapted')
      call remark(' to the new ion. To do this FOCUS transforms the actual parameter scheme to')
      call remark(' the so-called Akq-parameter which should be independent of the rare earth')
      call remark(' ion by definition. FOCUS sets then the rare earth ion to the new one and')
      call remark(' transforms the Akq-parameter back to the actual parameter scheme. E.g. type')
      call remark(' @FOCUS SET ION CE (or @ S I CE) to set the rare earth ion to Ce3+.')
      call prompt(' Hit a key to continue ')
      line_len = getlin(line)

      write(*,*)
      call section (out, ' WHAT IS NEW COMPARED TO FOCUS VERSION 1.0 ')
      write(*,*)
      call remark(' Switching to another CEF symmetry with the CEF parameters will be adapted ')
      call remark(' to the new symmetry. To do this FOCUS transforms the actual parameter scheme')
      call remark(' to the so-called Akq-parameter, applies the symmetry selection rules to    ')
      call remark(' this parameters and transforms the selected Akq-parameter back to the    ')
      call remark(' actual parameter scheme.                                                 ')
      write(*,*)
      call remark(' A lot of informations are now implemented. You can have an information')
      call remark(' about the used scattering function (@ I S), the implemented lineshapes') 
      call remark(' (@ I L), the resolution function (@ I R), how the CEF transition intensity')
      call remark(' is defined (@ I I), who contributed to the used FOCUS library (@ I A) or how')
      call remark(' the constant mu*t is defined (@ I M).                                    ')
      write(*,*)
      call remark(' To unlimit parameters use the new MULTI-FRILLS command ZERO. E.g. if you ')
      call remark(' want to clear the pmin and pmax values of the parameter 14 18 19 20, type')
      call remark(' ZERO 14 18-20.                                                           ')
      write(*,*)
      call remark(' You do not have to give a GO before if you want to calculate CEF properties')
      call remark(' only. FOCUS has now a build in calculation mode which is automatically')
      call remark(' called by FOCUS itself. E.g. @ D E will pick up the actual parameters of')
      call remark(' the parameter table and display the calculated excitations.              ')
      call prompt(' Hit a key to continue ')
      line_len = getlin(line)
      return
      end
!---------------------------------------------------------------------------------------------
! INFO about mu*t
!---------------------------------------------------------------------------------------------
      subroutine i_mut
      write(*,*)
      call remark(' ### The definition of mu*t ############################################# ')
      write(*,*)
      call remark(' Chemical formula of the sample : X1   X2   ... XN     (e.g. Ce  Fe   B ) ')
      call remark('                                    n1   n2       nN           2   14     ')
      write(*,*)
      call remark('                                                         N                ')
      call remark('                   density of sample         1          ---               ')
      call remark(' mu*t = t*N_A * ---------------------- * ------------ * >    ni * sigma(Xi)')
      call remark('                n1*A_X1 + .. + nN*A_XN   n1 + .. + nN   ---               ')
      call remark('                                                        i=1               ')
      call remark('                     t : effective thickness of the sample in cm          ')
      call remark(' Avogardo constant N_A = 0.622045 in particles per mol                    ')
      call remark(' density of the sample : in gr per cubic cm                               ')
      call remark('                  A_Xi : relative mass in gr per mol for chemical element i')
      call remark('             sigma(Xi) : cross section in barn for chemical element i     ')
      call remark('      mu*t (scattering): use the total scattering cross section for sigma ')
      call remark('      mu*t (absorption): use the absorption cross section for sigma       ')
      call remark('                         (for an incident neutron energy of 25.299 meV)   ')
      call remark('          unit of mu*t : 1                                                ')
      return
      end
!---------------------------------------------------------------------------------------------
! INFO about the used scattering function
!---------------------------------------------------------------------------------------------
      subroutine i_scattering_function
      write(*,*)
      call remark(' ### The used scattering function S(Q,E) for the elastic line ########### ')
      write(*,*)
      call remark('  S(QE) = intensity_factor/4pi * elastic_int * lineshape(E - elastic_pos) ')
      write(*,*)
      call remark('      elastic_int : elastic intensity. See D P.                           ')
      call remark('      elastic_pos : elastic position.  See D P.                           ')
      write(*,*)
      call remark(' ### The used scattering function S(Q,E) for the CEF excitations ######## ')
      write(*,*)
      call remark('                                 ----                                     ')
      call remark('  S(QE) = intensity_factor/4pi * >     I(i->k) * lineshape(E - [Ek-Ei])   ')
      call remark('                                 ----                                     ')
      call remark('                                  i k                                     ')
      write(*,*)
      call remark(' units of S(Q,E)  :  millibarn/meV/sr/R3+                                 ')
      call remark(' intensity_factor :  overall scale parameter. See D P.                    ')
      call remark('          Ek - Ei :  CEF transition energy for a transition i->k from     ')
      call remark('                     energy level Ei to energy level Ek. See @ D E.       ')
      call remark('          I(i->k) :  CEF transition intensity. See @ I IN.                ')
      call remark('     lineshape(E) :  FOCUS implemented lineshapes. See @ I L.             ')
      return
      end
!---------------------------------------------------------------------------------------------
! INFO about the implemented lineshapes
!---------------------------------------------------------------------------------------------
      subroutine i_lineshapes
      write(*,*)
      call remark(' ### FOCUS implemented lineshapes ####################################### ')
      call remark(' All implemented lineshapes have an intensity normalised to 1, i.e.       ')
      call remark(' integral(L(E),E=-infinity..infinity) = 1. Implemented functions are:     ')
      write(*,*)
      call remark(' 1 : Gaussian                                                             ')
      call remark(' 2 : Lorentzian                                                           ')
      call remark(' 3 : Lorentzian * Bose                                                    ')
      write(*,*)
      call remark(' Gaussian(E)   = 1/sqrt(2pi)/sigma * exp( - 0.5*[E/sigma]**2 )            ')
      call remark(' Lorentzian(E) = gamma/pi / ( gamma**2 + E**2 )                           ')
      call remark(' Bose(E,T)     = E/T / ( 1 - exp( -E/T ))                                 ')
      write(*,*)
      call remark('            E  : energy transfer in Kelvin                                ')
      call remark('            T  : sample temperature in Kelvin                             ')
      call remark('        sigma  = FWHM/sqrt(ln(256))                                       ')
      call remark('        gamma  = FWHM/2                                                   ')
      write(*,*)
      call remark(' Use @ D S to display what lineshapes the excitations have.               ')
      call remark(' Use @ s s to set the lineshapes.                                         ')
      return
      end
!---------------------------------------------------------------------------------------------
! INFO about the resolution function
!---------------------------------------------------------------------------------------------
      subroutine i_resolution
      write(*,*)
      call remark(' ### Resolution function definition ######################################')
      call remark(' The implemented resolution function R(E) is normalised to 1, i.e.     ')
      call remark(' the integral(R(E),E=-infinity..infinity) = 1.                            ')
      write(*,*)
      call remark(' S(Q,E)_R  =  Integral( S(Q,E-EE)*R(EE), EE=-infinity..infinity )         ')
      write(*,*)
      call remark(' The scattering function S(Q,E)_R is fitted to the data. If no resolution ')
      call remark(' function is calculated S(Q,E)_R is set to S(Q,E).                        ')
      write(*,*)
      call remark(' Type @ I S to print out how S(Q,E) is defined.                           ')
      call remark(' Type @ A R to change the parameters for the resolution function.         ')
      return
      end                                                                                                         
!---------------------------------------------------------------------------------------------
! INFO about the implemented intensity
!---------------------------------------------------------------------------------------------
      subroutine i_intensities
      write(*,*)
      call remark(' ### The magnetic CEF transition intensities #############################')
      call remark(' The magnetic CEF transition intensities for a powder                     ')
      write(*,*)
      call remark('                                  2   ----   exp(-Ei/T)             2     ')
      call remark('  I(i->k) = 4*pi*[0.5*r0*g_J*f(Q)] *  >    ------------ * |<i|jt|k>|    ')
      call remark('                                      ----     Z(T)                     ')
      call remark('                                       i k                              ')
      write(*,*)
      call remark('           I(i->k) : CEF transition intensity i->k in barn.               ')
      call remark('               r0  : magnetic scattering length (r0**2 = 0.290 barn).     ')
      call remark('               g_J : Lande factor of R3+ see @ D IO.                      ')
      call remark('              f(Q) : magnetic form factor of R3+.                         ')
      call remark('                Ei : energy of CEF level i.                               ')
      call remark('              Z(T) : occupation factor for sample temperature T.          ')
      call remark('     |<i|jt|k>|**2 : powder matrix element for a CEF transition i->k :    ')
      call remark('     |<i|jt|k>|**2 = 2/3 * (|<i|jx|k>|**2 + |<i|jy|k>|**2 +  |<i|jz|k>|**2)')
      write(*,*)
      call remark(' Type @ D IN to print out what the CEF transition intensities are. ')
      call remark(' Type @ D MA to print out what the CEF transition matrixelements are.')
      return
      end
!---------------------------------------------------------------------------------------------
! INFO about the implemented intensity
!---------------------------------------------------------------------------------------------
      subroutine i_authors
      write(*,*)
      call remark(' ### The origin of some FOCUS libraries ##################################')
      call remark(' FRILLS: an Interactive Least-Squares Fitting Package                     ')
      call remark('        Author: Ray Osborn                                                ')
      call remark('        Manual: Rutherford Appleton Laboratory report RAL-91-011          ')
      call remark('                January 1991                                              ')
      call remark('        Note  : Fits one genie workspace by a user defined function       ')
      call remark(' MULTI_FRILLS: an Interactive Least-Squares Fitting Package               ')
      call remark('        Author: Ray Osborn and Toby G Perring                            ')
      call remark('        Manual: Rutherford Appleton Laboratory preprint                  ')
      call remark('                March 1995                                               ')
      call remark('        Note  : Fits multiple genie workspaces by a user defined function')
      call remark(' RESOLUTION FUNCTION:                                                    ')
      call remark('        Author: Ray Osborn and Eugen Goremychkin                         ')
      call remark(' SELFSHIELDING AND ABSORPTION:                                           ')
      call remark('        Author: Ray Osborn and Eugen Goremychkin                         ')
      call remark(' FOCUS: an Interactive Crystal Electric Field Parameter Fitting          ')
      call remark('        Package Using Neutron Scattering Data                            ')
      call remark('        Author: Peter Fabi                                               ')
      call remark('        Manual: Council for the Central Laboratory of the Research       ')
      call remark('                Councils, Technical Report RAL-TR-95-023, June 1995      ')
      call remark('        Note  : The manual mentioned above, descripes FOCUS Version 1.0  ')
      return
      end
