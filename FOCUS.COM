!
!       -------------------------------------------------- 
!      |             Start up command file for            |
!      |         *****   ***    ***   *   *   ***         |     
!      |         *      *   *  *      *   *  *            |   
!      |         ***    *   *  *      *   *   ***         |   
!      |         *      *   *  *      *   *      *        |      
!      |         *       ***    ***    ***    ***         |
!      |      An Interactive Crystal Electric Field       |
!      |         Parameter Fitting Package using          |
!      |             Neutron Scattering Data              |
!      |                   Peter Fabi                     |
!      |  ISIS-Facility, Rutherford Appleton Laboratory   |
!      |  telephone: +44 1235 445428                      |
!      |  faxsimile: +44 1235 445383                      |
!      |                                                  |
!       --------------------------------------------------              
! ---------------------------------------------------------------
!| FOCUS Version for single workspace fit only                   |
! ---------------------------------------------------------------
! Please give in  first  the workspace number where the crystal 
! electric field data are you want to fit (e.g. enter 1 for w1). 
! Then enter a workspace number for the FOCUS output. 
$ INQUIRE n_in    " Enter the workspace number of neutron data "
$ INQUIRE n_out   " Enter the workspace number of FOCUS output "
!
$! ---------------------------------------------------------------
$!| Write data to GENIE intermediate file                         |
$!| (to be read by multi_frills program FOCUS)                    |
$! ---------------------------------------------------------------
> write w'n_in' ss:frills.in
> write/close
$! ---------------------------------------------------------------
$!| invoke the MULTI_FRILLS program FOCUS into GENIE              |
$! ---------------------------------------------------------------
> transform w'n_in' focus_ref_exe:focus  w'n_out'
$! ---------------------------------------------------------------
$!| read output from MULTI_FRILLS program into output workspaces  |
$! ---------------------------------------------------------------
> read w'n_out' ss:frills.out
$ exit
