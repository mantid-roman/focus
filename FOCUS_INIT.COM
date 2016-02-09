$ set noverify
!
!       -------------------------------------------------- 
!      |          Initialisation command file for         |
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
$  IF f$getsyi("arch_type") .EQ. 1
$  THEN
$       define/nolog focus_ref_exe   het$disk:[hetmgr.library.focus.ref.vax]
$  ELSE
$       define/nolog focus_ref_exe   het$disk:[hetmgr.library.focus.ref.alpha]
$  ENDIF   
$  write sys$output "FOCUS initialisation done." 
$  exit
