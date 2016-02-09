$ set verify
$  define/nolog focus_devel    het$disk:[hetmgr.library.focus.devel]
$  define/nolog focus_ref      het$disk:[hetmgr.library.focus.ref]
$  define/nolog focus_start    het$disk:[hetmgr.genie]
$!
$!-------------------------------------------
$! Release a new version of focus to public                    
$!-------------------------------------------
$ IF f$getsyi("arch_type") .EQ. 1
$ THEN
$!
$!    releasing on vax
$!
$!    change automatically the version number 
$!
$ for  focus_ref:version_change.for /obj=sys$scratch:version_change.obj
$ link sys$scratch:version_change,-
           het_library:MULTI_FRILLS/LIB,-
           het_library:PEAKS/LIB,-
           het_library:G_V2/LIB,-
           het_library:IO/LIB,-
           GENIE_LIB/LIB,-
           NAGF/OPT,-
           SYS$LIBRARY:GKSFORBND/LIB
$ run   version_change
$ purge focus_ref:focus.ver
$ copy  focus_ref:focus.ver focus_devel:
$ purge focus_devel:focus.ver
$ delete version_change.exe;*
$ delete sys$scratch:version_change.obj;*
$ purge   focus_ref:*.*
$ copy  focus_devel:*.com focus_ref:
$ purge   focus_ref:*.com
$ copy  focus_devel:*.for focus_ref:
$ purge   focus_ref:*.for
$ copy  focus_devel:*.inc focus_ref:
$ purge   focus_ref:*.inc
$ copy  focus_devel:*.hlp focus_ref:
$ purge   focus_ref:*.hlp
$ copy  focus_devel:*.par focus_ref:
$ purge   focus_ref:*.par
$ delete  focus_ref_exe:focus.exe;*
$ copy  focus_devel_exe:focus.exe focus_ref_exe:
$ copy    focus_ref:focus.com focus_start:
$ purge focus_start:focus.com
$ copy    focus_ref:focus_init.com focus_start:
$ purge focus_start:focus_init.com
$ else
$!
$! releasing on alpha
$! 
$ copy    focus_devel_exe:focus.exe ss:
$ delete  focus_ref_exe:focus.exe;*
$ delete  focus_devel_exe:focus.exe;*
$ copy    ss:focus.exe focus_ref_exe:
$ delete  ss:focus.exe;*
$ endif
$!*****************************************************
$! NEW VERSION OF FOCUS RELEASED 
$! ALL RIGHTS RESERVED BY 
$! P O FABI,  R OSBORN AND T G PERRING
$!*****************************************************
$ typ focus_ref:focus.ver
$!*****************************************************
$ set noverify
$ exit
