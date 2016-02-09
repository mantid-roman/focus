$ set verify
$ define/nolog frill_sources het$disk:[hetmgr.library.multi_frills.ref]
$ define/nolog genie_sources het$disk:[hetmgr.gsource]
$ define/nolog focus_devel   het$disk:[hetmgr.library.focus.devel] 
$ define/nolog focus_ref     het$disk:[hetmgr.library.focus.ref]
$  IF f$getsyi("arch_type") .EQ. 1
$  THEN
$       define/nolog focus_ref_exe   het$disk:[hetmgr.library.focus.ref.vax]
$       define/nolog focus_devel_exe het$disk:[hetmgr.library.focus.devel.vax]
$  ELSE
$       define/nolog focus_ref_exe   het$disk:[hetmgr.library.focus.ref.alpha]
$       define/nolog focus_devel_exe het$disk:[hetmgr.library.focus.devel.alpha]
$  ENDIF                                                                     
$ define/nolog focus_start   het$disk:[hetmgr.genie] 
$ define/nolog focus_sources focus_devel:
$ define/nolog cf_sources    focus_sources:
$!
$! - compile the source code on focus_devel:
$! - the object files are stored at sys$scratch:
$!
$ for focus_sources:3j6j_8.for         /obj=sys$scratch:3j6j_8.obj
$ for focus_sources:bndpar.for         /obj=sys$scratch:bndpar.obj
$ for focus_sources:clepar.for         /obj=sys$scratch:clepar.obj
$ for focus_sources:cf_fabi.for        /obj=sys$scratch:cf_fabi.obj
$ for focus_sources:cf_formfactor.for  /obj=sys$scratch:cf_formfactor.obj
$ for focus_sources:dermat.for         /obj=sys$scratch:dermat.obj
$ for focus_sources:fitdis.for         /obj=sys$scratch:fitdis.obj
$ for focus_sources:focus.for          /obj=sys$scratch:focus.obj  
$ for focus_sources:focus_calc.for     /obj=sys$scratch:focus_calc.obj  
$ for focus_sources:function.for       /obj=sys$scratch:function.obj
$ for focus_sources:help.for           /obj=sys$scratch:help.obj
$ for focus_sources:info.for           /obj=sys$scratch:info.obj
$ for focus_sources:limpar.for         /obj=sys$scratch:limpar.obj
$ for focus_sources:monte_carlo.for    /obj=sys$scratch:monte_carlo.obj
$ for focus_sources:multi_frills.for   /obj=sys$scratch:multi_frills.obj
$ for focus_sources:parin.for          /obj=sys$scratch:parin.obj
$ for focus_sources:parout.for         /obj=sys$scratch:parout.obj
$ for focus_sources:pardis.for         /obj=sys$scratch:pardis.obj
$ for focus_sources:pback.for          /obj=sys$scratch:pback.obj
$ for focus_sources:pmfit.for          /obj=sys$scratch:pmfit.obj
$ for focus_sources:ptrans.for         /obj=sys$scratch:ptrans.obj
$ for focus_sources:resolution.for     /obj=sys$scratch:resolution.obj   
$ for focus_sources:self.for           /obj=sys$scratch:self.obj          
$ for focus_sources:setpar.for         /obj=sys$scratch:setpar.obj
$ for focus_sources:unlimpar.for       /obj=sys$scratch:unlimpar.obj
$ for focus_sources:user_command.for   /obj=sys$scratch:user_command.obj
$ for focus_sources:user_help.for      /obj=sys$scratch:user_help.obj
$ for focus_sources:version_number.for /obj=sys$scratch:version_number.obj
$ for focus_sources:show.for           /obj=sys$scratch:show.obj
$!
$ purge sys$scratch:*.obj
$!
$ set noverify
$ exit
