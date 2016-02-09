$ set verify
$ IF f$getsyi("arch_type") .EQ. 1
$ THEN
$!linking on vax
$ link  sys$scratch:focus,-
        sys$scratch:bndpar,-
        sys$scratch:cf_fabi,-
        sys$scratch:cf_formfactor,-
        sys$scratch:clepar,-
        sys$scratch:dermat,-
        sys$scratch:fitdis,-
        sys$scratch:focus_calc,-
        sys$scratch:function,-
        sys$scratch:help,-
        sys$scratch:info,-
        sys$scratch:limpar,-
        sys$scratch:monte_carlo,-
        sys$scratch:multi_frills,-
        sys$scratch:parin,-
        sys$scratch:parout,-
        sys$scratch:pardis,-
        sys$scratch:pback,-
        sys$scratch:pmfit,-
        sys$scratch:ptrans,-
        sys$scratch:resolution,-
        sys$scratch:self,-
        sys$scratch:setpar,-
        sys$scratch:unlimpar,-
        sys$scratch:user_command,-
        sys$scratch:user_help,-
        sys$scratch:version_number,-
        sys$scratch:show,-
        sys$scratch:3j6j_8,-
        het_library:multi_frills/lib,-
        het_library:peaks/lib,-
        het_library:g_v2/lib,-
        het_library:io/lib,-
        genie_lib/lib,-
        nagf/opt,-
        sys$library:gksforbnd/lib,-
        sys$input/opt
        ISD_MAX=128 
$ else
$!linking on alpha
$ link  sys$scratch:focus,-
        sys$scratch:bndpar,-
        sys$scratch:cf_fabi,-
        sys$scratch:cf_formfactor,-
        sys$scratch:clepar,-
        sys$scratch:dermat,-
        sys$scratch:fitdis,-
        sys$scratch:focus_calc,-
        sys$scratch:function,-
        sys$scratch:help,-
        sys$scratch:info,-
        sys$scratch:limpar,-
        sys$scratch:monte_carlo,-
        sys$scratch:multi_frills,-
        sys$scratch:parin,-
        sys$scratch:parout,-
        sys$scratch:pardis,-
        sys$scratch:pback,-
        sys$scratch:pmfit,-
        sys$scratch:ptrans,-
        sys$scratch:resolution,-
        sys$scratch:self,-
        sys$scratch:setpar,-
        sys$scratch:unlimpar,-
        sys$scratch:user_command,-
        sys$scratch:user_help,-
        sys$scratch:version_number,-
        sys$scratch:show,-
        sys$scratch:3j6j_8,-
        het_library:multi_frills/lib,-
        het_library:peaks/lib,-
        het_library:g_v2/lib,-
        het_library:io/lib,-
        genie_lib/lib,-
        nagf/opt,-
        sys$library:gksforbnd/lib
$ endif
$ delete focus_devel_exe:focus.exe;*
$ copy   focus.exe  focus_devel_exe:focus.exe
$ delete focus.exe;*
$ copy   focus_devel:focus_devel.com focus_start:
$ purge  focus_start:focus_devel.com
$ copy   focus_devel:focus_init.com focus_start:
$ purge  focus_start:focus_init.com
$ delete sys$scratch:*.obj;*
$ set noverify
$ exit
