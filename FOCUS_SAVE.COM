$ set verify
$ define/nolog focus_devel    het$disk:[hetmgr.library.focus.devel] 
$ define/nolog focus_ref      het$disk:[hetmgr.library.focus.ref]
$ define/nolog focus_start    het$disk:[hetmgr.genie] 
$
$ delete/noconfirm focus_devel:*.for;*
$ copy user$disk:[pf01.focus]*.com focus_devel:
$ purge focus_devel:*.com
$ copy user$disk:[pf01.focus]*.for focus_devel:
$ purge focus_devel:*.for
$ copy user$disk:[pf01.focus]*.inc focus_devel:
$ purge focus_devel:*.inc
$ copy user$disk:[pf01.focus]*.hlp focus_devel:
$ purge focus_devel:*.hlp
$ copy user$disk:[pf01.focus]*.par focus_devel:
$ purge focus_devel:*.par
$
$ set noverify
$ exit
