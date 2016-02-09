	include 'FRILLS.PAR'

	character nam*20, message*50, fit_time*20, xlab*5, ylab*5, xcap*40, ycap*40, title*80, subtitle*60
	integer*4 nptot, np, inp, inr, nf, inf, inrf, ns, is, ins, nd, nv, inv, nc, status, max_iter, max_div, output
	integer*4 mk_type, mk_coli, ln_type, ln_coli, n_plt, smooth, hard_type, emu_type, wk_type
	real*4 mk_size, ln_size
	real*8 p, pmin, pmax, rp, v, vsig, x, y, ysig, yc, pf, pfsig, pfmin, pfmax, rf, chisq, dc, dv, dd, H, dp, acc

	common /fr_nams/ nptot, nam(MAX_PAR)
	common /fr_pars/ np, inp(MAX_PAR),  p(MAX_PAR),  pmin(MAX_PAR), pmax(MAX_PAR), inr(MAX_PAR), rp(MAX_PAR)
	common /fr_fits/ nf, inf(MAX_PAR), pf(MAX_PAR), pfsig(MAX_PAR), pfmin(MAX_PAR), pfmax(MAX_PAR), inrf(MAX_PAR), rf(MAX_PAR)
	common /fr_sets/ ns, is, ins(MAX_PAR)
	common /fr_dats/ nd(MAX_SET),   x(MAX_DAT,MAX_SET), y(MAX_DAT,MAX_SET), ysig(MAX_DAT,MAX_SET)
	common /fr_vals/ nv(MAX_SET), inv(MAX_DAT,MAX_SET), v(MAX_DAT,MAX_SET), vsig(MAX_DAT,MAX_SET)
	common /fr_cmps/ nc, yc(MAX_DAT,MAX_CMP)
	common /fr_chis/ chisq, status, message, fit_time
	common /fr_mats/ dc(MAX_SET*MAX_DAT), dv(MAX_PAR,MAX_SET*MAX_DAT), dd(MAX_PAR), H(MAX_PAR,MAX_PAR)
	common /fr_labs/ xlab(MAX_SET), ylab(MAX_SET), xcap(MAX_SET), ycap(MAX_SET), title(MAX_SET), subtitle(3,MAX_SET)
	common /fr_cons/ dp, acc, max_iter, max_div, output
	common /fr_plts/ mk_type, mk_size, mk_coli, ln_type, ln_size, ln_coli, n_plt, smooth
	common /fr_wrks/ hard_type, emu_type, wk_type
