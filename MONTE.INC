c file  MONTE.INC
c       --------------------------------------------------------------------        
c       FOR MONTE CARLO
        common /MONTE/ s_chisq(MAX_CHI2), s_ratio(MAX_CHI2), s_p(MAX_CHI2,MAX_PAR),
     *                 s_pmin(MAX_CHI2,MAX_PAR),s_pmax(MAX_CHI2,MAX_PAR),
     *                 if_monte,n_mc_par, if_plot, 
     *                 imc_array(MAX_CHI2), nset
        real*8 s_chisq, s_ratio, s_p, s_pmin, s_pmax
        integer if_monte, n_mc_par, if_plot
        integer*4  imc_array, nset
c       ---------------------------------------------------------------------        
        common /BOUNDARY/ min_x, max_x, min_y, max_y, if_first
        real*4  min_x, max_x, min_y, max_y
        integer if_first
c       --------------------------------------------------------------------        
