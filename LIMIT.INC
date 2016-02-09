! FILE LIMIT.INC

      common /LIMITS/ pmikq(  0:6,0:6),pmakq( 0:6,0:6),
     *                pmiokq( 0:6,0:6),pmaokq(0:6,0:6),
     *                indokq( 0:6,0:6),fixokq(0:6,0:6),
     *                knrkq(  0:6,0:6),kkq(   0:6,0:6),
     *                knpkq(  0:6,0:6),rpkq(  0:6,0:6),
     *             cf_pmin(MAX_PAR)    , cf_pmax(MAX_PAR),
     *             cf_pminold(MAX_PAR) , cf_pmaxold(MAX_PAR),
     *             cf_indvold(MAX_PAR) , cf_fixold(MAX_PAR),
     *             cf_knr(MAX_PAR), cf_rp(MAX_PAR), cf_k(MAX_PAR),
     *             mf_knr(MAX_PAR), mf_rp(MAX_PAR), mf_k(MAX_PAR),
     *             cf_knp(MAX_PAR)     , mf_knp(MAX_PAR),
     *             mf_pminold(MAX_PAR) , mf_pmaxold(MAX_PAR),
     *             mf_indvold(MAX_PAR) , mf_fixold(MAX_PAR)
      complex*16   pmikq,pmakq,pmiokq,pmaokq,indokq,fixokq,knpkq,knrkq,rpkq,kkq
      real*8       cf_pmin    ,cf_pmax 
      real*8       cf_pminold ,cf_pmaxold
      integer      cf_indvold ,cf_fixold 
      real*8       mf_pminold ,mf_pmaxold
      integer      mf_indvold ,mf_fixold 
      integer*4    cf_knr  
      integer*4    mf_knr 
      integer*4    cf_knp     ,mf_knp 
      integer*4    cf_k       ,mf_k
      real*8       cf_rp      ,mf_rp 
