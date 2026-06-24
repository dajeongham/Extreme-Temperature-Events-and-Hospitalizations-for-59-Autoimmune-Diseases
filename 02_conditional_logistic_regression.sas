/*==============================================================================
  02. Main analysis: time-stratified case-crossover, conditional logistic
  - Conditional logistic regression (PROC LOGISTIC with STRATA), one model per
    disease x exposure definition.
  - Daily mean relative humidity adjusted with a cubic B-spline (df = 3).
  - Case and matched control days share a stratum (residential district,
    calendar year, month, day of week).
  Input:  b.<dx>_flag_all  (case-crossover days with exposure indicators from 01,
          `case` (1/0), matching `stratum`, humidity `hum`)
  Output: <dx>_OR  (odds ratios across all exposure definitions)
==============================================================================*/

%macro run_regression(dxlist);
  %local i thisdx;
  %let n = %sysfunc(countw(&dxlist));

  %let heatvars = heat95 hw95_2 hw95_3 heat975 hw97_2 hw97_3 heat99 hw99_2 hw99_3;
  %let coldvars = cold1  cw1_2  cw1_3  cold2p5 cw2_2  cw2_3  cold5  cw5_2  cw5_3;

  %do i = 1 %to &n;
    %let thisdx  = %scan(&dxlist, &i);
    %let dataset = b.&thisdx._flag_all;

    /* Heat-wave models */
    %do h = 1 %to %sysfunc(countw(&heatvars));
      %let hv = %scan(&heatvars, &h);
      ods output OddsRatios = or_&hv;
      proc logistic data=&dataset;
        class &hv (ref='0') / param=reference;
        effect spl_hum = spline(hum / basis=bspline degree=3);
        model case(event='1') = &hv spl_hum;
        strata stratum;
      run;
      data or_&hv; set or_&hv; Disease="&thisdx"; Cutoff="&hv"; run;
    %end;

    /* Cold-spell models */
    %do c = 1 %to %sysfunc(countw(&coldvars));
      %let cv = %scan(&coldvars, &c);
      ods output OddsRatios = or_&cv;
      proc logistic data=&dataset;
        class &cv (ref='0') / param=reference;
        effect spl_hum = spline(hum / basis=bspline degree=3);
        model case(event='1') = &cv spl_hum;
        strata stratum;
      run;
      data or_&cv; set or_&cv; Disease="&thisdx"; Cutoff="&cv"; run;
    %end;

    data &thisdx._OR; set or_:; run;
    proc datasets lib=work nolist; delete or_:; quit;
  %end;
%mend run_regression;

/* Example:
   %let dxlist = IBD CD UC RA SLE ... ;   * all 59 autoimmune diseases ;
   %run_regression(&dxlist);
*/
