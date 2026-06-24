/*==============================================================================
  03. Subgroup analyses (sex, age, income, urbanicity)
  - The main model (02) is re-fitted on case-crossover datasets restricted to
    each subgroup level (e.g., Male / Female; <65 / >=65; income; metro / non-metro).
  - Effect modification was assessed with Wald tests comparing subgroup-specific
    log-odds-ratios (Altman & Bland, BMJ 2003):
        z = (b1 - b2) / sqrt(se1^2 + se2^2)
  Input:  b.<dx>_flag_<subgroup>   (case-crossover data restricted to one stratum)
==============================================================================*/

%macro run_subgroup(dxlist, subgroup);
  %local i thisdx;
  %let n = %sysfunc(countw(&dxlist));

  %let expovars = heat95 hw95_2 hw95_3 heat975 hw97_2 hw97_3 heat99 hw99_2 hw99_3
                  cold1  cw1_2  cw1_3  cold2p5 cw2_2  cw2_3  cold5  cw5_2  cw5_3;

  %do i = 1 %to &n;
    %let thisdx = %scan(&dxlist, &i);
    %do e = 1 %to %sysfunc(countw(&expovars));
      %let ev = %scan(&expovars, &e);
      ods output OddsRatios = or_&ev;
      proc logistic data=b.&thisdx._flag_&subgroup;
        class &ev (ref='0') / param=reference;
        effect spl_hum = spline(hum / basis=bspline degree=3);
        model case(event='1') = &ev spl_hum;
        strata stratum;
      run;
      data or_&ev; set or_&ev; Disease="&thisdx"; Subgroup="&subgroup"; Cutoff="&ev"; run;
    %end;
    data &thisdx._OR_&subgroup; set or_:; run;
    proc datasets lib=work nolist; delete or_:; quit;
  %end;
%mend run_subgroup;

/* Example:
   %run_subgroup(&dxlist, Male);
   %run_subgroup(&dxlist, Female);
*/
