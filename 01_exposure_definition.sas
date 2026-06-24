/*==============================================================================
  01. Exposure definition: heat-wave and cold-spell indicators
  - District- and year-specific percentile thresholds from daily mean
    temperature (tem), computed separately for each district x year:
      warm season (May-Sep): 95th / 97.5th / 99th
      cold season (Nov-Feb): 1st / 2.5th / 5th
  - Single-day and 2-/3-consecutive-day indicators, computed within district.
  Input:  b.<dx>_cco_merged  (case-crossover days with temperature `tem`,
          district `RVSN_ADDR_CD`, calendar `date`, `year`)
  Output: b.<dx>_summer_flags, b.<dx>_winter_flags
  NOTE: NHIS data not included; server paths and identifiers removed.
==============================================================================*/

%macro make_hw_cw(dxlist);
  %local i thisdx;
  %let n = %sysfunc(countw(&dxlist));

  %do i = 1 %to &n;
    %let thisdx = %scan(&dxlist, &i);

    /* Split warm (May-Sep) and cold (Nov-Feb) seasons */
    data work.&thisdx._summer work.&thisdx._winter;
      set b.&thisdx._cco_merged;
      month = month(date);
      if month in (5,6,7,8,9)  then output work.&thisdx._summer;
      if month in (11,12,1,2)  then output work.&thisdx._winter;
    run;

    /* District x year warm-season percentiles (95 / 97.5 / 99th) */
    proc sort data=work.&thisdx._summer; by RVSN_ADDR_CD year; run;
    proc univariate data=work.&thisdx._summer noprint;
      by RVSN_ADDR_CD year;
      var tem;
      output out=&thisdx._quant_summer pctlpre=P_ pctlpts=95 97.5 99;
    run;

    /* District x year cold-season percentiles (1 / 2.5 / 5th) */
    proc sort data=work.&thisdx._winter; by RVSN_ADDR_CD year; run;
    proc univariate data=work.&thisdx._winter noprint;
      by RVSN_ADDR_CD year;
      var tem;
      output out=&thisdx._quant_winter pctlpre=P_ pctlpts=1 2.5 5;
    run;

    /* Attach district x year thresholds to each warm-season day */
    data work.&thisdx._summer_thr;
      merge work.&thisdx._summer(in=a) &thisdx._quant_summer;
      by RVSN_ADDR_CD year;
      if a;
      heat95  = (tem >= P_95);
      heat975 = (tem >= P_97_5);
      heat99  = (tem >= P_99);
    run;

    /* Heat-wave indicators: single day + 2/3 consecutive days, within district */
    proc sort data=work.&thisdx._summer_thr; by RVSN_ADDR_CD date; run;
    data b.&thisdx._summer_flags;
      set work.&thisdx._summer_thr;
      by RVSN_ADDR_CD;
      retain hw95_2 hw95_3 hw97_2 hw97_3 hw99_2 hw99_3;
      if first.RVSN_ADDR_CD then do;
        hw95_2=0; hw95_3=0; hw97_2=0; hw97_3=0; hw99_2=0; hw99_3=0;
      end;
      hw95_2 = (lag1(heat95)=1  and heat95=1);
      hw95_3 = (lag1(heat95)=1  and lag2(heat95)=1  and heat95=1);
      hw97_2 = (lag1(heat975)=1 and heat975=1);
      hw97_3 = (lag1(heat975)=1 and lag2(heat975)=1 and heat975=1);
      hw99_2 = (lag1(heat99)=1  and heat99=1);
      hw99_3 = (lag1(heat99)=1  and lag2(heat99)=1  and heat99=1);
    run;

    /* Attach district x year thresholds to each cold-season day */
    data work.&thisdx._winter_thr;
      merge work.&thisdx._winter(in=a) &thisdx._quant_winter;
      by RVSN_ADDR_CD year;
      if a;
      cold1   = (tem <= P_1);
      cold2p5 = (tem <= P_2_5);
      cold5   = (tem <= P_5);
    run;

    /* Cold-spell indicators: single day + 2/3 consecutive days, within district */
    proc sort data=work.&thisdx._winter_thr; by RVSN_ADDR_CD date; run;
    data b.&thisdx._winter_flags;
      set work.&thisdx._winter_thr;
      by RVSN_ADDR_CD;
      retain cw1_2 cw1_3 cw2_2 cw2_3 cw5_2 cw5_3;
      if first.RVSN_ADDR_CD then do;
        cw1_2=0; cw1_3=0; cw2_2=0; cw2_3=0; cw5_2=0; cw5_3=0;
      end;
      cw1_2 = (lag1(cold1)=1   and cold1=1);
      cw1_3 = (lag1(cold1)=1   and lag2(cold1)=1   and cold1=1);
      cw2_2 = (lag1(cold2p5)=1 and cold2p5=1);
      cw2_3 = (lag1(cold2p5)=1 and lag2(cold2p5)=1 and cold2p5=1);
      cw5_2 = (lag1(cold5)=1   and cold5=1);
      cw5_3 = (lag1(cold5)=1   and lag2(cold5)=1   and cold5=1);
    run;
  %end;
%mend make_hw_cw;
