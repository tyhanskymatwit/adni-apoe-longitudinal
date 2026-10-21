/*==============================================================================
  Input : adni_wide.csv one row per participant; 
  Output: nominal years 0..6 in columns, adas13_0..adas13_6 (also mmse_*, cdrsb_*), visit time t_0..t_6
==============================================================================*/

%let root    = /home/<your-sas-user>/ADNI;   /* EDIT: project folder */
%let outcome = adas13;   /* adas13 (higher=worse) | mmse (higher=better) | cdrsb */
%let window  = 0.35;

libname adni "&root/output";

/* save the full log and all printed results + figures to output*/
proc printto log="&root/output/sas_adni_longitudinal.log" new; run;
ods pdf file="&root/output/sas_adni_results.pdf" style=journal;
options nofmterr ls=140 nocenter mprint;
ods graphics on / reset imagename="sas_fig" imagefmt=png width=7in;
ods listing gpath="&root/output";

proc format;
  value dxf     1='CN' 2='MCI' 3='Dementia';
  value apof    0='0 e4' 1='1 e4' 2='2 e4';
  value sexf    0='Male' 1='Female';
  value dgf     0='Retained/admin-censored' 1='Dropout by yr 2' 2='Dropout yr 3-5';
run;

/* IMPORT */
proc import datafile="&root/data/adni_wide.csv" out=work.raw dbms=csv replace;
  guessingrows=max;
  getnames=yes;
run;

data work.wide;
  set work.raw;
  /* baseline diagnosis (numeric so class order = CN, MCI, Dementia) */
  if      dx_bl = 'CN'  then dx = 1;
  else if dx_bl = 'MCI' then dx = 2;
  else if dx_bl = 'DEM' then dx = 3;
  if not missing(apoe4) then apoe4_carrier = (apoe4 >= 1);
  age_c  = age_bl - 72;          /* centred at 72 y */
  educ_c = educ - 16;            /* centred at 16 y */
  if nmiss(age_bl, female, educ, dx, apoe4, &outcome._0) = 0;   /* sample */
  format dx dxf. apoe4 apof. female sexf.;
  label age_c='Baseline age (centred 72)' educ_c='Education (centred 16)';
run;

title 'Analysis sample';
proc freq data=work.wide;
  tables dx*apoe4 / nocol nopercent;
  tables protocol female;
run;
proc means data=work.wide n mean std maxdec=2;
  class apoe4;
  var age_bl educ &outcome._0;
run;

/* MISSING-DATA PATTERNS AND DROPOUT DEFINITION */
data work.wide;
  set work.wide;
  array yy{0:6} &outcome._0 - &outcome._6;
  array el{0:6} elig_0 - elig_6;
  lastobs = 0; lastelig = 0; nobs = 0;
  do k = 0 to 6;
    el{k} = (fu_potential >= k + &window) or not missing(yy{k});
    if not missing(yy{k}) then do; lastobs = k; nobs + 1; end;
    if el{k} then lastelig = k;
  end;
  dropout = (lastobs < lastelig);
  if      dropout = 0  then dropgrp = 0;
  else if lastobs <= 2 then dropgrp = 1;
  else                      dropgrp = 2;
  format dropgrp dgf.;
  drop k;
run;

title 'Missing-data patterns of the outcome (wide format)';
proc mi data=work.wide nimpute=0 simple;
  var &outcome._0 - &outcome._6;
run;

title 'Retention among administratively eligible participants';
data work.retention;
  set work.wide;
  array yy{0:6} &outcome._0 - &outcome._6;
  array el{0:6} elig_0 - elig_6;
  do year = 0 to 6;
    if el{year} then do; observed = not missing(yy{year}); output; end;
  end;
  keep rid year observed;
run;
proc means data=work.retention n mean maxdec=3;
  class year;
  var observed;
run;

title 'Who drops out? (differences => not MCAR)';
proc freq data=work.wide;
  tables dropgrp*died / nocol nopercent;
  tables (dx apoe4)*dropout / chisq nocol nopercent;
run;
proc ttest data=work.wide;
  class dropout;
  var age_bl &outcome._0;
run;
proc logistic data=work.wide;
  class apoe4(ref='0 e4') dx(ref='CN') female(ref='Male') protocol / param=ref;
  model dropout(event='1') = age_c female educ_c apoe4 dx &outcome._0 protocol;
run;

/* WIDE TO LONG WITH A DATA STEP OUTPUT LOOP */
data adni.long_all;
  set work.wide;
  array tt{0:6} t_0      - t_6;
  array yy{0:6} &outcome._0 - &outcome._6;
  array mm{0:6} mmse_0   - mmse_6;
  array el{0:6} elig_0   - elig_6;
  do year = 0 to 6;
    time   = tt{year};
    y      = yy{year};
    mmse   = mm{year};
    elig   = el{year};
    output;
  end;
  keep rid year time y mmse elig age_bl age_c female educ educ_c dx apoe4
       apoe4_carrier protocol died fu_potential dropout dropgrp lastobs;
  label y="&outcome" time='Years since baseline (actual)'
        year='Nominal annual visit';
run;

proc sort data=adni.long_all; by rid year; run;

/* history variables: previous observed score, last change, baseline score */
data adni.long_all;
  set adni.long_all;
  by rid;
  retain y0;
  yl1 = lag1(y);
  yl2 = lag2(y);
  if first.rid then y0 = y;
  y_prev  = ifn(year >= 1, yl1, .);
  dy_prev = ifn(year >= 2, yl1 - yl2, .);
  miss    = missing(y);
  atrisk  = (year >= 1 and not missing(y_prev) and elig = 1);
  drop yl1 yl2;
run;

/* analysis file: observed outcomes only */
data adni.long_obs;
  set adni.long_all;
  where not missing(y) and not missing(time);
run;

proc sql;
  title 'Long file';
  select count(distinct rid) as subjects, count(*) as observations
  from adni.long_obs;
quit;

/* MISSING-DATA MECHANISM: MCAR vs MAR */
title 'Missingness model A: baseline covariates';
proc logistic data=adni.long_all(where=(atrisk=1));
  class apoe4(ref='0 e4') dx(ref='CN') female(ref='Male') protocol year
        / param=ref;
  model miss(event='1') = age_c female educ_c apoe4 dx y0 protocol year;
  ods output FitStatistics=work.missA_fit;
run;

title 'Missingness model B: + most recent observed score';
proc logistic data=adni.long_all(where=(atrisk=1));
  class apoe4(ref='0 e4') dx(ref='CN') female(ref='Male') protocol year
        / param=ref;
  model miss(event='1') = y_prev age_c female educ_c apoe4 dx y0 protocol year;
  ods output FitStatistics=work.missB_fit;
run;

/* likelihood-ratio test A vs B (same at-risk rows) */
data work.lrt_miss;
  merge work.missA_fit(where=(criterion='-2 Log L') keep=criterion interceptandcovariates
                       rename=(interceptandcovariates=m2ll_A))
        work.missB_fit(where=(criterion='-2 Log L') keep=criterion interceptandcovariates
                       rename=(interceptandcovariates=m2ll_B));
  chisq = m2ll_A - m2ll_B;
  p     = 1 - probchi(chisq, 1);
run;
title 'LR test: does the last observed score predict next-visit missingness?';
proc print data=work.lrt_miss noobs; run;

/* PROC GLIMMIX: same model with a random intercept for person, so repeated
   at-risk rows of the same participant are not treated as independent */
title 'Missingness GLMM (random intercept), Laplace';
proc glimmix data=adni.long_all(where=(atrisk=1)) method=laplace;
  class rid apoe4 dx female protocol year;
  model miss(event='1') = y_prev age_c female educ_c apoe4 dx y0 protocol year
        / dist=binary link=logit solution oddsratio;
  random intercept / subject=rid;
  ods output ParameterEstimates=adni.glmm_missingness;
run;

/* LINEAR MIXED-EFFECTS MODELS (PROC MIXED) */
/* macro: observed-mean weights for ESTIMATE statements (= Stata margins) */
%macro covmeans(data=, where=1);
  %global m_agec m_educc p_dx1 p_dx2 p_dx3 p_f0 p_f1 p_g0 p_g1 p_g2;
  proc sql noprint;
    select mean(age_c)      format=best16., mean(educ_c)     format=best16.,
           mean(dx=1)       format=best16., mean(dx=2)       format=best16.,
           mean(dx=3)       format=best16.,
           mean(female=0)   format=best16., mean(female=1)   format=best16.,
           mean(dropgrp=0)  format=best16., mean(dropgrp=1)  format=best16.,
           mean(dropgrp=2)  format=best16.
      into :m_agec trimmed, :m_educc trimmed, :p_dx1 trimmed, :p_dx2 trimmed,
           :p_dx3 trimmed, :p_f0 trimmed, :p_f1 trimmed,
           :p_g0 trimmed, :p_g1 trimmed, :p_g2 trimmed
      from &data where &where;
  quit;
  %put NOTE: weights age_c=&m_agec dx=&p_dx1/&p_dx2/&p_dx3 female=&p_f0/&p_f1;
%mend covmeans;

/* macro: ESTIMATE statements for the annual slope by APOE4 dose.
   extra= lets the pattern-mixture model add its time*dropgrp weights. */
%macro slope_estimates(extra=);
  %local common;
  %let common = time*age_c &m_agec time*dx &p_dx1 &p_dx2 &p_dx3
                time*female &p_f0 &p_f1 &extra;
  estimate 'Slope: 0 e4' time 1 time*apoe4 1 0 0 &common / cl;
  estimate 'Slope: 1 e4' time 1 time*apoe4 0 1 0 &common / cl;
  estimate 'Slope: 2 e4' time 1 time*apoe4 0 0 1 &common / cl;
  estimate 'Slope diff: 1 vs 0 e4' time*apoe4 -1 1 0 / cl;
  estimate 'Slope diff: 2 vs 0 e4' time*apoe4 -1 0 1 / cl;
  estimate 'Slope diff: 2 vs 1 e4' time*apoe4  0 -1 1 / cl;
  contrast 'APOE4 x time (2 df)'   time*apoe4 -1 1 0, time*apoe4 -1 0 1;
%mend slope_estimates;

/* macro: population-average trajectory at t = 0..6 for each APOE4 dose */
%macro traj_estimates;
  %local a t c0 c1 c2;
  %do a = 0 %to 2;
    %let c0 = 0; %let c1 = 0; %let c2 = 0; %let c&a = 1;
    %do t = 0 %to 6;
      estimate "Traj e4=&a t=&t"
        intercept 1 time &t
        apoe4 &c0 &c1 &c2
        time*apoe4 %sysevalf(&t*&c0) %sysevalf(&t*&c1) %sysevalf(&t*&c2)
        age_c &m_agec time*age_c %sysevalf(&t*&m_agec)
        dx &p_dx1 &p_dx2 &p_dx3
        time*dx %sysevalf(&t*&p_dx1) %sysevalf(&t*&p_dx2) %sysevalf(&t*&p_dx3)
        female &p_f0 &p_f1
        time*female %sysevalf(&t*&p_f0) %sysevalf(&t*&p_f1)
        educ_c &m_educc / cl;
    %end;
  %end;
%mend traj_estimates;

/* 5a. Random-effects structure (ML, so LR tests are valid) */
%macro fitre(label, random, out);
  title "Unconditional growth model: &label";
  proc mixed data=adni.long_obs method=ml covtest plots=none;
    class rid;
    model y = time / solution;
    &random;
    ods output FitStatistics=work.&out;
  run;
%mend fitre;
%fitre(random intercept,        %str(random intercept / subject=rid),                 fit1);
%fitre(random slope independent,%str(random intercept time / subject=rid type=vc),    fit2);
%fitre(random slope unstructured,%str(random intercept time / subject=rid type=un),   fit3);

data work.lrt_re;
  merge work.fit1(where=(descr='-2 Log Likelihood') rename=(value=m2ll_ri))
        work.fit2(where=(descr='-2 Log Likelihood') rename=(value=m2ll_rs_vc))
        work.fit3(where=(descr='-2 Log Likelihood') rename=(value=m2ll_rs_un));
  /* RI vs RS(vc): slope variance on boundary -> 50:50 mixture chi2(0)/chi2(1) */
  chi_slope = m2ll_ri - m2ll_rs_vc;
  p_slope   = 0.5*(1 - probchi(chi_slope, 1));
  chi_cov   = m2ll_rs_vc - m2ll_rs_un;
  p_cov     = 1 - probchi(chi_cov, 1);
  drop descr;
run;
title 'LR tests for random-effects structure';
proc print data=work.lrt_re noobs; run;

/* 5b. Linear vs quadratic time (ML) */
%macro fitquad(extra, out);
  proc mixed data=adni.long_obs method=ml order=internal plots=none;
    class rid apoe4 dx female;
    model y = time apoe4 time*apoe4 age_c time*age_c dx time*dx female
              time*female educ_c &extra / solution;   /* sign of time*time */
    random intercept time / subject=rid type=un;
    ods output FitStatistics=work.&out;
  run;
%mend fitquad;
title 'Linear time (ML)';     %fitquad(, fq1);
title 'Quadratic time (ML)';  %fitquad(time*time, fq2);
data work.lrt_quad;
  merge work.fq1(where=(descr='-2 Log Likelihood') rename=(value=m2ll_lin))
        work.fq2(where=(descr='-2 Log Likelihood') rename=(value=m2ll_quad));
  chisq = m2ll_lin - m2ll_quad;  p = 1 - probchi(chisq, 1);
  drop descr;
run;
title 'LR test: quadratic time term';
proc print data=work.lrt_quad noobs; run;

/* 5c. FINAL MODEL (REML)
   y_ij = b0 + b1*time + APOE4 + age + dx + sex + educ
          + time x (APOE4 + age + dx + sex) + u0_i + u1_i*time + e_ij */
data work.score;
  set adni.long_obs;
  by rid;
  if first.rid;
  y = .;
  score = 1;
  do time = 0, 1; output; end;
run;
data work.final_in;
  set adni.long_obs(in=a) work.score;
  by rid;
  if a then score = 0;
run;

%covmeans(data=adni.long_obs);

title "FINAL MODEL: &outcome ~ time x (APOE4 + age + dx + sex) + educ, REML, UN G";
proc mixed data=work.final_in method=reml covtest order=internal plots=none;
  class rid apoe4 dx female;
  model y = time apoe4 time*apoe4 age_c time*age_c dx time*dx female
            time*female educ_c
            / solution ddfm=bw residual outpm=work.pm outp=work.pc;
  random intercept time / subject=rid type=un g gcorr solution;
  %slope_estimates();
  %traj_estimates;
  ods output SolutionF=adni.fixed_main CovParms=adni.covparms_main
             Tests3=adni.tests3_main Estimates=work.est_main
             Contrasts=adni.contrast_main SolutionR=work.blup;
  ods exclude SolutionR Estimates;
run;

data adni.slopes_main adni.traj_main;
  set work.est_main;
  if index(label, 'Traj') then do;
    apoe4 = input(scan(scan(label, 2, ' '), 2, '='), 8.);
    time  = input(scan(scan(label, 3, ' '), 2, '='), 8.);
    output adni.traj_main;
  end;
  else output adni.slopes_main;
run;

title 'Annual rate of change by APOE4 dose (averaged over sample covariates)';
proc print data=adni.slopes_main noobs label;
  var label estimate stderr df tvalue probt lower upper;
run;
proc print data=adni.contrast_main noobs; run;

title 'Predicted trajectory by APOE4 dose';
proc sgplot data=adni.traj_main;
  band   x=time lower=lower upper=upper / group=apoe4 transparency=0.8;
  series x=time y=estimate / group=apoe4 lineattrs=(thickness=2) markers;
  xaxis label='Years since baseline' values=(0 to 6);
  yaxis label="Predicted &outcome";
  keylegend / title='APOE e4 alleles';
run;

/* 5d. Exploratory: APOE4 effect on slope by baseline diagnosis */
title 'Exploratory: APOE4 carrier x diagnosis x time';
proc mixed data=adni.long_obs method=reml order=internal plots=none;
  class rid apoe4_carrier dx female;
  model y = time|apoe4_carrier|dx age_c time*age_c female time*female educ_c
            / solution ddfm=bw;
  random intercept time / subject=rid type=un;
  /* carrier minus non-carrier difference in annual slope within each
     diagnosis; time*apoe4_carrier*dx levels: (0,CN)(0,MCI)(0,DEM)(1,CN)(1,MCI)(1,DEM) */
  estimate 'Carrier slope diff: CN'       time*apoe4_carrier -1 1
           time*apoe4_carrier*dx -1 0 0 1 0 0 / cl;
  estimate 'Carrier slope diff: MCI'      time*apoe4_carrier -1 1
           time*apoe4_carrier*dx 0 -1 0 0 1 0 / cl;
  estimate 'Carrier slope diff: Dementia' time*apoe4_carrier -1 1
           time*apoe4_carrier*dx 0 0 -1 0 0 1 / cl;
  ods select Tests3 Estimates;
run;

/* PREDICTIONS, BLUPs, INDIVIDUAL SLOPES, RESIDUAL DIAGNOSTICS */
/* individual annual change from the scoring rows */
proc transpose data=work.pc(where=(score=1)) out=work.cslope(drop=_name_) prefix=pc_t;
  by rid; id time; var pred;
run;
proc transpose data=work.pm(where=(score=1)) out=work.mslope(drop=_name_) prefix=pm_t;
  by rid; id time; var pred;
run;

/* BLUPs: one row per subject */
data work.blup;
  set work.blup;
  /* the subject column may be stored as character or numeric: vvalue()
     handles both */
  rid_n = input(strip(vvalue(rid)), 12.);
  drop rid;
  rename rid_n = rid;
run;
proc sort data=work.blup; by rid; run;
proc transpose data=work.blup out=work.blupw(drop=_name_) prefix=re_;
  by rid; id effect; var estimate;
run;

/* per-person OLS slopes (>= 3 visits) to illustrate shrinkage */
proc sql;
  create table work.ge3 as
  select * from adni.long_obs group by rid having count(*) >= 3
  order by rid, year;
quit;
proc reg data=work.ge3 outest=work.ols(keep=rid time rename=(time=ols_slope))
         noprint;
  by rid;
  model y = time;
run;

data adni.individual_slopes;
  merge work.cslope(in=a) work.mslope work.blupw work.ols;
  by rid;
  if a;
  indiv_slope = pc_t1 - pc_t0;     /* fixed slope + random-slope BLUP */
  fixed_slope = pm_t1 - pm_t0;     /* covariate-specific population slope */
  label indiv_slope='Predicted individual annual change'
        re_Intercept='Random intercept BLUP' re_time='Random slope BLUP';
  drop pc_t0 pc_t1 pm_t0 pm_t1;
run;
proc sort data=work.wide; by rid; run;
data adni.individual_slopes;
  merge adni.individual_slopes(in=a)
        work.wide(keep=rid apoe4 dx age_bl female educ);
  by rid;
  if a;
run;

title 'Individual annual rate of change (BLUP-based) by APOE4 and diagnosis';
proc means data=adni.individual_slopes n mean std p25 median p75 maxdec=2;
  class apoe4;
  var indiv_slope;
run;
proc means data=adni.individual_slopes n mean std p25 median p75 maxdec=2;
  class dx;
  var indiv_slope;
run;
title 'Shrinkage: BLUP slopes vs per-person OLS slopes';
proc corr data=adni.individual_slopes;
  var indiv_slope ols_slope;
run;

title 'Random-effects BLUPs';
proc sgplot data=adni.individual_slopes;
  scatter x=re_Intercept y=re_time / markerattrs=(size=4) transparency=0.6;
  xaxis label='Random intercept (BLUP)'; yaxis label='Random slope (BLUP)';
run;

/* residual diagnostics (conditional residuals, observed rows only) */
data work.pc_obs;
  set work.pc(where=(score=0));
run;
title 'Studentized conditional residuals';
proc univariate data=work.pc_obs noprint;
  var studentresid;
  histogram studentresid / normal;
  qqplot    studentresid / normal(mu=est sigma=est);
run;
proc sgplot data=work.pc_obs;
  scatter x=pred y=studentresid / markerattrs=(size=3) transparency=0.7;
  refline 0 / axis=y;
  xaxis label='Conditional prediction'; yaxis label='Studentized residual';
run;

/* observed vs fitted for 20 random participants with >= 4 visits */
proc sql;
  create table work.ids4 as
  select distinct rid from adni.long_obs group by rid having count(*) >= 4;
quit;
proc surveyselect data=work.ids4 out=work.ids20 method=srs sampsize=20
                  seed=20261002 noprint;
run;
data work.spag;
  merge work.pc_obs(in=a keep=rid time y pred rename=(pred=fit_cond))
        work.pm(where=(score=0) keep=rid time score pred rename=(pred=fit_marg))
        work.ids20(in=b);
  by rid;
  if a and b;
run;
title 'Observed (dots), subject-specific fit (solid), population fit (dashed)';
proc sgpanel data=work.spag;
  panelby rid / columns=5 rows=4 novarname;
  scatter x=time y=y        / markerattrs=(symbol=circlefilled size=5);
  series  x=time y=fit_cond / lineattrs=(thickness=2);
  series  x=time y=fit_marg / lineattrs=(pattern=dash);
  colaxis label='Years'; rowaxis label="&outcome";
run;

/* PROC GLIMMIX: BINOMIAL GLMM FOR MMSE ERRORS */
data work.mmse_bin;
  set adni.long_all;
  where not missing(mmse) and not missing(time);
  errors = 30 - mmse;
  items  = 30;
run;

title 'Binomial GLMM: MMSE errors out of 30 ~ time x APOE4 (Laplace)';
proc glimmix data=work.mmse_bin method=laplace order=internal;
  class rid apoe4 dx female;
  model errors/items = time apoe4 time*apoe4 age_c time*age_c dx time*dx
                       female educ_c / dist=binomial link=logit solution;
  random intercept time / subject=rid type=un;
  estimate 'OR/yr: 0 e4' time 1 time*apoe4 1 0 0 time*age_c &m_agec
           time*dx &p_dx1 &p_dx2 &p_dx3 / exp cl;
  estimate 'OR/yr: 1 e4' time 1 time*apoe4 0 1 0 time*age_c &m_agec
           time*dx &p_dx1 &p_dx2 &p_dx3 / exp cl;
  estimate 'OR/yr: 2 e4' time 1 time*apoe4 0 0 1 time*age_c &m_agec
           time*dx &p_dx1 &p_dx2 &p_dx3 / exp cl;
  estimate 'Ratio of ORs/yr: 1 vs 0 e4' time*apoe4 -1 1 0 / exp cl;
  estimate 'Ratio of ORs/yr: 2 vs 0 e4' time*apoe4 -1 0 1 / exp cl;
  covtest 'random slope variance = 0' . 0 0;   /* UN(1,1) . ; UN(2,1) 0; UN(2,2) 0 */
  ods output ParameterEstimates=adni.glmm_mmse Estimates=adni.glmm_mmse_or
             FitStatistics=work.glmm_fit;
run;
/* "Pearson Chi-Square / DF" in the conditional fit table >> 1 would indicate
   residual over-dispersion (then consider an observation-level random effect). */

/* SENSITIVITY ANALYSES FOR MISSING DATA (MNAR) */
/* 8a. Complete-case: participants retained to their last eligible visit */
%covmeans(data=adni.long_obs, where=%str(dropgrp=0));
title 'Sensitivity: complete-case (retained participants only)';
proc mixed data=adni.long_obs(where=(dropgrp=0)) method=reml order=internal plots=none;
  class rid apoe4 dx female;
  model y = time apoe4 time*apoe4 age_c time*age_c dx time*dx female
            time*female educ_c / solution ddfm=bw;
  random intercept time / subject=rid type=un;
  %slope_estimates();
  ods output Estimates=work.est_cc;
run;

/* 8b. Pattern-mixture: intercept and slope vary by dropout pattern; the
       APOE4 slopes average over the observed mix of patterns */
%covmeans(data=adni.long_obs);
title 'Sensitivity: pattern-mixture model';
proc mixed data=adni.long_obs method=reml order=internal plots=none;
  class rid apoe4 dx female dropgrp;
  model y = time apoe4 time*apoe4 age_c time*age_c dx time*dx female
            time*female educ_c dropgrp time*dropgrp / solution ddfm=bw;
  random intercept time / subject=rid type=un;
  %slope_estimates(extra=time*dropgrp &p_g0 &p_g1 &p_g2);
  estimate 'Slope diff: dropout by yr2 vs retained' time*dropgrp -1 1 0 / cl;
  estimate 'Slope diff: dropout yr3-5 vs retained'  time*dropgrp -1 0 1 / cl;
  ods output Estimates=work.est_pmm Tests3=work.tests3_pmm;
run;

data adni.slopes_compare;
  length analysis $16 label $60;
  set adni.slopes_main(in=a) work.est_cc(in=b) work.est_pmm(in=c);
  if a then analysis = 'MAR (main)';
  else if b then analysis = 'Complete-case';
  else analysis = 'Pattern-mixture';
  keep analysis label estimate stderr probt lower upper;
run;
title 'APOE4 slopes: main model vs MNAR sensitivity analyses';
proc print data=adni.slopes_compare noobs; run;

/* EXPORT */
%macro exp(ds, file);
  proc export data=&ds outfile="&root/output/&file..csv" dbms=csv replace; run;
%mend exp;
%exp(adni.fixed_main,        sas_fixed_effects_main);
%exp(adni.covparms_main,     sas_covparms_main);
%exp(adni.slopes_compare,    sas_slopes_by_apoe4);
%exp(adni.traj_main,         sas_trajectories_apoe4);
%exp(adni.individual_slopes, sas_individual_slopes);
%exp(adni.glmm_mmse_or,      sas_glimmix_mmse_or);
%exp(adni.glmm_missingness,  sas_glimmix_missingness);

title;
ods pdf close;
ods graphics off;
proc printto; run;