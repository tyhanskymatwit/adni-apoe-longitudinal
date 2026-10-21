/*==============================================================================
  Program : adni_tlf.sas
  Purpose : Tables, listings and figures (TLFs) for the ADNI APOE4 longitudinal
            analysis of ADAS-Cog 13, in submission style (PROC REPORT, ODS RTF).
  Inputs  : &root/data/adni_wide.csv           (one row per participant)
            adni.slopes_compare, adni.traj_main (from adni_longitudinal.sas)
  Outputs : &root/output/TLF/adni_tlfs.rtf      all TLFs, one per page
            &root/output/TLF/adni_tlfs.pdf      same content as PDF
            &root/output/TLF/qc_t*.csv          table cells for independent
                                                QC (tlf/adni_tlf_qc.R)
  Run     : after adni_longitudinal.sas (needs the adni library datasets)

  Outputs produced
    Table   14.1.1  Demographics and baseline characteristics
    Table   14.1.2  Participant disposition and visit retention
    Table   14.2.1  ADAS-Cog 13 observed values and change from baseline by visit
    Table   14.2.2  Annual rate of change by APOE e4 dose: main and sensitivity models
    Figure  14.2.1  Mean (+/- SE) change from baseline in ADAS-Cog 13 by visit
    Figure  14.2.2  Model-predicted ADAS-Cog 13 trajectory with 95% CI
    Figure  14.2.3  Forest plot of slope differences across missing-data analyses
    Listing 16.2.1  Discontinued participants with two APOE e4 alleles
==============================================================================*/

%let root   = /home/<your-sas-user>/ADNI;   /* EDIT: project folder */
%let tlfdir = &root/output/TLF;
%let window = 0.35;     /* visit k is "eligible" if follow-up potential >= k + window */
%let pgm    = adni_tlf.sas;

options dlcreatedir;
libname _mk "&tlfdir"; libname _mk clear;      /* creates the TLF folder */
libname adni "&root/output";

options nodate nonumber orientation=landscape missing=' ' nofmterr validvarname=v7
        leftmargin=0.6in rightmargin=0.6in topmargin=0.6in bottommargin=0.6in;
/* landscape Letter with 0.6in margins leaves 9.8in: every report below is sized in
   inches to fit inside it, so columns never wrap onto a second page */
ods escapechar='^';
%let e4 = APOE ^{unicode 03B5}4;

proc format;
  value trtf   0='0 e4' 1='1 e4' 2='2 e4' 3='Total';
  value dxf    1='Cognitively normal (CN)' 2='Mild cognitive impairment (MCI)' 3='Dementia';
  value sexf   0='Male' 1='Female';
  value amyf   0='Not elevated' 1='Elevated';
  value protf  1='ADNI1' 2='ADNIGO' 3='ADNI2' 4='ADNI3' 5='ADNI4';
  value dgf    0='Retained to last eligible visit' 1='Discontinued by year 2' 2='Discontinued in years 3-5';
  value ynf    0='No' 1='Yes';
run;

/*------------------------------------------------------------------------------
  1. ANALYSIS DATASETS
     adsl : one row per participant in the analysis set (ADaM ADSL-like)
     adqs : one row per participant x nominal year with ADAS-Cog 13 (BDS-like)
------------------------------------------------------------------------------*/
proc import datafile="&root/data/adni_wide.csv" out=work.raw dbms=csv replace;
  guessingrows=max;
  getnames=yes;
run;

data work.adsl;
  set work.raw;
  if      dx_bl = 'CN'  then dxn = 1;
  else if dx_bl = 'MCI' then dxn = 2;
  else if dx_bl = 'DEM' then dxn = 3;
  /* analysis set: same inclusion rule as adni_longitudinal.sas */
  if nmiss(age_bl, female, educ, dxn, apoe4, adas13_0) = 0;
  trtn = apoe4;
  select (protocol);
    when ('ADNI1')  protn = 1;
    when ('ADNIGO') protn = 2;
    when ('ADNI2')  protn = 3;
    when ('ADNI3')  protn = 4;
    when ('ADNI4')  protn = 5;
    otherwise       protn = .;
  end;

  /* dropout definition (separates administrative censoring from attrition) */
  array yy{0:6} adas13_0 - adas13_6;
  array el{0:6} elig_0 - elig_6;
  lastobs = 0; lastelig = 0; lastval = adas13_0;
  do k = 0 to 6;
    el{k} = (fu_potential >= k + &window) or not missing(yy{k});
    if not missing(yy{k}) then do; lastobs = k; lastval = yy{k}; end;
    if el{k} then lastelig = k;
  end;
  dropout = (lastobs < lastelig);
  if      dropout = 0  then dropgrp = 0;
  else if lastobs <= 2 then dropgrp = 1;
  else                      dropgrp = 2;
  drop k;
run;

/* add a "Total" column (trtn = 3) by duplicating every participant */
data work.adsl_t;
  set work.adsl;
  output;
  trtn = 3; output;
run;

/* big N per column -> macro variables n0..n3 for the column headers */
proc sql noprint;
  create table work.bign as
    select trtn, count(*) as bign from work.adsl_t group by trtn order by trtn;
  select bign into :n0 - :n3 from work.bign;
quit;
%put NOTE: big N = &n0 / &n1 / &n2 / &n3;

/* BDS-like: one record per participant x visit year with an observed score */
data work.adqs;
  set work.adsl;
  array yy{0:6} adas13_0 - adas13_6;
  base = adas13_0;
  do avisitn = 0 to 6;
    aval = yy{avisitn};
    if not missing(aval) then do;
      chg = ifn(avisitn > 0, aval - base, .);
      output;
    end;
  end;
  keep rid trtn avisitn aval base chg;
run;

/*------------------------------------------------------------------------------
  2. SUMMARY MACROS
     Every table row is (ord, sub, rowlbl, stat, trt0..trt3), so all tables share
     one PROC REPORT layout and one QC export format.
------------------------------------------------------------------------------*/

/* continuous variable: n, Mean (SD), Median, Min, Max. dec = raw decimals */
%macro cont(var=, label=, ord=, dec=1, data=work.adsl_t, out=work.tab);
  %local d1 d2;
  %let d1 = %eval(&dec + 1);
  %let d2 = %eval(&dec + 2);
  proc means data=&data noprint nway;
    class trtn;
    var &var;
    output out=work._s n=n mean=mean std=sd median=med min=min max=max;
  run;
  data work._r;
    set work._s;
    length rowlbl $80 stat $60 val $40;
    rowlbl = "&label"; ord = &ord;
    sub = 1; stat = 'n';         val = strip(put(n, 8.));                                   output;
    sub = 2; stat = 'Mean (SD)'; val = catx(' ', put(mean, 12.&d1), cats('(', put(sd, 12.&d2), ')')); output;
    sub = 3; stat = 'Median';    val = strip(put(med, 12.&d1));                             output;
    sub = 4; stat = 'Min, Max';  val = catx(', ', put(min, 12.&dec), put(max, 12.&dec));    output;
    keep ord sub rowlbl stat trtn val;
  run;
  %_stack(out=&out);
%mend cont;

/* categorical variable: n (%) with the column big N as denominator */
%macro cat(var=, label=, ord=, fmt=, data=work.adsl_t, out=work.tab);
  proc freq data=&data noprint;
    tables trtn * &var / out=work._f(drop=percent) missing sparse;
  run;
  proc sql;
    create table work._f2 as
      select a.*, b.bign from work._f a inner join work.bign b on a.trtn = b.trtn;
  quit;
  data work._r;
    set work._f2;
    length rowlbl $80 stat $60 val $40;
    rowlbl = "&label"; ord = &ord;
    if missing(&var) then do; sub = 99; stat = 'Missing'; end;
    else do; sub = &var; stat = put(&var, &fmt); end;
    val = catx(' ', put(count, 8.), cats('(', put(100 * count / bign, 8.1), ')'));
    keep ord sub rowlbl stat trtn val;
  run;
  %_stack(out=&out);
%mend cat;

/* transpose long cells (trtn) to columns trt0..trt3 and append to &out */
%macro _stack(out=);
  proc sort data=work._r; by ord sub rowlbl stat; run;
  proc transpose data=work._r out=work._t(drop=_name_) prefix=trt;
    by ord sub rowlbl stat;
    id trtn;
    var val;
  run;
  proc append base=&out data=work._t force; run;
%mend _stack;

/* standard titles and footnotes */
%macro tf(num=, ttl=, pop=Analysis set, fn1=, fn2=);
  title1 j=l 'ADNI: APOE e4 and longitudinal cognitive decline' j=r 'Page ^{thispage} of ^{lastpage}';
  title2 j=c "&num";
  title3 j=c "&ttl";
  title4 j=c "&pop";
  footnote1 j=l "&fn1";
  footnote2 j=l "&fn2";
  footnote3 j=l "Program: &pgm   Source: adni_wide.csv (data cut 2026-07-06)   Run: %sysfunc(today(), date9.) %sysfunc(time(), time5.)";
%mend tf;

/* export a finished table for QC */
%macro qcexp(ds, file);
  proc export data=&ds outfile="&tlfdir/&file..csv" dbms=csv replace; run;
%mend qcexp;

/* shared PROC REPORT layout for trt0..trt3 tables */
%macro rpt(ds);
  proc report data=&ds nowd split='|' missing
       style(report)=[cellpadding=3pt] style(header)=[just=c font_weight=bold]
       style(column)=[just=c];
    column ord rowlbl sub stat trt0 trt1 trt2 trt3;
    define ord    / order order=internal noprint;
    define rowlbl / order order=data ' ' style(column)=[just=l cellwidth=2.0in font_weight=bold];
    define sub    / order order=internal noprint;
    define stat   / display ' ' style(column)=[just=l cellwidth=2.2in];
    define trt0   / display "&e4 = 0|(N=&n0)" style(column)=[cellwidth=1.1in];
    define trt1   / display "&e4 = 1|(N=&n1)" style(column)=[cellwidth=1.1in];
    define trt2   / display "&e4 = 2|(N=&n2)" style(column)=[cellwidth=1.1in];
    define trt3   / display "Total|(N=&n3)"   style(column)=[cellwidth=1.1in];
    compute after ord; line ' '; endcomp;
  run;
%mend rpt;

/*------------------------------------------------------------------------------
  3. BUILD TABLE DATA
------------------------------------------------------------------------------*/

/* Table 14.1.1 Demographics and baseline characteristics */
proc datasets lib=work nolist; delete t1411 t1412 t1421; quit;
%cont(var=age_bl,   label=Age at baseline (years),   ord=1, dec=1, out=work.t1411);
%cat (var=female,   label=Sex n (%),                 ord=2, fmt=sexf.,  out=work.t1411);
%cont(var=educ,     label=Education (years),         ord=3, dec=0, out=work.t1411);
%cat (var=dxn,      label=Baseline diagnosis n (%),  ord=4, fmt=dxf.,   out=work.t1411);
%cat (var=amyloid,  label=Amyloid PET status n (%),  ord=5, fmt=amyf.,  out=work.t1411);
%cat (var=protn,    label=Enrolment protocol n (%),  ord=6, fmt=protf., out=work.t1411);
%cont(var=adas13_0, label=ADAS-Cog 13 at baseline,   ord=7, dec=2, out=work.t1411);
%cont(var=mmse_0,   label=MMSE at baseline,          ord=8, dec=0, out=work.t1411);
%cont(var=cdrsb_0,  label=CDR-SB at baseline,        ord=9, dec=1, out=work.t1411);

/* Table 14.1.2 Disposition (part 1) and visit retention (part 2) */
%cat(var=dropgrp, label=Disposition n (%),     ord=1, fmt=dgf., out=work.t1412);
%cat(var=died,    label=Death recorded n (%),  ord=2, fmt=ynf., out=work.t1412);

data work.ret;
  set work.adsl_t;
  array yy{0:6} adas13_0 - adas13_6;
  array el{0:6} elig_0 - elig_6;
  do year = 1 to 6;
    if el{year} then do; observed = not missing(yy{year}); output; end;
  end;
  keep trtn year observed;
run;
proc means data=work.ret noprint nway;
  class trtn year;
  var observed;
  output out=work._rs n=elig sum=obs;
run;
data work._r;
  set work._rs;
  length rowlbl $80 stat $60 val $40;
  rowlbl = 'Visit retention: observed / eligible (%)';
  ord = 3; sub = year; stat = catx(' ', 'Year', year);
  val = catx(' ', cats(obs, '/', elig), cats('(', put(100 * obs / elig, 8.1), ')'));
  keep ord sub rowlbl stat trtn val;
run;
%_stack(out=work.t1412);

/* Table 14.2.1 ADAS-Cog 13 by visit: observed value and change from baseline
   columns c<trtn><param>: param 1 = observed, 2 = change from baseline */
data work.adqs_l;
  set work.adqs;
  param = 1; v = aval; output;
  if avisitn > 0 then do; param = 2; v = chg; output; end;
  keep trtn avisitn param v;
run;
proc means data=work.adqs_l noprint nway;
  class avisitn trtn param;
  var v;
  output out=work._s n=n mean=mean std=sd median=med min=min max=max;
run;
data work._r;
  set work._s;
  length rowlbl $80 stat $60 val $40;
  ord = avisitn;
  rowlbl = ifc(avisitn = 0, 'Baseline', catx(' ', 'Year', avisitn));
  sub = 1; stat = 'n';         val = strip(put(n, 8.));                                  output;
  sub = 2; stat = 'Mean (SD)'; val = catx(' ', put(mean, 12.3), cats('(', put(sd, 12.4), ')')); output;
  sub = 3; stat = 'Median';    val = strip(put(med, 12.3));                              output;
  sub = 4; stat = 'Min, Max';  val = catx(', ', put(min, 12.2), put(max, 12.2));         output;
  keep ord sub rowlbl stat trtn param val;
run;
proc sort data=work._r; by ord sub rowlbl stat; run;
proc transpose data=work._r out=work.t1421(drop=_name_) prefix=c;
  by ord sub rowlbl stat;
  id trtn param;
  var val;
run;

/* Table 14.2.2 Annual slopes by APOE e4 dose (from adni_longitudinal.sas) */
data work.t1422;
  set adni.slopes_compare;
  where index(label, 'e4') > 0;
  length analysis_l $40 est_se ci $30 pval $10;
  select (analysis);
    when ('MAR (main)')      do; aord = 1; analysis_l = 'Main model (MAR)';           end;
    when ('Complete-case')   do; aord = 2; analysis_l = 'Complete-case (retained)';   end;
    when ('Pattern-mixture') do; aord = 3; analysis_l = 'Pattern-mixture (MNAR)';     end;
    otherwise aord = 9;
  end;
  lord   = _n_;
  est_se = catx(' ', put(estimate, 8.2), cats('(', put(stderr, 8.3), ')'));
  ci     = '(' || strip(put(lower, 8.2)) || ', ' || strip(put(upper, 8.2)) || ')';
  pval   = ifc(probt < 0.001, '<0.001', strip(put(probt, 6.3)));
run;

/*------------------------------------------------------------------------------
  4. RENDER: RTF + PDF
------------------------------------------------------------------------------*/
ods listing close;
ods graphics on / reset width=9in height=5in imagefmt=png noborder;
ods rtf file="&tlfdir/adni_tlfs.rtf" style=journal;
ods pdf file="&tlfdir/adni_tlfs.pdf" style=journal notoc;

/* ---- Table 14.1.1 ---- */
%tf(num=Table 14.1.1, ttl=Demographics and Baseline Characteristics by &e4 Allele Count,
    fn1=Percentages use the column N as denominator. SD = standard deviation.,
    fn2=Analysis set: non-missing age%str(,) sex%str(,) education%str(,) baseline diagnosis%str(,) APOE genotype and baseline ADAS-Cog 13.);
%rpt(work.t1411);

/* ---- Table 14.1.2 ---- */
%tf(num=Table 14.1.2, ttl=Participant Disposition and Visit Retention by &e4 Allele Count,
    fn1=Discontinued = last observed ADAS-Cog 13 visit precedes the last visit the participant was administratively eligible for.,
    fn2=Eligible for year k = follow-up potential at data cut >= k + &window years%str(,) or the visit was observed.);
%rpt(work.t1412);

/* ---- Table 14.2.1 ---- */
%tf(num=Table 14.2.1, ttl=ADAS-Cog 13: Observed Values and Change from Baseline by Nominal Visit Year,
    fn1=ADAS-Cog 13 range 0-85 (higher scores indicate worse cognition). Change from baseline = visit value - baseline value.,
    fn2=Visits mapped to nominal annual years from actual visit dates (window +/- &window years). Observed data only%str(,) no imputation.);
proc report data=work.t1421 nowd split='|' missing
     style(report)=[cellpadding=2pt fontsize=9pt] style(header)=[just=c font_weight=bold fontsize=9pt]
     style(column)=[just=c fontsize=9pt];
  column ord rowlbl sub stat
         ("&e4 = 0 (N=&n0)" c01 c02)
         ("&e4 = 1 (N=&n1)" c11 c12)
         ("&e4 = 2 (N=&n2)" c21 c22);
  define ord    / order order=internal noprint;
  define rowlbl / order order=data 'Visit' style(column)=[just=l cellwidth=0.75in font_weight=bold];
  define sub    / order order=internal noprint;
  define stat   / display 'Statistic' style(column)=[just=l cellwidth=0.85in];
  define c01    / display 'Observed'          style(column)=[cellwidth=1.15in];
  define c02    / display 'Change from|baseline' style(column)=[cellwidth=1.15in];
  define c11    / display 'Observed'          style(column)=[cellwidth=1.15in];
  define c12    / display 'Change from|baseline' style(column)=[cellwidth=1.15in];
  define c21    / display 'Observed'          style(column)=[cellwidth=1.15in];
  define c22    / display 'Change from|baseline' style(column)=[cellwidth=1.15in];
  compute after ord; line ' '; endcomp;
run;

/* ---- Table 14.2.2 ---- */
%tf(num=Table 14.2.2, ttl=Annual Rate of Change in ADAS-Cog 13 by &e4 Allele Count: Main and Missing-Data Sensitivity Analyses,
    fn1=Linear mixed-effects model (REML%str(,) random intercept and slope%str(,) unstructured G)%str(,) adjusted for age%str(,) sex%str(,) education and baseline diagnosis%str(,) each interacted with time.,
    fn2=Slopes are points per year averaged over the observed covariate distribution. Positive values indicate decline.);
proc report data=work.t1422 nowd split='|' missing
     style(report)=[cellpadding=3pt] style(header)=[just=c font_weight=bold] style(column)=[just=c];
  column aord analysis_l lord label est_se ci pval;
  define aord       / order order=internal noprint;
  define analysis_l / order 'Analysis' style(column)=[just=l cellwidth=2.0in font_weight=bold];
  define lord       / order order=internal noprint;
  define label      / display 'Parameter' style(column)=[just=l cellwidth=2.2in];
  define est_se     / display 'Estimate (SE)|points/year' style(column)=[cellwidth=1.6in];
  define ci         / display '95% CI' style(column)=[cellwidth=1.6in];
  define pval       / display 'p-value' style(column)=[cellwidth=0.9in];
  compute after aord; line ' '; endcomp;
run;

/* ---- Figure 14.2.1: observed mean change from baseline ---- */
proc means data=work.adqs noprint nway;
  class trtn avisitn;
  var chg;
  output out=work.f1 mean=mean stderr=se n=n;
run;
data work.f1;
  set work.f1;
  if avisitn = 0 then do; mean = 0; se = 0; end;
  lo = mean - se; hi = mean + se;
  /* small horizontal offset so error bars do not overlap */
  xpos = avisitn + (trtn - 1) * 0.08;
  format trtn trtf.;
run;
%tf(num=Figure 14.2.1, ttl=Mean (+/- SE) Change from Baseline in ADAS-Cog 13 by Visit Year and &e4 Allele Count,
    fn1=Observed data only. Error bars show +/- 1 standard error of the mean. Positive change indicates worsening.,
    fn2=);
proc sgplot data=work.f1 noautolegend;
  styleattrs datacontrastcolors=(cx3987E5 cxD95926 cx199E70) datasymbols=(circlefilled squarefilled trianglefilled);
  highlow x=xpos low=lo high=hi / group=trtn lowcap=serif highcap=serif;
  series  x=xpos y=mean / group=trtn markers lineattrs=(thickness=2) name='s';
  refline 0 / axis=y lineattrs=(pattern=dot);
  xaxis label='Years since baseline' values=(0 to 6 by 1);
  yaxis label='Mean change from baseline in ADAS-Cog 13';
  keylegend 's' / title='APOE e4 alleles' location=inside position=topleft;
run;

/* ---- Figure 14.2.2: model-predicted trajectories ---- */
data work.f2;
  set adni.traj_main;
  format apoe4 trtf.;
run;
%tf(num=Figure 14.2.2, ttl=Model-Predicted ADAS-Cog 13 Trajectory with 95% Confidence Band by &e4 Allele Count,
    fn1=Population-average predictions from the main mixed-effects model (Table 14.2.2)%str(,) averaged over the observed covariate distribution.,
    fn2=);
proc sgplot data=work.f2;
  styleattrs datacolors=(cx3987E5 cxD95926 cx199E70) datacontrastcolors=(cx3987E5 cxD95926 cx199E70);
  band   x=time lower=lower upper=upper / group=apoe4 transparency=0.8;
  series x=time y=estimate / group=apoe4 markers lineattrs=(thickness=2) name='s';
  xaxis label='Years since baseline' values=(0 to 6 by 1);
  yaxis label='Predicted ADAS-Cog 13';
  keylegend 's' / title='APOE e4 alleles' location=inside position=topleft;
run;

/* ---- Figure 14.2.3: forest plot of slope differences ---- */
data work.f3;
  set work.t1422;
  where index(label, 'diff') > 0;
  length contrast $30;
  contrast = tranwrd(label, 'Slope diff: ', '');
run;
%tf(num=Figure 14.2.3, ttl=Difference in Annual ADAS-Cog 13 Slope Between &e4 Groups Across Missing-Data Assumptions,
    fn1=Points show estimated slope differences (points/year). Horizontal bars show 95% confidence intervals.,
    fn2=);
proc sgplot data=work.f3;
  styleattrs datacontrastcolors=(cx3987E5 cxD95926 cx199E70) datasymbols=(circlefilled squarefilled trianglefilled);
  scatter y=contrast x=estimate / group=analysis_l groupdisplay=cluster clusterwidth=0.5
          xerrorlower=lower xerrorupper=upper markerattrs=(size=9) name='p';
  refline 0 / axis=x lineattrs=(pattern=dash);
  yaxis display=(nolabel) discreteorder=data reverse;
  xaxis label='Difference in annual slope (points/year), 95% CI';
  keylegend 'p' / title='Analysis' position=bottom;
run;

/* ---- Listing 16.2.1 ---- */
data work.l1621;
  set work.adsl;
  where apoe4 = 2 and dropout = 1;
  format female sexf. dxn dxf. died ynf.;
run;
proc sort data=work.l1621; by protn rid; run;
%tf(num=Listing 16.2.1, ttl=Discontinued Participants with Two &e4 Alleles,
    fn1=Last visit = last nominal year with an observed ADAS-Cog 13. Last eligible = last year the participant could have been seen before the data cut.,
    fn2=);
proc report data=work.l1621 nowd split='|' missing
     style(report)=[cellpadding=3pt] style(header)=[just=c font_weight=bold] style(column)=[just=c];
  column protocol rid age_bl female dxn apoe_geno adas13_0 lastobs lastelig lastval died;
  define protocol  / order 'Protocol' style(column)=[just=l];
  define rid       / display 'RID';
  define age_bl    / display 'Age|(years)' format=5.1;
  define female    / display 'Sex';
  define dxn       / display 'Baseline|diagnosis' style(column)=[just=l cellwidth=2.0in];
  define apoe_geno / display 'APOE|genotype';
  define adas13_0  / display 'Baseline|ADAS-Cog 13' format=6.2;
  define lastobs   / display 'Last visit|(year)';
  define lastelig  / display 'Last eligible|(year)';
  define lastval   / display 'ADAS-Cog 13 at|last visit' format=6.2;
  define died      / display 'Death|recorded';
run;

ods pdf close;
ods rtf close;
ods listing;
title; footnote;

/*------------------------------------------------------------------------------
  5. QC EXPORTS (compared cell by cell by tlf/adni_tlf_qc.R)
------------------------------------------------------------------------------*/
%qcexp(work.t1411, qc_t1411);
%qcexp(work.t1412, qc_t1412);
%qcexp(work.t1421, qc_t1421);
%qcexp(work.l1621(keep=rid protocol age_bl female dxn adas13_0 lastobs lastelig lastval died), qc_l1621);
