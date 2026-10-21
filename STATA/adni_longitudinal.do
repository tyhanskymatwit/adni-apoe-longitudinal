*==============================================================================
* 02_adni_longitudinal.do
* Longitudinal cognitive decline in ADNI: ADAS-Cog13 ~ time x APOE4,
* adjusted for baseline age, sex, education and baseline diagnosis.
*
* Input : data/adni_wide.csv  (created by 01_extract_adni.R from ADNIMERGE2)
*         one row per participant; outcome at nominal years 0..6 in columns
*         adas13_0..adas13_6 (also mmse_*, cdrsb_*), actual visit time t_0..t_6
*
* Sections
*   1  Import, recode, analysis sample            (wide)
*   2  Missing-data patterns & dropout definition (wide)
*   3  reshape long                               (wide -> long)
*   4  Missing-data mechanism: MCAR vs MAR         (long, discrete-time model)
*   5  Linear mixed models: random-effects structure, final model
*   6  predict: marginal/conditional fits, BLUPs, individual slopes, residuals
*   7  Sensitivity to MNAR: complete-case & pattern-mixture models
*   8  Export results
*
* Built-in commands only (Stata 16+). Outputs go to output/.
*==============================================================================
version 16
clear all
set more off
set linesize 120

* ---- settings (EDIT root if the project folder moves) ------------------------
global root    "."
global outcome "adas13"      // adas13 (higher=worse) | mmse (higher=better) | cdrsb
global window  0.35          // visit window used in 01_extract_adni.R (years)

cd "$root"
capture mkdir "output"
capture log close
log using "output/stata_adni_longitudinal.log", replace text

*==============================================================================
* 1. IMPORT, RECODE, ANALYSIS SAMPLE (WIDE: one row per participant)
*==============================================================================
import delimited using "data/adni_wide.csv", clear varnames(1) ///
    stringcols(2 6 7 10)
isid rid

* baseline diagnosis
gen byte dx = .
replace dx = 1 if dx_bl == "CN"
replace dx = 2 if dx_bl == "MCI"
replace dx = 3 if dx_bl == "DEM"
label define dxl 1 "CN" 2 "MCI" 3 "Dementia"
label values dx dxl

* APOE e4 dose (0/1/2) and carrier status
label define apol 0 "0 e4" 1 "1 e4" 2 "2 e4"
label values apoe4 apol
gen byte apoe4_carrier = apoe4 >= 1 if !missing(apoe4)

label define fem 0 "Male" 1 "Female"
label values female fem
encode protocol, gen(prot)

* centred covariates (age 72 and 16 y education ~ sample medians)
gen age_c  = age_bl - 72
gen educ_c = educ - 16
label variable age_c  "Baseline age (centred at 72)"
label variable educ_c "Education (centred at 16 y)"

* analysis sample: complete baseline covariates and a baseline outcome
count
gen byte insample = !missing(age_bl, female, educ, dx, apoe4, ${outcome}_0)
tab insample
keep if insample
drop insample

di as result _n "Analysis sample: " _N " participants"
tab dx apoe4, row
tabstat age_bl educ ${outcome}_0, by(apoe4) stat(n mean sd) format(%6.2f)

*==============================================================================
* 2. MISSING-DATA PATTERNS AND DROPOUT DEFINITION (WIDE)
*==============================================================================
* Administrative eligibility: a visit at year k could only have happened if
* the visit window had closed before the data cut-off. Missing values at
* years a person was not yet eligible for are administrative censoring
* (e.g. recent ADNI4 enrolees), not dropout.
local yv ""
forvalues k = 0/6 {
    gen byte elig_`k' = (fu_potential >= `k' + $window) | !missing(${outcome}_`k')
    local yv "`yv' ${outcome}_`k'"
}

misstable summarize `yv'
misstable patterns  `yv', frequency

* observed outcome by year, among those eligible
forvalues k = 0/6 {
    quietly count if elig_`k'
    local ne = r(N)
    quietly count if elig_`k' & !missing(${outcome}_`k')
    di as text "Year `k': eligible = " %5.0f `ne' "   observed = " %5.0f r(N) ///
       "   retention = " %5.1f 100*r(N)/`ne' "%"
}

* last observed year, last eligible year, dropout pattern
gen byte lastobs  = 0
gen byte lastelig = 0
gen byte nobs     = 0
forvalues k = 0/6 {
    replace lastobs  = `k' if !missing(${outcome}_`k')
    replace lastelig = `k' if elig_`k'
    replace nobs     = nobs + !missing(${outcome}_`k')
}
gen byte dropout = lastobs < lastelig
label variable dropout "Not seen at last eligible annual visit"

* pattern groups for the pattern-mixture model (Section 7)
gen byte dropgrp = 0
replace  dropgrp = 1 if dropout & lastobs <= 2
replace  dropgrp = 2 if dropout & lastobs >= 3
label define dgl 0 "Retained/admin-censored" 1 "Dropout by yr 2" 2 "Dropout yr 3-5"
label values dropgrp dgl
tab dropgrp
tab dropgrp died, row

* who drops out? (baseline comparison; differences => not MCAR)
tab dx dropout, row chi2
tab apoe4 dropout, row chi2
ttest age_bl,          by(dropout)
ttest ${outcome}_0,    by(dropout)
logit dropout age_c i.female educ_c i.apoe4 i.dx ${outcome}_0 i.prot, or

* Little's MCAR test (optional, user-written: ssc install mcartest)
capture which mcartest
if !_rc {
    mcartest ${outcome}_0 ${outcome}_1 ${outcome}_2 ${outcome}_3
}
else di as text "mcartest not installed - skipping Little's test"

save "output/adni_wide_analysis.dta", replace

*==============================================================================
* 3. RESHAPE WIDE -> LONG (one row per participant x nominal year)
*==============================================================================
reshape long t_ adas13_ mmse_ cdrsb_ elig_, i(rid) j(year)
rename (t_ adas13_ mmse_ cdrsb_ elig_) (time adas13 mmse cdrsb elig)

gen y = ${outcome}
label variable y    "Outcome: ${outcome}"
label variable time "Years since baseline (actual visit time)"
label variable year "Nominal annual visit"
sort rid year
xtset rid year

* rows: 7 per person, including missing visits (needed for Section 4)
count
count if !missing(y)
save "output/adni_long_all.dta", replace

*==============================================================================
* 4. MISSING-DATA MECHANISM (LONG): MCAR vs MAR
*    Discrete-time missingness model: among participants observed at year k-1
*    and eligible at year k, is the visit at k missing?
*    Model A: baseline covariates (incl. baseline outcome)
*    Model B: + most recent OBSERVED outcome (y at k-1)
*    - If A's covariates predict missingness  -> not MCAR
*    - If y_prev adds beyond baseline          -> missingness depends on the
*      evolving observed trajectory: MAR is the appropriate working assumption.
*      Likelihood-based mixed models (Section 5) are valid under MAR without
*      modelling the missingness, because y_prev is part of the observed data.
*    - MNAR (dependence on the unobserved value) cannot be tested: Section 7.
*==============================================================================
by rid: gen y_prev = y[_n-1] if year == year[_n-1] + 1
gen byte miss   = missing(y)
gen byte atrisk = year >= 1 & !missing(y_prev) & elig == 1

tab year miss if atrisk, row

bysort rid (year): gen y0 = y[1]           // baseline outcome (year 0)

logit miss age_c i.female educ_c i.apoe4 i.dx y0 i.prot i.year ///
      if atrisk, vce(cluster rid) or
estimates store missA

logit miss y_prev age_c i.female educ_c i.apoe4 i.dx y0 i.prot i.year ///
      if atrisk, vce(cluster rid) or
estimates store missB
test y_prev
di as text "Effect of most recent observed score on next-visit missingness: " ///
   "OR per point = " as result %6.3f exp(_b[y_prev])

* compare: change between the last two OBSERVED values predicts dropout?
by rid: gen dy_prev = y[_n-1] - y[_n-2] if year == year[_n-2] + 2
logit miss dy_prev y_prev age_c i.female educ_c i.apoe4 i.dx i.prot i.year ///
      if atrisk, vce(cluster rid) or

*==============================================================================
* 5. LINEAR MIXED-EFFECTS MODELS
*==============================================================================
use "output/adni_long_all.dta", clear
keep if !missing(y, time)
bysort rid: gen byte first = _n == 1

* 5a. Random-effects structure (ML so LR tests are valid)
mixed y c.time || rid: , mle
estimates store m_ri
estat icc

mixed y c.time || rid: time, covariance(independent) mle
estimates store m_rs_ind

mixed y c.time || rid: time, covariance(unstructured) mle
estimates store m_rs_un

* random slope variance is on the boundary => LR p-value is conservative
lrtest m_rs_ind m_ri
lrtest m_rs_un  m_rs_ind
estimates stats m_ri m_rs_ind m_rs_un

* 5b. Is linear time adequate? (fixed quadratic term, ML)
mixed y c.time##(i.apoe4 c.age_c i.dx i.female) c.educ_c ///
      || rid: time, covariance(unstructured) mle
estimates store m_lin_ml
mixed y c.time##(i.apoe4 c.age_c i.dx i.female) c.educ_c c.time#c.time ///
      || rid: time, covariance(unstructured) mle
estimates store m_quad_ml
lrtest m_quad_ml m_lin_ml

* 5c. FINAL MODEL (REML)
*     y_ij = b0 + b1*time + APOE4 + age + dx + sex + educ
*            + time x (APOE4 + age + dx + sex)
*            + u0_i + u1_i*time + e_ij,   (u0,u1) ~ N(0, G unstructured)
mixed y c.time##i.apoe4 c.time##c.age_c c.time##i.dx c.time##i.female c.educ_c ///
      || rid: time, covariance(unstructured) reml
estimates store m_main

testparm i.apoe4#c.time                   // does APOE4 change the rate of change?
testparm i.apoe4                          // baseline differences by APOE4
estat recovariance, correlation           // intercept-slope correlation
estat ic

* model-based annual rate of change by APOE4 dose, and pairwise differences
margins apoe4, dydx(time) saving("output/stata_slopes_apoe4", replace)
margins apoe4, dydx(time) pwcompare(effects)
* ... by APOE4 within each baseline diagnosis (population-averaged over age etc.)
margins apoe4, dydx(time) over(dx)

* population-average trajectories
margins apoe4, at(time = (0(1)6)) saving("output/stata_traj_apoe4", replace)
marginsplot, recast(line) recastci(rarea) ciopts(color(%20))          ///
    title("Predicted ${outcome} by APOE e4 dose")                      ///
    subtitle("Adjusted for age, sex, education, baseline diagnosis")   ///
    xtitle("Years since baseline") ytitle("${outcome}") legend(rows(1))
graph export "output/stata_fig_traj_apoe4.png", replace width(1600)

* 5d. Exploratory: does the APOE4 effect on slope differ by diagnosis?
mixed y c.time##i.apoe4_carrier##i.dx c.time##c.age_c c.time##i.female c.educ_c ///
      || rid: time, covariance(unstructured) reml
estimates store m_apoe_dx
testparm i.apoe4_carrier#i.dx#c.time
margins apoe4_carrier, dydx(time) over(dx)

*==============================================================================
* 6. PREDICT
*==============================================================================
* per-person OLS slopes (>= 3 visits), used to illustrate BLUP shrinkage
preserve
    bysort rid: gen nvis = _N
    keep if nvis >= 3
    statsby ols_slope = _b[time], by(rid) clear nodots: regress y time
    save "output/tmp_ols.dta", replace
restore
merge m:1 rid using "output/tmp_ols.dta", nogenerate keep(master match)
capture erase "output/tmp_ols.dta"
sort rid year

estimates restore m_main

predict xb_fixed,  xb                   // marginal (population) prediction
predict yhat_fit,  fitted                // conditional prediction incl. BLUPs
predict resid_c,   residuals             // conditional residuals
predict resid_std, rstandard             // standardised conditional residuals
predict re_slope re_int, reffects        // BLUPs: order = time, _cons
predict se_slope se_int, reses

* individual annual rate of change = fixed slope for the person's covariates
* + their random-slope BLUP. The fixed part is obtained from predict, xb at
* time+1 minus time (exact, because the model is linear in time).
predict xb_t0, xb
replace time = time + 1
predict xb_t1, xb
replace time = time - 1
gen fixed_slope = xb_t1 - xb_t0
gen indiv_slope = fixed_slope + re_slope
label variable indiv_slope "Predicted individual annual change (points/yr)"
drop xb_t0 xb_t1

summarize re_int re_slope indiv_slope resid_std
tabstat indiv_slope if first, by(apoe4) stat(n mean sd p25 p50 p75) format(%6.2f)
tabstat indiv_slope if first, by(dx)    stat(n mean sd p25 p50 p75) format(%6.2f)

* shrinkage check: BLUP-based slopes vs per-person OLS slopes
correlate indiv_slope ols_slope if first
summarize indiv_slope ols_slope if first & !missing(ols_slope)
di as text "(BLUP slopes are shrunk toward the covariate-specific mean: smaller SD)"

* 6a. Residual diagnostics
qnorm resid_std, title("Normal Q-Q: standardised conditional residuals")
graph export "output/stata_fig_qq.png", replace width(1200)
scatter resid_std yhat_fit, msize(vsmall) mcolor(%30) yline(0) ///
    xtitle("Fitted value (conditional)") ytitle("Standardised residual")
graph export "output/stata_fig_resid_fitted.png", replace width(1200)

* 6b. Random-effects distribution
scatter re_slope re_int if first, msize(vsmall) mcolor(%40)          ///
    xtitle("Random intercept BLUP") ytitle("Random slope BLUP")      ///
    title("Empirical Bayes estimates of random effects")
graph export "output/stata_fig_blups.png", replace width(1200)

* 6c. Observed vs fitted for 20 random participants with >= 4 visits
preserve
    bysort rid: gen nv = _N
    keep if nv >= 4
    set seed 20261002
    gen double u = runiform() if first
    bysort rid (u): replace u = u[1]
    egen rk = group(u)
    keep if rk <= 20
    sort rid time
    twoway (scatter y time, mcolor(gs8) msize(small))                    ///
           (line yhat_fit time, lcolor(navy))                           ///
           (line xb_fixed time, lcolor(cranberry) lpattern(dash)),      ///
           by(rid, cols(5) note("") legend(off)) ytitle("${outcome}")   ///
           xtitle("Years since baseline")
    graph export "output/stata_fig_individual_fits.png", replace width(1800)
restore

* 6d. Save subject-level predictions
preserve
    keep if first
    keep rid apoe4 dx age_bl female educ re_int re_slope indiv_slope ols_slope
    export delimited using "output/stata_individual_slopes.csv", replace
restore
save "output/adni_long_pred.dta", replace

*==============================================================================
* 7. SENSITIVITY ANALYSES FOR MISSING DATA (MNAR)
*==============================================================================
use "output/adni_long_all.dta", clear
keep if !missing(y, time)

* 7a. Complete-case: only participants retained to their last eligible visit
mixed y c.time##i.apoe4 c.time##c.age_c c.time##i.dx c.time##i.female c.educ_c ///
      if dropgrp == 0 || rid: time, covariance(unstructured) reml
estimates store m_cc
margins apoe4, dydx(time) saving("output/stata_slopes_cc", replace)

* 7b. Pattern-mixture model: intercept and slope differ by dropout pattern.
*     The marginal APOE4 slopes average over the observed pattern mix.
mixed y c.time##i.apoe4 c.time##c.age_c c.time##i.dx c.time##i.female c.educ_c ///
      c.time##i.dropgrp || rid: time, covariance(unstructured) reml
estimates store m_pmm
testparm i.dropgrp#c.time                 // do dropouts decline faster?
margins dropgrp, dydx(time)
margins apoe4,  dydx(time) saving("output/stata_slopes_pmm", replace)

* 7c. Side-by-side
estimates table m_main m_cc m_pmm, b(%8.3f) se(%8.3f) stats(N ll aic) ///
    title("Final model vs MNAR sensitivity models")

*==============================================================================
* 8. EXPORT RESULTS (CSV)
*==============================================================================
* slopes by APOE4 across the three analyses
tempfile acc
local first 1
foreach f in slopes_apoe4 slopes_cc slopes_pmm {
    use "output/stata_`f'.dta", clear
    gen str20 analysis = "`f'"
    if `first' == 0 append using `acc'
    save `acc', replace
    local first 0
}
use `acc', clear
rename (_m1 _margin _se _ci_lb _ci_ub _pvalue) (apoe4 slope se ci_lb ci_ub p)
keep analysis apoe4 slope se ci_lb ci_ub p
replace analysis = "MAR (main)"      if analysis == "slopes_apoe4"
replace analysis = "Complete-case"   if analysis == "slopes_cc"
replace analysis = "Pattern-mixture" if analysis == "slopes_pmm"
list, sepby(analysis) noobs
export delimited using "output/stata_slopes_by_apoe4.csv", replace

use "output/stata_traj_apoe4.dta", clear
export delimited using "output/stata_trajectories_apoe4.csv", replace

log close
*==============================================================================
* END
*==============================================================================
