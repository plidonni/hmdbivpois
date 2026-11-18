*! version 1.0.0  18nov2026
*! hmdbivpois: hidden Markov dropout bivariate Poisson regression for panel counts
*! Model: Carrieri, Gyrd-Hansen & Li Donni (2026, JRSS-A, doi:10.1093/jrsssa/qnag108)
*! Author: Paolo Li Donni, https://github.com/plidonni/hmdbivpois

/*
    Model
    -----
    Two latent variables:
      C_i  in {1..kc}   latent class, constant over time (persistent frailty),
                        prior P(C_i = u) multinomial logit on (classcov, y1_i0, y2_i0)
      V_it in {1..k}    latent state, first-order Markov chain (time-varying
                        heterogeneity); V_i1 multinomial logit on (initcov, y1_i0,
                        y2_i0); transitions: each row a logit with base =
                        staying, homogeneous by default; xtrans(varlist) adds
                        covariates measured at the period of arrival, tvtrans
                        period dummies (a free matrix per period); none with
                        notransition
    C and V are independent a priori.

    Measurement: y_it = (y1_it, y2_it) | C_i = u, V_it = v ~ bivariate Poisson in the
    Famoye (Lakshminarayana et al. 1999) form
        p(y1,y2) = Po(y1;mu1) Po(y2;mu2) [1 + lambda (exp(-y1)-exp(-q mu1)) (exp(-y2)-exp(-q mu2))]
        log mu_jit = x_jit b_j + g_j1 y1_i,t-1 + g_j2 y2_i,t-1 + alpha_jv + delta_ju,
        q = 1 - exp(-1), delta_j1 = 0
    Dropout (optional third equation): D_it = 1 if the individual leaves the panel
    after period t (e.g. death in the following interval)
        logit P(D_it = 1 | u, v) = xd_it bd + gd1 y1_i,t-1 + gd2 y2_i,t-1 + alpha_dv + delta_du
    Rows at risk of dropout: every estimation period except the last period of an
    individual with D = 0 there (end of follow-up, or censoring): in that period
    the dropout is not observed and the period does not contribute to the
    dropout likelihood (Bartolucci & Farcomeni 2015).
    Given (C_i, V_i1..V_iT) the periods are independent; y1 and y2 depend on each
    other within the period through lambda; the dropout depends on the counts
    only through the latent variables and the lags (shared-parameter MNAR).

    Without the dropout equation and with classes(1) the model is the hidden
    Markov bivariate Poisson model of Laudicella & Li Donni (2022) (hmbivpois).

    Lags are built by the command: the first observed period of each individual
    is the presample (y_i0) and is not an estimation period.
    lagtransform() applies g() to both lags, in the three equations, and to y_i0
    in the initial probabilities of the states and in the prior probabilities
    of the classes: linear (default), log1p = ln(1+y), logc(#) = ln(max(y,#)),
    logzero = ln(y) for y>0 plus a dummy for y=0, indicator = 1(y>0),
    topcode(#) = min(y,#).

    Computation: the pair (C, V) is a single hidden Markov chain on K = kc*k
    states, s = (u-1)*k + v, with transition matrix I(kc) # P and initial
    probabilities P(C=u) P(V_1=v); forward-backward with scaled recursions,
    vectorised over individuals by period.
    Estimation: EM (closed form for the transitions, Newton steps for the
    measurement and dropout equations and for the two multinomial logits), then
    direct maximisation of the log likelihood with the analytic score (Fisher
    identity) if the EM does not converge within emiterate() or if technique()
    is specified. Variance: inverse of the observed information, from central
    differences of the analytic score; vce(robust): sandwich clustered on
    individuals. Lambda is kept inside the interval that makes the pmf
    nonnegative at the fitted means (all states and classes).
    Labels: states in increasing order of alpha_1v (sortby(eq1), default) or
    alpha_2v (sortby(eq2)); classes in increasing order of the dropout shift
    delta_du (of delta_1u, or delta_2u with sortby(eq2), without the dropout
    equation). sortby(none) keeps the order of the estimation.

    starts(#): EM from # starting values (the default or from(), and #-1
    random ones), for startiterate() iterations each; the EM goes on from the
    start with the largest log likelihood (Carrieri et al. 2026 use multiple
    starts against local maxima).

    estat lrtrans tests notransition (states constant over time) against the
    model with homogeneous transitions: LR statistic, chi-bar-squared null
    distribution with weights estimated by simulation (Bartolucci 2006).
    statevary(varlist) and lagvary(depvars), inside an equation, make the
    coefficients of those regressors or lagged counts specific to the state
    (equations <depvar>_s1, _s2, ...), as in hmbivpois; the classes keep
    shifting only the intercepts.

    dropcounts(current|both) puts the counts of the period (observed before
    the dropout) in the dropout equation, in place of or with the lagged ones;
    classeq(dropout|counts) removes the class shifts from the count equations
    or from the dropout equation.

    estat stability checks the dynamics of the lagged counts (coefficients,
    spectral radius or index by latent type, deterministic and simulated
    paths), as in hmbivpois.
*/

program define hmdbivpois, eclass
    version 16
    if replay() {
        if "`e(cmd)'" != "hmdbivpois" error 301
        Display `0'
        exit
    }
    // entry point of predict (hmdbivpois_p.ado)
    gettoken sub rest : 0
    if `"`sub'"' == "__predict" {
        if "`e(cmd)'" != "hmdbivpois" error 301
        Predict `rest'
        exit
    }
    // entry point of estat (hmdbivpois_estat.ado)
    if `"`sub'"' == "__estat" {
        if "`e(cmd)'" != "hmdbivpois" error 301
        Estat `rest'
        exit
    }
    Estimate `0'
end

program define Estimate, eclass sortpreserve
    version 16
    local cmdline `"hmdbivpois `0'"'

    // ---------------------------------------------------------------- parsing
    gettoken eq1 0 : 0, match(paren)
    if "`paren'" != "(" {
        di as err "equation 1 must be enclosed in parentheses: (depvar1 = indepvars1)"
        exit 198
    }
    gettoken eq2 0 : 0, match(paren)
    if "`paren'" != "(" {
        di as err "equation 2 must be enclosed in parentheses: (depvar2 = indepvars2)"
        exit 198
    }
    // optional third equation: dropout
    local eq3
    gettoken tok rest : 0, match(paren)
    if "`paren'" == "(" {
        local eq3 `"`tok'"'
        local 0 `"`rest'"'
    }
    // inside each equation: statevary() lists the regressors and lagvary()
    // the lagged outcomes whose coefficients differ across the states
    ParseEq `eq1'
    local y1 `s(dep)'
    local x1 `s(indep)'
    local svspec1 `s(statevary)'
    local lvspec1 `s(lagvary)'
    ParseEq `eq2'
    local y2 `s(dep)'
    local x2 `s(indep)'
    local svspec2 `s(statevary)'
    local lvspec2 `s(lagvary)'
    local yd
    local xd
    local svspecd
    local lvspecd
    if `"`eq3'"' != "" {
        ParseEq `eq3'
        local yd `s(dep)'
        local xd `s(indep)'
        local svspecd `s(statevary)'
        local lvspecd `s(lagvary)'
    }
    local hasD = ("`yd'" != "")
    if "`y1'" == "`y2'" {
        di as err "the two dependent variables must be different"
        exit 198
    }
    if `hasD' {
        if "`yd'" == "`y1'" | "`yd'" == "`y2'" {
            di as err "the dropout indicator must differ from the two counts"
            exit 198
        }
    }

    syntax [if] [in], ID(varname numeric) TIME(varname numeric)           ///
        [ K(integer 0) CLASSES(integer 1) INITcov(varlist numeric)        ///
          CLASSCOV(varlist numeric) EMITERate(integer 500)                ///
          EMTOLerance(real 1e-8) ITERate(integer 300) TECHnique(string)   ///
          FROM(name) SORTby(string) VCE(string) Level(cilevel)            ///
          NOLOg EMONLY NOTRansition XTRans(varlist numeric fv) TVTRans    ///
          LAGTransform(string) STARTS(integer 1) STARTITerate(integer -1) ///
          STARTSEED(integer 20261117) DROPCounts(string) CLASSEQ(string) ]

    if `k' == 0 {
        di as err "option k() required"
        exit 198
    }
    if `k' < 1 | `k' > 20 {
        di as err "k() must be between 1 and 20"
        exit 198
    }
    if `classes' < 1 | `classes' > 10 {
        di as err "classes() must be between 1 and 10"
        exit 198
    }
    if `k'*`classes' > 60 {
        di as err "k() times classes() may not exceed 60"
        exit 198
    }
    if `classes' == 1 & "`classcov'" != "" {
        di as txt "note: classcov() ignored with classes(1)"
        local classcov
    }
    if `k' == 1 & "`initcov'" != "" {
        di as txt "note: initcov() ignored with k(1)"
        local initcov
    }
    local notr = ("`notransition'" != "")
    if `k' == 1 & `notr' {
        di as txt "note: notransition ignored with k(1)"
        local notr = 0
    }
    // transitions by transition row: covariates measured at the period of
    // arrival (xtrans), period dummies (tvtrans)
    local het = ("`xtrans'" != "" | "`tvtrans'" != "")
    if `k' == 1 & `het' {
        di as txt "note: xtrans() and tvtrans ignored with k(1)"
        local het = 0
        local xtrans
        local tvtrans
    }
    local tvt = ("`tvtrans'" != "")
    if `notr' & `het' {
        di as err "notransition may not be combined with xtrans() or tvtrans"
        exit 198
    }
    // direct maximisation: 0 = never, 1 = only if EM did not converge, 2 = always
    local doml = 1
    if "`technique'" != "" local doml = 2
    if "`emonly'" != ""    local doml = 0
    if "`technique'" == "" local technique "bfgs"
    if !inlist("`technique'", "bfgs", "dfp", "nr") {
        di as err "technique() must be bfgs, dfp or nr"
        exit 198
    }
    if "`sortby'" == "" local sortby "eq1"
    if !inlist("`sortby'", "eq1", "eq2", "none") {
        di as err "sortby() must be eq1, eq2 or none"
        exit 198
    }
    if "`vce'" == "" local vce "oim"
    if !inlist("`vce'", "oim", "robust") {
        di as err "vce() must be oim or robust"
        exit 198
    }
    if `emiterate' < 0 | `iterate' < 0 {
        di as err "emiterate() and iterate() must be nonnegative"
        exit 198
    }
    if "`from'" != "" confirm matrix `from'
    local dolog = ("`nolog'" == "")
    // multiple starting values: EM for startiterate() iterations from each
    // (default: emiterate()), then the EM goes on from the best
    if `starts' < 1 | `starts' > 1000 {
        di as err "starts() must be between 1 and 1000"
        exit 198
    }
    if `startiterate' < 0 local startiterate = `emiterate'
    if `starts' > 1 & `startiterate' < 1 {
        di as err "starts() requires startiterate() and emiterate() of at least 1"
        exit 198
    }

    // transformation g() of the lagged counts, in the three equations, and of
    // the presample counts, in the initial probabilities of the states and in
    // the prior probabilities of the classes
    local lagt : subinstr local lagtransform " " "", all
    if "`lagt'" == "" local lagt "linear"
    local lagfull "`lagt'"
    local lagarg
    if regexm("`lagt'", "^(logc|topcode)\(([0-9.]+)\)$") {
        local lagarg = regexs(2)
        local lagt   = regexs(1)
    }
    if !inlist("`lagt'", "linear", "log1p", "logc", "logzero", "indicator", "topcode") {
        di as err "lagtransform() must be linear, log1p, logc(#), logzero, indicator or topcode(#)"
        exit 198
    }
    if inlist("`lagt'", "logc", "topcode") & "`lagarg'" == "" {
        di as err "lagtransform(`lagt'(#)) requires a number"
        exit 198
    }
    if "`lagt'" == "logc" {
        if !(`lagarg' > 0 & `lagarg' <= 1) {
            di as err "lagtransform(logc(#)): # must be in (0, 1]"
            exit 198
        }
    }
    if "`lagt'" == "topcode" {
        if !(`lagarg' >= 1) {
            di as err "lagtransform(topcode(#)): # must be at least 1"
            exit 198
        }
    }

    // counts in the dropout equation: lagged (default), current (those of the
    // period after which the individual may leave) or both; equations where
    // the classes enter: all (default), dropout only, counts only
    if "`dropcounts'" == "" | "`dropcounts'" == "lagged" local dropcounts "lag"
    if !inlist("`dropcounts'", "lag", "current", "both") {
        di as err "dropcounts() must be lag, current or both"
        exit 198
    }
    if "`dropcounts'" != "lag" & !`hasD' {
        di as err "dropcounts() requires the dropout equation"
        exit 198
    }
    if "`classeq'" == "" local classeq "all"
    if !inlist("`classeq'", "all", "dropout", "counts") {
        di as err "classeq() must be all, dropout or counts"
        exit 198
    }
    if "`classeq'" != "all" & `classes' == 1 {
        di as txt "note: classeq() ignored with classes(1)"
        local classeq "all"
    }
    if "`classeq'" != "all" & !`hasD' {
        di as err "classeq() requires the dropout equation"
        exit 198
    }
    // read by the Mata setup
    local clseq   "`classeq'"
    local dropcnt "`dropcounts'"

    // state-specific coefficients
    if `k' == 1 & "`svspec1'`svspec2'`svspecd'`lvspec1'`lvspec2'`lvspecd'" != "" {
        di as txt "note: statevary() and lagvary() ignored with k(1)"
        foreach j in 1 2 d {
            local svspec`j'
            local lvspec`j'
        }
    }
    if "`lvspecd'" != "" & "`dropcounts'" == "current" {
        di as err "lagvary() in the dropout equation requires the lagged counts there (dropcounts(lag) or dropcounts(both))"
        exit 198
    }
    foreach j in 1 2 d {
        // lagS`j': lags (1 = first count, 2 = second) whose coefficients are
        // state specific in equation `j'
        local lagS`j'
        foreach v of local lvspec`j' {
            capture unab vv : `v', max(1)
            if _rc local vv `v'
            if "`vv'" == "`y1'" {
                local lagS`j' `lagS`j'' 1
            }
            else if "`vv'" == "`y2'" {
                local lagS`j' `lagS`j'' 2
            }
            else {
                di as err "lagvary(): `v' is not one of the two counts"
                exit 198
            }
        }
        local lagS`j' : list uniq lagS`j'
    }

    // ---------------------------------------------------------------- sample
    // touse: rows usable as observed periods (outcomes, dropout, id, time)
    // xok  : rows where all regressors are also observed (estimation periods)
    marksample touse, novarlist
    markout `touse' `y1' `y2' `id' `time'
    if `hasD' markout `touse' `yd'

    tempvar xok first last gap gap2 L1 L2 es Tn negT risk
    qui gen byte `xok' = `touse'
    if "`x1'`x2'`xd'" != "" {
        fvrevar `x1' `x2' `xd', list
        markout `xok' `r(varlist)'
    }
    if "`initcov'`classcov'" != "" markout `xok' `initcov' `classcov'
    if "`xtrans'" != "" {
        fvrevar `xtrans', list
        markout `xok' `r(varlist)'
    }

    capture assert `y1' >= 0 & `y1' == floor(`y1') & ///
                   `y2' >= 0 & `y2' == floor(`y2') if `touse'
    if _rc {
        di as err "dependent variables must be nonnegative integers"
        exit 459
    }
    if `hasD' {
        capture assert inlist(`yd', 0, 1) if `touse'
        if _rc {
            di as err "the dropout indicator `yd' must be 0 or 1"
            exit 459
        }
    }
    capture assert `time' == floor(`time') if `touse'
    if _rc {
        di as err "`time' must be integer-valued"
        exit 459
    }

    sort `touse' `id' `time'
    qui by `touse' `id': gen byte `gap' = (`time' != `time'[_n-1] + 1) ///
        if `touse' & _n > 1
    qui count if `gap' == 1
    if r(N) > 0 {
        di as err "repeated or non-consecutive values of `time' within `id';"
        di as err "gaps in the observed periods are not allowed"
        exit 459
    }
    qui by `touse' `id': gen byte `first' = (_n == 1) if `touse'
    qui by `touse' `id': gen byte `last'  = (_n == _N) if `touse'
    qui by `touse' `id': gen double `L1' = `y1'[_n-1] if `touse'
    qui by `touse' `id': gen double `L2' = `y2'[_n-1] if `touse'
    tempvar G1 G2 Z1 Z2
    forvalues j = 1/2 {
        if "`lagt'" == "linear"    qui gen double `G`j'' = `L`j''
        if "`lagt'" == "log1p"     qui gen double `G`j'' = ln(1 + `L`j'')
        if "`lagt'" == "logc"      qui gen double `G`j'' = ln(max(`L`j'', `lagarg')) if `L`j'' < .
        if "`lagt'" == "indicator" qui gen double `G`j'' = (`L`j'' > 0) if `L`j'' < .
        if "`lagt'" == "topcode"   qui gen double `G`j'' = min(`L`j'', `lagarg') if `L`j'' < .
        if "`lagt'" == "logzero" {
            qui gen double `G`j'' = cond(`L`j'' > 0, ln(`L`j''), 0) if `L`j'' < .
            qui gen double `Z`j'' = (`L`j'' == 0) if `L`j'' < .
        }
    }
    local lagvars  `G1' `G2'
    local lagnames `=substr("lag_`y1'", 1, 32)' `=substr("lag_`y2'", 1, 32)'
    if "`lagt'" == "logzero" {
        local lagvars  `lagvars' `Z1' `Z2'
        local lagnames `lagnames' `=substr("lag0_`y1'", 1, 32)' `=substr("lag0_`y2'", 1, 32)'
    }
    // the same transformation of the current counts, for dropcounts()
    tempvar C1 C2 CZ1 CZ2
    local curvars
    local curnames
    if "`dropcounts'" != "lag" {
        forvalues j = 1/2 {
            if "`lagt'" == "linear"    qui gen double `C`j'' = `y`j'' if `touse'
            if "`lagt'" == "log1p"     qui gen double `C`j'' = ln(1 + `y`j'') if `touse'
            if "`lagt'" == "logc"      qui gen double `C`j'' = ln(max(`y`j'', `lagarg')) if `touse'
            if "`lagt'" == "indicator" qui gen double `C`j'' = (`y`j'' > 0) if `touse'
            if "`lagt'" == "topcode"   qui gen double `C`j'' = min(`y`j'', `lagarg') if `touse'
            if "`lagt'" == "logzero" {
                qui gen double `C`j'' = cond(`y`j'' > 0, ln(`y`j''), 0) if `touse'
                qui gen double `CZ`j'' = (`y`j'' == 0) if `touse'
            }
        }
        local curvars  `C1' `C2'
        local curnames `=substr("cur_`y1'", 1, 32)' `=substr("cur_`y2'", 1, 32)'
        if "`lagt'" == "logzero" {
            local curvars  `curvars' `CZ1' `CZ2'
            local curnames `curnames' `=substr("cur0_`y1'", 1, 32)' `=substr("cur0_`y2'", 1, 32)'
        }
    }
    qui gen byte `es' = `touse' & `first' == 0 & `xok'

    // estimation periods must be consecutive within individual
    sort `es' `id' `time'
    qui by `es' `id': gen byte `gap2' = (`time' != `time'[_n-1] + 1) ///
        if `es' & _n > 1
    qui count if `gap2' == 1
    if r(N) > 0 {
        di as err "missing regressors inside an individual's sequence of periods;"
        di as err "estimation periods must be consecutive"
        exit 459
    }

    // dropout: D = 1 only in the last observed period of an individual. Rows
    // at risk: estimation periods, except the last one of an individual when
    // D = 0 there (dropout not observed)
    local ndrop = 0
    local nlost = 0
    if `hasD' {
        qui count if `touse' & `yd' == 1 & `last' == 0
        if r(N) > 0 {
            di as err "`yd' = 1 in a period that is not the last observed period of the individual;"
            di as err "the dropout indicator may be 1 only in the period after which the individual leaves the panel"
            exit 459
        }
        qui count if `touse' & `yd' == 1 & `es' == 0
        local nlost = r(N)
        sort `es' `id' `time'
        qui by `es' `id': gen byte `risk' = `es' & (_n < _N | `yd' == 1)
        qui count if `es' & `yd' == 1
        local ndrop = r(N)
        if `ndrop' == 0 {
            di as err "no dropouts (`yd' = 1) in the estimation periods"
            exit 2000
        }
    }
    else qui gen byte `risk' = 0

    sort `es' `id' `time'
    qui by `es' `id': gen long `Tn' = _N if `es'
    // transition rows: estimation periods after the first one of each individual
    tempvar tr
    qui by `es' `id': gen byte `tr' = `es' & _n > 1
    if `het' {
        qui count if `tr'
        if r(N) == 0 {
            di as err "no transitions observed: xtrans() and tvtrans need T > 1"
            exit 459
        }
    }
    qui count if `es'
    local N = r(N)
    if `N' == 0 error 2000
    if inlist("`lagt'", "indicator", "logzero") {
        // the lag must take both zero and positive values in the sample
        forvalues j = 1/2 {
            qui count if `es' & `L`j'' == 0
            if r(N) == 0 | r(N) == `N' {
                di as err "lagtransform(`lagt'): lagged `y`j'' is never zero or always zero in the sample"
                exit 459
            }
        }
        if "`dropcounts'" != "lag" {
            forvalues j = 1/2 {
                qui count if `risk' & `y`j'' == 0
                local nz = r(N)
                qui count if `risk'
                if `nz' == 0 | `nz' == r(N) {
                    di as err "lagtransform(`lagt'): `y`j'' is never zero or always zero in the rows at risk"
                    exit 459
                }
            }
        }
    }

    // ------------------------------------------------- regressors (fv aware)
    // cnam`j': names of the terms (regressors, then the two lags);
    // cvar`j': one variable per term
    foreach j in 1 2 d {
        local xnames`j'
        local xvars`j'
        if "`x`j''" != "" {
            qui _rmcoll `x`j'' if `es', expand
            local full `r(varlist)'
            foreach t of local full {
                if !ustrregexm("`t'", "(^|#)[0-9]*[bo]\.") {
                    local xnames`j' `xnames`j'' `t'
                }
            }
            foreach t of local xnames`j' {
                if !ustrregexm("`t'", "[.#]") {
                    local xvars`j' `xvars`j'' `t'
                    continue
                }
                tempvar tv
                TermVar `tv' `t'
                local xvars`j' `xvars`j'' `tv'
            }
        }
        // regressors with state-specific coefficients: the terms of
        // statevary(), which must be among the regressors of the equation
        local svterm`j'
        if "`svspec`j''" != "" {
            fvexpand `svspec`j'' if `es'
            local svfull `r(varlist)'
            foreach t of local svfull {
                if !ustrregexm("`t'", "(^|#)[0-9]*[bo]\.") {
                    local isin : list t in xnames`j'
                    if !`isin' {
                        di as err "statevary(): `t' is not among the regressors of the equation of `y`j''"
                        exit 198
                    }
                    local svterm`j' `svterm`j'' `t'
                }
            }
        }
        // common terms (cnam, cvar) and state-specific terms (snam, svar):
        // regressors, then the lags (in the dropout equation, those of
        // dropcounts()), then the current counts (common)
        local cnam`j'
        local cvar`j'
        local snam`j'
        local svar`j'
        local i = 0
        foreach t of local xnames`j' {
            local ++i
            local v : word `i' of `xvars`j''
            local isin : list t in svterm`j'
            if `isin' {
                local snam`j' `snam`j'' `t'
                local svar`j' `svar`j'' `v'
            }
            else {
                local cnam`j' `cnam`j'' `t'
                local cvar`j' `cvar`j'' `v'
            }
        }
        local withlag = ("`j'" != "d" | inlist("`dropcounts'", "lag", "both"))
        if `withlag' {
            local i = 0
            foreach t of local lagnames {
                local ++i
                local v : word `i' of `lagvars'
                // features 1 and 3 belong to the lag of the first count, 2
                // and 4 to the second
                local li = 2 - mod(`i', 2)
                local isin : list li in lagS`j'
                if `isin' {
                    local snam`j' `snam`j'' `t'
                    local svar`j' `svar`j'' `v'
                }
                else {
                    local cnam`j' `cnam`j'' `t'
                    local cvar`j' `cvar`j'' `v'
                }
            }
        }
        if "`j'" == "d" & inlist("`dropcounts'", "current", "both") {
            local cnamd `cnamd' `curnames'
            local cvard `cvard' `curvars'
        }
    }
    if !`hasD' {
        local cnamd
        local cvard
        local snamd
        local svard
    }

    // covariates of the transition probabilities (measured at the period of
    // arrival), and the periods of arrival for tvtrans
    local xnamesT
    local xvarsT
    if "`xtrans'" != "" {
        qui _rmcoll `xtrans' if `tr', expand
        local full `r(varlist)'
        foreach t of local full {
            if !ustrregexm("`t'", "(^|#)[0-9]*[bo]\.") {
                local xnamesT `xnamesT' `t'
            }
        }
        foreach t of local xnamesT {
            if !ustrregexm("`t'", "[.#]") {
                local xvarsT `xvarsT' `t'
                continue
            }
            tempvar tv
            TermVar `tv' `t'
            local xvarsT `xvarsT' `tv'
        }
    }
    local tvals
    local tvals0
    if `het' qui levelsof `time' if `tr', local(tvals0)
    if `tvt' local tvals `tvals0'

    // ------------------------------------------------- starting values (k=1)
    tempname bs1 bs2 bsd
    if "`from'" == "" {
        qui poisson `y1' `cvar1' `svar1' if `es'
        matrix `bs1' = e(b)
        qui poisson `y2' `cvar2' `svar2' if `es'
        matrix `bs2' = e(b)
        if `hasD' {
            capture qui logit `yd' `cvard' `svard' if `risk'
            if _rc {
                di as err "starting values: logit of `yd' on the rows at risk failed"
                exit 430
            }
            matrix `bsd' = e(b)
        }
    }

    // ------------------------------------------------------------ estimation
    qui gen long `negT' = -`Tn'
    sort `es' `negT' `id' `time'

    mata: hmdbp_main("`y1'", "`y2'", "`cvar1'", "`cvar2'", "`yd'",      ///
        "`cvard'", "`risk'", "`initcov' `lagvars'",                         ///
        "`classcov' `lagvars'", "`id'", "`es'", `k', `classes',             ///
        "`bs1'", "`bs2'", "`bsd'", "`from'", `emiterate', `emtolerance',    ///
        `iterate', "`technique'", `dolog', "`sortby'", "`vce'", `doml',     ///
        `notr', "`xvarsT'", "`tr'", "`time'", `het', `tvt', `starts',      ///
        `startiterate', `startseed')

    // --------------------------------------------------------------- results
    tempname b V Pi om omt al de pi0 pcl wcl lrng Pit
    if `het' {
        matrix `Pit' = __hmdbp_Pit
        capture matrix drop __hmdbp_Pit __hmdbp_tvals
    }
    matrix `b'   = __hmdbp_b
    matrix `V'   = __hmdbp_V
    matrix `Pi'  = __hmdbp_Pi
    matrix `om'  = __hmdbp_om
    matrix `omt' = __hmdbp_omt
    matrix `al'  = __hmdbp_al
    matrix `de'  = __hmdbp_de
    matrix `pi0' = __hmdbp_pi0
    matrix `pcl' = __hmdbp_pcl
    matrix `wcl' = __hmdbp_wcl
    matrix `lrng' = __hmdbp_lrng
    capture matrix drop __hmdbp_b __hmdbp_V __hmdbp_Pi __hmdbp_om ///
        __hmdbp_omt __hmdbp_al __hmdbp_de __hmdbp_pi0 __hmdbp_pcl   ///
        __hmdbp_wcl __hmdbp_lrng

    // names: in each equation the regressors, the lags, the state intercepts
    // and the class shifts (class 1 is the base)
    local eqs1 `y1'
    local eqs2 `y2'
    local eqsd `yd'
    local names
    foreach j in 1 2 d {
        if "`j'" == "d" & !`hasD' continue
        foreach v of local cnam`j' {
            local names `names' `eqs`j'':`v'
        }
        forvalues c = 1/`k' {
            local names `names' `eqs`j'':state`c'
        }
        local hc = (`classes' > 1)
        if "`j'" != "d" & "`classeq'" == "dropout" local hc = 0
        if "`j'" == "d" & "`classeq'" == "counts"  local hc = 0
        if `hc' {
            forvalues u = 2/`classes' {
                local names `names' `eqs`j'':class`u'
            }
        }
        // state-specific coefficients: equations <depvar>_s1, _s2, ...
        if "`snam`j''" != "" {
            forvalues c = 1/`k' {
                local eqs = substr("`eqs`j''", 1, 32 - strlen("_s`c'")) + "_s`c'"
                foreach v of local snam`j' {
                    local names `names' `eqs':`v'
                }
            }
        }
        if "`j'" == "2" local names `names' lambda:_cons
    }
    if `k' > 1 {
        forvalues c = 2/`k' {
            foreach v of local initcov {
                local names `names' init`c':`v'
            }
            foreach v of local lagnames {
                local names `names' init`c':`v'
            }
            local names `names' init`c':_cons
        }
        // transitions: covariates of xtrans(), then the constant or, with
        // tvtrans, one dummy per period of arrival
        local trn `xnamesT'
        if `tvt' {
            foreach tv of local tvals {
                local trn `trn' t_`=subinstr("`tv'", "-", "m", .)'
            }
        }
        else local trn `trn' _cons
        if !`notr' {
            forvalues c = 1/`k' {
                forvalues d = 1/`k' {
                    if `c' != `d' {
                        foreach v of local trn {
                            local names `names' tr`c'_`d':`v'
                        }
                    }
                }
            }
        }
    }
    if `classes' > 1 {
        forvalues u = 2/`classes' {
            foreach v of local classcov {
                local names `names' class`u':`v'
            }
            foreach v of local lagnames {
                local names `names' class`u':`v'
            }
            local names `names' class`u':_cons
        }
    }
    matrix colnames `b' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'

    local snames
    forvalues c = 1/`k' {
        local snames `snames' state`c'
    }
    local cnames
    forvalues u = 1/`classes' {
        local cnames `cnames' class`u'
    }
    local eqnames `y1' `y2'
    if `hasD' local eqnames `eqnames' `yd'
    matrix rownames `Pi'  = `snames'
    matrix colnames `Pi'  = `snames'
    if `het' {
        // one row per period of arrival: the k*k probabilities, by row of the
        // transition matrix
        local pn
        forvalues c = 1/`k' {
            forvalues d = 1/`k' {
                local pn `pn' p`c'_`d'
            }
        }
        local tn
        foreach tv of local tvals0 {
            local tn `tn' t_`=subinstr("`tv'", "-", "m", .)'
        }
        matrix colnames `Pit' = `pn'
        matrix rownames `Pit' = `tn'
    }
    matrix colnames `om'  = `snames'
    matrix rownames `om'  = share
    matrix colnames `omt' = `snames'
    matrix colnames `pi0' = `snames'
    matrix rownames `pi0' = prob
    matrix rownames `al'  = `snames'
    matrix colnames `al'  = `eqnames'
    matrix rownames `de'  = `cnames'
    matrix colnames `de'  = `eqnames'
    matrix colnames `pcl' = `cnames'
    matrix rownames `pcl' = prior
    matrix colnames `wcl' = `cnames'
    matrix rownames `wcl' = share
    matrix colnames `lrng' = lower upper
    matrix rownames `lrng' = lambda
    local tnames
    forvalues t = 1/`=rowsof(`omt')' {
        local tnames `tnames' t`t'
    }
    matrix rownames `omt' = `tnames'

    ereturn post `b' `V', esample(`es') obs(`N')
    ereturn scalar N_g       = __hmdbp_n
    ereturn scalar k         = `k'
    ereturn scalar classes   = `classes'
    ereturn scalar ll        = __hmdbp_ll
    ereturn scalar npar      = __hmdbp_npar
    ereturn scalar rank      = __hmdbp_rank
    ereturn scalar bic       = -2*__hmdbp_ll + __hmdbp_npar*ln(__hmdbp_n)
    ereturn scalar bic_N     = -2*__hmdbp_ll + __hmdbp_npar*ln(`N')
    ereturn scalar aic       = -2*__hmdbp_ll + 2*__hmdbp_npar
    ereturn scalar converged = __hmdbp_conv
    ereturn scalar em_iter   = __hmdbp_emit
    ereturn scalar eval_only = (`emiterate' == 0 & `doml' == 0)
    ereturn scalar starts    = `starts'
    if `starts' > 1 {
        local sb = __hmdbp_sbest
        ereturn scalar start_best = `sb'
        ereturn scalar start_iter = `startiterate'
        ereturn scalar start_seed = `startseed'
        // starts that reach the best log likelihood (within 0.1)
        tempname lls
        matrix `lls' = __hmdbp_lls
        local nb = 0
        forvalues r = 1/`starts' {
            if `lls'[`r', 1] < . {
                if `lls'[`r', 1] >= `lls'[`sb', 1] - 0.1 local ++nb
            }
        }
        ereturn scalar n_best = `nb'
        matrix colnames `lls' = ll iterations
        local rn
        forvalues r = 1/`starts' {
            local rn `rn' start`r'
        }
        matrix rownames `lls' = `rn'
        ereturn matrix ll_starts = `lls'
        capture matrix drop __hmdbp_lls
        capture scalar drop __hmdbp_sbest
    }
    ereturn scalar g_min     = __hmdbp_gmin
    ereturn scalar g_avg     = __hmdbp_gavg
    ereturn scalar g_max     = __hmdbp_gmax
    ereturn scalar lambda_lo = __hmdbp_lamlo
    ereturn scalar lambda_hi = __hmdbp_lamhi
    ereturn scalar n_bound   = __hmdbp_nbnd
    ereturn scalar N_drop    = `ndrop'
    ereturn scalar N_risk    = __hmdbp_nrisk
    ereturn scalar N_droplost = `nlost'
    if "`vce'" == "robust" {
        ereturn scalar N_clust = __hmdbp_n
        ereturn local clustvar "`id'"
        ereturn local vcetype  "Robust"
    }
    capture scalar drop __hmdbp_n __hmdbp_ll __hmdbp_npar __hmdbp_rank  ///
        __hmdbp_conv __hmdbp_emit __hmdbp_gmin __hmdbp_gavg __hmdbp_gmax ///
        __hmdbp_lamlo __hmdbp_lamhi __hmdbp_nbnd __hmdbp_nrisk
    ereturn matrix Pi           = `Pi'
    ereturn matrix omega        = `om'
    ereturn matrix omega_t      = `omt'
    ereturn matrix alpha        = `al'
    ereturn matrix delta        = `de'
    ereturn matrix pi0          = `pi0'
    ereturn matrix pclass       = `pcl'
    ereturn matrix wclass       = `wcl'
    ereturn matrix lambda_range = `lrng'
    ereturn local vce          "`vce'"
    ereturn local technique    "`technique'"
    if `k' == 1        ereturn local transitions "none (k = 1)"
    else if `notr'     ereturn local transitions "none"
    else if `het'      ereturn local transitions "heterogeneous"
    else               ereturn local transitions "homogeneous"
    ereturn local xtrans       "`xtrans'"
    ereturn local xtrans_terms "`xnamesT'"
    ereturn local tvtrans      "`tvtrans'"
    if `het' ereturn matrix Pi_t = `Pit'

    ereturn local sortby       "`sortby'"
    ereturn local initcov      "`initcov'"
    ereturn local classcov     "`classcov'"
    ereturn local dropout      "`yd'"
    ereturn local tvar         "`time'"
    ereturn local ivar         "`id'"
    ereturn local depvar       "`y1' `y2'"
    if `hasD' ereturn local title "Hidden Markov dropout bivariate Poisson regression"
    else      ereturn local title "Hidden Markov bivariate Poisson regression with latent classes"
    ereturn local common1      "`cnam1'"
    ereturn local common2      "`cnam2'"
    ereturn local commond      "`cnamd'"
    ereturn local lagnames     "`lagnames'"
    ereturn local lagtransform "`lagfull'"
    ereturn local dropcounts   "`dropcounts'"
    ereturn local statevary1   "`snam1'"
    ereturn local statevary2   "`snam2'"
    ereturn local statevaryd   "`snamd'"
    ereturn local classeq      "`classeq'"
    ereturn local curnames     "`curnames'"
    ereturn local marginsnotok "_ALL"
    ereturn local predict      "hmdbivpois_p"
    ereturn local estat_cmd    "hmdbivpois_estat"
    ereturn local cmdline      `"`cmdline'"'
    ereturn local cmd          "hmdbivpois"

    Display, level(`level')
end

// Variable of one term in factor-variable notation, generated as `name':
// level.var gives the indicator of the level (missing where var is missing),
// c.var and var the variable, # the product of the parts. Other operators are
// left to fvrevar.
program define TermVar
    version 16
    args name term
    mata: st_local("tex", hmdbp_termexpr(`"`term'"'))
    if `"`tex'"' != "" {
        qui gen double `name' = `tex'
    }
    else {
        fvrevar `term'
        qui gen double `name' = `r(varlist)'
    }
end

program define ParseEq, sclass
    version 16
    sreturn clear
    local eq : subinstr local 0 "=" " ", all
    gettoken dep rest : eq
    capture unab dep : `dep', max(1)
    if _rc {
        di as err "`dep': dependent variable not found"
        exit 111
    }
    confirm numeric variable `dep'
    local 0 `rest'
    syntax [varlist(numeric fv default=none)]                             ///
        [, STATEVary(varlist numeric fv) LAGVary(namelist) ]
    sreturn local dep       "`dep'"
    sreturn local indep     "`varlist'"
    sreturn local statevary "`statevary'"
    sreturn local lagvary   "`lagvary'"
end

program define Display
    version 16
    syntax [, Level(cilevel)]
    di _n as txt "`e(title)'"
    di _n as txt "Dependence: " as res "Famoye (Sarmanov)" _col(49) as txt ///
        "Number of obs" _col(67) "=" _col(69) as res %10.0fc e(N)
    di as txt "Group variable: " as res abbrev("`e(ivar)'", 16) _col(49) as txt ///
        "Number of groups" _col(67) "=" _col(69) as res %10.0fc e(N_g)
    di as txt "Latent states:  " as res e(k) _col(49) as txt "Obs per group:"
    di as txt "Latent classes: " as res e(classes) ///
        _col(63) as txt "min" _col(67) "=" _col(69) as res %10.0fc e(g_min)
    if "`e(dropout)'" != "" {
        di as txt "Dropout: " as res abbrev("`e(dropout)'", 16) as txt ///
            " (" as res e(N_drop) as txt " events)" ///
            _col(63) as txt "avg" _col(67) "=" _col(69) as res %10.1f e(g_avg)
    }
    else {
        di as txt "Dropout: " as res "ignorable (no dropout equation)" ///
            _col(63) as txt "avg" _col(67) "=" _col(69) as res %10.1f e(g_avg)
    }
    di as txt "Transitions: " as res "`e(transitions)'" ///
        _col(63) as txt "max" _col(67) "=" _col(69) as res %10.0fc e(g_max)
    di as txt "Lagged outcomes: " as res "`e(lagtransform)'"
    if "`e(dropcounts)'" == "current" {
        di as txt "Counts in the dropout equation: " as res "current (period of the dropout)"
    }
    if "`e(dropcounts)'" == "both" {
        di as txt "Counts in the dropout equation: " as res "lagged and current"
    }
    local dv `e(depvar)' `e(dropout)'
    local i = 0
    foreach j in 1 2 d {
        if "`j'" != "d" local ++i
        else            local i = 3
        if "`e(statevary`j')'" != "" {
            local dn : word `i' of `dv'
            di as txt "State-specific coefficients (`dn'): " as res "`e(statevary`j')'"
        }
    }
    if "`e(classeq)'" == "dropout" {
        di as txt "Classes: " as res "in the dropout equation only"
    }
    if "`e(classeq)'" == "counts" {
        di as txt "Classes: " as res "in the count equations only"
    }
    if e(starts) > 1 & e(starts) < . {
        di as txt "Starting values: " as res e(starts) as txt " (best: start " ///
            as res e(start_best) as txt "; " as res e(n_best) as txt ///
            " within 0.1 of its log likelihood)"
    }
    di as txt "Log likelihood = " as res %12.0g e(ll) ///
        _col(49) as txt "BIC" _col(67) "=" _col(69) as res %10.0g e(bic)
    di
    ereturn display, level(`level')
    if e(k) > 1 {
        if "`e(transitions)'" == "none" {
            di as txt _n "No transitions between the states (latent class model for the states)"
        }
        else {
            if "`e(transitions)'" == "heterogeneous" {
                di as txt _n "Mean transition probabilities over the transition rows (rows: state at t-1;"
                di as txt "columns: state at t); by period of arrival in e(Pi_t)"
            }
            else di as txt _n "Transition probabilities (rows: state at t-1; columns: state at t)"
            matlist e(Pi), format(%9.4f) border(rows) rowtitle("from")
        }
        di as txt _n "Average posterior state shares"
        matlist e(omega), format(%9.4f) border(rows)
    }
    if e(classes) > 1 {
        di as txt _n "Latent classes: average prior and posterior probabilities"
        tempname cl
        matrix `cl' = e(pclass) \ e(wclass)
        matlist `cl', format(%9.4f) border(rows)
    }
    di as txt _n "Admissible range of lambda at the fitted means: [" ///
        as res %7.4f e(lambda_lo) as txt ", " as res %7.4f e(lambda_hi) as txt "]"
    if e(n_bound) > 0 & e(n_bound) < . {
        di as txt "note: lambda on a bound of its admissible range. It is held fixed in the"
        di as txt "      variance estimate (no standard error); the other standard errors are"
        di as txt "      conditional on it."
    }
    if e(N_droplost) > 0 & e(N_droplost) < . {
        di as txt "note: " as res e(N_droplost) as txt " dropouts occur in periods that are not estimation periods"
        di as txt "      (presample period or missing regressors) and are not used"
    }
    if e(eval_only) == 1 {
        di as txt "note: evaluation at the starting values, no iterations performed"
    }
    else if e(converged) == 0 {
        di as err "warning: convergence not achieved"
    }
end

// ===========================================================================
// predict (called by hmdbivpois_p.ado)
//   predict [type] newvar(s) [if] [in] [, statistic weights]
//   statistic: mu1 (default), mu2, pdrop, postprob, classprob, state, class,
//              viterbi, rho, pr(a,b)
//   weights:   posterior (default), filtered, onestep, prior;
//              instate(#), inclass(#) condition on a state, on a class
// ===========================================================================
program define Predict
    version 16
    syntax anything(name=vlist id="newvarlist") [if] [in]                  ///
        [, MU1 MU2 PDrop POSTprob CLASSPRob STate CLASS VITerbi RHO        ///
           PR(string) PRIor ONEstep FILTered INSTate(integer 0)            ///
           INCLass(integer 0) ]

    local k  = e(k)
    local kc = e(classes)
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'

    // statistic
    local nstat = ("`mu1'" != "") + ("`mu2'" != "") + ("`pdrop'" != "") +  ///
        ("`postprob'" != "") + ("`classprob'" != "") + ("`state'" != "") + ///
        ("`class'" != "") + ("`viterbi'" != "") + ("`rho'" != "") +        ///
        (`"`pr'"' != "")
    if `nstat' > 1 {
        di as err "only one statistic may be specified"
        exit 198
    }
    local stat "`mu1'`mu2'`pdrop'`postprob'`classprob'`state'`class'`viterbi'`rho'"
    local ya = 0
    local yb = 0
    if `"`pr'"' != "" {
        local stat "pr"
        local prs : subinstr local pr "," " ", all
        capture numlist "`prs'", min(2) max(2) integer range(>=0)
        if _rc {
            di as err "pr(a,b): a and b must be two nonnegative integers"
            exit 198
        }
        local prs `r(numlist)'
        local ya : word 1 of `prs'
        local yb : word 2 of `prs'
    }
    if `nstat' == 0 {
        local stat "mu1"
        di as txt "(option mu1 assumed; expected value of `y1')"
    }
    if "`stat'" == "pdrop" & "`e(dropout)'" == "" {
        di as err "pdrop requires a model with the dropout equation"
        exit 198
    }

    // probabilities of the latent variables used as weights
    local nw = ("`prior'" != "") + ("`onestep'" != "") + ("`filtered'" != "")
    if `nw' > 1 {
        di as err "only one of prior, onestep and filtered may be specified"
        exit 198
    }
    local ptype = 0
    if "`filtered'" != "" local ptype = 1
    if "`onestep'"  != "" local ptype = 2
    if "`prior'"    != "" local ptype = 3
    if `instate' != 0 {
        if `instate' < 1 | `instate' > `k' {
            di as err "instate() must be an integer between 1 and `k'"
            exit 198
        }
    }
    if `inclass' != 0 {
        if `inclass' < 1 | `inclass' > `kc' {
            di as err "inclass() must be an integer between 1 and `kc'"
            exit 198
        }
    }
    if inlist("`stat'", "postprob", "classprob", "state", "class", "viterbi") & ///
        (`instate' > 0 | `inclass' > 0) {
        di as err "instate() and inclass() may not be combined with `stat'"
        exit 198
    }
    if "`stat'" == "viterbi" & `nw' > 0 {
        di as err "viterbi may not be combined with prior, onestep or filtered"
        exit 198
    }

    // new variables: k for postprob, kc for classprob (stub* or names), one otherwise
    local need = 1
    if "`stat'" == "postprob"  local need = `k'
    if "`stat'" == "classprob" local need = `kc'
    gettoken first rest : vlist
    local type "`c(type)'"
    if inlist("`first'", "byte", "int", "long", "float", "double") {
        local type "`first'"
        local vlist `rest'
    }
    local nv : word count `vlist'
    if `nv' == 1 {
        if substr("`vlist'", -1, 1) == "*" {
            local stub = substr("`vlist'", 1, strlen("`vlist'") - 1)
            local vlist
            forvalues c = 1/`need' {
                local vlist `vlist' `stub'`c'
            }
            local nv = `need'
        }
    }
    if `nv' != `need' {
        if `need' == 1 di as err "one new variable name required"
        else           di as err "`stat' requires `need' new variable names, or stub*"
        exit 198
    }
    confirm new variable `vlist'

    // the results are computed in temporary variables on the estimation sample
    marksample touse, novarlist
    local rv
    forvalues i = 1/`need' {
        tempvar r`i'
        qui gen double `r`i'' = .
        local rv `rv' `r`i''
    }
    PostRun, task(predict) out(`rv') stat(`stat') ptype(`ptype')           ///
        instate(`instate') inclass(`inclass') ya(`ya') yb(`yb')

    local wlab "posterior"
    if `ptype' == 1 local wlab "filtered"
    if `ptype' == 2 local wlab "one-step-ahead"
    if `ptype' == 3 local wlab "prior"
    local cond
    if `instate' > 0 local cond " in state `instate'"
    if `inclass' > 0 local cond "`cond' in class `inclass'"
    local i = 0
    foreach v of local vlist {
        local ++i
        qui gen `type' `v' = `r`i'' if `touse'
        if "`stat'" == "postprob"  label variable `v' "Pr(state `i'), `wlab'"
        if "`stat'" == "classprob" label variable `v' "Pr(class `i'), `wlab'"
        if "`stat'" == "state"     label variable `v' "Most likely state, `wlab' probabilities"
        if "`stat'" == "class"     label variable `v' "Most likely class, `wlab' probabilities"
        if "`stat'" == "viterbi"   label variable `v' "Most likely sequence of states (Viterbi)"
        if inlist("`stat'", "mu1", "mu2") {
            local dv = cond("`stat'" == "mu1", "`y1'", "`y2'")
            label variable `v' "Expected `dv'`cond', `wlab' probabilities"
        }
        if "`stat'" == "pdrop" label variable `v' "Pr(dropout after the period)`cond', `wlab' probabilities"
        if "`stat'" == "rho"   label variable `v' "Correlation of the counts`cond', `wlab' probabilities"
        if "`stat'" == "pr"    label variable `v' "Pr(`y1'=`ya', `y2'=`yb')`cond', `wlab' probabilities"
    }
    local v1 : word 1 of `vlist'
    qui count if `v1' >= .
    if r(N) > 0 di as txt "(" r(N) " missing values generated)"
end

// ---------------------------------------------------------------------------
// Postestimation engine: rebuilds the data of the estimation from e() and the
// variables in memory, and runs a task in Mata. The estimation sample is
// e(sample); the lags of the first estimation period of each individual come
// from the observation of the period before, which must be identified by the
// individual and the period.
// ---------------------------------------------------------------------------
program define PostRun, sortpreserve
    version 16
    syntax , TASK(string) [ OUT(string) STAT(string) PTYPE(integer 0)      ///
        INSTATE(integer 0) INCLASS(integer 0) YA(integer 0) YB(integer 0)  ///
        OVER(string) MAXC(integer 5) BYCL(integer 0) PER(string)            ///
        REPS(integer 0) SEED(integer 0) HORizon(integer 50)                 ///
        LIMit(real 1000) NSIM(integer 1) ]

    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'
    local yd   "`e(dropout)'"
    local id   "`e(ivar)'"
    local time "`e(tvar)'"
    local k  = e(k)
    local kc = e(classes)
    local notr = ("`e(transitions)'" == "none")
    local het  = ("`e(transitions)'" == "heterogeneous")
    local tvt  = ("`e(tvtrans)'" != "")
    local initcov  "`e(initcov)'"
    local classcov "`e(classcov)'"
    // read by the Mata setup
    local clseq    "`e(classeq)'"
    local dropcnt  "`e(dropcounts)'"
    foreach v in `y1' `y2' `yd' `id' `time' `initcov' `classcov' {
        capture confirm numeric variable `v'
        if _rc {
            di as err "variable `v' not found; the data must be those used by hmdbivpois"
            exit 111
        }
    }

    // ---------------------------------------------------------------- sample
    tempvar es t0 pre npre touse L1 L2 Tn negT risk
    qui gen byte `es' = e(sample)
    qui count if `es'
    if r(N) == 0 | r(N) != e(N) {
        di as err "estimation sample not available; the data must be those used by hmdbivpois"
        exit 459
    }
    // the observation that precedes the first estimation period of each individual
    sort `id' `es' `time'
    qui by `id' `es': gen double `t0' = `time'[1] if `es'
    qui by `id': replace `t0' = `t0'[_N]
    qui gen byte `pre' = (`es' == 0 & `t0' < . & `time' == `t0' - 1 &      ///
        `y1' < . & `y2' < .)
    qui by `id': gen long `npre' = sum(`pre')
    qui by `id': replace `npre' = `npre'[_N]
    capture assert `npre' == 1 if `es'
    if _rc {
        di as err "cannot identify the period before the first estimation period of each"
        di as err "individual: `id' and `time' must identify the observations, and the data"
        di as err "must be those used by hmdbivpois"
        exit 459
    }
    qui gen byte `touse' = (`es' | `pre')
    sort `touse' `id' `time'
    qui by `touse' `id': gen double `L1' = `y1'[_n-1] if `touse'
    qui by `touse' `id': gen double `L2' = `y2'[_n-1] if `touse'
    // lags and their transformation, as in the estimation
    local lagt "`e(lagtransform)'"
    if "`lagt'" == "" local lagt "linear"
    local lagarg
    if regexm("`lagt'", "^(logc|topcode)\(([0-9.]+)\)$") {
        local lagarg = regexs(2)
        local lagt   = regexs(1)
    }
    tempvar G1 G2 Z1 Z2
    forvalues j = 1/2 {
        if "`lagt'" == "linear"    qui gen double `G`j'' = `L`j''
        if "`lagt'" == "log1p"     qui gen double `G`j'' = ln(1 + `L`j'')
        if "`lagt'" == "logc"      qui gen double `G`j'' = ln(max(`L`j'', `lagarg')) if `L`j'' < .
        if "`lagt'" == "indicator" qui gen double `G`j'' = (`L`j'' > 0) if `L`j'' < .
        if "`lagt'" == "topcode"   qui gen double `G`j'' = min(`L`j'', `lagarg') if `L`j'' < .
        if "`lagt'" == "logzero" {
            qui gen double `G`j'' = cond(`L`j'' > 0, ln(`L`j''), 0) if `L`j'' < .
            qui gen double `Z`j'' = (`L`j'' == 0) if `L`j'' < .
        }
    }
    local lagvars `G1' `G2'
    if "`lagt'" == "logzero" local lagvars `lagvars' `Z1' `Z2'
    local lagnames `e(lagnames)'
    // untransformed lagged counts, read by estat ame
    local rawlags `L1' `L2'
    // the same transformation of the current counts (dropcounts())
    local curnames `e(curnames)'
    local curvars
    if "`curnames'" != "" {
        tempvar C1 C2 CZ1 CZ2
        forvalues j = 1/2 {
            if "`lagt'" == "linear"    qui gen double `C`j'' = `y`j'' if `touse'
            if "`lagt'" == "log1p"     qui gen double `C`j'' = ln(1 + `y`j'') if `touse'
            if "`lagt'" == "logc"      qui gen double `C`j'' = ln(max(`y`j'', `lagarg')) if `touse'
            if "`lagt'" == "indicator" qui gen double `C`j'' = (`y`j'' > 0) if `touse'
            if "`lagt'" == "topcode"   qui gen double `C`j'' = min(`y`j'', `lagarg') if `touse'
            if "`lagt'" == "logzero" {
                qui gen double `C`j'' = cond(`y`j'' > 0, ln(`y`j''), 0) if `touse'
                qui gen double `CZ`j'' = (`y`j'' == 0) if `touse'
            }
        }
        local curvars `C1' `C2'
        if "`lagt'" == "logzero" local curvars `curvars' `CZ1' `CZ2'
    }

    // rows at risk of dropout, as in the estimation
    sort `es' `id' `time'
    if "`yd'" != "" qui by `es' `id': gen byte `risk' = `es' & (_n < _N | `yd' == 1)
    else            qui gen byte `risk' = 0

    // design variables, from the names of the terms stored in e(): common
    // terms (cvar) and terms with state-specific coefficients (svar); the
    // names (cnam, snam) are read by the Mata setup
    foreach j in 1 2 d {
        local cvar`j'
        local svar`j'
        local cnam`j'
        local snam`j'
        if "`j'" == "d" & "`yd'" == "" continue
        local cnam`j' `e(common`j')'
        local snam`j' `e(statevary`j')'
        foreach kd in c s {
            foreach t of local `kd'nam`j' {
                local pos : list posof "`t'" in lagnames
                local pc  : list posof "`t'" in curnames
                if `pos' > 0 {
                    local v : word `pos' of `lagvars'
                }
                else if `pc' > 0 {
                    local v : word `pc' of `curvars'
                }
                else if !ustrregexm("`t'", "[.#]") local v `t'
                else {
                    tempvar tv
                    TermVar `tv' `t'
                    local v `tv'
                }
                local `kd'var`j' ``kd'var`j'' `v'
            }
        }
    }

    // estimation periods of each individual, sorted as in the estimation
    sort `es' `id' `time'
    qui by `es' `id': gen long `Tn' = _N if `es'
    tempvar tr
    qui by `es' `id': gen byte `tr' = (`es' & _n > 1)
    qui gen long `negT' = -`Tn'
    sort `es' `negT' `id' `time'

    // covariates of the transitions
    local xvarsT
    foreach t in `e(xtrans_terms)' {
        if !ustrregexm("`t'", "[.#]") {
            local xvarsT `xvarsT' `t'
            continue
        }
        tempvar tv
        TermVar `tv' `t'
        local xvarsT `xvarsT' `tv'
    }

    // estat ame: counterfactual designs. For each effect, two versions (a and
    // b) of the design variables of the three equations: the variable shifted
    // by +h and -h (continuous), or set to a level and to the base level
    // (factor). The lags are shifted in Mata. The names, the equations where
    // each effect enters and its type are passed back to the caller.
    if "`task'" == "ame" {
        local eqs 1 2
        if "`yd'" != "" local eqs 1 2 d
        mata: hmdbp_amevars("`e(common1)' `e(common2)' `e(commond)' `e(statevary1)' `e(statevary2)' `e(statevaryd)'", "`lagnames' `curnames'")
        local ne = 0
        local i  = 0
        foreach v of local amev {
            local ++i
            local isf : word `i' of `amef'
            if `isf' {
                qui levelsof `v' if `es', local(alllev)
                local lvs `amelv`i''
                local base
                foreach l of local alllev {
                    local inl : list l in lvs
                    if !`inl' & "`base'" == "" local base `l'
                }
                if "`base'" == "" {
                    di as err "estat ame: base level of `v' not found in the estimation sample"
                    exit 459
                }
                foreach l of local lvs {
                    local ++ne
                    tempvar sa`ne' sb`ne'
                    qui gen double `sa`ne'' = `l'
                    qui gen double `sb`ne'' = `base'
                    local ame_name`ne' "`l'.`v'"
                    local ame_t`ne' "x"
                    local ame_div`ne' = 1
                    local ame_v`ne' "`v'"
                }
            }
            else {
                qui summarize `v' if `es'
                local h = 1e-3*cond(r(sd) > 0 & r(sd) < ., r(sd), 1)
                local ++ne
                tempvar sa`ne' sb`ne'
                qui gen double `sa`ne'' = `v' + `h'
                qui gen double `sb`ne'' = `v' - `h'
                local ame_name`ne' "`v'"
                local ame_t`ne' "x"
                local ame_div`ne' = 2*`h'
                local ame_v`ne' "`v'"
            }
        }
        forvalues e = 1/`ne' {
            local v "`ame_v`e''"
            foreach j in 1 2 d {
                local ame_f`j'_`e' = 0
            }
            // designs a and b of the common terms (ame_a#_e) and of the
            // state-specific terms (ame_as#_e)
            foreach sd in a b {
                local rep "`s`sd'`e''"
                foreach j of local eqs {
                    foreach kd in c s {
                        local lst
                        local p = 0
                        foreach t of local `kd'nam`j' {
                            local ++p
                            local w : word `p' of ``kd'var`j''
                            mata: st_local("tt", hmdbp_subst(`"`t'"', "`v'", "`rep'"))
                            if `"`tt'"' != `"`t'"' {
                                tempvar tv
                                TermVar `tv' `tt'
                                local w `tv'
                                local ame_f`j'_`e' = 1
                            }
                            local lst `lst' `w'
                        }
                        if "`kd'" == "c" local ame_`sd'`j'_`e' `lst'
                        else             local ame_`sd's`j'_`e' `lst'
                    }
                }
            }
        }
        // one more count of each lagged outcome (all the equations; the
        // dropout equation only if the lags enter it)
        forvalues l = 1/2 {
            local ++ne
            local ame_name`ne' : word `l' of `lagnames'
            local ame_t`ne' "lag`l'"
            local ame_div`ne' = 1
            foreach j in 1 2 d {
                local ame_f`j'_`ne' = 1
            }
            if "`yd'" == "" | "`e(dropcounts)'" == "current" local ame_fd_`ne' = 0
        }
        // one more current count (dropout equation only, with dropcounts())
        if "`curnames'" != "" {
            forvalues l = 1/2 {
                local ++ne
                local ame_name`ne' : word `l' of `curnames'
                local ame_t`ne' "cur`l'"
                local ame_div`ne' = 1
                local ame_f1_`ne' = 0
                local ame_f2_`ne' = 0
                local ame_fd_`ne' = 1
            }
        }
        local ame_n = `ne'
        c_local ame_n `ne'
        forvalues e = 1/`ne' {
            c_local ame_name`e' "`ame_name`e''"
            foreach j in 1 2 d {
                c_local ame_f`j'_`e' `ame_f`j'_`e''
            }
        }
    }

    // estat survival: groups (numeric codes 1, 2, ... of the variable over())
    local gvar
    if "`over'" != "" {
        tempvar gv
        qui egen long `gv' = group(`over') if `es'
        local gvar `gv'
    }

    mata: hmdbp_post("`task'", "`y1'", "`y2'", "`cvar1'", "`cvar2'",        ///
        "`yd'", "`cvard'", "`risk'", "`initcov' `lagvars'",                 ///
        "`classcov' `lagvars'", "`id'", "`es'", `k', `kc', `notr',         ///
        "`out'", "`stat'", `ptype', `instate', `inclass', `ya', `yb',       ///
        "`gvar'", `maxc', `bycl', "`xvarsT'", "`tr'", "`time'", `het', `tvt',  ///
        `reps', `seed')
end

// ---------------------------------------------------------------------------
// estat (called by hmdbivpois_estat.ado)
//   estat lrtrans [, reps(#) seed(#)]   test of no transitions between states
//   estat stability [, horizon(#) limit(#) nsim(#) seed(#)]   explosive lags
//   estat ame [, prior onestep filtered instate(#) inclass(#) level(#)]
//   estat survival [, over(varname) onestep prior filtered posterior]
//   estat fit [, max(#) classes posterior onestep prior filtered]
//   estat transition [, period(#) level(#)]
//   estat ic [, n(#)]
//   other subcommands: those of estat for every estimation command (vce,
//   summarize)
// ---------------------------------------------------------------------------
program define Estat
    version 16
    gettoken sub rest : 0, parse(" ,")
    if `"`sub'"' == "lrtrans" {
        EstatLR `rest'
    }
    else if `"`sub'"' == "ic" {
        EstatIC `rest'
    }
    else if inlist(`"`sub'"', "transition", "trans") {
        EstatTrans `rest'
    }
    else if `"`sub'"' == "ame" {
        EstatAME `rest'
    }
    else if inlist(`"`sub'"', "survival", "surv") {
        EstatSurv `rest'
    }
    else if `"`sub'"' == "fit" {
        EstatFit `rest'
    }
    else if inlist(`"`sub'"', "stability", "stab") {
        EstatStab `rest'
    }
    else estat_default `0'
end

// Likelihood-ratio test of H0: no transitions between the Markov states (all
// the off-diagonal transition probabilities are zero; the states, like the
// classes, are constant over time). H0 is on the boundary of the parameter
// space, so the LR statistic has a chi-bar-squared null distribution: a
// mixture of chi2(j), j = 0..k(k-1), whose weights are the probabilities that
// the projection of a N(0, inv(I)) vector onto the nonnegative orthant has j
// positive components, with I the information of the off-diagonal transition
// probabilities net of the other parameters, at the restricted estimates
// (Bartolucci 2006). The information is estimated by the outer product of the
// individual scores; the weights by simulation.
// The restricted model is fitted twice, from the default starting values and
// from the estimates of the model with transitions, and the larger log
// likelihood is kept.
program define EstatLR, rclass
    version 16
    syntax [, REPS(integer 10000) SEED(integer 20261117) ]
    if e(k) < 2 {
        di as err "estat lrtrans requires a model with two or more states, k() > 1"
        exit 198
    }
    if "`e(transitions)'" == "none" {
        di as err "estat lrtrans: the model in memory has no transitions (notransition);"
        di as err "fit the model with transitions and run estat lrtrans after it"
        exit 198
    }
    if "`e(transitions)'" != "homogeneous" {
        di as err "estat lrtrans requires homogeneous transitions (no xtrans() or tvtrans)"
        exit 198
    }
    if `reps' < 100 {
        di as err "reps() must be at least 100"
        exit 198
    }
    local k  = e(k)
    local nr = `k'*(`k'-1)
    tempname ll1 ll0 lla llb b1 bfrom ba W full
    scalar `ll1' = e(ll)
    local N1 = e(N)
    matrix `b1' = e(b)
    // starting values of the restricted model: e(b) without the equations of
    // the transitions (tr#_#)
    local eqn : coleq `b1'
    local keep
    local j = 0
    foreach q of local eqn {
        local ++j
        if !regexm("`q'", "^tr[0-9]+_[0-9]+$") local keep `keep' `j'
    }
    local nk : word count `keep'
    if `nk' != colsof(`b1') - `nr' {
        di as err "estat lrtrans: cannot find the `nr' transition parameters in e(b)"
        exit 498
    }
    matrix `bfrom' = J(1, `nk', 0)
    local j = 0
    foreach c of local keep {
        local ++j
        matrix `bfrom'[1, `j'] = `b1'[1, `c']
    }

    // command line of the restricted model: the same equations, sample and
    // options, without from(), emiterate(), emonly and sortby()
    local cmd `"`e(cmdline)'"'
    gettoken c0 cmd : cmd
    gettoken eq1 cmd : cmd, match(par)
    gettoken eq2 cmd : cmd, match(par)
    local eqs `"(`eq1') (`eq2')"'
    gettoken tok rest : cmd, match(par)
    if "`par'" == "(" {
        local eqs `"`eqs' (`tok')"'
        local cmd `"`rest'"'
    }
    local 0 `"`cmd'"'
    syntax [if] [in] [, *]
    local ifin `"`if' `in'"'
    local 0 `", `options'"'
    syntax [, FROM(string) NOLOg EMITERate(string) EMONLY SORTby(string) *]
    local opts0 `"`options'"'
    local opts  `"`options'"'
    if "`emiterate'" != "" local opts `"`opts' emiterate(`emiterate')"'
    if "`emonly'" != ""    local opts `"`opts' emonly"'

    _estimates hold `full', copy restore
    di as txt "fitting the model without transitions..."
    capture quietly hmdbivpois `eqs' `ifin', `opts' notransition nolog
    local rca = _rc
    if `rca' == 1 exit 1
    scalar `lla' = .
    if `rca' == 0 {
        scalar `lla' = e(ll)
        matrix `ba' = e(b)
    }
    capture quietly hmdbivpois `eqs' `ifin', `opts' notransition nolog ///
        from(`bfrom')
    local rcb = _rc
    if `rcb' == 1 exit 1
    scalar `llb' = .
    if `rcb' == 0 scalar `llb' = e(ll)
    if `rca' & `rcb' {
        di as err "estat lrtrans: the model without transitions could not be fitted"
        exit 498
    }
    // the first fit is better: its estimates are made current again
    if `rca' == 0 & (`rcb' != 0 | `lla' > `llb' + 1e-8) {
        quietly hmdbivpois `eqs' `ifin', `opts0' notransition nolog ///
            from(`ba') emiterate(0) emonly sortby(none)
    }
    scalar `ll0' = e(ll)
    if e(N) != `N1' {
        di as err "estat lrtrans: the model without transitions has a different sample"
        exit 498
    }

    // weights of the chi-bar-squared distribution, at the restricted estimates
    PostRun, task(lrtrans) out(`W') reps(`reps') seed(`seed')
    local cn
    forvalues j = 0/`nr' {
        local cn `cn' chi2_`j'
    }
    matrix colnames `W' = `cn'
    matrix rownames `W' = weight

    tempname lr p pc
    scalar `lr' = 2*(`ll1' - `ll0')
    local neg = (`lr' < 0)
    if `neg' scalar `lr' = 0
    scalar `p' = 0
    if `lr' <= 0 scalar `p' = 1
    else {
        forvalues j = 1/`nr' {
            scalar `p' = `p' + `W'[1,`j'+1]*chi2tail(`j', `lr')
        }
    }
    scalar `pc' = chi2tail(`nr', `lr')

    di _n as txt "Likelihood-ratio test of no transitions between the states"
    di as txt "H0: the " as res `nr' as txt " transition probabilities off the diagonal are zero"
    di as txt "    (the states, like the classes, are constant over time)"
    di _n as txt "  log likelihood, model with transitions    = " as res %14.4f `ll1'
    di as txt "  log likelihood, model without transitions = " as res %14.4f `ll0'
    di as txt "  LR statistic                              = " as res %14.4f `lr'
    di _n as txt "  p-value, chi-bar-squared distribution     = " as res %10.4f `p'
    di as txt "  p-value, chi2(" as res `nr' as txt ") (conservative)" _col(45) "= " as res %10.4f `pc'
    di _n as txt "Weights of the chi-bar-squared distribution, mixture of chi2(j) (" ///
        as res `reps' as txt " simulations)"
    matlist `W', format(%7.4f) border(rows) names(columns)
    if `neg' {
        di as txt "note: the log likelihood of the model with transitions is below that of the"
        di as txt "      model without transitions; the model with transitions has not reached"
        di as txt "      its maximum, and LR is set to zero"
    }
    di as txt "note: under H0 the maximum of the model with transitions is on the boundary,"
    di as txt "      where the EM converges slowly and may stop early: LR is then too small and"
    di as txt "      the test conservative. The chi2(" as res `nr' as txt ") p-value is an upper bound."

    return scalar lr     = `lr'
    return scalar p      = `p'
    return scalar p_chi2 = `pc'
    return scalar df     = `nr'
    return scalar ll_1   = `ll1'
    return scalar ll_0   = `ll0'
    return scalar reps   = `reps'
    return matrix weights = `W'
end

// estat stability [, horizon(#) limit(#) nsim(#) seed(#)]: diagnostics of the
// explosive behaviour of the lagged counts, as in hmbivpois: coefficients of
// the lags (common to the latent types); with log-type transformations the
// spectral radius of their 2 x 2 matrix, with linear lags the index of the
// univariate map y -> exp(a + g y) by type; expected counts iterated for
// horizon() periods in each type from every observation; simulated paths of
// the model over horizon() periods from the last observed period, with the
// probability of still being in the panel along them
program define EstatStab, rclass
    version 16
    syntax [, HORizon(integer 50) LIMit(real 1000) NSIM(integer 1)        ///
        SEED(integer 20261117) ]
    if `horizon' < 1 | `nsim' < 1 | `limit' <= 0 {
        di as err "horizon() and nsim() must be positive integers, limit() a positive number"
        exit 198
    }
    local k  = e(k)
    local kc = e(classes)
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'
    local yd "`e(dropout)'"
    local lagt "`e(lagtransform)'"
    if "`lagt'" == "" local lagt "linear"
    local lagt0 = regexr("`lagt'", "\(.*\)$", "")
    tempname R
    local stub = substr("`R'", 1, 25)
    PostRun, task(stability) out(`stub') horizon(`horizon') limit(`limit') ///
        nsim(`nsim') seed(`seed')
    tempname G RR U D S Q
    foreach m in G R U D S Q {
        local mm = cond("`m'" == "R", "RR", "`m'")
        matrix ``mm'' = `stub'`m'
        matrix drop `stub'`m'
    }
    // latent types: states within classes
    local tn
    forvalues u = 1/`kc' {
        forvalues v = 1/`k' {
            if `kc' > 1 local tn `tn' class`u':state`v'
            else        local tn `tn' state`v'
        }
    }
    // coefficients of the lag features by type, equations y1, y2 and (when the
    // lags enter it) the dropout
    local nlf : word count `e(lagnames)'
    local gn
    local eqg `y1' `y2'
    if colsof(`G') > 2*`nlf' local eqg `eqg' `yd'
    foreach q of local eqg {
        foreach v in `e(lagnames)' {
            local gn `gn' `q':`v'
        }
    }
    matrix rownames `G'  = `tn'
    matrix colnames `G'  = `gn'
    matrix rownames `RR' = `tn'
    matrix colnames `RR' = radius norm
    matrix rownames `U'  = `tn'
    matrix colnames `U'  = `y1' `y2'
    matrix rownames `D'  = `tn'
    matrix colnames `D'  = share
    matrix colnames `S'  = share paths horizon limit stay
    matrix rownames `Q'  = sim_`y1' sim_`y2' fit_`y1' fit_`y2'
    matrix colnames `Q'  = p50 p90 p99 max

    di _n as txt "Stability of the dynamics of the counts (lagged outcomes: " as res "`lagt'" as txt ")"
    di _n as txt "1. Coefficients of the transformed lags, by latent type (they differ across"
    di as txt "   the states only with lagvary())"
    matlist `G', format(%9.4f) border(rows) showcoleq(combined)
    if inlist("`lagt0'", "log1p", "logc", "logzero") {
        di _n as txt "   Spectral radius and norm of the 2 x 2 matrix of the coefficients in the"
        di as txt "   count equations. In logs the dynamics is close to linear for large counts:"
        di as txt "   a radius below 1 rules out explosive paths, a norm below 1 is a sufficient"
        di as txt "   condition"
        matlist `RR', format(%9.4f) border(rows)
    }
    if "`lagt0'" == "linear" {
        di _n as txt "   Index g_jj exp(a), with a the rest of the linear predictor at its 99th"
        di as txt "   percentile in the type: the map y -> exp(a + g_jj y) has a stable fixed point"
        di as txt "   only if the index is at most exp(-1) = 0.3679 (cross-lags ignored)"
        matlist `U', format(%9.4f) border(rows) showcoleq(combined)
    }
    if inlist("`lagt0'", "indicator", "topcode") {
        di _n as txt "   The transformed lags are bounded: the expected counts cannot explode"
    }
    di _n as txt "2. Deterministic paths: from each observation, with its covariates, the counts"
    di as txt "   are replaced by their expected values for " as res `horizon' as txt " periods in each type."
    di as txt "   Share of paths whose expected count exceeds " as res `limit'
    matlist `D', format(%9.4f) border(rows)
    di _n as txt "3. Simulated paths of the model: " as res `S'[1,2] as txt " paths of " as res `horizon' as txt " periods, from"
    di as txt "   the last observed period of each individual (its covariates; type drawn"
    di as txt "   from its posterior probabilities; the class stays, the state moves with the"
    di as txt "   estimated chain; the counts are drawn as if the individual stayed)"
    di as txt "   Share of paths whose expected count exceeds " as res `limit' as txt ": " as res %6.4f `S'[1,1]
    if "`yd'" != "" {
        di as txt "   Mean probability of still being in the panel after " as res `horizon' ///
            as txt " periods: " as res %6.4f `S'[1,5]
    }
    di as txt "   Expected counts at the horizon (paths below the limit) and fitted in the sample"
    matlist `Q', format(%10.3f) border(rows)

    return matrix G      = `G'
    return matrix radius = `RR'
    return matrix index  = `U'
    return matrix det    = `D'
    return matrix sim    = `S'
    return matrix q      = `Q'
end

// estat ic [, n(#)]: information criteria and relative entropy of the
// classification into the states (periods) and into the classes (individuals)
program define EstatIC, rclass
    version 16
    syntax [, N(real 0) ]
    local nn = e(N_g)
    if `n' > 0 local nn = `n'
    tempname E S
    PostRun, task(ic) out(`E')
    local df = e(npar)
    local aic = -2*e(ll) + 2*`df'
    local bic = -2*e(ll) + `df'*ln(`nn')
    matrix `S' = (`nn', ., e(ll), `df', `aic', `bic')
    matrix colnames `S' = N ll0 ll df AIC BIC
    matrix rownames `S' = .

    di _n as txt "Akaike's and Bayesian information criteria"
    di as txt "{hline 13}{c TT}{hline 63}"
    di as txt "       Model {c |}  N indiv.    N obs.      ll(model)    df          AIC          BIC"
    di as txt "{hline 13}{c +}{hline 63}"
    di as txt "           . {c |}" as res %10.0fc `nn' %10.0fc e(N) "  " %13.3f e(ll) %6.0f `df' ///
        "  " %11.3f `aic' "  " %11.3f `bic'
    di as txt "{hline 13}{c BT}{hline 63}"
    di as txt "Note: BIC uses N = number of individuals" _c
    if `n' > 0 di as txt " (n() specified)."
    else       di as txt "."
    di as txt "BIC with N = number of person-periods (as in Carrieri et al. 2026): " ///
        as res %11.3f -2*e(ll) + `df'*ln(e(N))
    di as txt "Relative entropy of the classification (1 = known with certainty):"
    di as txt "  states, by period:      " as res %6.4f `E'[1, 1]
    di as txt "  classes, by individual: " as res %6.4f `E'[1, 2]

    return scalar entropy_state = `E'[1, 1]
    return scalar entropy_class = `E'[1, 2]
    return scalar N_obs = e(N)
    return matrix S     = `S'
end

// estat transition [, level(#)]: transition probabilities, mean initial
// probabilities of the states and mean prior probabilities of the classes,
// with delta-method standard errors
program define EstatTrans, rclass
    version 16
    syntax [, PERiod(string) Level(cilevel) ]
    local k  = e(k)
    local kc = e(classes)
    if `k' < 2 & `kc' < 2 {
        di as err "estat transition requires two or more states or classes"
        exit 198
    }
    local het = ("`e(transitions)'" == "heterogeneous")
    if "`period'" != "" {
        capture confirm number `period'
        if _rc {
            di as err "period() must be a number"
            exit 198
        }
        if !`het' {
            di as txt "note: period() ignored with homogeneous transitions"
            local period
        }
    }
    tempname R P SP P0 SP0 C0 SC0
    PostRun, task(trans) out(`R') per(`period')
    local kk = `k'*`k'
    matrix `P'   = J(`k', `k', .)
    matrix `SP'  = J(`k', `k', .)
    local j0 = `kk' + 1
    local j1 = `kk' + `k'
    matrix `P0'  = `R'[1, `j0'..`j1']
    matrix `SP0' = `R'[2, `j0'..`j1']
    local j0 = `j1' + 1
    local j1 = `j1' + `kc'
    matrix `C0'  = `R'[1, `j0'..`j1']
    matrix `SC0' = `R'[2, `j0'..`j1']
    forvalues c = 1/`k' {
        forvalues d = 1/`k' {
            local j = (`c'-1)*`k' + `d'
            matrix `P'[`c', `d']  = `R'[1, `j']
            matrix `SP'[`c', `d'] = `R'[2, `j']
        }
    }
    local sn
    forvalues c = 1/`k' {
        local sn `sn' state`c'
    }
    local cn
    forvalues u = 1/`kc' {
        local cn `cn' class`u'
    }
    matrix rownames `P'  = `sn'
    matrix colnames `P'  = `sn'
    matrix rownames `SP' = `sn'
    matrix colnames `SP' = `sn'
    matrix colnames `P0'  = `sn'
    matrix colnames `SP0' = `sn'
    matrix rownames `P0'  = prob
    matrix rownames `SP0' = se
    matrix colnames `C0'  = `cn'
    matrix colnames `SC0' = `cn'
    matrix rownames `C0'  = prob
    matrix rownames `SC0' = se

    local z = invnormal(1 - (1 - `level'/100)/2)
    if `k' > 1 & "`e(transitions)'" != "none" {
        local ttl "Transition probabilities"
        if `het' & "`period'" == "" local ttl "Mean transition probabilities over the transition rows"
        if `het' & "`period'" != "" local ttl "Mean transition probabilities into period `period'"
        di _n as txt "`ttl' (delta-method standard errors;"
        di as txt "confidence intervals computed on the logit scale)"
        di as txt "{hline 13}{c TT}{hline 50}"
        di as txt "  from -> to {c |}      Prob.   Std. err.     [`level'% conf. interval]"
        di as txt "{hline 13}{c +}{hline 50}"
        forvalues c = 1/`k' {
            forvalues d = 1/`k' {
                local p  = `P'[`c', `d']
                local se = `SP'[`c', `d']
                local lo = .
                local hi = .
                if `p' > 1e-8 & `p' < 1 - 1e-8 & `se' > 0 {
                    local lo = invlogit(logit(`p') - `z'*`se'/(`p'*(1-`p')))
                    local hi = invlogit(logit(`p') + `z'*`se'/(`p'*(1-`p')))
                }
                di as txt %12s "`c' -> `d'" " {c |}" as res %11.4f `p' %12.4f `se' ///
                    %12.4f `lo' %11.4f `hi'
            }
            if `c' < `k' di as txt "{hline 13}{c +}{hline 50}"
        }
        di as txt "{hline 13}{c BT}{hline 50}"

        di _n as txt "Expected duration of a stay in each state, 1/(1 - P(c -> c)), in periods"
        di as txt "{hline 13}{c TT}{hline 24}"
        di as txt "       state {c |}   Duration   Std. err."
        di as txt "{hline 13}{c +}{hline 24}"
        tempname DU
        matrix `DU' = J(2, `k', .)
        forvalues c = 1/`k' {
            local p  = `P'[`c', `c']
            local du = 1/(1 - `p')
            local sd = `SP'[`c', `c']/(1 - `p')^2
            matrix `DU'[1, `c'] = `du'
            matrix `DU'[2, `c'] = `sd'
            di as txt %12s "`c'" " {c |}" as res %11.3f `du' %12.3f `sd'
        }
        di as txt "{hline 13}{c BT}{hline 24}"
        matrix colnames `DU' = `sn'
        matrix rownames `DU' = duration se
        return matrix duration = `DU'
    }
    if `k' > 1 {
        di _n as txt "Initial probabilities of the states, mean over the individuals"
        di as txt "{hline 13}{c TT}{hline 24}"
        di as txt "       state {c |}      Prob.   Std. err."
        di as txt "{hline 13}{c +}{hline 24}"
        forvalues c = 1/`k' {
            di as txt %12s "`c'" " {c |}" as res %11.4f `P0'[1, `c'] %12.4f `SP0'[1, `c']
        }
        di as txt "{hline 13}{c BT}{hline 24}"
    }
    if `kc' > 1 {
        di _n as txt "Prior probabilities of the classes, mean over the individuals"
        di as txt "{hline 13}{c TT}{hline 24}"
        di as txt "       class {c |}      Prob.   Std. err."
        di as txt "{hline 13}{c +}{hline 24}"
        forvalues u = 1/`kc' {
            di as txt %12s "`u'" " {c |}" as res %11.4f `C0'[1, `u'] %12.4f `SC0'[1, `u']
        }
        di as txt "{hline 13}{c BT}{hline 24}"
    }

    return matrix P        = `P'
    return matrix se_P     = `SP'
    return matrix pi0      = `P0'
    return matrix se_pi0   = `SP0'
    return matrix pclass   = `C0'
    return matrix se_pclass = `SC0'
end

// estat fit [, max(#) classes posterior onestep prior filtered]: observed
// distribution of each count (0, 1, ..., max-1, max or more) against the
// distributions predicted for each latent state (or class, with classes) and
// their mixture, the predicted marginal distribution (Table 5-6 and Fig. 7 of
// Carrieri et al. 2026). The distribution of a type is the mean over the
// periods of the Poisson probabilities in that type, weighted by the
// probabilities of the latent variables (posterior by default); the mixture
// weights the types by their mean probabilities, so it is the mean over the
// periods of the predicted probabilities. Pearson chi2 over the person-periods,
// dissimilarity index (half the sum of the absolute differences between the
// observed and predicted shares) and largest difference of a cell.
program define EstatFit, rclass
    version 16
    syntax [, MAX(integer 5) CLASSes POSTerior ONEstep PRIor FILTered ]
    if `max' < 1 | `max' > 100 {
        di as err "max() must be between 1 and 100"
        exit 198
    }
    local nw = ("`posterior'" != "") + ("`onestep'" != "") + ("`prior'" != "") + ("`filtered'" != "")
    if `nw' > 1 {
        di as err "only one of posterior, onestep, prior and filtered may be specified"
        exit 198
    }
    local ptype = 0
    if "`filtered'" != "" local ptype = 1
    if "`onestep'"  != "" local ptype = 2
    if "`prior'"    != "" local ptype = 3
    local wlab "posterior"
    if `ptype' == 1 local wlab "filtered"
    if `ptype' == 2 local wlab "one-step-ahead"
    if `ptype' == 3 local wlab "prior"
    local bycl = ("`classes'" != "")
    local tname = cond(`bycl', "class", "state")
    local nt = cond(`bycl', e(classes), e(k))
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'

    tempname R F1 F2 G
    local stub = substr("`R'", 1, 25)
    PostRun, task(fit) out(`stub') ptype(`ptype') maxc(`max') bycl(`bycl')
    local m1 = `max' + 1
    local m2 = 2*`m1'
    matrix `F1' = `stub'[1..`m1', 1...]
    matrix `F2' = `stub'[`=`m1'+1'..`m2', 1...]
    matrix `G'  = `stub'G
    matrix drop `stub' `stub'G

    local cn obs
    forvalues c = 1/`nt' {
        local cn `cn' `tname'`c'
    }
    local cn `cn' mixture
    local rn
    forvalues c = 0/`=`max'-1' {
        local rn `rn' `c'
    }
    local rn `rn' `max'_more
    foreach j in 1 2 {
        matrix colnames `F`j'' = `cn'
        matrix rownames `F`j'' = `rn'
    }
    matrix colnames `G' = chi2 df dissimilarity maxdev N
    matrix rownames `G' = `y1' `y2'

    di _n as txt "Observed and predicted distributions of the counts, by latent `tname'"
    di as txt "Weights: " as res "`wlab'" as txt "; mixture = types weighted by their mean probabilities"
    foreach j in 1 2 {
        di _n as txt "`y`j''"
        matlist `F`j'', format(%9.4f) border(rows) rowtitle("count")
    }
    di _n as txt "Goodness of fit of the mixture (Pearson chi2 over the " ///
        as res e(N) as txt " person-periods)"
    di as txt "{hline 13}{c TT}{hline 50}"
    di as txt "    response {c |}      chi2   df   dissimilarity   max. cell dev."
    di as txt "{hline 13}{c +}{hline 50}"
    foreach j in 1 2 {
        di as txt %12s abbrev("`y`j''", 12) " {c |}" as res %10.1f `G'[`j', 1] %5.0f `G'[`j', 2] ///
            %15.2f 100*`G'[`j', 3] "%" %15.2f 100*`G'[`j', 4] " pp"
    }
    di as txt "{hline 13}{c BT}{hline 50}"
    di as txt "With many person-periods the chi2 rejects exact fit for negligible differences:"
    di as txt "read it, with the dissimilarity index, as a descriptive measure."

    return matrix F1  = `F1'
    return matrix F2  = `F2'
    return matrix gof = `G'
end

// estat ame [, prior onestep filtered instate(#) inclass(#) level(#)]:
// average marginal effects on E(y1), E(y2) and on the probability of dropout,
// with the probabilities of the latent variables held fixed
program define EstatAME, rclass
    version 16
    syntax [, PRIor ONEstep FILTered INSTate(integer 0) INCLass(integer 0) ///
        Level(cilevel) ]
    local k  = e(k)
    local kc = e(classes)
    local nw = ("`prior'" != "") + ("`onestep'" != "") + ("`filtered'" != "")
    if `nw' > 1 {
        di as err "only one of prior, onestep and filtered may be specified"
        exit 198
    }
    local ptype = 0
    if "`filtered'" != "" local ptype = 1
    if "`onestep'"  != "" local ptype = 2
    if "`prior'"    != "" local ptype = 3
    if `instate' != 0 {
        if `instate' < 1 | `instate' > `k' {
            di as err "instate() must be an integer between 1 and `k'"
            exit 198
        }
    }
    if `inclass' != 0 {
        if `inclass' < 1 | `inclass' > `kc' {
            di as err "inclass() must be an integer between 1 and `kc'"
            exit 198
        }
    }
    local y1 : word 1 of `e(depvar)'
    local y2 : word 2 of `e(depvar)'
    local yd "`e(dropout)'"
    local eqs 1 2
    if "`yd'" != "" local eqs 1 2 d
    local lab1 "E(`y1')"
    local lab2 "E(`y2')"
    local labd "Pr(`yd')"

    tempname R b V
    local stub = substr("`R'", 1, 25)
    PostRun, task(ame) out(`stub') ptype(`ptype') instate(`instate') inclass(`inclass')
    matrix `b' = `stub'
    matrix `V' = `stub'V
    matrix drop `stub' `stub'V
    // only the effects of the variables that enter each equation
    local names
    local keep
    local q = 0
    foreach j of local eqs {
        local ++q
        forvalues e = 1/`ame_n' {
            if `ame_f`j'_`e'' {
                if "`j'" == "d" local en `yd'
                else            local en `y`j''
                local names `names' `en':`ame_name`e''
                local keep `keep' `=(`q'-1)*`ame_n' + `e''
            }
        }
    }
    mata: hmdbp_msel("`b'", "`V'", "`keep'")
    matrix colnames `b' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'

    local wlab "posterior (given all the observed outcomes of the individual)"
    if `ptype' == 1 local wlab "filtered (given the outcomes up to the period)"
    if `ptype' == 2 local wlab "one-step-ahead (given the outcomes up to the period before)"
    if `ptype' == 3 local wlab "prior (given no outcome)"
    if `instate' > 0 local wlab "`wlab'; in state `instate'"
    if `inclass' > 0 local wlab "`wlab'; in class `inclass'"
    di _n as txt "Average marginal effects on the expected counts and on the probability of"
    di as txt "dropout, with the probabilities of the latent variables held fixed."
    di as txt "Weights: " as res "`wlab'"
    di as txt "dy/dx for continuous variables; discrete change from the base level for factor"
    local ln1 : word 1 of `e(lagnames)'
    local ln2 : word 2 of `e(lagnames)'
    di as txt "levels; effect of one more count of the lagged outcome for " as res "`ln1' `ln2'"
    if "`e(curnames)'" != "" {
        local cn1 : word 1 of `e(curnames)'
        local cn2 : word 2 of `e(curnames)'
        di as txt "and of one more count of the current period, on the dropout, for " as res "`cn1' `cn2'"
    }
    local z = invnormal(1 - (1 - `level'/100)/2)
    di as txt "{hline 13}{c TT}{hline 64}"
    di as txt "             {c |}            Delta-method"
    di as txt "             {c |}     Effect   std. err.      z    P>|z|     [`level'% conf. interval]"
    di as txt "{hline 13}{c +}{hline 64}"
    local c = 0
    local q = 0
    local neq : word count `eqs'
    foreach j of local eqs {
        local ++q
        di as res %-12s abbrev("`lab`j''", 12) as txt " {c |}"
        forvalues e = 1/`ame_n' {
            if `ame_f`j'_`e'' {
                local ++c
                local ef = `b'[1, `c']
                local se = sqrt(`V'[`c', `c'])
                local zz = `ef'/`se'
                di as txt %12s abbrev("`ame_name`e''", 12) " {c |}" as res "  " %9.0g `ef' "  " %9.0g `se' ///
                    %9.2f `zz' %8.3f 2*normal(-abs(`zz')) "    " %9.0g `ef' - `z'*`se' "   " %9.0g `ef' + `z'*`se'
            }
        }
        if `q' < `neq' di as txt "{hline 13}{c +}{hline 64}"
    }
    di as txt "{hline 13}{c BT}{hline 64}"

    return matrix b = `b'
    return matrix V = `V'
end

// estat survival [, over(varname) onestep prior filtered posterior]:
// Kaplan-Meier survival by period of follow-up (estimation periods of each
// individual: 1, 2, ...) against the survival predicted by the model, the
// product over the periods of 1 - the mean predicted probability of dropout
// among the individuals at risk. Weights of the latent variables: one-step-
// ahead by default (given the outcomes, and the survival, up to the period
// before)
program define EstatSurv, rclass
    version 16
    syntax [, OVER(varname) ONEstep PRIor FILTered POSTerior ]
    if "`e(dropout)'" == "" {
        di as err "estat survival requires a model with the dropout equation"
        exit 198
    }
    local nw = ("`prior'" != "") + ("`onestep'" != "") + ("`filtered'" != "") + ("`posterior'" != "")
    if `nw' > 1 {
        di as err "only one of onestep, prior, filtered and posterior may be specified"
        exit 198
    }
    local ptype = 2
    if "`posterior'" != "" local ptype = 0
    if "`filtered'"  != "" local ptype = 1
    if "`prior'"     != "" local ptype = 3
    local wlab "one-step-ahead"
    if `ptype' == 0 local wlab "posterior"
    if `ptype' == 1 local wlab "filtered"
    if `ptype' == 3 local wlab "prior"

    tempname S
    PostRun, task(surv) out(`S') ptype(`ptype') over(`over')
    // S: one row per group and period: group, period, at risk, dropouts,
    // Kaplan-Meier, model
    local glab
    if "`over'" != "" {
        qui levelsof `over' if e(sample), local(glev)
        local glab `glev'
    }
    di _n as txt "Survival by period of follow-up: Kaplan-Meier and model (" as res "`wlab'" as txt " weights)"
    di as txt "{hline 9}{c TT}{hline 58}"
    di as txt "  period {c |}    at risk   dropouts   Kaplan-Meier      Model   Difference"
    local gold = .
    forvalues r = 1/`=rowsof(`S')' {
        local g = `S'[`r', 1]
        if `g' != `gold' {
            di as txt "{hline 9}{c +}{hline 58}"
            if "`over'" != "" {
                local gl : word `g' of `glab'
                di as txt "`over' = " as res "`gl'"
            }
            local gold = `g'
        }
        di as txt %8.0f `S'[`r', 2] " {c |}" as res %11.0fc `S'[`r', 3] %11.0fc `S'[`r', 4] ///
            %15.4f `S'[`r', 5] %11.4f `S'[`r', 6] %13.4f `S'[`r', 6] - `S'[`r', 5]
    }
    di as txt "{hline 9}{c BT}{hline 58}"
    di as txt "Model: product over the periods of 1 - the mean predicted probability of"
    di as txt "dropout among the individuals at risk."
    matrix colnames `S' = group period atrisk dropouts km model
    return matrix S = `S'
end

// ===========================================================================
// Mata
// ===========================================================================
mata:

// Data and dimensions. Expanded latent state s = (u-1)*k2 + v, K = k1*k2,
// with u the class (k1 = classes()) and v the Markov state (k2 = k()).
// C = (Ev, Eu) maps the reduced intercepts (alpha_1..alpha_k2, delta_2..delta_k1)
// of one equation into the K intercepts of the expanded states.
struct hmdbp_dat {
    real scalar     n, N, k1, k2, K, p1, p2, pd, nc, q, qu, Tmax
    real scalar     m, md, nD, ntau, nU, npar, hasD, notr, chk, het, qw
    real scalar     emsoft, cm, cdd, nc1, ncd, nl, s1, s2, sd, nt1, nt2, ntd
    real rowvector  cc1, sc1, cc2, sc2, ccd, scd
    real colvector  y1, y2, ey1, ey2, lnfy, D, risk, Ts, ns, trow, tt
    real matrix     X1, X2, Xd, Z, ZU, info, Ev, Euf, C, Jm, Jd, W, C1, Cd
    real matrix     XS1, XS2, XSd
    pointer(real colvector) colvector R, RT
}

// measurement model (N x K): means, Famoye terms, log densities; dropout
struct hmdbp_mm {
    real matrix     lm0, lm, Mu1, Mu2, D1, D2, U1, U2, G, ld, Pd
}

// ------------------------------------------------------------ parameters
// Layout of the parameter vector b:
//   equation 1: slopes (p1: regressors, then the two lags), state intercepts
//               (k2), class shifts (k1-1)
//   equation 2: the same (p2, k2, k1-1)
//   lambda                                     -> end of the measurement block (m)
//   dropout: slopes (pd), state intercepts (k2), class shifts (k1-1)   (md)
//   initial probabilities of the states: q coefficients for each state 2..k2
//   transitions: one coefficient for each pair (from, to), from != to
//   prior probabilities of the classes: qu coefficients for each class 2..k1
real matrix hmdbp_P(real rowvector tau, real scalar k)
{
    real matrix    P
    real rowvector e
    real scalar    c, cb, j

    P = J(k, k, 0)
    j = 0
    for (c=1; c<=k; c++) {
        e = J(1, k, 0)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                e[cb] = tau[j]
            }
        }
        e = exp(e :- max(e))
        P[c,.] = e :/ sum(e)
    }
    return(P)
}

real rowvector hmdbp_tau(real matrix P)
{
    real rowvector tau
    real scalar    k, c, cb, j

    k = rows(P)
    tau = J(1, k*(k-1), .)
    j = 0
    for (c=1; c<=k; c++) {
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                tau[j] = ln(P[c,cb]) - ln(P[c,c])
            }
        }
    }
    return(tau)
}

// pieces of the measurement block bm (1 x m)
void hmdbp_mparts(struct hmdbp_dat scalar M, real rowvector bm,
    real rowvector b1, real rowvector t1, real rowvector b2,
    real rowvector t2, real scalar lam)
{
    real scalar o

    o = 0
    if (M.p1 > 0) {
        b1 = bm[|o+1 \ o+M.p1|]
        o  = o + M.p1
    }
    else b1 = J(1, 0, .)
    t1 = bm[|o+1 \ o+M.nt1|]
    o  = o + M.nt1
    if (M.p2 > 0) {
        b2 = bm[|o+1 \ o+M.p2|]
        o  = o + M.p2
    }
    else b2 = J(1, 0, .)
    t2  = bm[|o+1 \ o+M.nt2|]
    o   = o + M.nt2
    lam = bm[o+1]
}

// pieces of the dropout block bdv (1 x md)
void hmdbp_dparts(struct hmdbp_dat scalar M, real rowvector bdv,
    real rowvector bd, real rowvector td)
{
    if (M.pd > 0) {
        bd = bdv[|1 \ M.pd|]
        td = bdv[|M.pd+1 \ M.pd+M.ntd|]
    }
    else {
        bd = J(1, 0, .)
        td = bdv[|1 \ M.ntd|]
    }
}

void hmdbp_unpack(struct hmdbp_dat scalar M, real rowvector b,
    real rowvector b1, real rowvector t1, real rowvector b2,
    real rowvector t2, real scalar lam, real rowvector bd,
    real rowvector td, real matrix Dv, real rowvector tau, real matrix Gu)
{
    real scalar o

    hmdbp_mparts(M, b[|1 \ M.m|], b1, t1, b2, t2, lam)
    o = M.m
    if (M.hasD) {
        hmdbp_dparts(M, b[|o+1 \ o+M.md|], bd, td)
        o = o + M.md
    }
    else {
        bd = J(1, 0, .)
        td = J(1, 0, .)
    }
    if (M.k2 > 1) {
        Dv = colshape(b[|o+1 \ o+M.nD|], M.q)'          // q x (k2-1)
        o  = o + M.nD
    }
    else Dv = J(M.q, 0, .)
    if (M.ntau > 0) {
        tau = b[|o+1 \ o+M.ntau|]
        o   = o + M.ntau
    }
    else tau = J(1, 0, .)
    if (M.k1 > 1) Gu = colshape(b[|o+1 \ o+M.nU|], M.qu)'   // qu x (k1-1)
    else          Gu = J(M.qu, 0, .)
}

real rowvector hmdbp_pack(struct hmdbp_dat scalar M, real rowvector b1,
    real rowvector t1, real rowvector b2, real rowvector t2, real scalar lam,
    real rowvector bd, real rowvector td, real matrix Dv, real rowvector tau,
    real matrix Gu)
{
    real rowvector b

    b = (b1, t1, b2, t2, lam)
    if (M.hasD)     b = (b, bd, td)
    if (M.k2 > 1)   b = (b, vec(Dv)')
    if (M.ntau > 0) b = (b, tau)
    if (M.k1 > 1)   b = (b, vec(Gu)')
    return(b)
}

// multinomial logit probabilities, base = first category
real matrix hmdbp_mlogit(real matrix Z, real matrix B, real scalar k)
{
    real matrix E

    if (k == 1) return(J(rows(Z), 1, 1))
    E = (J(rows(Z), 1, 0), Z*B)
    E = exp(E :- rowmax(E))
    return(E :/ rowsum(E))
}

// ln(invlogit(x)), without overflow
real matrix hmdbp_lnsig(real matrix X)
{
    return((X :< 0):*X - ln(1 :+ exp(-abs(X))))
}

// ------------------------------------------------------- measurement model
// linear predictor of one equation (N x K, one column per expanded state):
// regressors with common coefficients b, intercepts of the states and shifts
// of the classes (the first nc elements of t, mapped by C), and regressors XS
// with state-specific coefficients (the rest of t: s coefficients for state
// 1, then for state 2, ...)
real matrix hmdbp_eta(struct hmdbp_dat scalar M, real matrix X,
    real matrix XS, real rowvector b, real rowvector t, real matrix C,
    real scalar nc, real scalar s)
{
    real matrix E

    E = J(rows(X), 1, (C*t[|1 \ nc|]')')
    if (cols(b) > 0) E = E :+ (X*b')
    if (s > 0) E = E + (XS*colshape(t[|nc+1 \ cols(t)|], s)')*M.Ev'
    return(E)
}

// means and everything else that does not depend on lambda
struct hmdbp_mm scalar hmdbp_means(struct hmdbp_dat scalar M,
    real rowvector b1, real rowvector t1, real rowvector b2, real rowvector t2)
{
    struct hmdbp_mm scalar mm
    real matrix eta1, eta2
    real scalar q

    q = 1 - exp(-1)
    eta1 = hmdbp_eta(M, M.X1, M.XS1, b1, t1, M.C1, M.nc1, M.s1)
    eta2 = hmdbp_eta(M, M.X2, M.XS2, b2, t2, M.C1, M.nc1, M.s2)
    mm.Mu1 = exp(eta1)
    mm.Mu2 = exp(eta2)
    mm.D1  = exp(-q*mm.Mu1)
    mm.D2  = exp(-q*mm.Mu2)
    mm.U1  = M.ey1 :- mm.D1
    mm.U2  = M.ey2 :- mm.D2
    mm.lm0 = M.y1:*eta1 - mm.Mu1 + M.y2:*eta2 - mm.Mu2
    return(mm)
}

void hmdbp_dep(struct hmdbp_dat scalar M, struct hmdbp_mm scalar mm,
    real scalar lam)
{
    mm.G  = 1 :+ lam*(mm.U1:*mm.U2)
    mm.lm = mm.lm0 + ln(mm.G) :- M.lnfy
}

// dropout: log contributions (N x K, zero on the rows not at risk) and
// probabilities
void hmdbp_dropf(struct hmdbp_dat scalar M, real rowvector bd,
    real rowvector td, struct hmdbp_mm scalar mm)
{
    real matrix E

    if (!M.hasD) {
        mm.ld = J(M.N, M.K, 0)
        mm.Pd = J(0, 0, .)
        return
    }
    E = hmdbp_eta(M, M.Xd, M.XSd, bd, td, M.Cd, M.ncd, M.sd)
    mm.Pd = invlogit(E)
    mm.ld = M.risk :* (M.D :* hmdbp_lnsig(E) + (1 :- M.D) :* hmdbp_lnsig(-E))
}

// admissible interval of lambda (pmf nonnegative for every pair of counts at
// the current means, in every expanded state): (lower \ upper)
real colvector hmdbp_lambnd(struct hmdbp_mm scalar mm)
{
    real matrix A, B, X
    real scalar lo, hi

    A  = mm.D1:*mm.D2
    B  = (1:-mm.D1):*(1:-mm.D2)
    X  = A:*(A:>=B) + B:*(A:<B)
    lo = max(-1:/X)
    A  = (1:-mm.D1):*mm.D2
    B  = mm.D1:*(1:-mm.D2)
    X  = A:*(A:>=B) + B:*(A:<B)
    hi = min(1:/X)
    return(lo \ hi)
}

// ------------------------------------------------ forward / backward
// individuals are sorted by decreasing T_i: those active at period t are 1..ns[t]
real colvector hmdbp_fwd(struct hmdbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix P, real matrix A, real matrix F)
{
    real colvector m, ll, s, r
    real matrix    a
    real scalar    t, nt

    m = rowmax(lf)
    F = exp(lf :- m)
    A = J(M.N, M.K, 0)
    r = *M.R[1]
    a = Pi0 :* F[r,.]
    s = rowsum(a)
    a = a :/ s
    A[r,.] = a
    ll = ln(s) + m[r]
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        a  = (a[|1,1 \ nt,M.K|]*P) :* F[r,.]
        s  = rowsum(a)
        a  = a :/ s
        A[r,.] = a
        ll[|1 \ nt|] = ll[|1 \ nt|] + ln(s) + m[r]
    }
    return(ll)
}

// posteriors W (N x K), expected transition counts Zs (K x K) and, if wantxi,
// their individual contributions XiI (n x K*K, column (s-1)*K+s' = s -> s')
void hmdbp_bwd(struct hmdbp_dat scalar M, real matrix F, real matrix A,
    real matrix P, real matrix W, real matrix Zs, real scalar wantxi,
    real matrix XiI)
{
    real matrix    B, G, Ap, bb
    real colvector r, rp, nrm
    real scalar    t, nt, K, c

    K  = M.K
    W  = A
    Zs = J(K, K, 0)
    if (wantxi) XiI = J(M.n, K*K, 0)
    if (M.Tmax >= 2) {
        B = J(M.ns[M.Tmax], K, 1)
        for (t=M.Tmax; t>=2; t--) {
            nt  = M.ns[t]
            r   = *M.R[t]
            rp  = (*M.R[t-1])[|1 \ nt|]
            G   = F[r,.] :* B
            Ap  = A[rp,.]
            nrm = rowsum((Ap*P) :* G)
            Ap  = Ap :/ nrm
            Zs  = Zs + P :* cross(Ap, G)
            if (wantxi) {
                for (c=1; c<=K; c++) {
                    XiI[|1,(c-1)*K+1 \ nt,c*K|] = XiI[|1,(c-1)*K+1 \ nt,c*K|] + (Ap[.,c] :* G) :* P[c,.]
                }
            }
            bb = G*P'
            bb = bb :/ rowsum(bb)
            B  = J(M.ns[t-1], K, 1)
            B[|1,1 \ nt,K|] = bb
            rp = *M.R[t-1]
            W[rp,.] = A[rp,.] :* B
        }
    }
    W = W :/ rowsum(W)
}

// expected transition counts of the Markov states: sum over the classes of the
// diagonal blocks of Zs
real matrix hmdbp_zv(struct hmdbp_dat scalar M, real matrix Zs)
{
    real matrix    Z
    real rowvector r
    real scalar    u

    Z = J(M.k2, M.k2, 0)
    for (u=1; u<=M.k1; u++) {
        r = ((u-1)*M.k2+1)..(u*M.k2)
        Z = Z + Zs[r, r]
    }
    return(Z)
}

// ------------------------------------------------ heterogeneous transitions
// Transition probabilities of the Markov states by transition row (xtrans,
// tvtrans): rows(M.W) x k2*k2, columns (c-1)*k2+d = P(state d | state c).
// tau holds, for each pair (c, d), d != c, qw coefficients of the design M.W;
// the base of each row is "staying".
real matrix hmdbp_PT(struct hmdbp_dat scalar M, real rowvector tau)
{
    real matrix PT, TT, G, E
    real scalar c, cb, j, k

    k  = M.k2
    TT = colshape(tau, M.qw)
    PT = J(rows(M.W), k*k, .)
    j  = 0
    for (c=1; c<=k; c++) {
        G = J(M.qw, k, 0)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                G[.,cb] = TT[j,.]'
            }
        }
        E = M.W*G
        E = exp(E :- rowmax(E))
        PT[|1,(c-1)*k+1 \ .,c*k|] = E :/ rowsum(E)
    }
    return(PT)
}

// one step of the expanded chain with the transition rows rt: the class does
// not change, the state moves with PT
real matrix hmdbp_proph(struct hmdbp_dat scalar M, real matrix ap,
    real matrix PT, real colvector rt)
{
    real matrix    a
    real rowvector dest, cv
    real scalar    u, v, s, k2

    k2 = M.k2
    a  = J(rows(ap), M.K, 0)
    for (u=1; u<=M.k1; u++) {
        dest = ((u-1)*k2+1)..(u*k2)
        for (v=1; v<=k2; v++) {
            s  = (u-1)*k2 + v
            cv = ((v-1)*k2+1)..(v*k2)
            a[., dest] = a[., dest] + ap[., s] :* PT[rt, cv]
        }
    }
    return(a)
}

// forward recursion with transition probabilities by transition row
real colvector hmdbp_fwdh(struct hmdbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix PT, real matrix A, real matrix F)
{
    real colvector m, ll, s, r, rt
    real matrix    a
    real scalar    t, nt

    m = rowmax(lf)
    F = exp(lf :- m)
    A = J(M.N, M.K, 0)
    r = *M.R[1]
    a = Pi0 :* F[r,.]
    s = rowsum(a)
    a = a :/ s
    A[r,.] = a
    ll = ln(s) + m[r]
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        rt = *M.RT[t]
        a  = hmdbp_proph(M, a[|1,1 \ nt,M.K|], PT, rt) :* F[r,.]
        s  = rowsum(a)
        a  = a :/ s
        A[r,.] = a
        ll[|1 \ nt|] = ll[|1 \ nt|] + ln(s) + m[r]
    }
    return(ll)
}

// backward recursion with transition probabilities by transition row. W:
// posteriors of the expanded states; XI (rows(M.W) x k2*k2): joint posteriors
// of the states of the two periods of each transition row, summed over the
// classes
void hmdbp_bwdh(struct hmdbp_dat scalar M, real matrix F, real matrix A,
    real matrix PT, real matrix W, real matrix XI)
{
    real matrix    B, G, Ap, bb
    real colvector r, rt, rp, nrm
    real rowvector dest, cv
    real scalar    t, nt, K, k2, u, v, s

    K  = M.K
    k2 = M.k2
    W  = A
    XI = J(rows(M.W), k2*k2, 0)
    if (M.Tmax >= 2) {
        B = J(M.ns[M.Tmax], K, 1)
        for (t=M.Tmax; t>=2; t--) {
            nt  = M.ns[t]
            r   = *M.R[t]
            rt  = *M.RT[t]
            rp  = (*M.R[t-1])[|1 \ nt|]
            G   = F[r,.] :* B
            Ap  = A[rp,.]
            // bb[., s] = sum over the arrival states of P(s -> s') G[., s']
            bb  = J(nt, K, 0)
            for (u=1; u<=M.k1; u++) {
                dest = ((u-1)*k2+1)..(u*k2)
                for (v=1; v<=k2; v++) {
                    s  = (u-1)*k2 + v
                    cv = ((v-1)*k2+1)..(v*k2)
                    bb[., s] = rowsum(PT[rt, cv] :* G[., dest])
                }
            }
            nrm = rowsum(Ap :* bb)
            for (u=1; u<=M.k1; u++) {
                dest = ((u-1)*k2+1)..(u*k2)
                for (v=1; v<=k2; v++) {
                    s  = (u-1)*k2 + v
                    cv = ((v-1)*k2+1)..(v*k2)
                    XI[rt, cv] = XI[rt, cv] + ((Ap[., s] :/ nrm) :* PT[rt, cv]) :* G[., dest]
                }
            }
            bb = bb :/ rowsum(bb)
            B  = J(M.ns[t-1], K, 1)
            B[|1,1 \ nt,K|] = bb
            rp = *M.R[t-1]
            W[rp,.] = A[rp,.] :* B
        }
    }
    W = W :/ rowsum(W)
}

real matrix hmdbp_lnsoft(real matrix E)
{
    real colvector m

    m = rowmax(E)
    return(E :- m :- ln(rowsum(exp(E :- m))))
}

// weighted (fractional) multinomial logit by Newton-Raphson with step-halving.
// Y: n x k nonnegative responses (row sums are the weights); B0: q x k
// coefficients with a zero column for the base category. Returns B (q x k).
real matrix hmdbp_fmlogit(real matrix X, real matrix Y, real matrix B0,
    real scalar base, real scalar nsteps)
{
    real matrix    B, Bn, lnP, lnPn, Pr, H, blk
    real colvector w, g, step, wg
    real rowvector free
    real scalar    k, q, it, a, e, j, st, Q0, Q1

    k    = cols(Y)
    q    = cols(X)
    w    = rowsum(Y)
    free = select(1..k, (1..k) :!= base)
    B    = B0
    lnP  = hmdbp_lnsoft(X*B)
    Q0   = quadsum(Y :* lnP)
    for (it=1; it<=nsteps; it++) {
        Pr = exp(lnP)
        g  = vec(cross(X, Y[., free] - Pr[., free] :* w))
        H  = J(q*(k-1), q*(k-1), 0)
        for (a=1; a<=k-1; a++) {
            for (e=a; e<=k-1; e++) {
                wg  = w :* Pr[., free[a]] :* ((a==e) :- Pr[., free[e]])
                blk = -cross(X, wg, X)
                H[|(a-1)*q+1,(e-1)*q+1 \ a*q,e*q|] = blk
                H[|(e-1)*q+1,(a-1)*q+1 \ e*q,a*q|] = blk'
            }
        }
        step = cholsolve(-H, g)
        if (hasmissing(step)) step = invsym(-H)*g
        if (hasmissing(step)) break
        st = 1
        for (j=1; j<=30; j++) {
            Bn = B
            Bn[., free] = B[., free] + st*colshape(step', q)'
            lnPn = hmdbp_lnsoft(X*Bn)
            Q1   = quadsum(Y :* lnPn)
            if (!hasmissing(lnPn) & Q1 >= Q0 - 1e-12*abs(Q0)) break
            st = st/2
        }
        if (j > 30) break
        B   = Bn
        lnP = lnPn
        Q0  = Q1
        if (max(abs(st*step)) < 1e-10) break
    }
    return(B)
}

// M-step for the transitions by transition row: one fractional multinomial
// logit per origin state, on the joint posteriors XI
real rowvector hmdbp_mtrans(struct hmdbp_dat scalar M, real rowvector tau,
    real matrix XI)
{
    real matrix    TT, B0, Bn
    real rowvector cc, idx
    real scalar    c, cb, j, i, k

    k  = M.k2
    TT = colshape(tau, M.qw)
    j  = 0
    for (c=1; c<=k; c++) {
        cc  = ((c-1)*k+1)..(c*k)
        B0  = J(M.qw, k, 0)
        idx = J(1, k, 0)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                j++
                B0[.,cb] = TT[j,.]'
                idx[cb]  = j
            }
        }
        Bn = hmdbp_fmlogit(M.W, XI[., cc], B0, c, 5)
        for (cb=1; cb<=k; cb++) {
            if (cb != c) {
                i = idx[cb]
                TT[i,.] = Bn[.,cb]'
            }
        }
    }
    return(vec(TT')')
}

// log likelihood (and posteriors if post != 0); missing if b is not admissible
real scalar hmdbp_estep(struct hmdbp_dat scalar M, real rowvector b,
    real matrix W, real matrix Zs, real matrix PiV, real matrix PiU,
    real matrix Pv, struct hmdbp_mm scalar mm, real scalar post,
    real scalar wantxi, real matrix XiI)
{
    real rowvector b1, t1, b2, t2, bd, td, tau
    real matrix    Dv, Gu, Pi0, P, A, F
    real colvector bnd, ll
    real scalar    lam

    hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    mm = hmdbp_means(M, b1, t1, b2, t2)
    if (hasmissing(mm.lm0)) return(.)
    if (M.chk) {
        // lambda must be inside its admissible interval
        bnd = hmdbp_lambnd(mm)
        if (lam <= bnd[1] | lam >= bnd[2]) return(.)
    }
    hmdbp_dep(M, mm, lam)
    if (hasmissing(mm.lm)) return(.)
    hmdbp_dropf(M, bd, td, mm)
    if (hasmissing(mm.ld)) return(.)
    PiV = hmdbp_mlogit(M.Z,  Dv, M.k2)
    PiU = hmdbp_mlogit(M.ZU, Gu, M.k1)
    Pi0 = (PiU*M.Euf') :* (PiV*M.Ev')
    if (M.het) {
        // transitions by transition row: Pv is rows(M.W) x k2*k2, XiI holds
        // the joint posteriors by transition row (summed over the classes) and
        // Zs their sum (k2 x k2)
        Pv = hmdbp_PT(M, tau)
        if (hasmissing(Pv)) return(.)
        ll = hmdbp_fwdh(M, mm.lm + mm.ld, Pi0, Pv, A, F)
        if (hasmissing(ll)) return(.)
        if (post) {
            hmdbp_bwdh(M, F, A, Pv, W, XiI)
            Zs = colshape(colsum(XiI), M.k2)
        }
        return(quadsum(ll))
    }
    if (M.ntau > 0) Pv = hmdbp_P(tau, M.k2)
    else            Pv = I(M.k2)
    P  = I(M.k1) # Pv
    ll = hmdbp_fwd(M, mm.lm + mm.ld, Pi0, P, A, F)
    if (hasmissing(ll)) return(.)
    if (post) hmdbp_bwd(M, F, A, P, W, Zs, wantxi, XiI)
    return(quadsum(ll))
}

// ------------------------------------------------- scores and derivatives
// observation-level scores of log f with respect to eta1, eta2, lambda
void hmdbp_mscore(struct hmdbp_dat scalar M, struct hmdbp_mm scalar mm,
    real scalar lam, real matrix S1, real matrix S2, real matrix SL)
{
    real scalar q

    q  = 1 - exp(-1)
    S1 = M.y1 :- mm.Mu1 + lam*(mm.U2:*(q*mm.Mu1:*mm.D1)):/mm.G
    S2 = M.y2 :- mm.Mu2 + lam*(mm.U1:*(q*mm.Mu2:*mm.D2)):/mm.G
    SL = (mm.U1:*mm.U2):/mm.G
}

// analytic score of the log likelihood (Fisher identity)
real rowvector hmdbp_grad(struct hmdbp_dat scalar M, real rowvector b,
    real matrix W, real matrix Zs, real matrix PiV, real matrix PiU,
    real matrix Pv, struct hmdbp_mm scalar mm, real matrix XiI)
{
    real rowvector g, gt, cc
    real matrix    S1, S2, SL, wS1, wS2, Rd, W1, Zv, Xc, Gm
    real scalar    c, cb, j

    hmdbp_mscore(M, mm, b[M.m], S1, S2, SL)
    wS1 = W:*S1
    wS2 = W:*S2
    g = J(1, 0, .)
    if (M.p1 > 0) g = (g, cross(M.X1, rowsum(wS1))')
    g = (g, colsum(wS1)*M.C1)
    if (M.s1 > 0) g = (g, vec(cross(M.XS1, wS1*M.Ev))')
    if (M.p2 > 0) g = (g, cross(M.X2, rowsum(wS2))')
    g = (g, colsum(wS2)*M.C1)
    if (M.s2 > 0) g = (g, vec(cross(M.XS2, wS2*M.Ev))')
    g = (g, quadsum(W:*SL))
    if (M.hasD) {
        Rd = W :* (M.risk :* (M.D :- mm.Pd))
        if (M.pd > 0) g = (g, cross(M.Xd, rowsum(Rd))')
        g = (g, colsum(Rd)*M.Cd)
        if (M.sd > 0) g = (g, vec(cross(M.XSd, Rd*M.Ev))')
    }
    W1 = W[M.info[.,1], .]
    if (M.k2 > 1) {
        g = (g, vec(cross(M.Z, (W1*M.Ev)[|1,2 \ .,.|] - PiV[|1,2 \ .,.|]))')
        if (M.het) {
            // one logit per origin state: sum over the transition rows of
            // w (xi(c, cb) - P(c, cb) sum_d xi(c, d))
            gt = J(1, 0, .)
            for (c=1; c<=M.k2; c++) {
                cc = ((c-1)*M.k2+1)..(c*M.k2)
                Xc = XiI[., cc]
                Gm = cross(M.W, Xc - Pv[., cc] :* rowsum(Xc))
                for (cb=1; cb<=M.k2; cb++) {
                    if (cb != c) gt = (gt, Gm[.,cb]')
                }
            }
            g = (g, gt)
        }
        else if (M.ntau > 0) {
            Zv = hmdbp_zv(M, Zs)
            gt = J(1, M.ntau, .)
            j  = 0
            for (c=1; c<=M.k2; c++) {
                for (cb=1; cb<=M.k2; cb++) {
                    if (cb != c) {
                        j++
                        gt[j] = Zv[c,cb] - Pv[c,cb]*sum(Zv[c,.])
                    }
                }
            }
            g = (g, gt)
        }
    }
    if (M.k1 > 1) {
        g = (g, vec(cross(M.ZU, (W1*M.Euf)[|1,2 \ .,.|] - PiU[|1,2 \ .,.|]))')
    }
    return(g)
}

real rowvector hmdbp_gradonly(struct hmdbp_dat scalar M, real rowvector b)
{
    struct hmdbp_mm scalar mm
    real matrix W, Zs, PiV, PiU, Pv, XiI
    real scalar ll

    ll = hmdbp_estep(M, b, W, Zs, PiV, PiU, Pv, mm, 1, 0, XiI)
    if (ll >= .) return(J(1, cols(b), .))
    return(hmdbp_grad(M, b, W, Zs, PiV, PiU, Pv, mm, XiI))
}

// individual scores of the state-specific coefficients (n x k2*s): for each
// state, the weighted residuals summed over its expanded states times XS
real matrix hmdbp_sscore(struct hmdbp_dat scalar M, real matrix XS,
    real matrix wS)
{
    real matrix R, Sv
    real scalar v

    Sv = wS*M.Ev
    R  = J(M.n, 0, .)
    for (v=1; v<=M.k2; v++) R = (R, panelsum(Sv[., v]:*XS, M.info))
    return(R)
}

// individual contributions to the score (n x npar), for vce(robust)
real matrix hmdbp_scores(struct hmdbp_dat scalar M, real rowvector b)
{
    struct hmdbp_mm scalar mm
    real matrix  W, Zs, PiV, PiU, Pv, XiI, S1, S2, SL, wS1, wS2, Rd, E, Sc, St, Xv, Xc, Rc, Ef
    real rowvector r, cc
    real scalar  ll, c, cb, j, u, k2, K

    k2 = M.k2
    K  = M.K
    ll = hmdbp_estep(M, b, W, Zs, PiV, PiU, Pv, mm, 1, 1, XiI)
    hmdbp_mscore(M, mm, b[M.m], S1, S2, SL)
    wS1 = W:*S1
    wS2 = W:*S2
    Sc = J(M.n, 0, .)
    if (M.p1 > 0) Sc = (Sc, panelsum(rowsum(wS1):*M.X1, M.info))
    Sc = (Sc, panelsum(wS1, M.info)*M.C1)
    if (M.s1 > 0) Sc = (Sc, hmdbp_sscore(M, M.XS1, wS1))
    if (M.p2 > 0) Sc = (Sc, panelsum(rowsum(wS2):*M.X2, M.info))
    Sc = (Sc, panelsum(wS2, M.info)*M.C1)
    if (M.s2 > 0) Sc = (Sc, hmdbp_sscore(M, M.XS2, wS2))
    Sc = (Sc, panelsum(rowsum(W:*SL), M.info))
    if (M.hasD) {
        Rd = W :* (M.risk :* (M.D :- mm.Pd))
        if (M.pd > 0) Sc = (Sc, panelsum(rowsum(Rd):*M.Xd, M.info))
        Sc = (Sc, panelsum(Rd, M.info)*M.Cd)
        if (M.sd > 0) Sc = (Sc, hmdbp_sscore(M, M.XSd, Rd))
    }
    if (k2 > 1) {
        E = W[M.info[.,1], .]*M.Ev - PiV
        for (c=2; c<=k2; c++) Sc = (Sc, M.Z :* E[.,c])
        if (M.het) {
            // by transition row, then summed within the individual
            for (c=1; c<=k2; c++) {
                cc = ((c-1)*k2+1)..(c*k2)
                Xc = XiI[., cc]
                Rc = Xc - Pv[., cc] :* rowsum(Xc)
                for (cb=1; cb<=k2; cb++) {
                    if (cb != c) {
                        Ef = J(M.N, M.qw, 0)
                        Ef[M.trow, .] = M.W :* Rc[.,cb]
                        Sc = (Sc, panelsum(Ef, M.info))
                    }
                }
            }
        }
        else if (M.ntau > 0) {
            // expected transitions of each individual, from c to cb (summed
            // over the classes)
            St = J(M.n, M.ntau, .)
            j  = 0
            for (c=1; c<=k2; c++) {
                Xv = J(M.n, k2, 0)
                for (u=1; u<=M.k1; u++) {
                    r  = ((u-1)*k2 + c - 1)*K :+ ((u-1)*k2 :+ (1..k2))
                    Xv = Xv + XiI[., r]
                }
                for (cb=1; cb<=k2; cb++) {
                    if (cb != c) {
                        j++
                        St[.,j] = Xv[.,cb] - Pv[c,cb]*rowsum(Xv)
                    }
                }
            }
            Sc = (Sc, St)
        }
    }
    if (M.k1 > 1) {
        E = W[M.info[.,1], .]*M.Euf - PiU
        for (u=2; u<=M.k1; u++) Sc = (Sc, M.ZU :* E[.,u])
    }
    return(Sc)
}

// H[r0.., c0..] = B (nothing to do when B is empty)
void hmdbp_blk(real matrix H, real scalar r0, real scalar c0, real matrix B)
{
    if (rows(B) > 0 & cols(B) > 0) {
        H[|r0,c0 \ r0+rows(B)-1,c0+cols(B)-1|] = B
    }
}

// gradient and Hessian of the expected complete log likelihood, measurement
// block. They are computed for the expanded parameters of each equation
// (common coefficients, K intercepts, K x s state-specific coefficients) and
// for lambda, and mapped to the parameters of the model with Jm. In expanded
// state s the linear predictor of equation j is Zc_j b_j(s), with
// Zc_j = (X_j, 1, XS_j): the contributions are accumulated state by state on
// the elements idx_j(s) of the expanded vector.
void hmdbp_mderiv(struct hmdbp_dat scalar M, real matrix W,
    struct hmdbp_mm scalar mm, real scalar lam, real rowvector g,
    real matrix H)
{
    real matrix    S1, S2, SL, h1, h2, G2, r1, r2, wS1, wS2
    real matrix    w11, w22, w12, w1L, w2L, He, Z1, Z2, B
    real rowvector ge, i1, i2
    real scalar    q, K, n1, n2, ne, s

    q  = 1 - exp(-1)
    h1 = q*mm.Mu1:*mm.D1
    h2 = q*mm.Mu2:*mm.D2
    G2 = mm.G:^2
    r1 = lam*(mm.U2:*h1):/mm.G
    r2 = lam*(mm.U1:*h2):/mm.G
    S1 = M.y1 :- mm.Mu1 + r1
    S2 = M.y2 :- mm.Mu2 + r2
    SL = (mm.U1:*mm.U2):/mm.G
    // weighted second derivatives with respect to (eta1, eta2, lambda)
    w11 = W:*(-mm.Mu1 + r1:*(1:-q*mm.Mu1) - r1:^2)
    w22 = W:*(-mm.Mu2 + r2:*(1:-q*mm.Mu2) - r2:^2)
    w12 = W:*((lam*h1:*h2):/G2)
    w1L = W:*((mm.U2:*h1):/G2)
    w2L = W:*((mm.U1:*h2):/G2)
    wS1 = W:*S1
    wS2 = W:*S2

    K  = M.K
    n1 = M.p1 + K + K*M.s1
    n2 = M.p2 + K + K*M.s2
    ne = n1 + n2 + 1
    ge = J(1, ne, 0)
    He = J(ne, ne, 0)
    Z1 = (M.X1, J(M.N, 1, 1), M.XS1)
    Z2 = (M.X2, J(M.N, 1, 1), M.XS2)
    for (s=1; s<=K; s++) {
        i1 = hmdbp_eidx(M.p1, K, M.s1, s, 0)
        i2 = hmdbp_eidx(M.p2, K, M.s2, s, n1)
        ge[i1] = ge[i1] + cross(Z1, wS1[., s])'
        ge[i2] = ge[i2] + cross(Z2, wS2[., s])'
        He[i1, i1] = He[i1, i1] + cross(Z1, w11[., s], Z1)
        He[i2, i2] = He[i2, i2] + cross(Z2, w22[., s], Z2)
        B = cross(Z1, w12[., s], Z2)
        He[i1, i2] = He[i1, i2] + B
        He[i2, i1] = He[i2, i1] + B'
        B = cross(Z1, w1L[., s])
        He[i1, ne] = He[i1, ne] + B
        He[ne, i1] = He[ne, i1] + B'
        B = cross(Z2, w2L[., s])
        He[i2, ne] = He[i2, ne] + B
        He[ne, i2] = He[ne, i2] + B'
    }
    ge[ne]    = quadsum(W:*SL)
    He[ne,ne] = -quadsum(W:*(SL:^2))

    g = ge*M.Jm
    H = M.Jm'*He*M.Jm
}

// elements of the expanded vector of one equation (p common coefficients, K
// intercepts, K blocks of s state-specific coefficients) that enter the
// linear predictor of expanded state st, after an offset o
real rowvector hmdbp_eidx(real scalar p, real scalar K, real scalar s,
    real scalar st, real scalar o)
{
    real rowvector i

    i = J(1, 0, .)
    if (p > 0) i = 1..p
    i = (i, p + st)
    if (s > 0) i = (i, (p + K + (st-1)*s + 1)..(p + K + st*s))
    return(o :+ i)
}

// ---------------------------------------------------------------- M-step
// Measurement block: Newton steps on the expected complete log likelihood,
// starting from bm0, where the measurement model is mm0 (from the E-step). The
// step is halved until the objective does not decrease. Every trial point is
// made admissible by projecting lambda onto its admissible interval at the
// trial means. A lambda on a bound of its interval with the gradient pointing
// outwards is left out of the Newton system (active set).
real rowvector hmdbp_mmeas(struct hmdbp_dat scalar M, real matrix W,
    real rowvector bm0, real scalar nsteps, struct hmdbp_mm scalar mm0)
{
    struct hmdbp_mm scalar mm, mn
    real rowvector bm, bn, g, step, s, b1, t1, b2, t2, fr
    real matrix    H
    real colvector bnd
    real scalar    it, j, st, Q0, Q1, eps, lam, w, lv, blk

    eps = 1e-6
    bm  = bm0
    mm  = mm0
    for (it=1; it<=nsteps; it++) {
        lam = bm[M.m]
        Q0  = quadsum(W:*mm.lm)
        hmdbp_mderiv(M, W, mm, lam, g, H)
        // lambda on a bound, with the gradient pointing outwards
        bnd = hmdbp_lambnd(mm)
        w   = bnd[2] - bnd[1]
        blk = ((lam >= bnd[2] - 2*eps*w) & (g[M.m] > 0)) | ((lam <= bnd[1] + 2*eps*w) & (g[M.m] < 0))
        if (blk) {
            fr = 1..(M.m-1)
            s  = cholsolve(-H[fr, fr], g[fr]')'
            if (hasmissing(s)) s = (invsym(-H[fr, fr])*g[fr]')'
            if (hasmissing(s) | (g[fr]*s') <= 0) s = g[fr] :/ max(abs(diagonal(H)))
            step = (s, 0)
        }
        else {
            step = cholsolve(-H, g')'
            if (hasmissing(step)) step = (invsym(-H)*g')'
            if (hasmissing(step) | (g*step') <= 0) step = g :/ max(abs(diagonal(H)))
        }
        st = 1
        for (j=1; j<=40; j++) {
            bn = bm + st*step
            hmdbp_mparts(M, bn, b1, t1, b2, t2, lam)
            mn = hmdbp_means(M, b1, t1, b2, t2)
            if (!hasmissing(mn.lm0)) {
                // projection of lambda onto the admissible interval at the trial means
                bnd = hmdbp_lambnd(mn)
                w   = bnd[2] - bnd[1]
                lv  = max((min((bn[M.m], bnd[2] - eps*w)), bnd[1] + eps*w))
                bn[M.m] = lv
                hmdbp_dep(M, mn, lv)
                if (!hasmissing(mn.lm)) {
                    Q1 = quadsum(W:*mn.lm)
                    if (Q1 >= Q0 - 1e-12*abs(Q0)) break
                }
            }
            st = st/2
        }
        if (j > 40) break
        bm = bn
        mm = mn
        if (max(abs(st*step)) < 1e-10) break
    }
    return(bm)
}

// Dropout block: Newton steps with step-halving on the weighted logit of the
// rows at risk (expanded parameters mapped with Jd, as in hmdbp_mderiv)
real rowvector hmdbp_mdrop(struct hmdbp_dat scalar M, real matrix W,
    real rowvector bdv0, real scalar nsteps)
{
    real rowvector bdv, bn, bd, td, g, step, ge, i
    real matrix    Wr, E, En, Pr, Rr, Ww, He, H, Zd
    real scalar    it, j, st, Q0, Q1, K, nd, s

    K   = M.K
    nd  = M.pd + K + K*M.sd
    Wr  = W :* M.risk
    Zd  = (M.Xd, J(M.N, 1, 1), M.XSd)
    bdv = bdv0
    hmdbp_dparts(M, bdv, bd, td)
    E   = hmdbp_eta(M, M.Xd, M.XSd, bd, td, M.Cd, M.ncd, M.sd)
    Q0  = quadsum(Wr :* (M.D :* hmdbp_lnsig(E) + (1 :- M.D) :* hmdbp_lnsig(-E)))
    for (it=1; it<=nsteps; it++) {
        Pr = invlogit(E)
        Rr = Wr :* (M.D :- Pr)
        Ww = Wr :* Pr :* (1 :- Pr)
        ge = J(1, nd, 0)
        He = J(nd, nd, 0)
        for (s=1; s<=K; s++) {
            i = hmdbp_eidx(M.pd, K, M.sd, s, 0)
            ge[i]    = ge[i] + cross(Zd, Rr[., s])'
            He[i, i] = He[i, i] - cross(Zd, Ww[., s], Zd)
        }
        g = ge*M.Jd
        H = M.Jd'*He*M.Jd
        step = cholsolve(-H, g')'
        if (hasmissing(step)) step = (invsym(-H)*g')'
        if (hasmissing(step)) break
        st = 1
        for (j=1; j<=30; j++) {
            bn = bdv + st*step
            hmdbp_dparts(M, bn, bd, td)
            En = hmdbp_eta(M, M.Xd, M.XSd, bd, td, M.Cd, M.ncd, M.sd)
            Q1 = quadsum(Wr :* (M.D :* hmdbp_lnsig(En) + (1 :- M.D) :* hmdbp_lnsig(-En)))
            if (Q1 >= Q0 - 1e-12*abs(Q0)) break
            st = st/2
        }
        if (j > 30) break
        bdv = bn
        E   = En
        Q0  = Q1
        if (max(abs(st*step)) < 1e-10) break
    }
    return(bdv)
}

// fractional multinomial logit (Newton-Raphson with step-halving): responses
// W1 (n x k, posterior probabilities), coefficients D0 (q x (k-1)), base = 1
real matrix hmdbp_minit(real matrix Z, real matrix W1, real matrix D0,
    real scalar k, real scalar nsteps)
{
    real matrix    D, Dn, Pn, Pi, Hm, blk
    real colvector g, step, wg
    real scalar    q, it, c, e, j, st, Q0, Q1

    q  = cols(Z)
    D  = D0
    Pi = hmdbp_mlogit(Z, D, k)
    Q0 = quadsum(W1:*ln(Pi))
    for (it=1; it<=nsteps; it++) {
        g  = vec(cross(Z, W1[|1,2 \ .,.|] - Pi[|1,2 \ .,.|]))
        Hm = J(q*(k-1), q*(k-1), 0)
        for (c=2; c<=k; c++) {
            for (e=c; e<=k; e++) {
                wg  = Pi[.,c] :* ((c==e) :- Pi[.,e])
                blk = -cross(Z, wg, Z)
                Hm[|(c-2)*q+1,(e-2)*q+1 \ (c-1)*q,(e-1)*q|] = blk
                Hm[|(e-2)*q+1,(c-2)*q+1 \ (e-1)*q,(c-1)*q|] = blk'
            }
        }
        step = cholsolve(-Hm, g)
        if (hasmissing(step)) step = invsym(-Hm)*g
        if (hasmissing(step)) break
        st = 1
        for (j=1; j<=30; j++) {
            Dn = D + st*colshape(step', q)'
            Pn = hmdbp_mlogit(Z, Dn, k)
            Q1 = quadsum(W1:*ln(Pn))
            if (Q1 >= Q0 - 1e-12*abs(Q0)) break
            st = st/2
        }
        if (j > 30) break
        D  = Dn
        Pi = Pn
        Q0 = Q1
        if (max(abs(st*step)) < 1e-10) break
    }
    return(D)
}

// ------------------------------------------------------------------ EM
real rowvector hmdbp_em(struct hmdbp_dat scalar M, real rowvector b0,
    real scalar maxit, real scalar tol, real scalar dolog,
    real scalar itout, real scalar llout)
{
    struct hmdbp_mm scalar mm
    real rowvector b, bm, bdv, b1, t1, b2, t2, bd, td, tau
    real matrix    W, Zs, PiV, PiU, Pv, XiI, Dv, Gu, W1, Zv, Pn
    real colvector rs
    real scalar    it, ll, llold, c, lam

    b = b0
    llold = .
    ll = .
    for (it=1; it<=maxit; it++) {
        ll = hmdbp_estep(M, b, W, Zs, PiV, PiU, Pv, mm, 1, 0, XiI)
        if (ll >= .) {
            // with multiple starting values a failed start is discarded
            if (M.emsoft == 1) {
                itout = it
                llout = .
                return(b)
            }
            errprintf("EM: log likelihood not finite at iteration %5.0f\n", it)
            exit(430)
        }
        if (dolog) {
            printf("{txt}EM iteration %5.0f:{col 22}log likelihood = {res}%12.0g\n", it, ll)
            displayflush()
        }
        if (it > 1) {
            if (abs(ll - llold) <= tol*abs(llold)) break
        }
        llold = ll
        hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
        // measurement and dropout blocks: Newton steps with step-halving
        bm = hmdbp_mmeas(M, W, b[|1 \ M.m|], 2, mm)
        hmdbp_mparts(M, bm, b1, t1, b2, t2, lam)
        if (M.hasD) {
            bdv = hmdbp_mdrop(M, W, (bd, td), 2)
            hmdbp_dparts(M, bdv, bd, td)
        }
        W1 = W[M.info[.,1], .]
        if (M.k2 > 1) {
            // initial probabilities of the states: fractional multinomial logit
            Dv = hmdbp_minit(M.Z, W1*M.Ev, Dv, M.k2, 5)
            if (M.het) {
                // transitions by transition row: one logit per origin state
                tau = hmdbp_mtrans(M, tau, XiI)
            }
            else if (M.ntau > 0) {
                // transitions: closed form
                Zv = hmdbp_zv(M, Zs)
                rs = rowsum(Zv)
                Pn = Pv
                for (c=1; c<=M.k2; c++) {
                    if (rs[c] > 1e-300) Pn[c,.] = Zv[c,.] :/ rs[c]
                }
                Pn  = Pn :+ 1e-10
                Pn  = Pn :/ rowsum(Pn)
                tau = hmdbp_tau(Pn)
            }
        }
        // prior probabilities of the classes: fractional multinomial logit
        if (M.k1 > 1) Gu = hmdbp_minit(M.ZU, W1*M.Euf, Gu, M.k1, 5)
        b = hmdbp_pack(M, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    }
    itout = (it > maxit ? maxit : it)
    llout = ll
    return(b)
}

// evaluator for optimize(), type d1
void hmdbp_eval(real scalar todo, real rowvector b,
    pointer(struct hmdbp_dat scalar) scalar pM,
    real scalar v, real rowvector g, real matrix H)
{
    struct hmdbp_mm scalar mm
    real matrix W, Zs, PiV, PiU, Pv, XiI

    v = hmdbp_estep(*pM, b, W, Zs, PiV, PiU, Pv, mm, (todo >= 1), 0, XiI)
    if (v >= .) return
    if (todo >= 1) g = hmdbp_grad(*pM, b, W, Zs, PiV, PiU, Pv, mm, XiI)
}

// Observed information: central differences of the analytic score with respect
// to the free parameters fr (a lambda on a bound is held fixed). The
// admissibility check of the E-step is switched off here: the log likelihood is
// smooth across the bounds.
real matrix hmdbp_hess(struct hmdbp_dat scalar M, real rowvector b,
    real rowvector fr)
{
    real matrix    H
    real rowvector bp, bq, gp, gq
    real scalar    i, j, h, nf

    nf = cols(fr)
    H  = J(nf, nf, .)
    M.chk = 0
    for (i=1; i<=nf; i++) {
        j  = fr[i]
        h  = 1e-5*(abs(b[j]) + 1)
        bp = b; bp[j] = b[j] + h
        bq = b; bq[j] = b[j] - h
        gp = hmdbp_gradonly(M, bp)
        gq = hmdbp_gradonly(M, bq)
        H[.,i] = (gp[fr] - gq[fr])' :/ (2*h)
    }
    M.chk = 1
    return((H + H')/2)
}

// Relabel the states in increasing order of alpha_1v (how = "eq1") or alpha_2v
// (how = "eq2"), and the classes in increasing order of the dropout shift
// (delta_1u or delta_2u without the dropout equation). The log likelihood
// does not change.
real rowvector hmdbp_sort(struct hmdbp_dat scalar M, real rowvector b,
    string scalar how)
{
    real rowvector b1, t1, b2, t2, bd, td, tau, key, df, v1, v2, vd
    real matrix    Dv, Gu, Df, TT, TN, pos
    real colvector perm
    real scalar    c, cb, j, k, k1, s0, lam

    if (how == "none") return(b)
    k  = M.k2
    k1 = M.k1
    hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    // state-specific coefficients apart (they follow the intercepts)
    v1 = hmdbp_sub(t1, M.nc1+1, M.nt1)
    v2 = hmdbp_sub(t2, M.nc1+1, M.nt2)
    t1 = t1[|1 \ M.nc1|]
    t2 = t2[|1 \ M.nc1|]
    if (M.hasD) {
        vd = hmdbp_sub(td, M.ncd+1, M.ntd)
        td = td[|1 \ M.ncd|]
    }
    // ------------------------------------------------------------ states
    if (k > 1) {
        if (how == "eq2") key = t2[|1 \ k|]
        else              key = t1[|1 \ k|]
        perm = order(key', 1)
        t1[|1 \ k|] = t1[perm']
        t2[|1 \ k|] = t2[perm']
        if (M.hasD) td[|1 \ k|] = td[perm']
        v1 = hmdbp_sperm(v1, M.s1, perm)
        v2 = hmdbp_sperm(v2, M.s2, perm)
        if (M.hasD) vd = hmdbp_sperm(vd, M.sd, perm)
        Df = (J(M.q, 1, 0), Dv)[., perm']
        Df = Df :- Df[., 1]
        Dv = Df[|1,2 \ .,.|]
        if (M.ntau > 0) {
            // the base of each row is "staying", which relabelling preserves
            pos = J(k, k, 0)
            j = 0
            for (c=1; c<=k; c++) {
                for (cb=1; cb<=k; cb++) {
                    if (cb != c) {
                        j++
                        pos[c,cb] = j
                    }
                }
            }
            // rows of TT: pairs (from, to); columns: the qw coefficients
            TT = colshape(tau, M.qw)
            TN = TT
            for (c=1; c<=k; c++) {
                for (cb=1; cb<=k; cb++) {
                    if (cb != c) TN[pos[c,cb], .] = TT[pos[perm[c],perm[cb]], .]
                }
            }
            tau = vec(TN')'
        }
    }
    // ----------------------------------------------------------- classes
    if (k1 > 1) {
        // key: the dropout shifts when the classes enter the dropout
        // equation, otherwise those of the count equation of sortby()
        if (M.cdd)             key = (0, td[|k+1 \ M.ncd|])
        else if (how == "eq2") key = (0, t2[|k+1 \ M.nc1|])
        else                  key = (0, t1[|k+1 \ M.nc1|])
        perm = order(key', 1)
        // shifts relative to the new base class; the state intercepts absorb
        // the shift of the new base
        if (M.cm) {
            df = (0, t1[|k+1 \ M.nc1|])[perm']
            s0 = df[1]
            t1 = (t1[|1 \ k|] :+ s0, df[|2 \ k1|] :- s0)
            df = (0, t2[|k+1 \ M.nc1|])[perm']
            s0 = df[1]
            t2 = (t2[|1 \ k|] :+ s0, df[|2 \ k1|] :- s0)
        }
        if (M.cdd) {
            df = (0, td[|k+1 \ M.ncd|])[perm']
            s0 = df[1]
            td = (td[|1 \ k|] :+ s0, df[|2 \ k1|] :- s0)
        }
        Df = (J(M.qu, 1, 0), Gu)[., perm']
        Df = Df :- Df[., 1]
        Gu = Df[|1,2 \ .,.|]
    }
    t1 = (t1, v1)
    t2 = (t2, v2)
    if (M.hasD) td = (td, vd)
    return(hmdbp_pack(M, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu))
}

// ------------------------------------------------------------------ driver
// data matrix for a list of variable names (N x 0 when the list is empty)
real matrix hmdbp_getx(string scalar names, string scalar touse, real scalar N)
{
    if (strtrim(names) == "") return(J(N, 0, .))
    return(st_data(., tokens(names), touse))
}

// Data and dimensions of the model, from the variables in memory: the rows
// selected by touse, sorted by decreasing T, individual and period
void hmdbp_setup(struct hmdbp_dat scalar M, string scalar y1n,
    string scalar y2n, string scalar x1n, string scalar x2n,
    string scalar ydn, string scalar xdn, string scalar riskn,
    string scalar zn, string scalar zun, string scalar idn,
    string scalar touse, real scalar k2, real scalar k1, real scalar notr,
    string scalar wn, string scalar trn, string scalar timen,
    real scalar het, real scalar tvt, real colvector u, real scalar nwx)
{
    real colvector id, sel
    real matrix    zr, Wx
    real scalar    t, o, oo, ne, K, nc, nt
    string scalar  ce

    M.k1   = k1
    M.k2   = k2
    M.K    = k1*k2
    M.hasD = (strtrim(ydn) != "")
    M.notr = notr
    M.chk  = 1
    M.y1   = st_data(., y1n, touse)
    M.y2   = st_data(., y2n, touse)
    M.ey1  = exp(-M.y1)
    M.ey2  = exp(-M.y2)
    M.lnfy = lngamma(M.y1 :+ 1) + lngamma(M.y2 :+ 1)
    M.N    = rows(M.y1)
    M.X1   = hmdbp_getx(x1n, touse, M.N)
    M.X2   = hmdbp_getx(x2n, touse, M.N)
    if (M.hasD) {
        M.Xd   = hmdbp_getx(xdn, touse, M.N)
        M.D    = st_data(., ydn, touse)
        M.risk = st_data(., riskn, touse)
    }
    else {
        M.Xd   = J(M.N, 0, .)
        M.D    = J(M.N, 1, 0)
        M.risk = J(M.N, 1, 0)
    }
    // regressors with state-specific coefficients (locals svar1, svar2, svard
    // of the caller) and the codes of the terms of each equation (0
    // regressor, m = 1..4 lag feature m, 10 + m feature m of the current
    // counts), from the names in the locals cnam* (common terms) and snam*
    // (state-specific terms)
    M.XS1 = hmdbp_getx(st_local("svar1"), touse, M.N)
    M.XS2 = hmdbp_getx(st_local("svar2"), touse, M.N)
    if (M.hasD) M.XSd = hmdbp_getx(st_local("svard"), touse, M.N)
    else        M.XSd = J(M.N, 0, .)
    M.cc1 = hmdbp_codes(st_local("cnam1"))
    M.sc1 = hmdbp_codes(st_local("snam1"))
    M.cc2 = hmdbp_codes(st_local("cnam2"))
    M.sc2 = hmdbp_codes(st_local("snam2"))
    M.ccd = hmdbp_codes(st_local("cnamd"))
    M.scd = hmdbp_codes(st_local("snamd"))
    if (cols(M.cc1) != cols(M.X1) | cols(M.sc1) != cols(M.XS1) |
        cols(M.cc2) != cols(M.X2) | cols(M.sc2) != cols(M.XS2) |
        cols(M.ccd) != cols(M.Xd) | cols(M.scd) != cols(M.XSd)) {
        errprintf("hmdbp_setup(): names and designs do not match\n")
        exit(498)
    }
    id     = st_data(., idn, touse)
    M.info = panelsetup(id, 1)
    M.n    = rows(M.info)
    M.Ts   = M.info[.,2] - M.info[.,1] :+ 1
    M.Tmax = max(M.Ts)
    zr     = st_data(., tokens(zn), touse)
    M.Z    = (zr[M.info[.,1], .], J(M.n, 1, 1))
    zr     = st_data(., tokens(zun), touse)
    M.ZU   = (zr[M.info[.,1], .], J(M.n, 1, 1))
    if (M.n > 1) {
        if (any(M.Ts[|2 \ M.n|] :> M.Ts[|1 \ M.n-1|])) {
            errprintf("hmdbp_setup(): individuals must be sorted by decreasing T\n")
            exit(498)
        }
    }
    M.ns = J(M.Tmax, 1, .)
    M.R  = J(M.Tmax, 1, NULL)
    for (t=1; t<=M.Tmax; t++) {
        M.ns[t] = sum(M.Ts :>= t)
        M.R[t]  = &(M.info[|1,1 \ M.ns[t],1|] :+ (t-1))
    }

    // transition design (xtrans, tvtrans): one row per transition row (the
    // estimation rows after the first of each individual), in the same order
    // as the data; tt is the period of arrival, u its distinct values
    M.het = het
    M.qw  = 1
    nwx   = 0
    u     = J(0, 1, .)
    M.tt  = J(0, 1, .)
    if (het) {
        M.tt = st_data(., timen, trn)
        if (strtrim(wn) != "") Wx = st_data(., tokens(wn), trn)
        else                   Wx = J(rows(M.tt), 0, .)
        nwx = cols(Wx)
        u   = uniqrows(M.tt)
        if (tvt) M.W = (Wx, (J(1, rows(u), M.tt) :== J(rows(M.tt), 1, u')))
        else     M.W = (Wx, J(rows(M.tt), 1, 1))
        M.qw = cols(M.W)
        if (rows(M.W) != M.N - M.n) {
            errprintf("hmdbp_setup(): transition rows do not match the panel structure\n")
            exit(498)
        }
        sel = J(M.N, 1, 1)
        sel[M.info[.,1]] = J(M.n, 1, 0)
        M.trow = selectindex(sel)
        M.RT = J(M.Tmax, 1, NULL)
        for (t=2; t<=M.Tmax; t++) {
            M.RT[t] = &((*M.R[t]) - (1::M.ns[t]))
        }
    }

    // dimensions
    K    = M.K
    M.p1 = cols(M.X1)
    M.p2 = cols(M.X2)
    M.pd = cols(M.Xd)
    M.nc = k2 + k1 - 1
    nc   = M.nc
    // classes in the count equations (cm) and in the dropout equation (cdd):
    // local clseq of the caller (all, dropout, counts)
    ce    = st_local("clseq")
    M.cm  = (k1 > 1 & ce != "dropout")
    M.cdd = (k1 > 1 & ce != "counts" & M.hasD)
    M.nc1 = k2 + M.cm*(k1-1)
    M.ncd = k2 + M.cdd*(k1-1)
    // state-specific coefficients: s per state in each equation; the vector
    // t of an equation holds the intercepts, the shifts of the classes and
    // the k2*s state-specific coefficients
    M.s1  = cols(M.XS1)
    M.s2  = cols(M.XS2)
    M.sd  = cols(M.XSd)
    M.nt1 = M.nc1 + k2*M.s1
    M.nt2 = M.nc1 + k2*M.s2
    M.ntd = M.ncd + k2*M.sd
    M.q  = cols(M.Z)
    M.qu = cols(M.ZU)
    M.m  = M.p1 + M.p2 + M.nt1 + M.nt2 + 1
    if (M.hasD) M.md = M.pd + M.ntd
    else        M.md = 0
    // lag features per set (4 with logzero)
    M.nl = (st_local("lagt") == "logzero" ? 4 : 2)
    M.nD   = M.q*(k2-1)
    if (notr | k2 == 1) M.ntau = 0
    else                M.ntau = k2*(k2-1)*M.qw
    M.nU   = M.qu*(k1-1)
    M.npar = M.m + M.md + M.nD + M.ntau + M.nU

    // expanded states: Ev (K x k2) and Euf (K x k1) are the indicators of the
    // state and of the class of each expanded state
    M.Ev  = J(k1, 1, 1) # I(k2)
    M.Euf = I(k1) # J(k2, 1, 1)
    if (k1 > 1) M.C = (M.Ev, M.Euf[|1,2 \ .,.|])
    else        M.C = M.Ev
    // maps of the intercepts of the count equations (C1) and of the dropout
    // equation (Cd): with the class shifts, or the states only
    M.C1 = (M.cm  ? M.C : M.Ev)
    M.Cd = (M.cdd ? M.C : M.Ev)
    // Jm: expanded measurement parameters = Jm * measurement block. Expanded:
    // for each count equation the common coefficients, K intercepts and K x s
    // state-specific coefficients (expanded state 1 first); then lambda
    ne   = M.p1 + K + K*M.s1 + M.p2 + K + K*M.s2 + 1
    M.Jm = J(ne, M.m, 0)
    hmdbp_jblk(M, M.Jm, 0, 0, M.p1, M.C1, M.nc1, M.s1)
    hmdbp_jblk(M, M.Jm, M.p1 + K + K*M.s1, M.p1 + M.nt1, M.p2, M.C1, M.nc1, M.s2)
    M.Jm[ne, M.m] = 1
    // Jd: expanded dropout parameters = Jd * dropout block
    M.Jd = J(M.pd + K + K*M.sd, M.pd + M.ntd, 0)
    hmdbp_jblk(M, M.Jd, 0, 0, M.pd, M.Cd, M.ncd, M.sd)
}

// state-specific coefficients (s per state, state 1 first) of the states in
// the order perm
real rowvector hmdbp_sperm(real rowvector v, real scalar s,
    real colvector perm)
{
    real matrix S

    if (s == 0) return(v)
    S = colshape(v, s)
    S = S[perm, .]
    return(vec(S')')
}

// elements i0..i1 of a row vector (1 x 0 when i1 < i0)
real rowvector hmdbp_sub(real rowvector x, real scalar i0, real scalar i1)
{
    if (i1 < i0) return(J(1, 0, .))
    return(x[|i0 \ i1|])
}

// block of one equation in Jm or Jd, at row offset ro and column offset co:
// identity for the p common coefficients, C for the intercepts, Ev # I(s)
// for the state-specific coefficients
void hmdbp_jblk(struct hmdbp_dat scalar M, real matrix Jb, real scalar ro,
    real scalar co, real scalar p, real matrix C, real scalar nc,
    real scalar s)
{
    real scalar K

    K = M.K
    if (p > 0) Jb[|ro+1,co+1 \ ro+p,co+p|] = I(p)
    Jb[|ro+p+1,co+p+1 \ ro+p+K,co+p+nc|] = C
    if (s > 0) Jb[|ro+p+K+1,co+p+nc+1 \ ro+p+K+K*s,co+p+nc+M.k2*s|] = M.Ev # I(s)
}

// codes of a list of term names: 0 regressor, m lag feature m (local
// lagnames of the caller), 10 + m feature m of the current counts (local
// curnames)
real rowvector hmdbp_codes(string scalar names)
{
    string rowvector T, L, C
    real rowvector   cd, w
    real scalar      i

    T  = tokens(names)
    L  = tokens(st_local("lagnames"))
    C  = tokens(st_local("curnames"))
    cd = J(1, cols(T), 0)
    for (i=1; i<=cols(T); i++) {
        w = selectindex(L :== T[i])
        if (cols(w) > 0) cd[i] = w[1]
        w = selectindex(C :== T[i])
        if (cols(w) > 0) cd[i] = 10 + w[1]
    }
    return(cd)
}

// starting intercepts of one equation: alpha_v = cons + sv*offv - mean(delta),
// delta_u = su*(u-1)/(k1-1)
real rowvector hmdbp_t0(struct hmdbp_dat scalar M, real scalar cons,
    real scalar sv, real scalar su, real scalar hc)
{
    real rowvector offv, du

    if (M.k2 > 1) offv = 0.8*((1..M.k2) :- (M.k2+1)/2)/((M.k2-1)/2)
    else          offv = 0
    if (M.k1 > 1 & hc) {
        du = su*(0..(M.k1-1))/(M.k1-1)
        return((cons :+ sv*offv :- sum(du)/M.k1, du[|2 \ M.k1|]))
    }
    return(cons :+ sv*offv)
}

// Random starting values (starts() > 1), around the starting vector b0: the
// slopes are those of b0 and lambda is zero; the intercepts of the states are
// the mean intercept of b0 plus sorted normal offsets (s.d. 0.7) in the count
// equations and unsorted offsets (s.d. 0.5) in the dropout equation; the
// shifts of the classes are normal (s.d. 0.3 for the counts, 1.5 for the
// dropout); the constants of the initial and prior probabilities are N(0,1);
// each row of the transition matrix has staying probability U(0.6, 0.95) and
// the rest split at random (with xtrans() the covariates start at zero).
real rowvector hmdbp_rstart(struct hmdbp_dat scalar M, real rowvector b0,
    real scalar nwx)
{
    real rowvector b1, t1, b2, t2, bd, td, tau, th
    real matrix    Dv, Gu, P
    real scalar    lam, j, c

    hmdbp_unpack(M, b0, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    t1 = hmdbp_rint(M, t1, 0.7, 1, 0.3, M.cm)
    t2 = hmdbp_rint(M, t2, 0.7, 1, 0.3, M.cm)
    if (M.hasD) td = hmdbp_rint(M, td, 0.5, 0, 1.5, M.cdd)
    if (M.k2 > 1) {
        Dv = J(M.q, M.k2-1, 0)
        Dv[M.q, .] = rnormal(1, M.k2-1, 0, 1)
        if (M.ntau > 0) {
            P = J(M.k2, M.k2, 0)
            for (c=1; c<=M.k2; c++) {
                P[c,.] = -ln(runiform(1, M.k2))
                P[c,c] = 0
                P[c,.] = (1 - (0.6 + 0.35*runiform(1,1)))*P[c,.]/sum(P[c,.])
                P[c,c] = 1 - sum(P[c,.])
            }
            th = hmdbp_tau(P)
            if (M.het) {
                tau = J(1, 0, .)
                for (j=1; j<=cols(th); j++) {
                    tau = (tau, J(1, nwx, 0), J(1, M.qw - nwx, th[j]))
                }
            }
            else tau = th
        }
    }
    if (M.k1 > 1) {
        Gu = J(M.qu, M.k1-1, 0)
        Gu[M.qu, .] = rnormal(1, M.k1-1, 0, 1)
    }
    return(hmdbp_pack(M, b1, t1, b2, t2, 0, bd, td, Dv, tau, Gu))
}

// random intercepts of one equation (states, then shifts of classes 2..k1),
// with the same mean intercept as t0
real rowvector hmdbp_rint(struct hmdbp_dat scalar M, real rowvector t0,
    real scalar sv, real scalar srt, real scalar su, real scalar hc)
{
    real rowvector o, du, sl
    real scalar    cons, nci

    // the state-specific coefficients after the intercepts are kept
    nci = M.k2 + (M.k1 > 1 & hc)*(M.k1-1)
    sl  = hmdbp_sub(t0, nci+1, cols(t0))
    cons = mean(t0[|1 \ M.k2|]')
    if (M.k1 > 1 & hc) cons = cons + sum(t0[|M.k2+1 \ M.nc|])/M.k1
    o = rnormal(1, M.k2, 0, sv)
    if (srt) o = sort(o', 1)'
    o = o :- mean(o')
    if (M.k1 > 1 & hc) {
        du = (0, rnormal(1, M.k1-1, 0, su))
        return((cons :+ o :- sum(du)/M.k1, du[|2 \ M.k1|], sl))
    }
    return((cons :+ o, sl))
}

void hmdbp_main(string scalar y1n, string scalar y2n, string scalar x1n,
    string scalar x2n, string scalar ydn, string scalar xdn,
    string scalar riskn, string scalar zn, string scalar zun,
    string scalar idn, string scalar touse, real scalar k2,
    real scalar k1, string scalar s1n, string scalar s2n,
    string scalar sdn, string scalar fromn, real scalar emit,
    real scalar emtol, real scalar maxit, string scalar tech,
    real scalar dolog, string scalar sortby, string scalar vce,
    real scalar doml, real scalar notr, string scalar wn, string scalar trn,
    string scalar timen, real scalar het, real scalar tvt, real scalar nst,
    real scalar sit, real scalar sseed)
{
    struct hmdbp_dat scalar M
    struct hmdbp_mm scalar  mm
    real matrix    LLs
    real rowvector b0, br, bb
    real scalar    r, its, llr, best, bestll
    string scalar  rs
    real matrix    W, Zs, PiV, PiU, Pv, XiI, H, Vo, V, Sc, omt, Dv, Gu, al, de, Pm, Pit
    real rowvector b, bs1, bs2, bsd, b1, t1, b2, t2, bd, td, tau, frv, fr
    real colvector bnd, u
    real scalar    t, emits, ll, conv, rc, rank, atb, wl, lam, nwx
    transmorphic   S

    // ---------------------------------------------------------- data
    hmdbp_setup(M, y1n, y2n, x1n, x2n, ydn, xdn, riskn, zn, zun, idn,
        touse, k2, k1, notr, wn, trn, timen, het, tvt, u, nwx)

    // ---------------------------------------------------------- start
    if (fromn != "") {
        b = st_matrix(fromn)
        if (cols(b) != M.npar) {
            errprintf("from(): vector has %5.0f elements; %5.0f required\n", cols(b), M.npar)
            exit(503)
        }
    }
    else {
        // Poisson regressions and logit (regressors, lags, _cons): the state
        // intercepts of the counts are spread around the constant, the classes
        // start with increasing dropout risk
        bs1 = st_matrix(s1n)
        bs2 = st_matrix(s2n)
        // the state-specific coefficients start equal in all the states, at
        // the coefficient of the regression (which lists the common terms,
        // then the state-specific ones, then the constant)
        b1  = hmdbp_sub(bs1, 1, M.p1)
        b2  = hmdbp_sub(bs2, 1, M.p2)
        t1  = (hmdbp_t0(M, bs1[cols(bs1)], 1, 0.2, M.cm), J(1, M.k2, 1) # hmdbp_sub(bs1, M.p1+1, M.p1+M.s1))
        t2  = (hmdbp_t0(M, bs2[cols(bs2)], 1, 0.2, M.cm), J(1, M.k2, 1) # hmdbp_sub(bs2, M.p2+1, M.p2+M.s2))
        if (M.hasD) {
            bsd = st_matrix(sdn)
            bd  = hmdbp_sub(bsd, 1, M.pd)
            td  = (hmdbp_t0(M, bsd[cols(bsd)], 0, 1.5, M.cdd), J(1, M.k2, 1) # hmdbp_sub(bsd, M.pd+1, M.pd+M.sd))
        }
        else {
            bd = J(1, 0, .)
            td = J(1, 0, .)
        }
        Dv = J(M.q, M.k2-1, 0)
        // transitions start from a persistent chain (staying probability 0.9);
        // with xtrans the covariates start at zero, with tvtrans every period
        // dummy has the constant
        if (M.ntau > 0) {
            tau = J(1, M.k2*(M.k2-1), 1) # (J(1, nwx, 0), J(1, M.qw - nwx, ln(0.1/(M.k2-1)/0.9)))
        }
        else tau = J(1, 0, .)
        Gu = J(M.qu, M.k1-1, 0)
        b  = hmdbp_pack(M, b1, t1, b2, t2, 0, bd, td, Dv, tau, Gu)
    }
    ll = hmdbp_estep(M, b, W, Zs, PiV, PiU, Pv, mm, 0, 0, XiI)
    if (ll >= .) {
        errprintf("log likelihood not finite at the starting values\n")
        exit(430)
    }

    // ------------------------------------------- multiple starting values
    // EM for sit iterations from each start (the first is the default or
    // from()), then the EM below continues from the start with the largest
    // log likelihood. Starts whose log likelihood is not finite are dropped.
    if (nst > 1) {
        b0  = b
        bb  = b
        LLs = J(nst, 2, .)
        best   = .
        bestll = .
        rs = rngstate()
        rseed(sseed)
        M.emsoft = 1
        if (dolog) printf("{txt}\nEM from %g starting values (up to %g iterations each)\n", nst, sit)
        for (r=1; r<=nst; r++) {
            if (r == 1) br = b0
            else        br = hmdbp_rstart(M, b0, nwx)
            its = 0
            llr = .
            br  = hmdbp_em(M, br, sit, emtol, 0, its, llr)
            // log likelihood at the returned vector (after the last M step)
            if (llr < .) llr = hmdbp_estep(M, br, W, Zs, PiV, PiU, Pv, mm, 0, 0, XiI)
            LLs[r,.] = (llr, its)
            if (dolog) {
                if (llr < .) printf("{txt}start %4.0f:{col 15}log likelihood = {res}%14.0g{txt}  (%g iterations)\n", r, llr, its)
                else         printf("{txt}start %4.0f:{col 15}log likelihood not finite (dropped)\n", r)
                displayflush()
            }
            if (llr < .) {
                if (bestll >= . | llr > bestll) {
                    bestll = llr
                    bb     = br
                    best   = r
                }
            }
        }
        M.emsoft = 0
        rngstate(rs)
        if (best >= .) {
            errprintf("no starting value gave a finite log likelihood\n")
            exit(430)
        }
        if (dolog) printf("{txt}best: start %g\n", best)
        b = bb
        st_matrix("__hmdbp_lls", LLs)
        st_numscalar("__hmdbp_sbest", best)
    }

    // ---------------------------------------------------------- EM
    emits = 0
    if (emit > 0) {
        if (dolog) printf("{txt}\nEM algorithm\n")
        b = hmdbp_em(M, b, emit, emtol, dolog, emits, ll)
    }
    conv = (emits < emit)

    // ---------------------------------------------------------- direct ML
    if (doml == 2 | (doml == 1 & conv == 0)) {
        if (dolog) printf("{txt}\nDirect maximisation of the log likelihood\n")
        S = optimize_init()
        optimize_init_evaluator(S, &hmdbp_eval())
        optimize_init_evaluatortype(S, "d1")
        optimize_init_argument(S, 1, &M)
        optimize_init_params(S, b)
        optimize_init_technique(S, tech)
        optimize_init_conv_maxiter(S, maxit)
        optimize_init_tracelevel(S, (dolog ? "value" : "none"))
        optimize_init_conv_warning(S, "off")
        optimize_init_singularHmethod(S, "hybrid")
        rc = _optimize(S)
        if (rc) {
            printf("{err}warning: direct maximisation failed (%s); EM estimates kept\n",
                optimize_result_errortext(S))
        }
        else {
            b    = optimize_result_params(S)
            conv = optimize_result_converged(S)
        }
    }

    // ---------------------------------------------------------- final
    b  = hmdbp_sort(M, b, sortby)
    ll = hmdbp_estep(M, b, W, Zs, PiV, PiU, Pv, mm, 1, 0, XiI)
    if (ll >= .) {
        errprintf("log likelihood not finite at the final estimates\n")
        exit(430)
    }
    // a lambda on a bound of its admissible interval is held fixed in the
    // variance estimate: its row and column of V are zero
    bnd = hmdbp_lambnd(mm)
    wl  = bnd[2] - bnd[1]
    atb = (b[M.m] >= bnd[2] - 1e-4*wl) | (b[M.m] <= bnd[1] + 1e-4*wl)
    frv = J(1, M.npar, 1)
    frv[M.m] = (atb == 0)
    fr  = selectindex(frv)
    if (dolog) printf("{txt}\nComputing the observed information matrix...\n")
    displayflush()
    H  = hmdbp_hess(M, b, fr)
    Vo = J(M.npar, M.npar, 0)
    Vo[fr, fr] = invsym(-H)
    rank = sum(diagonal(Vo) :!= 0)
    if (vce == "robust") {
        Sc = hmdbp_scores(M, b)
        V  = Vo*((M.n/(M.n-1))*cross(Sc, Sc))*Vo
    }
    else V = Vo
    V = (V + V')/2

    // ---------------------------------------------------------- summaries
    omt = J(M.Tmax, M.k2, .)
    for (t=1; t<=M.Tmax; t++) omt[t,.] = mean(W[*M.R[t], .]*M.Ev)
    hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    // intercepts of the states (k2 x 2 or 3) and shifts of the classes (k1 x 2
    // or 3, zero for class 1), by equation
    al = (t1[|1 \ M.k2|]', t2[|1 \ M.k2|]')
    if (M.hasD) al = (al, td[|1 \ M.k2|]')
    de = J(M.k1, cols(al), 0)
    if (M.k1 > 1) {
        if (M.cm) {
            de[|2,1 \ M.k1,1|] = t1[|M.k2+1 \ M.nc1|]'
            de[|2,2 \ M.k1,2|] = t2[|M.k2+1 \ M.nc1|]'
        }
        if (M.cdd) de[|2,3 \ M.k1,3|] = td[|M.k2+1 \ M.ncd|]'
    }

    st_matrix("__hmdbp_b",   b)
    st_matrix("__hmdbp_V",   V)
    // transition matrix: with transitions by transition row, the mean over
    // the transition rows, and (Pit) the mean by period of arrival, one row
    // per period with the k2*k2 probabilities by row of the matrix
    if (M.het) {
        Pm  = colshape(mean(Pv), M.k2)
        Pit = J(rows(u), M.k2*M.k2, .)
        for (t=1; t<=rows(u); t++) {
            Pit[t,.] = mean(Pv[selectindex(M.tt :== u[t]), .])
        }
        st_matrix("__hmdbp_Pit", Pit)
        st_matrix("__hmdbp_tvals", u')
    }
    else Pm = Pv
    st_matrix("__hmdbp_Pi",  Pm)
    st_matrix("__hmdbp_om",  mean(panelsum(W*M.Ev, M.info) :/ M.Ts))
    st_matrix("__hmdbp_omt", omt)
    st_matrix("__hmdbp_al",  al)
    st_matrix("__hmdbp_de",  de)
    st_matrix("__hmdbp_pi0", mean(PiV))
    st_matrix("__hmdbp_pcl", mean(PiU))
    st_matrix("__hmdbp_wcl", mean(W[M.info[.,1], .]*M.Euf))
    st_matrix("__hmdbp_lrng", bnd')
    st_numscalar("__hmdbp_n",     M.n)
    st_numscalar("__hmdbp_ll",    ll)
    st_numscalar("__hmdbp_npar",  M.npar)
    st_numscalar("__hmdbp_rank",  rank)
    st_numscalar("__hmdbp_conv",  conv)
    st_numscalar("__hmdbp_emit",  emits)
    st_numscalar("__hmdbp_gmin",  min(M.Ts))
    st_numscalar("__hmdbp_gavg",  mean(M.Ts))
    st_numscalar("__hmdbp_gmax",  M.Tmax)
    st_numscalar("__hmdbp_lamlo", bnd[1])
    st_numscalar("__hmdbp_lamhi", bnd[2])
    st_numscalar("__hmdbp_nbnd",  atb)
    st_numscalar("__hmdbp_nrisk", sum(M.risk))
}

// ------------------------------------------------------------ fv helpers
string rowvector hmdbp_split(string scalar t)
{
    string rowvector P
    string scalar    rest
    real scalar      pos

    P    = J(1, 0, "")
    rest = t
    pos  = strpos(rest, "#")
    while (pos > 0) {
        P    = (P, substr(rest, 1, pos-1))
        rest = substr(rest, pos+1, .)
        pos  = strpos(rest, "#")
    }
    return((P, rest))
}

// Stata expression of a term (see TermVar); "" for operators other than
// level indicators and c.
string scalar hmdbp_termexpr(string scalar t)
{
    string rowvector P
    string scalar    ex, p, pre, nm, f
    real scalar      i

    P  = hmdbp_split(t)
    ex = ""
    for (i=1; i<=cols(P); i++) {
        p = P[i]
        if (regexm(p, "^(.*)\.([^.]+)$")) {
            pre = regexs(1)
            nm  = regexs(2)
            if (pre == "c") f = nm
            else if (regexm(pre, "^[0-9]+$")) f = "cond(missing(" + nm + "),.," + nm + "==" + pre + ")"
            else return("")
        }
        else f = p
        ex = ex + (i > 1 ? "*" : "") + "(" + f + ")"
    }
    return(ex)
}


// ===========================================================================
// Postestimation
// ===========================================================================

// column of the largest element of each row
real colvector hmdbp_argmax(real matrix X)
{
    real colvector best, bv, upd
    real scalar    c

    best = J(rows(X), 1, 1)
    bv   = X[.,1]
    for (c=2; c<=cols(X); c++) {
        upd  = X[.,c] :> bv
        best = upd:*c + (1 :- upd):*best
        bv   = upd:*X[.,c] + (1 :- upd):*bv
    }
    return(best)
}

// Probabilities of the expanded states by observation: filtered (ptype = 1:
// given the outcomes up to the period), one-step-ahead (2: given the outcomes
// up to the period before) or prior (3: given no outcome). The outcomes of a
// period are the two counts and, in the periods at risk, the dropout.
real matrix hmdbp_sprob(struct hmdbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix P, real scalar ptype)
{
    real matrix    A, F, Q, a, ap
    real colvector ll, r, rp
    real scalar    t, nt, K

    K = M.K
    if (ptype != 3) {
        if (M.het) ll = hmdbp_fwdh(M, lf, Pi0, P, A, F)
        else       ll = hmdbp_fwd(M, lf, Pi0, P, A, F)
        if (ptype == 1) return(A)
    }
    Q = J(M.N, K, 0)
    r = *M.R[1]
    Q[r,.] = Pi0
    a = Pi0
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        // probabilities of the period before: filtered or prior
        if (ptype == 2) {
            rp = (*M.R[t-1])[|1 \ nt|]
            ap = A[rp, .]
        }
        else ap = a[|1,1 \ nt,K|]
        if (M.het) a = hmdbp_proph(M, ap, P, *M.RT[t])
        else       a = ap*P
        Q[r,.] = a
    }
    return(Q)
}

// Most likely sequence of expanded states of each individual (Viterbi
// algorithm, in logs). Dl holds the log probability of the best sequence that
// ends in each state, Psi the state of the period before along that sequence.
real colvector hmdbp_viterbi(struct hmdbp_dat scalar M, real matrix lf,
    real matrix Pi0, real matrix P)
{
    real matrix    lnP, Dl, Psi, d, dn, pb, cand
    real colvector r, r1, rt, S, cur, nxt, bv, arg, upd
    real scalar    t, nt, n1, K, c, e, ue, ve, s, k2

    K   = M.K
    k2  = M.k2
    Dl  = J(M.N, K, .)
    Psi = J(M.N, K, 1)
    // log transition probabilities; an impossible transition (from one class
    // to another) gets -1e300, so that it is never chosen
    lnP = ln(P :+ (P :<= 0)) :- 1e300*(P :<= 0)
    r   = *M.R[1]
    d   = ln(Pi0 :+ 1e-300) + lf[r,.]
    Dl[r,.] = d
    for (t=2; t<=M.Tmax; t++) {
        nt = M.ns[t]
        r  = *M.R[t]
        d  = d[|1,1 \ nt,K|]
        dn = J(nt, K, .)
        pb = J(nt, K, 1)
        if (M.het) rt = *M.RT[t]
        for (e=1; e<=K; e++) {
            // best state of the period before, for arrival state e
            if (M.het) {
                // only the states of the same class: P(state v -> state ve)
                ue   = floor((e-1)/k2) + 1
                ve   = e - (ue-1)*k2
                cand = J(nt, K, -1e300)
                for (s=(ue-1)*k2+1; s<=ue*k2; s++) {
                    cand[., s] = d[., s] + lnP[rt, (s-(ue-1)*k2-1)*k2 + ve]
                }
            }
            else cand = d :+ lnP[., e]'
            // running maximum by selection (an update by differences,
            // bv + upd*(cand - bv), loses the value when bv is -1e300)
            bv   = cand[.,1]
            arg  = J(nt, 1, 1)
            for (c=2; c<=K; c++) {
                upd = cand[.,c] :> bv
                arg = upd:*c + (1 :- upd):*arg
                bv  = upd:*cand[.,c] + (1 :- upd):*bv
            }
            dn[.,e] = bv
            pb[.,e] = arg
        }
        d = dn + lf[r,.]
        Dl[r,.]  = d
        Psi[r,.] = pb
    }
    // backtracking: individuals are sorted by decreasing T, so those observed
    // at t+1 are the first ns[t+1] of those observed at t
    S = J(M.N, 1, .)
    for (t=M.Tmax; t>=1; t--) {
        r   = *M.R[t]
        cur = hmdbp_argmax(Dl[r,.])
        if (t < M.Tmax) {
            n1 = M.ns[t+1]
            r1 = *M.R[t+1]
            cur[|1 \ n1|] = rowsum(Psi[r1, .] :* (J(1, K, nxt) :== (1..K)))
        }
        S[r] = cur
        nxt  = cur
    }
    return(S)
}

// Weights of the expanded states (N x K): posterior (ptype = 0) or as in
// hmdbp_sprob; ist > 0 and icl > 0 condition on state ist and on class icl
// (the other latent variable keeps its weights, renormalised)
real matrix hmdbp_pr(struct hmdbp_dat scalar M, real matrix W,
    real matrix lf, real matrix Pi0, real matrix P, real scalar ptype,
    real scalar ist, real scalar icl)
{
    real matrix    Pr
    real rowvector msk
    real colvector s

    if (ptype == 0) Pr = W
    else            Pr = hmdbp_sprob(M, lf, Pi0, P, ptype)
    if (ist > 0 | icl > 0) {
        msk = J(1, M.K, 1)
        if (ist > 0) msk = msk :* (M.Ev[., ist]')
        if (icl > 0) msk = msk :* (M.Euf[., icl]')
        Pr = Pr :* msk
        s  = rowsum(Pr)
        // rows where the conditioning event has probability zero: equal weights
        Pr = Pr + (s :<= 1e-300)*msk
        Pr = Pr :/ rowsum(Pr)
    }
    return(Pr)
}

// E-step at b for the postestimation tasks, with the quantities they need
real scalar hmdbp_pestep(struct hmdbp_dat scalar M, real rowvector b,
    real matrix W, struct hmdbp_mm scalar mm, real matrix lf,
    real matrix Pi0, real matrix P, real matrix PiV, real matrix PiU,
    real matrix Pv)
{
    real matrix Zs, XiI
    real scalar ll

    ll = hmdbp_estep(M, b, W, Zs, PiV, PiU, Pv, mm, 1, 0, XiI)
    if (ll >= .) return(.)
    lf  = mm.lm + mm.ld
    Pi0 = (PiU*M.Euf') :* (PiV*M.Ev')
    // with transitions by transition row, P holds them (rows(M.W) x k2*k2)
    if (M.het) P = Pv
    else       P = I(M.k1) # Pv
    return(ll)
}

// Predictions for the estimation sample, one row per observation.
// stat: mu1, mu2, pdrop, postprob, classprob, state, class, viterbi, rho, pr.
// pr: probability of the pair (ya, yb) at the regressors and lags of each
// observation.
real matrix hmdbp_predict(struct hmdbp_dat scalar M, real rowvector b,
    string scalar stat, real scalar ptype, real scalar ist, real scalar icl,
    real scalar ya, real scalar yb)
{
    struct hmdbp_mm scalar mm
    real matrix    W, lf, Pi0, P, PiV, PiU, Pv, Pr, f
    real colvector S
    real scalar    ll, lam, q

    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    if (reldif(ll, st_numscalar("e(ll)")) > 1e-8) {
        printf("{txt}note: the log likelihood at e(b) with the data in memory, {res}%12.0g{txt}, differs from e(ll);\n", ll)
        printf("{txt}      the data are not those used in the estimation\n")
    }
    if (stat == "viterbi") {
        S = hmdbp_viterbi(M, lf, Pi0, P)
        return(mod(S :- 1, M.k2) :+ 1)
    }
    Pr = hmdbp_pr(M, W, lf, Pi0, P, ptype, ist, icl)
    if (stat == "postprob")  return(Pr*M.Ev)
    if (stat == "classprob") return(Pr*M.Euf)
    if (stat == "state")     return(hmdbp_argmax(Pr*M.Ev))
    if (stat == "class")     return(hmdbp_argmax(Pr*M.Euf))
    if (stat == "mu1")       return(rowsum(Pr:*mm.Mu1))
    if (stat == "mu2")       return(rowsum(Pr:*mm.Mu2))
    if (stat == "pdrop")     return(rowsum(Pr:*mm.Pd))
    lam = b[M.m]
    q   = 1 - exp(-1)
    if (stat == "rho") return(rowsum(Pr:*(lam*q^2*sqrt(mm.Mu1:*mm.Mu2):*mm.D1:*mm.D2)))
    // pr(a,b)
    f = exp(ya*ln(mm.Mu1) - mm.Mu1 :- lngamma(ya+1) + yb*ln(mm.Mu2) - mm.Mu2 :- lngamma(yb+1))
    f = f :* (1 :+ lam*(exp(-ya) :- mm.D1):*(exp(-yb) :- mm.D2))
    f = f:*(f :> 0)
    return(rowsum(Pr:*f))
}

// ------------------------------------------------------------- estat lrtrans
// Individual derivatives of ln L_i with respect to the transition
// probabilities of the Markov states Pv(c,d), at Pv = I (n x k2*k2, columns
// (c-1)*k2+d). In the expanded chain P = I(k1) # Pv, so the derivative is the
// sum over the classes u of the derivatives with respect to the expanded
// probabilities P((u-1)k2+c, (u-1)k2+d), each the sum over t of
// A_t-1(s) F_t(s') B_t(s') / norm_t, from the scaled forward and backward
// variables of the model without transitions
real matrix hmdbp_dpsc(struct hmdbp_dat scalar M, real matrix lf,
    real matrix Pi0)
{
    real matrix    P, A, F, Dn, B, G, Ap, bb, Dv
    real colvector ll, r, rp, nrm
    real scalar    t, nt, K, k2, c, u, s0

    K  = M.K
    k2 = M.k2
    P  = I(K)
    ll = hmdbp_fwd(M, lf, Pi0, P, A, F)
    Dv = J(M.n, k2*k2, 0)
    if (M.Tmax < 2) return(Dv)
    Dn = J(M.n, K*K, 0)
    B  = J(M.ns[M.Tmax], K, 1)
    for (t=M.Tmax; t>=2; t--) {
        nt  = M.ns[t]
        r   = *M.R[t]
        rp  = (*M.R[t-1])[|1 \ nt|]
        G   = F[r,.] :* B
        Ap  = A[rp,.]
        nrm = rowsum((Ap*P) :* G)
        Ap  = Ap :/ nrm
        for (c=1; c<=K; c++) {
            Dn[|1,(c-1)*K+1 \ nt,c*K|] = Dn[|1,(c-1)*K+1 \ nt,c*K|] + Ap[.,c] :* G
        }
        bb = G*P'
        bb = bb :/ rowsum(bb)
        B  = J(M.ns[t-1], K, 1)
        B[|1,1 \ nt,K|] = bb
    }
    for (u=1; u<=M.k1; u++) {
        s0 = (u-1)*k2
        for (c=1; c<=k2; c++) {
            r = ((s0+c-1)*K+s0+1)::((s0+c-1)*K+s0+k2)
            Dv[|1,(c-1)*k2+1 \ .,c*k2|] = Dv[|1,(c-1)*k2+1 \ .,c*k2|] + Dn[., r]
        }
    }
    return(Dv)
}

// min 0.5 th'Q th - c'th subject to th >= 0 (Lawson-Hanson active set on the
// normal equations); Q positive definite
real colvector hmdbp_nnqp(real matrix Q, real colvector c)
{
    real colvector th, s, w, Ps, cand, idx, neg
    real scalar    r, it, inner, j, alpha

    r  = rows(c)
    th = J(r, 1, 0)
    Ps = J(r, 1, 0)
    w  = c
    for (it=1; it<=3*r+10; it++) {
        cand = selectindex((Ps :== 0) :& (w :> 1e-10))
        if (rows(cand) == 0) break
        j = cand[hmdbp_argmax(w[cand]')]
        Ps[j] = 1
        for (inner=1; inner<=3*r+10; inner++) {
            idx = selectindex(Ps)
            s   = J(r, 1, 0)
            s[idx] = lusolve(Q[idx, idx], c[idx])
            if (all(s[idx] :> 1e-12)) {
                th = s
                break
            }
            neg   = idx[selectindex(s[idx] :<= 1e-12)]
            alpha = min(th[neg] :/ (th[neg] - s[neg]))
            if (alpha >= .) alpha = 0
            th = th + alpha*(s - th)
            Ps = Ps :& (th :> 1e-12)
            th = th :* Ps
        }
        w = c - Q*th
    }
    return(th)
}

// weights of the chi-bar-squared distribution of the LR test of no
// transitions between the Markov states (j = 0, 1, ..., k2(k2-1) positive
// components), at the estimates b of the model without transitions. The parameters are the off-diagonal
// probabilities of Pv, with Pv(c,c) = 1 - the sum of the others; their
// information net of the other parameters (measurement, dropout, initial
// probabilities, classes) is estimated by the outer product of the individual
// scores. The state of the random-number generator is restored at the end.
real rowvector hmdbp_lrw(struct hmdbp_dat scalar M, real rowvector b,
    real scalar reps, real scalar seed)
{
    struct hmdbp_mm scalar mm
    real matrix    W, lf, Pi0, P, PiV, PiU, Pv, Sc, Dv, Sp, S, Jm, Q, L, iJ
    real colvector cv, th
    real rowvector w
    real scalar    ll, k2, nr, c, e, j, i, np
    string scalar  rs

    k2 = M.k2
    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) {
        errprintf("the log likelihood of the model without transitions cannot be evaluated\n")
        exit(459)
    }
    Sc = hmdbp_scores(M, b)
    Dv = hmdbp_dpsc(M, lf, Pi0)
    nr = k2*(k2-1)
    Sp = J(M.n, nr, .)
    j  = 0
    for (c=1; c<=k2; c++) {
        for (e=1; e<=k2; e++) {
            if (e != c) {
                j++
                Sp[.,j] = Dv[.,(c-1)*k2+e] - Dv[.,(c-1)*k2+c]
            }
        }
    }
    S  = (Sc, Sp)
    Jm = quadcross(S, S)
    np = cols(Sc)
    iJ = invsym(Jm[|1,1 \ np,np|])
    Q  = Jm[|np+1,np+1 \ .,.|] - Jm[|np+1,1 \ .,np|]*iJ*Jm[|1,np+1 \ np,.|]
    Q  = (Q + Q')/2
    L  = cholesky(Q)
    if (hasmissing(L)) {
        errprintf("the information of the transition probabilities is not positive definite\n")
        exit(498)
    }
    rs = rngstate()
    rseed(seed)
    w  = J(1, nr+1, 0)
    for (i=1; i<=reps; i++) {
        cv = L*rnormal(nr, 1, 0, 1)
        th = hmdbp_nnqp(Q, cv)
        j  = sum(th :> 1e-12)
        w[j+1] = w[j+1] + 1
    }
    rngstate(rs)
    // No score (Kuhn-Tucker) statistic: the individual scores of the
    // transition probabilities at Pv = I are sums of likelihood ratios
    // f(y | state d)/f(y | state c), whose distribution has very heavy tails;
    // their outer product in one sample is far from the information, and a
    // statistic s' inv(OPG) s is badly calibrated (hmdbivpois_sim_04: rejection
    // rate 0.01 at the 10% level; unprojected, at the true values, mean 44
    // instead of 6). The weights depend only on the shape of Q and the LR test
    // keeps its size.
    return(w/reps)
}

// relative entropy of the classification: states by period, classes by
// individual (1 = known with certainty)
real rowvector hmdbp_entropy(struct hmdbp_dat scalar M, real rowvector b)
{
    struct hmdbp_mm scalar mm
    real matrix W, lf, Pi0, P, PiV, PiU, Pv, Wv, Wc
    real scalar ll, es, ec

    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    es = 1
    ec = 1
    if (M.k2 > 1) {
        Wv = W*M.Ev
        es = 1 - quadsum(-Wv:*ln(Wv :+ (Wv :<= 0)))/(M.N*ln(M.k2))
    }
    if (M.k1 > 1) {
        Wc = W[M.info[.,1], .]*M.Euf
        ec = 1 - quadsum(-Wc:*ln(Wc :+ (Wc :<= 0)))/(M.n*ln(M.k1))
    }
    return((es, ec))
}

// transition probabilities (by rows: element (c-1)*k+d = P(state d | state c)),
// mean initial probabilities of the states, mean prior probabilities of the
// classes
real rowvector hmdbp_trf(struct hmdbp_dat scalar M, real rowvector b,
    real colvector sel)
{
    real rowvector b1, t1, b2, t2, bd, td, tau, Pm
    real matrix    Dv, Gu, Pv
    real scalar    lam

    hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    // transitions by transition row: mean over the rows sel
    if (M.het) {
        Pv = hmdbp_PT(M, tau)
        Pm = mean(Pv[sel, .])
    }
    else if (M.ntau > 0) Pm = vec(hmdbp_P(tau, M.k2)')'
    else                Pm = vec(I(M.k2))'
    return((Pm, mean(hmdbp_mlogit(M.Z, Dv, M.k2)), mean(hmdbp_mlogit(M.ZU, Gu, M.k1))))
}

// the same quantities (row 1) and their delta-method standard errors (row 2),
// from e(V), with the Jacobian by central differences
real matrix hmdbp_trse(struct hmdbp_dat scalar M, real rowvector b)
{
    real rowvector f, fp, fq, bp, bq, idx
    real colvector sel
    real matrix    G, V
    real scalar    j, i, h, per

    // transitions by transition row: period of arrival per (local of the
    // caller), all the transition rows when empty
    sel = J(0, 1, .)
    if (M.het) {
        per = strtoreal(st_local("per"))
        if (per < .) sel = selectindex(M.tt :== per)
        else         sel = (1::rows(M.tt))
        if (rows(sel) == 0) {
            errprintf("period(): no transition into period %g in the estimation sample\n", per)
            exit(498)
        }
    }
    f = hmdbp_trf(M, b, sel)
    if (M.m + M.md >= M.npar) return(f \ J(1, cols(f), 0))
    idx = (M.m+M.md+1)..M.npar
    G   = J(cols(f), cols(idx), 0)
    for (i=1; i<=cols(idx); i++) {
        j  = idx[i]
        h  = 1e-6*(abs(b[j]) + 1)
        bp = b; bp[j] = b[j] + h
        bq = b; bq[j] = b[j] - h
        fp = hmdbp_trf(M, bp, sel)
        fq = hmdbp_trf(M, bq, sel)
        G[., i] = (fp - fq)' :/ (2*h)
    }
    V = st_matrix("e(V)")
    V = V[idx, idx]
    return(f \ sqrt(diagonal(G*V*G'))')
}

// --------------------------------------------------------------- estat ame
// variables of the terms (lags excluded): names (local amev), 1 if they appear
// as factor variables (amef), their levels (amelv#)
void hmdbp_amevars(string scalar terms, string scalar lagn)
{
    string rowvector T, L, P, vars, lev
    real rowvector   fac, w
    string scalar    p, pre, nm, lv
    real scalar      i, j, wi, isf

    T    = tokens(terms)
    L    = tokens(lagn)
    vars = J(1, 0, "")
    lev  = J(1, 0, "")
    fac  = J(1, 0, .)
    for (i=1; i<=cols(T); i++) {
        if (anyof(L, T[i])) continue
        P = hmdbp_split(T[i])
        for (j=1; j<=cols(P); j++) {
            p   = P[j]
            pre = ""
            nm  = p
            if (regexm(p, "^(.*)\.([^.]+)$")) {
                pre = regexs(1)
                nm  = regexs(2)
            }
            isf = regexm(pre, "^[0-9]+")
            lv  = (isf ? regexs(0) : "")
            w   = selectindex(vars :== nm)
            if (cols(w) == 0) {
                vars = (vars, nm)
                fac  = (fac, isf)
                lev  = (lev, "")
                wi   = cols(vars)
            }
            else wi = w[1]
            if (isf) {
                fac[wi] = 1
                lev[wi] = strtrim(lev[wi] + " " + lv)
            }
        }
    }
    st_local("amev", invtokens(vars))
    st_local("amef", invtokens(strofreal(fac)))
    for (i=1; i<=cols(vars); i++) {
        st_local("amelv" + strofreal(i), invtokens(uniqrows(tokens(lev[i])')'))
    }
}

// a term with the variable v replaced by r (t is a term as in e(common1):
// var, c.var, #.var, interactions joined by #)
string scalar hmdbp_subst(string scalar t, string scalar v, string scalar r)
{
    string rowvector P
    string scalar    out, p
    real scalar      i

    P   = hmdbp_split(t)
    out = ""
    for (i=1; i<=cols(P); i++) {
        p = P[i]
        if (regexm(p, "^(.*)\.([^.]+)$")) {
            if (regexs(2) == v) p = regexs(1) + "." + r
        }
        else if (p == v) p = r
        out = out + (i > 1 ? "#" : "") + p
    }
    return(out)
}

// features of the lagged counts (N x 4): g(y1), g(y2) and, with logzero, the
// indicators of y1 = 0 and y2 = 0 (zero otherwise). lagt and larg: the
// transformation and its argument, as in lagtransform()
real matrix hmdbp_lagf(real colvector z1, real colvector z2,
    string scalar lagt, real scalar larg)
{
    real matrix F

    F = J(rows(z1), 4, 0)
    if (lagt == "linear") F[., 1..2] = (z1, z2)
    else if (lagt == "log1p") F[., 1..2] = ln(1 :+ (z1, z2))
    else if (lagt == "logc") F[., 1..2] = (ln(rowmax((z1, J(rows(z1), 1, larg)))), ln(rowmax((z2, J(rows(z2), 1, larg)))))
    else if (lagt == "indicator") F[., 1..2] = ((z1 :> 0), (z2 :> 0))
    else if (lagt == "topcode") F[., 1..2] = (rowmin((z1, J(rows(z1), 1, larg))), rowmin((z2, J(rows(z2), 1, larg))))
    else if (lagt == "logzero") {
        F[., 1] = (z1 :> 0):*ln(z1 + (z1 :<= 0))
        F[., 2] = (z2 :> 0):*ln(z2 + (z2 :<= 0))
        F[., 3] = (z1 :<= 0)
        F[., 4] = (z2 :<= 0)
    }
    return(F)
}

// Average marginal effects (1 x nq*nE: equation 1, then 2, then dropout): for
// each effect e and equation j, mean over the observations of
// sum_s Pr_s (m_s(A) - m_s(B)) / div, with m = exp(eta) for the counts and
// invlogit(eta) for the dropout. For the effects of the regressors PX[e, j]
// and PX[e, 3+j] point to the designs A and B of the common terms of
// equation j, PX[e, 6+j] and PX[e, 9+j] to those of the state-specific
// terms; for one more lagged count (or count of the period) PX[e, 1] points
// to the four features recomputed, which replace those of the design
// (codes of the terms: 1-4 lag features, 11-14 current counts).
real rowvector hmdbp_amef(struct hmdbp_dat scalar M, real rowvector b,
    pointer(real matrix) matrix PX, string rowvector typ, real rowvector div,
    real scalar ptype, real scalar ist, real scalar icl)
{
    struct hmdbp_mm scalar mm
    real rowvector b1, t1, b2, t2, bd, td, tau, bb, tt, cc, sc
    real matrix    W, lf, Pi0, P, PiV, PiU, Pv, Pr, Dv, Gu, out, X, XS, C
    real matrix    XA, XB, XSA, XSB, eA, eB, dm
    real scalar    ll, lam, e, j, nq, nE, nc, ns, base

    nE = cols(div)
    nq = 2 + M.hasD
    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) return(J(1, nq*nE, .))
    Pr = hmdbp_pr(M, W, lf, Pi0, P, ptype, ist, icl)
    hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    out = J(nq, nE, .)
    for (j=1; j<=nq; j++) {
        if (j == 1) {
            X = M.X1; XS = M.XS1; bb = b1; tt = t1; C = M.C1; nc = M.nc1
            ns = M.s1; cc = M.cc1; sc = M.sc1
        }
        else if (j == 2) {
            X = M.X2; XS = M.XS2; bb = b2; tt = t2; C = M.C1; nc = M.nc1
            ns = M.s2; cc = M.cc2; sc = M.sc2
        }
        else {
            X = M.Xd; XS = M.XSd; bb = bd; tt = td; C = M.Cd; nc = M.ncd
            ns = M.sd; cc = M.ccd; sc = M.scd
        }
        eB = hmdbp_eta(M, X, XS, bb, tt, C, nc, ns)
        for (e=1; e<=nE; e++) {
            if (typ[e] == "x") {
                XA  = *PX[e, j]
                XB  = *PX[e, 3+j]
                XSA = *PX[e, 6+j]
                XSB = *PX[e, 9+j]
                eA  = hmdbp_eta(M, XA, XSA, bb, tt, C, nc, ns)
                if (j < 3) dm = exp(eA) - exp(hmdbp_eta(M, XB, XSB, bb, tt, C, nc, ns))
                else       dm = invlogit(eA) - invlogit(hmdbp_eta(M, XB, XSB, bb, tt, C, nc, ns))
            }
            else {
                // one more count: its features recomputed wherever they enter
                base = (substr(typ[e], 1, 3) == "lag" ? 0 : 10)
                XA   = hmdbp_frep(X, cc, *PX[e, 1], base)
                XSA  = hmdbp_frep(XS, sc, *PX[e, 1], base)
                eA   = hmdbp_eta(M, XA, XSA, bb, tt, C, nc, ns)
                if (j < 3) dm = exp(eA) - exp(eB)
                else       dm = invlogit(eA) - invlogit(eB)
            }
            out[j, e] = mean(rowsum(Pr:*dm))/div[e]
        }
    }
    return(vec(out')')
}

// design X with the columns of code base + m (m = 1..4) replaced by column m
// of F
real matrix hmdbp_frep(real matrix X, real rowvector cd, real matrix F,
    real scalar base)
{
    real matrix XA
    real scalar i

    XA = X
    for (i=1; i<=cols(cd); i++) {
        if (cd[i] > base & cd[i] <= base + 4) XA[., i] = F[., cd[i] - base]
    }
    return(XA)
}

// estat ame: designs of the effects from the locals of the caller (ame_n,
// ame_t#, ame_div#, ame_[ab][12d]_#), estimates and their delta-method
// variance (Jacobian by central differences). Results: matrix stub (estimates)
// and stub + "V" (variance)
void hmdbp_ame(struct hmdbp_dat scalar M, real rowvector b, string scalar stub,
    real scalar ptype, real scalar ist, real scalar icl, string scalar touse)
{
    pointer(real matrix) matrix PX
    string rowvector typ, eqs
    real rowvector   div, f, fp, fq, bp, bq, fr
    real matrix      G, V, Lr, FA
    real scalar      nE, nq, e, s, j, h, i, nl, larg
    string scalar    sd, lagt

    nE  = strtoreal(st_local("ame_n"))
    nq  = 2 + M.hasD
    eqs = ("1", "2", "d")
    // untransformed lagged counts and transformation, for the effects of one
    // more lagged count; nl lag features (4 with logzero)
    lagt = st_local("lagt")
    larg = strtoreal(st_local("lagarg"))
    Lr   = st_data(., tokens(st_local("rawlags")), touse)
    nl   = (lagt == "logzero" ? 4 : 2)
    div = J(1, nE, .)
    typ = J(1, nE, "")
    PX  = J(nE, 12, NULL)
    for (e=1; e<=nE; e++) {
        typ[e] = st_local("ame_t" + strofreal(e))
        div[e] = strtoreal(st_local("ame_div" + strofreal(e)))
        if (typ[e] == "x") {
            // common terms (ame_a#_e, ame_b#_e), then state-specific terms
            // (ame_as#_e, ame_bs#_e)
            for (s=1; s<=2; s++) {
                sd = (s == 1 ? "a" : "b")
                for (j=1; j<=3; j++) {
                    PX[e, 3*(s-1)+j]   = &(hmdbp_getx(st_local("ame_" + sd + eqs[j] + "_" + strofreal(e)), touse, M.N))
                    PX[e, 6+3*(s-1)+j] = &(hmdbp_getx(st_local("ame_" + sd + "s" + eqs[j] + "_" + strofreal(e)), touse, M.N))
                }
            }
        }
        else {
            if (typ[e] == "lag1")      FA = hmdbp_lagf(Lr[., 1] :+ 1, Lr[., 2], lagt, larg)
            else if (typ[e] == "lag2") FA = hmdbp_lagf(Lr[., 1], Lr[., 2] :+ 1, lagt, larg)
            else if (typ[e] == "cur1") FA = hmdbp_lagf(M.y1 :+ 1, M.y2, lagt, larg)
            else                       FA = hmdbp_lagf(M.y1, M.y2 :+ 1, lagt, larg)
            PX[e, 1] = &(FA :+ 0)
        }
    }
    f = hmdbp_amef(M, b, PX, typ, div, ptype, ist, icl)
    if (hasmissing(f)) {
        errprintf("the marginal effects cannot be computed at e(b) with the data in memory\n")
        exit(459)
    }
    // Jacobian with respect to the parameters with nonzero variance
    V  = st_matrix("e(V)")
    fr = selectindex(diagonal(V)' :!= 0)
    G  = J(cols(f), cols(fr), 0)
    for (i=1; i<=cols(fr); i++) {
        j  = fr[i]
        h  = 1e-4*(abs(b[j]) + 1)
        bp = b; bp[j] = b[j] + h
        bq = b; bq[j] = b[j] - h
        fp = hmdbp_amef(M, bp, PX, typ, div, ptype, ist, icl)
        fq = hmdbp_amef(M, bq, PX, typ, div, ptype, ist, icl)
        G[., i] = (fp - fq)' :/ (2*h)
    }
    V = G*V[fr, fr]*G'
    st_matrix(stub, f)
    st_matrix(stub + "V", (V + V')/2)
}

// keeps the elements idx of the row vector bn and the rows and columns idx of
// the matrix vn (Stata matrices)
void hmdbp_msel(string scalar bn, string scalar vn, string scalar idx)
{
    real rowvector k, bb
    real matrix    V

    k  = strtoreal(tokens(idx))
    bb = st_matrix(bn)
    V  = st_matrix(vn)
    st_matrix(bn, bb[k])
    st_matrix(vn, V[k, k])
}

// ---------------------------------------------------------- estat survival
// one row per group and period of follow-up: group, period, individuals at
// risk, dropouts, Kaplan-Meier survival, survival predicted by the model
// (product of 1 - the mean predicted probability of dropout among those at
// risk)
real matrix hmdbp_surv(struct hmdbp_dat scalar M, real rowvector b,
    real scalar ptype, string scalar gvar, string scalar touse)
{
    struct hmdbp_mm scalar mm
    real matrix    W, lf, Pi0, P, PiV, PiU, Pv, Pr, out
    real colvector h, g, r, sel, rr
    real scalar    ll, gg, ng, t, nr, d, km, ms

    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    Pr = hmdbp_pr(M, W, lf, Pi0, P, ptype, 0, 0)
    h  = rowsum(Pr:*mm.Pd)
    if (gvar != "") g = st_data(., gvar, touse)
    else            g = J(M.N, 1, 1)
    ng  = max(g)
    out = J(0, 6, .)
    for (gg=1; gg<=ng; gg++) {
        km = 1
        ms = 1
        for (t=1; t<=M.Tmax; t++) {
            r   = *M.R[t]
            sel = selectindex((g[r] :== gg) :& (M.risk[r] :== 1))
            // nobody at risk in period t: nobody is observed later
            if (rows(sel) == 0) break
            rr  = r[sel]
            nr  = rows(rr)
            d   = sum(M.D[rr])
            km  = km*(1 - d/nr)
            ms  = ms*(1 - mean(h[rr]))
            out = out \ (gg, t, nr, d, km, ms)
        }
    }
    return(out)
}

// driver of the postestimation tasks: the data are rebuilt by the Stata
// program PostRun, the parameters are e(b)
void hmdbp_post(string scalar task, string scalar y1n, string scalar y2n,
    string scalar x1n, string scalar x2n, string scalar ydn,
    string scalar xdn, string scalar riskn, string scalar zn,
    string scalar zun, string scalar idn, string scalar touse,
    real scalar k2, real scalar k1, real scalar notr, string scalar outn,
    string scalar stat, real scalar ptype, real scalar ist, real scalar icl,
    real scalar ya, real scalar yb, string scalar gvar, real scalar maxc,
    real scalar bycl, string scalar wn, string scalar trn,
    string scalar timen, real scalar het, real scalar tvt, real scalar reps,
    real scalar seed)
{
    struct hmdbp_dat scalar M
    real rowvector b
    real matrix    R
    real colvector u
    real scalar    nwx

    hmdbp_setup(M, y1n, y2n, x1n, x2n, ydn, xdn, riskn, zn, zun, idn,
        touse, k2, k1, notr, wn, trn, timen, het, tvt, u, nwx)
    b = st_matrix("e(b)")
    if (cols(b) != M.npar) {
        errprintf("e(b) has %5.0f elements; the model rebuilt from the data has %5.0f\n", cols(b), M.npar)
        exit(459)
    }
    // the estimates are taken as they are, also with a lambda on a bound
    M.chk = 0
    if (task == "predict") {
        R = hmdbp_predict(M, b, stat, ptype, ist, icl, ya, yb)
        st_store(., tokens(outn), touse, R)
    }
    if (task == "lrtrans") st_matrix(outn, hmdbp_lrw(M, b, reps, seed))
    if (task == "ic")    st_matrix(outn, hmdbp_entropy(M, b))
    if (task == "trans") st_matrix(outn, hmdbp_trse(M, b))
    if (task == "ame")   hmdbp_ame(M, b, outn, ptype, ist, icl, touse)
    if (task == "surv")  st_matrix(outn, hmdbp_surv(M, b, ptype, gvar, touse))
    if (task == "fit")   hmdbp_fit(M, b, outn, ptype, maxc, bycl)
    if (task == "stability") hmdbp_stab(M, b, outn, seed, touse)
}

// ---------------------------------------------------------- estat stability
// p-th quantile of a column (the smallest value with at least a share p of
// the values below or at it)
real scalar hmdbp_pctile(real colvector x, real scalar p)
{
    real colvector s

    s = sort(x, 1)
    return(s[max((1, ceil(p*rows(s))))])
}

// one draw of the pair of counts for each row (Famoye/Sarmanov, one lambda):
// y2 is Poisson; y1 | y2 by inversion of its conditional pmf
void hmdbp_draw(real colvector mu1, real colvector mu2, real scalar lam,
    real colvector y1, real colvector y2)
{
    real colvector u2, d1, uu, cum, todo, hit, pr
    real scalar    q, y

    q    = 1 - exp(-1)
    y2   = rpoisson(1, 1, mu2)
    u2   = exp(-y2) - exp(-q*mu2)
    d1   = exp(-q*mu1)
    uu   = runiform(rows(mu1), 1)
    y1   = J(rows(mu1), 1, .)
    cum  = J(rows(mu1), 1, 0)
    todo = J(rows(mu1), 1, 1)
    pr   = exp(-mu1)
    for (y=0; y<=20000; y++) {
        cum  = cum + pr:*rowmax((1 :+ lam*(exp(-y) :- d1):*u2, J(rows(mu1), 1, 0)))
        hit  = todo :& (cum :>= uu)
        if (any(hit)) y1[selectindex(hit)] = J(sum(hit), 1, y)
        todo = todo :& (hit :== 0)
        if (!any(todo)) break
        pr   = pr:*mu1:/(y+1)
    }
    if (any(todo)) y1[selectindex(todo)] = J(sum(todo), 1, y)
}

// Pieces of the linear predictor of one equation, by expanded state: G (4 x K)
// and GC (4 x K), the coefficients of the four lag features and of the four
// features of the current counts (common or state-specific, from the codes
// of the terms); NL (N x K), the rest of the linear predictor, intercepts
// included. F and Fc: the features at the observed lags and current counts.
// chk: largest difference between the columns of the designs and the features
// rebuilt from the counts.
void hmdbp_lagparts(struct hmdbp_dat scalar M, real matrix X, real matrix XS,
    real rowvector bb, real rowvector tt, real matrix C, real scalar nc,
    real scalar ns, real rowvector cc, real rowvector sc, real matrix F,
    real matrix Fc, real matrix NL, real matrix G, real matrix GC,
    real scalar chk)
{
    real matrix    S
    real scalar    m, i

    G  = J(4, M.K, 0)
    GC = J(4, M.K, 0)
    if (ns > 0) S = colshape(tt[|nc+1 \ cols(tt)|], ns)    // k2 x ns
    for (i=1; i<=cols(cc); i++) {
        m = cc[i]
        if (m >= 1 & m <= 4) {
            G[m, .] = G[m, .] :+ bb[i]
            chk = max((chk, mreldif(X[., i], F[., m])))
        }
        if (m >= 11 & m <= 14) {
            GC[m-10, .] = GC[m-10, .] :+ bb[i]
            chk = max((chk, mreldif(X[., i], Fc[., m-10])))
        }
    }
    for (i=1; i<=cols(sc); i++) {
        m = sc[i]
        if (m >= 1 & m <= 4) {
            G[m, .] = G[m, .] + (M.Ev*S[., i])'
            chk = max((chk, mreldif(XS[., i], F[., m])))
        }
        if (m >= 11 & m <= 14) {
            GC[m-10, .] = GC[m-10, .] + (M.Ev*S[., i])'
            chk = max((chk, mreldif(XS[., i], Fc[., m-10])))
        }
    }
    NL = hmdbp_eta(M, X, XS, bb, tt, C, nc, ns) - F*G - Fc*GC
}

// Stability of the dynamics of the counts. The coefficients of the lags can be
// common or state-specific (lagvary()); the rest of the linear predictor
// varies by latent type (state x class). Results in Stata matrices whose
// names start with stub: stubG (K x 2nl or 3nl, coefficients of the nl lag
// features in the equations of y1, y2 and, when the lags enter it, of the
// dropout, by type), stubR (K x 2, spectral radius and norm of the 2 x 2
// matrix of the coefficients of g(y1), g(y2) in the count equations), stubU
// (K x 2, index gamma_jj exp(q99 of the rest of the linear predictor)), stubD
// (K x 1, share of deterministic paths that exceed the limit), stubS
// (simulated paths: share that exceeds the limit, paths, horizon, limit, mean
// probability of still being in the panel at the horizon), stubQ (quantiles
// of the expected counts at the horizon and in the sample). The simulated
// paths ignore the dropout when they draw the counts (paths of those who
// stay); the probability of staying is computed along them. Options from the
// locals of the caller: horizon, limit, nsim, lagt, lagarg, rawlags.
void hmdbp_stab(struct hmdbp_dat scalar M, real rowvector b, string scalar stub,
    real scalar seed, string scalar touse)
{
    struct hmdbp_mm scalar mm
    real rowvector b1, t1, b2, t2, bd, td, tau, qq, sv1, sv2
    real matrix    W, lf, Pi0, P, PiV, PiU, Pv, Dv, Gu, Lr, F, Fc, NL1, NL2, NLd
    real matrix    G1, G2, Gd, C1, C2, Cd, Gm, Q, U, R, Ri, Uk, PTi, F0, sv
    real colvector z1, z2, m1, m2, dv, ok, last, s, ex, mx1, mx2
    real colvector y1, y2, ind, fit1, fit2, rowof, cl, vv, e, surv, pdr
    real scalar    K, k2, H, lim, nsim, larg, nl, ll, c, h, i, tr, dt, chk, lam
    string scalar  lagt, rs

    K    = M.K
    k2   = M.k2
    H    = strtoreal(st_local("horizon"))
    lim  = strtoreal(st_local("limit"))
    nsim = strtoreal(st_local("nsim"))
    lagt = st_local("lagt")
    larg = strtoreal(st_local("lagarg"))
    Lr   = st_data(., tokens(st_local("rawlags")), touse)
    nl   = (lagt == "logzero" ? 4 : 2)

    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    hmdbp_unpack(M, b, b1, t1, b2, t2, lam, bd, td, Dv, tau, Gu)
    F   = hmdbp_lagf(Lr[., 1], Lr[., 2], lagt, larg)
    Fc  = hmdbp_lagf(M.y1, M.y2, lagt, larg)
    chk = 0
    hmdbp_lagparts(M, M.X1, M.XS1, b1, t1, M.C1, M.nc1, M.s1, M.cc1, M.sc1, F, Fc, NL1, G1, C1, chk)
    hmdbp_lagparts(M, M.X2, M.XS2, b2, t2, M.C1, M.nc1, M.s2, M.cc2, M.sc2, F, Fc, NL2, G2, C2, chk)
    if (M.hasD) {
        hmdbp_lagparts(M, M.Xd, M.XSd, bd, td, M.Cd, M.ncd, M.sd, M.ccd, M.scd, F, Fc, NLd, Gd, Cd, chk)
    }
    if (chk > 1e-8) {
        errprintf("hmdbp_stab(): the lag features are not rebuilt (%g)\n", chk)
        exit(498)
    }

    // coefficients of the lag features by type (the dropout only when the
    // lags enter its equation)
    Gm = (G1[|1,1 \ nl,.|]', G2[|1,1 \ nl,.|]')
    if (M.hasD) {
        if (any((M.ccd :>= 1) :& (M.ccd :<= 4)) | any((M.scd :>= 1) :& (M.scd :<= 4))) {
            Gm = (Gm, Gd[|1,1 \ nl,.|]')
        }
    }
    st_matrix(stub + "G", Gm)

    // spectral radius and norm of [g11 g12 \ g21 g22], by type
    R = J(K, 2, .)
    for (c=1; c<=K; c++) {
        Q  = (G1[1,c], G1[2,c] \ G2[1,c], G2[2,c])
        tr = Q[1,1] + Q[2,2]
        dt = Q[1,1]*Q[2,2] - Q[1,2]*Q[2,1]
        if (tr^2/4 - dt >= 0) R[c,1] = max(abs((tr/2 + sqrt(tr^2/4 - dt), tr/2 - sqrt(tr^2/4 - dt))))
        else                  R[c,1] = sqrt(dt)
        R[c,2] = sqrt(max(symeigenvalues(Q'*Q)))
    }
    st_matrix(stub + "R", R)

    // index of the univariate map y -> exp(a + g y): explosive above exp(-1)
    U = J(K, 2, .)
    for (c=1; c<=K; c++) {
        U[c, 1] = G1[1,c]*exp(hmdbp_pctile(NL1[., c], 0.99))
        U[c, 2] = G2[2,c]*exp(hmdbp_pctile(NL2[., c], 0.99))
    }
    st_matrix(stub + "U", U)

    // deterministic paths: in each type, from the lags and with the
    // covariates of each observation, the counts are replaced by their
    // expected values for H periods
    dv = J(K, 1, .)
    for (c=1; c<=K; c++) {
        z1 = Lr[., 1]
        z2 = Lr[., 2]
        e  = J(M.N, 1, 0)
        for (h=1; h<=H; h++) {
            F0 = hmdbp_lagf(z1, z2, lagt, larg)
            m1 = exp(NL1[., c] + F0*G1[., c])
            m2 = exp(NL2[., c] + F0*G2[., c])
            ok = (m1 :<= lim) :& (m2 :<= lim)
            e  = e :| (ok :== 0)
            z1 = ok:*m1 + (1:-ok)*lim
            z2 = ok:*m2 + (1:-ok)*lim
        }
        dv[c] = mean(e)
    }
    st_matrix(stub + "D", dv)

    // simulated paths: nsim per individual, from the last observed period,
    // with its covariates, the type drawn from its posterior probabilities;
    // the class stays, the state moves with the estimated chain (with
    // transitions by row, those of the last transition of the individual,
    // or their mean for those observed once)
    rs   = rngstate()
    rseed(seed)
    last = M.info[., 2]
    ind  = J(nsim, 1, 1) # (1::M.n)
    last = last[ind]
    Uk   = uppertriangle(J(K, K, 1))
    s    = 1 :+ rowsum(runiform(rows(ind), 1) :> (W[last, .]*Uk))
    s    = rowmin((s, J(rows(ind), 1, K)))
    if (M.het) {
        rowof = J(M.N, 1, 0)
        rowof[M.trow] = (1::rows(M.trow))
        PTi = J(rows(ind), k2*k2, .)
        for (i=1; i<=rows(ind); i++) {
            if (rowof[last[i]] > 0) PTi[i, .] = P[rowof[last[i]], .]
            else                    PTi[i, .] = mean(P)
        }
    }
    z1   = M.y1[last]
    z2   = M.y2[last]
    ex   = J(rows(ind), 1, 0)
    mx1  = J(rows(ind), 1, 0)
    mx2  = J(rows(ind), 1, 0)
    surv = J(rows(ind), 1, 1)
    Uk   = uppertriangle(J(k2, k2, 1))
    for (h=1; h<=H; h++) {
        F0 = hmdbp_lagf(z1, z2, lagt, larg)
        sv = (J(1, K, s) :== (1..K))
        m1 = exp(rowsum(NL1[last, .] :* sv) + rowsum(F0 :* G1[., s]'))
        m2 = exp(rowsum(NL2[last, .] :* sv) + rowsum(F0 :* G2[., s]'))
        ok = (m1 :<= lim) :& (m2 :<= lim)
        ex = ex :| (ok :== 0)
        // paths above the limit are not drawn (mean 1 in place of theirs)
        hmdbp_draw(ok:*m1 + (1:-ok), ok:*m2 + (1:-ok), lam, y1, y2)
        // probability of dropping out after this period, from the lags and
        // from the counts drawn (when they enter the dropout equation)
        if (M.hasD) {
            pdr = rowsum(NLd[last, .] :* sv) + rowsum(F0 :* Gd[., s]')
            F0  = hmdbp_lagf(ok:*y1 + (1:-ok)*lim, ok:*y2 + (1:-ok)*lim, lagt, larg)
            pdr = pdr + rowsum(F0 :* Cd[., s]')
            surv = surv :* (1 :- invlogit(pdr))
        }
        // paths that exceeded the limit stay at the limit
        z1 = ok:*y1 + (1:-ok)*lim
        z2 = ok:*y2 + (1:-ok)*lim
        if (h == H) {
            mx1 = m1
            mx2 = m2
        }
        // next state; the class does not change
        cl = ceil(s/k2)
        vv = s - (cl :- 1)*k2
        if (M.het) {
            Ri = J(rows(ind), k2, 0)
            for (c=1; c<=k2; c++) Ri = Ri + (vv :== c):*PTi[|1,(c-1)*k2+1 \ .,c*k2|]
        }
        else Ri = Pv[vv, .]
        vv = 1 :+ rowsum(runiform(rows(ind), 1) :> (Ri*Uk))
        vv = rowmin((vv, J(rows(ind), 1, k2)))
        s  = (cl :- 1)*k2 + vv
    }
    rngstate(rs)
    // quantiles of the expected counts at the horizon (paths below the limit)
    // and of the fitted expected counts in the sample (posterior weights)
    qq  = (0.5, 0.9, 0.99)
    ok  = selectindex(ex :== 0)
    sv1 = J(1, 4, .)
    sv2 = J(1, 4, .)
    if (rows(ok) > 0) {
        for (i=1; i<=3; i++) {
            sv1[i] = hmdbp_pctile(mx1[ok], qq[i])
            sv2[i] = hmdbp_pctile(mx2[ok], qq[i])
        }
        sv1[4] = max(mx1[ok])
        sv2[4] = max(mx2[ok])
    }
    fit1 = rowsum(W:*mm.Mu1)
    fit2 = rowsum(W:*mm.Mu2)
    st_matrix(stub + "S", (mean(ex), rows(ind), H, lim, (M.hasD ? mean(surv) : .)))
    st_matrix(stub + "Q", (sv1 \ sv2 \ (hmdbp_pctile(fit1, 0.5), hmdbp_pctile(fit1, 0.9), hmdbp_pctile(fit1, 0.99), max(fit1)) \ (hmdbp_pctile(fit2, 0.5), hmdbp_pctile(fit2, 0.9), hmdbp_pctile(fit2, 0.99), max(fit2))))
}

// ---------------------------------------------------------------- estat fit
// distributions of the counts (0..maxc-1, maxc or more): rows maxc+1 for y1,
// then for y2; columns: observed, one per type (states, or classes with
// bycl), mixture. Goodness of fit (matrix stub + "G"): chi2, df,
// dissimilarity, largest absolute difference, N, by count
void hmdbp_fit(struct hmdbp_dat scalar M, real rowvector b, string scalar stub,
    real scalar ptype, real scalar maxc, real scalar bycl)
{
    struct hmdbp_mm scalar mm
    real matrix    W, lf, Pi0, P, PiV, PiU, Pv, Pr, E, F, G, Mu, cum, Pc
    real colvector y, O, Ex
    real rowvector sw, num
    real scalar    ll, j, c, row, o, nt, r0

    ll = hmdbp_pestep(M, b, W, mm, lf, Pi0, P, PiV, PiU, Pv)
    if (ll >= .) {
        errprintf("the log likelihood cannot be evaluated at e(b) with the data in memory\n")
        exit(459)
    }
    Pr = hmdbp_pr(M, W, lf, Pi0, P, ptype, 0, 0)
    if (bycl) E = M.Euf
    else      E = M.Ev
    nt = cols(E)
    sw = colsum(Pr*E)
    F  = J(2*(maxc+1), nt+2, .)
    G  = J(2, 5, .)
    for (j=1; j<=2; j++) {
        if (j == 1) {
            y  = M.y1
            Mu = mm.Mu1
        }
        else {
            y  = M.y2
            Mu = mm.Mu2
        }
        cum = J(M.N, M.K, 0)
        for (c=0; c<=maxc; c++) {
            row = (j-1)*(maxc+1) + c + 1
            if (c < maxc) {
                Pc  = poissonp(Mu, c)
                cum = cum + Pc
                o   = mean(y :== c)
            }
            else {
                Pc = 1 :- cum
                Pc = Pc:*(Pc :> 0)
                o  = mean(y :>= maxc)
            }
            num = colsum((Pr:*Pc)*E)
            F[row, .] = (o, num:/sw, sum(num)/M.N)
        }
        r0 = (j-1)*(maxc+1)
        O  = F[|r0+1,1 \ r0+maxc+1,1|]
        Ex = F[|r0+1,nt+2 \ r0+maxc+1,nt+2|]
        G[j, .] = (M.N*sum((O - Ex):^2 :/ Ex), maxc, 0.5*sum(abs(O - Ex)), max(abs(O - Ex)), M.N)
    }
    st_matrix(stub, F)
    st_matrix(stub + "G", G)
}

end
