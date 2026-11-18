{smcl}
{* *! version 1.0.0  18nov2026}{...}
{vieweralsosee "hmbivpois" "help hmbivpois"}{...}
{vieweralsosee "poisson" "help poisson"}{...}
{vieweralsosee "xtset" "help xtset"}{...}
{viewerjumpto "Syntax" "hmdbivpois##syntax"}{...}
{viewerjumpto "Description" "hmdbivpois##description"}{...}
{viewerjumpto "Model" "hmdbivpois##model"}{...}
{viewerjumpto "Data" "hmdbivpois##data"}{...}
{viewerjumpto "Options" "hmdbivpois##options"}{...}
{viewerjumpto "Postestimation" "hmdbivpois##postestimation"}{...}
{viewerjumpto "Stored results" "hmdbivpois##results"}{...}
{viewerjumpto "Examples" "hmdbivpois##examples"}{...}
{viewerjumpto "References" "hmdbivpois##references"}{...}
{title:Title}

{phang}
{bf:hmdbivpois} {hline 2} Hidden Markov dropout bivariate Poisson regression for panel counts
with informative dropout


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:hmdbivpois}
{cmd:(}{it:depvar1} {cmd:=} [{it:indepvars1}] [{cmd:,} {it:eqopts}]{cmd:)}
{cmd:(}{it:depvar2} {cmd:=} [{it:indepvars2}] [{cmd:,} {it:eqopts}]{cmd:)}
[{cmd:(}{it:dropvar} {cmd:=} [{it:indepvarsd}] [{cmd:,} {it:eqopts}]{cmd:)}]
{ifin}{cmd:,} {opth id(varname)} {opth time(varname)} {opt k(#)} [{it:options}]

{synoptset 26 tabbed}{...}
{synopthdr:eqopts}
{synoptline}
{synopt:{opth statev:ary(varlist)}}regressors of the equation whose coefficients differ across the states{p_end}
{synopt:{opt lagv:ary(depvars)}}lagged counts ({it:depvar1}, {it:depvar2}) whose coefficients differ across the states{p_end}
{synoptline}

{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{p2coldent:* {opth id(varname)}}individual identifier{p_end}
{p2coldent:* {opth time(varname)}}integer time variable{p_end}
{p2coldent:* {opt k(#)}}number of latent states of the Markov chain, 1 to 20{p_end}
{synopt:{opt classes(#)}}number of latent classes (persistent heterogeneity), 1 to 10; default {cmd:1}{p_end}
{synopt:{opth initcov(varlist)}}covariates of the initial probabilities of the states{p_end}
{synopt:{opth classcov(varlist)}}covariates of the prior probabilities of the classes{p_end}
{synopt:{opt notr:ansition}}no transitions: the state of each individual is constant over time{p_end}
{synopt:{opt lagt:ransform(g)}}transformation of the lagged counts: {opt linear} (default), {opt log1p}, {opt logc(#)}, {opt logzero}, {opt indicator}, {opt topcode(#)}{p_end}
{synopt:{opt dropc:ounts(which)}}counts in the dropout equation: {opt lag} (default), {opt current} or {opt both}{p_end}
{synopt:{opt classeq(which)}}equations shifted by the classes: {opt all} (default), {opt dropout} or {opt counts}{p_end}

{syntab:Transitions}
{synopt:{opth xtr:ans(varlist)}}covariates of the transition probabilities, measured at the period of arrival{p_end}
{synopt:{opt tvtr:ans}}transition probabilities that vary freely by period (period dummies){p_end}

{syntab:Estimation}
{synopt:{opt emiter:ate(#)}}maximum number of EM iterations; default {cmd:500}{p_end}
{synopt:{opt emtol:erance(#)}}tolerance of the EM on the relative change of the log likelihood; default {cmd:1e-8}{p_end}
{synopt:{opt emonly}}no direct maximization after the EM{p_end}
{synopt:{opt tech:nique(name)}}{opt bfgs} (default), {opt dfp} or {opt nr} for the direct maximization; specifying it forces the direct maximization{p_end}
{synopt:{opt iter:ate(#)}}maximum number of iterations of the direct maximization; default {cmd:300}{p_end}
{synopt:{opt from(matname)}}starting values, a row vector in the order of {cmd:e(b)}{p_end}
{synopt:{opt starts(#)}}number of starting values; default {cmd:1}{p_end}
{synopt:{opt startit:erate(#)}}EM iterations from each starting value; default as {opt emiterate()}{p_end}
{synopt:{opt startseed(#)}}seed of the random starting values; default {cmd:20261117}{p_end}
{synopt:{opt sort:by(how)}}labelling of the states and classes: {opt eq1} (default), {opt eq2} or {opt none}{p_end}

{syntab:SE/Reporting}
{synopt:{opt vce(vcetype)}}{opt oim} (default) or {opt robust} (clustered on {cmd:id()}){p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}
{synopt:{opt nolo:g}}do not display the iteration log{p_end}
{synoptline}
{p 4 6 2}* {opt id()}, {opt time()} and {opt k()} are required.{p_end}
{p 4 6 2}{it:indepvars1}, {it:indepvars2}, {it:indepvarsd}, {opt statevary()} and {opt xtrans()} may contain factor variables; see {help fvvarlist}.{p_end}
{p 4 6 2}{cmd:hmdbivpois} is made of three files that must be in the same folder:
{cmd:hmdbivpois.ado} (the whole command), {cmd:hmdbivpois_p.ado} (for {cmd:predict}) and {cmd:hmdbivpois_estat.ado} (for {cmd:estat}).{p_end}


{marker description}{...}
{title:Description}

{pstd}
{cmd:hmdbivpois} fits the hidden Markov dropout bivariate Poisson (HMDBPo) model of Carrieri,
Gyrd-Hansen and Li Donni (2026) to two counts observed on the same individuals over time, when
individuals may leave the panel (for example, because they die) for reasons related to the
outcomes. The model combines

{p 8 11 2}- a bivariate Poisson distribution for the two counts in each period, with means that
depend on regressors and on the lagged counts of both outcomes (state dependence within and
between outcomes);{p_end}
{p 8 11 2}- a latent state that follows a first-order Markov chain (time-varying unobserved
heterogeneity);{p_end}
{p 8 11 2}- a latent class, constant over time (persistent heterogeneity, for example frailty);{p_end}
{p 8 11 2}- a discrete-time logit model for the dropout, which depends on the same latent
state and class (shared-parameter model, dropout missing not at random).{p_end}

{pstd}
Without the dropout equation and with {cmd:classes(1)}, the model is the hidden Markov bivariate
Poisson model of Laudicella and Li Donni (2022), fitted by {helpb hmbivpois}, and the two commands
return the same estimates.


{marker model}{...}
{title:Model}

{pstd}
Let C_i = 1..{it:kc} be the class of individual i and V_it = 1..{it:k} its state in period t.
Given C_i = u and V_it = v, the two counts follow the bivariate Poisson distribution of
Lakshminarayana et al. (1999) (Famoye 2010):

{p 8 8 2}
p(y1,y2) = Po(y1;mu1) Po(y2;mu2) [1 + lambda (e^-y1 - e^-q mu1)(e^-y2 - e^-q mu2)], q = 1 - e^-1,

{p 8 8 2}
ln mu_jit = x_jit b_j + g_j1 y1_i,t-1 + g_j2 y2_i,t-1 + alpha_jv + delta_ju,  j = 1, 2,

{pstd}
with delta_j1 = 0. The margins are exactly Poisson: lambda gives the sign of the dependence between
the two counts within the period, given the latent variables, and does not change the means.
lambda must lie in an interval that depends on the means; the command keeps it inside the interval
at the fitted means, which is reported after the estimates.

{pstd}
Dropout: D_it = 1 if the individual leaves the panel after period t. In each period at risk

{p 8 8 2}
logit P(D_it = 1 | u, v) = xd_it bd + gd1 y1_i,t-1 + gd2 y2_i,t-1 + alpha_dv + delta_du.

{pstd}
With {cmd:dropcounts(current)} the counts of period t, y1_it and y2_it, replace the lagged ones in the
dropout equation; with {cmd:dropcounts(both)} both enter. With {cmd:classeq(dropout)} delta_1u =
delta_2u = 0 (the classes shift only the dropout); with {cmd:classeq(counts)} delta_du = 0.

{pstd}
State-specific coefficients. With {cmd:statevary()} and {cmd:lagvary()} inside an equation, the
coefficients of the listed regressors and lagged counts depend on the state V_it = v (not on the
class): x_jit b_j becomes x_jit b_jv for those terms. They are reported in the equations
{it:depvar}{cmd:_s1}, ..., {it:depvar}{cmd:_s}{it:k}.

{pstd}
Latent variables: P(C_i = u) is a multinomial logit on {opt classcov()} and on the presample counts
y1_i0 and y2_i0 (class 1 is the base); P(V_i1 = v) is a multinomial logit on {opt initcov()} and
on y1_i0 and y2_i0 (state 1 is the base); the transitions P(V_it = w | V_i,t-1 = v) are constant
over time (homogeneous), and each row is a logit with "staying" as the base. With {opt xtrans()}
each row is a multinomial logit on covariates measured at the period of arrival (for the transition
from t-1 to t, their values in period t), with {opt tvtrans} on period dummies (a free transition
matrix for each period), or on both. C and V are independent a priori.
Given the class and the sequence of states, the periods are independent; within the period the
two counts depend on each other through lambda, and the dropout depends on the counts only through
the latent variables and the lags.

{pstd}
The class shifts the intercepts of all the equations: in the dropout equation it separates
persistent frailty from the time-varying health dynamics captured by the states.

{pstd}
Estimation. The pair (C_i, V_it) is treated as one hidden Markov chain on {it:kc}{c 215}{it:k}
states, whose transition matrix is block diagonal (the class never changes). The log likelihood
is computed with the scaled forward recursion; the EM algorithm uses the forward-backward
recursions in the E-step and, in the M-step, the closed form for the transitions, Newton steps
for the measurement and dropout equations and fractional multinomial logits for the initial
states and for the classes. Direct maximization with the analytic score follows if the EM does
not converge within {opt emiterate()}, or if {opt technique()} is specified. The variance is the
inverse of the observed information, from central differences of the analytic score.

{pstd}
Labels. The states are ordered by increasing alpha_1v ({cmd:sortby(eq1)}) or alpha_2v
({cmd:sortby(eq2)}); the classes by increasing dropout shift delta_du (class 1 has the lowest
dropout risk), or by delta_1u (delta_2u with {cmd:sortby(eq2)}) without the dropout equation or with
{cmd:classeq(counts)}.


{marker data}{...}
{title:Data}

{pstd}
Long format, one row per individual and period, with consecutive periods. The command builds the
lagged counts: the first observed period of each individual is the presample period and is not an
estimation period; its counts enter the initial probabilities of the states and the prior
probabilities of the classes, and the lags of the first estimation period.

{pstd}
{it:dropvar} is 0/1 and may be 1 only in the last observed period of an individual: 1 means that
the individual leaves the panel after that period. All the estimation periods of an individual are
at risk of dropout except the last one when {it:dropvar} = 0 there: in that period the follow-up
ends (or the individual is censored) and the dropout is not observed. This is the interval-censored
dropout of Carrieri et al. (2026), eq. (6). {it:dropvar} must not be missing in the presample period.
A dropout in the presample period, or in a period with missing regressors, is not used (a note
reports how many).


{marker options}{...}
{title:Options}

{dlgtab:Model}

{phang}
{opth id(varname)} and {opth time(varname)} identify individuals and periods.

{phang}
{opt k(#)} sets the number of latent states of the Markov chain.

{phang}
{opt classes(#)} sets the number of latent classes. With {cmd:classes(1)} the persistent
heterogeneity is absent and the dropout depends on the states only.

{phang}
{opth initcov(varlist)} and {opth classcov(varlist)} add covariates, taken in the first estimation
period of each individual, to the initial probabilities of the states and to the prior
probabilities of the classes.

{phang}
{opt notransition} fixes the state of each individual over time; {cmd:estat lrtrans} tests it.

{phang}
{opt lagtransform(g)} replaces the lagged counts y_j,t-1 with g(y_j,t-1) in the three equations, and the
presample counts y_j0 with g(y_j0) in the initial probabilities of the states and in the prior
probabilities of the classes. With linear lags, a positive coefficient and large counts the means can
grow fast; a concave g() keeps the dynamics bounded ({cmd:estat stability} checks it). {opt log1p}: ln(1 + y); {opt logc(#)}:
ln(max(y, #)), 0 < # <= 1; {opt logzero}: ln(y) for y > 0 and 0 for y = 0, plus a dummy for y = 0
(two more coefficients per equation, named {cmd:lag0_}{it:depvar}); {opt indicator}: 1(y > 0);
{opt topcode(#)}: min(y, #), # >= 1. With {opt indicator} and {opt logzero} each lagged count must
take both zero and positive values in the estimation sample. {cmd:estat ame} computes the effect of
one more lagged count through g().

{phang}
{opt statevary(varlist)} and {opt lagvary(depvars)} go inside the parentheses of an equation, after a
comma: {cmd:(y1 = age income, statevary(income) lagvary(y1))}. The terms of {opt statevary()} must be
among the regressors of the equation; {opt lagvary()} takes {it:depvar1} and/or {it:depvar2} (with
{opt lagtransform(logzero)} it moves both the log and the zero dummy of that lag). In the dropout
equation {opt lagvary()} needs the lagged counts ({opt dropcounts(lag)} or {opt dropcounts(both)});
the counts of the period keep common coefficients. They are ignored with {cmd:k(1)}. The state-specific
coefficients start from the common ones of the starting regressions and follow the relabelling of
the states. {cmd:estat ame} averages their effects with the weights of the states, and
{cmd:estat stability} reports the coefficients of the lags and their spectral radius by latent type.

{phang}
{opt dropcounts(which)} sets the counts in the dropout equation: {opt lag} (default) the lagged counts,
{opt current} the counts of the period after which the individual may leave, observed before the
dropout (for example the use of care in the quarter before death), {opt both} the two. The current
counts go through the same g() of {opt lagtransform()} and are named {cmd:cur_}{it:depvar}
({cmd:cur0_}{it:depvar} for the dummies of {opt logzero}). The dependence of the dropout on the
counts not observed after it still passes only through the latent variables. {cmd:estat ame} reports
the effect on the dropout of one more count of the current period.

{phang}
{opt classeq(which)} sets the equations whose intercepts are shifted by the classes: {opt all}
(default, as in Carrieri et al. 2026), {opt dropout} (the classes are a frailty that affects only the
risk of dropout; the counts depend on the states only) or {opt counts} (the classes do not affect
the dropout). The restricted models are nested in the default one, with the shifts at zero inside
the parameter space: an LR test against the default model has the usual chi-squared distribution
with 2({it:kc}-1) ({opt dropout}) or {it:kc}-1 ({opt counts}) degrees of freedom; give them to
{cmd:lrtest} with {opt df()}, because {cmd:lrtest} takes them from {cmd:e(rank)}, which drops when some
parameters are weakly identified. With {opt classeq(dropout)} the classes are identified only from
the dropout and its covariates, and the intercepts of the dropout equation can be imprecise: one
dropout per individual at most, over a short follow-up, gives little information on the mixture.

{phang}
{opth xtrans(varlist)} makes the transition probabilities depend on covariates. The row of the
transition matrix from state c is a multinomial logit with base "staying in c": P(c -> d) is
proportional to exp(w_it g_cd), d != c, where w_it holds the covariates in the period of arrival t
and a constant. The covariates must not be missing in the estimation periods.

{phang}
{opt tvtrans} replaces the constant of the transitions with one dummy for each period of arrival.
With both options the covariates shift the period-specific transitions.

{dlgtab:Estimation}

{phang}
{opt emiterate(#)}, {opt emtolerance(#)}, {opt emonly}, {opt technique(name)}, {opt iterate(#)}
control the EM and the direct maximization, as in {helpb hmbivpois}.

{phang}
{opt from(matname)} gives starting values in the order of {cmd:e(b)}. By default the starting values
come from Poisson regressions for the counts and a logit for the dropout on the rows at risk; the
state intercepts of the counts are spread around the constants, and the classes start with
increasing dropout risk.

{phang}
{opt starts(#)} reduces the risk of a local maximum, as in Carrieri et al. (2026): the EM is run from
{it:#} starting values, the default ones (or {opt from()}) and {it:#}-1 random ones, for
{opt startiterate()} iterations each (by default {opt emiterate()}, that is to convergence), and goes on
from the one with the largest log likelihood. The random starting values keep the slopes of the default
ones and set lambda to zero; they draw the intercepts of the states around the mean intercept (sorted
in the count equations), the shifts of the classes, the constants of the initial and prior
probabilities, and transition matrices with staying probabilities between 0.6 and 0.95. A smaller
{opt startiterate()} (for example 50) screens many starts at a lower cost. The log likelihood reached
from each start is in {cmd:e(ll_starts)}; {cmd:e(n_best)} counts the starts that reach the best one
(within 0.1): when it is 1, try more starts. {opt startseed(#)} sets the seed of the random starting
values; the state of the random-number generator is restored afterwards.

{phang}
{opt sortby(how)} sets the labelling; see {help hmdbivpois##model:Model}.

{dlgtab:SE/Reporting}

{phang}
{opt vce(robust)} gives the sandwich estimator clustered on individuals.


{marker postestimation}{...}
{title:Postestimation}

{pstd}
The postestimation commands rebuild the data of the estimation from the variables in memory; they must
be those used by {cmd:hmdbivpois}, including the presample period of each individual.

{p 8 17 2}
{cmd:predict} [{it:type}] {c -(}{it:newvar} | {it:stub}{cmd:*} | {it:newvarlist}{c )-} {ifin}
[{cmd:,} {it:statistic} {it:weights}]

{synoptset 20 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt:{opt mu1}}expected value of {it:depvar1}; the default{p_end}
{synopt:{opt mu2}}expected value of {it:depvar2}{p_end}
{synopt:{opt pd:rop}}probability of dropout after the period{p_end}
{synopt:{opt post:prob}}probabilities of the states ({opt k()} variables){p_end}
{synopt:{opt classpr:ob}}probabilities of the classes ({opt classes()} variables){p_end}
{synopt:{opt st:ate}}most likely state of the period{p_end}
{synopt:{opt class}}most likely class{p_end}
{synopt:{opt vit:erbi}}state of the period in the most likely sequence of class and states of the individual{p_end}
{synopt:{opt rho}}correlation of the two counts given the latent variables{p_end}
{synopt:{opt pr(a,b)}}probability that the two counts equal {it:a} and {it:b}{p_end}
{synoptline}

{synopthdr:weights}
{synoptline}
{synopt:(default)}posterior probabilities of the latent variables, given all the observed outcomes of the individual{p_end}
{synopt:{opt filt:ered}}given the outcomes up to the period{p_end}
{synopt:{opt one:step}}given the outcomes up to the period before (one-step-ahead){p_end}
{synopt:{opt pri:or}}given the covariates and the presample counts only{p_end}
{synopt:{opt inst:ate(#)}}in state #; the class keeps its weights{p_end}
{synopt:{opt incl:ass(#)}}in class #; the state keeps its weights{p_end}
{synoptline}
{p 4 6 2}The outcomes of a period are the two counts and, in the periods at risk, the dropout. The
posterior probabilities of the classes are constant within the individual.{p_end}

{p 8 17 2}
{cmd:estat ame} [{cmd:,} {opt pri:or} {opt one:step} {opt filt:ered} {opt inst:ate(#)} {opt incl:ass(#)} {opt l:evel(#)}]

{pstd}
Average marginal effects on E({it:depvar1}), E({it:depvar2}) and on the probability of dropout, over the
estimation periods, with the probabilities of the latent variables held fixed (weights as in
{cmd:predict}): dy/dx for continuous variables, discrete change from the base level for factor levels,
effect of one more count of each lagged outcome. Standard errors by the delta method.

{p 8 17 2}
{cmd:estat survival} [{cmd:,} {opth over(varname)} {opt one:step} {opt pri:or} {opt filt:ered} {opt post:erior}]

{pstd}
Kaplan-Meier survival by period of follow-up (the estimation periods of each individual: 1, 2, ...)
against the survival predicted by the model: the product over the periods of 1 minus the mean predicted
probability of dropout among the individuals at risk (Fig. 8 of Carrieri et al. 2026). The default
weights are one-step-ahead: given the counts and the survival up to the period before. Prior weights
ignore the selection of the survivors and usually overpredict the dropout in the later periods.

{p 8 17 2}
{cmd:estat fit} [{cmd:,} {opt max(#)} {opt classes} {opt post:erior} {opt one:step} {opt pri:or} {opt filt:ered}]

{pstd}
Observed distribution of each count over the categories 0, 1, ..., {it:#}-1 and {it:#} or more
(default {cmd:max(5)}) against the distribution predicted in each latent state (each class with
{opt classes}) and against their mixture, the predicted marginal distribution (Table 5, Table 6 and
Fig. 7 of Carrieri et al. 2026). The distribution of a type is the mean over the periods of the
Poisson probabilities in that type, weighted by the probabilities of the latent variables (posterior
by default); the mixture weights the types by their mean probabilities. Reported: Pearson chi2 over
the person-periods, dissimilarity index (the share of probability that should be moved to match the
observed distribution) and largest difference of a cell. With many person-periods the chi2 rejects
for negligible differences: it is a descriptive measure. Posterior weights use the observed outcomes;
{opt onestep} or {opt prior} give a check out of sample.

{p 8 17 2}
{cmd:estat transition} [{cmd:,} {opt per:iod(#)} {opt l:evel(#)}]

{pstd}
Transition probabilities, expected durations, mean initial probabilities of the states and mean prior
probabilities of the classes, with delta-method standard errors. With {opt xtrans()} or {opt tvtrans},
the means over the transition rows, or over the transitions into period {it:#} with {opt period(#)}.

{p 8 17 2}
{cmd:estat lrtrans} [{cmd:,} {opt reps(#)} {opt seed(#)}]

{pstd}
Tests the model without transitions between the states ({opt notransition}: the states, like the
classes, are constant over time) against the model with homogeneous transitions in memory. The
restricted model has the same equations, sample and options, and is fitted from two sets of starting
values (the default ones and the estimates in memory without the transitions); the larger log
likelihood is kept. The null hypothesis is on the boundary of the parameter space: the LR statistic
has a chi-bar-squared distribution, a mixture of chi2(0), ..., chi2(k(k-1)), whose weights are
estimated by {opt reps()} simulations (default 10000; {opt seed()} sets the seed, and the state of the
random-number generator is restored afterwards) from the information of the off-diagonal transition
probabilities net of the other parameters at the restricted estimates, estimated by the outer product
of the individual scores (Bartolucci 2006). The p-value of the chi2(k(k-1)) distribution, an upper
bound, is also reported. Under the null hypothesis the maximum of the model with transitions is on the
boundary, where the EM converges slowly and may stop early, which makes LR too small and the test
conservative. No score test is reported: the individual scores of the transition probabilities are
sums of likelihood ratios between states, with very heavy tails, and a score statistic based on their
outer product is badly calibrated. Not available with {opt xtrans()} or {opt tvtrans}.

{p 8 17 2}
{cmd:estat stability} [{cmd:,} {opt horizon(#)} {opt limit(#)} {opt nsim(#)} {opt seed(#)}]

{pstd}
Checks whether the dynamics of the lagged counts can explode, as in {helpb hmbivpois}. The coefficients
of the lags are common to the latent types (state and class) unless {opt lagvary()} makes them state
specific; the rest of the linear predictor varies by type. The command reports, by type, the
coefficients of the lags in the count equations and in the dropout equation; with log-type lags
({opt log1p}, {opt logc()}, {opt logzero}) the spectral radius and the norm of the 2 x 2 matrix of the
coefficients in the count equations, with linear lags the index g_jj exp(a),
which must not exceed exp(-1); the share of deterministic paths, from every observation and in every
type, whose expected count exceeds {opt limit()} (default 1000) within {opt horizon()} periods
(default 50); and the same share for {opt nsim()} simulated paths per individual (default 1) of the
estimated model, from the last observed period, with the type drawn from its posterior probabilities.
The simulated counts are drawn as if the individual stayed in the panel; the mean probability of
still being in the panel at the horizon is also reported. {opt seed()} sets the seed of the
simulation; the state of the random-number generator is restored afterwards.

{p 8 17 2}
{cmd:estat ic} [{cmd:,} {opt n(#)}]

{pstd}
AIC, BIC with the number of individuals (or {it:#}) and with the number of person-periods, relative
entropy of the classification into states (by period) and into classes (by individual).


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:hmdbivpois} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of estimation periods{p_end}
{synopt:{cmd:e(N_g)}}number of individuals{p_end}
{synopt:{cmd:e(k)}}number of states{p_end}
{synopt:{cmd:e(classes)}}number of classes{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(npar)}}number of parameters{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(bic)}}BIC with ln(number of individuals){p_end}
{synopt:{cmd:e(bic_N)}}BIC with ln(number of estimation periods), as in Table 3 of Carrieri et al. (2026){p_end}
{synopt:{cmd:e(aic)}}AIC{p_end}
{synopt:{cmd:e(N_drop)}}number of dropouts used{p_end}
{synopt:{cmd:e(N_risk)}}number of periods at risk of dropout{p_end}
{synopt:{cmd:e(N_droplost)}}number of dropouts not used{p_end}
{synopt:{cmd:e(converged)}}1 if converged{p_end}
{synopt:{cmd:e(em_iter)}}number of EM iterations (after the choice of the start, with {opt starts()}){p_end}
{synopt:{cmd:e(starts)}}number of starting values{p_end}
{synopt:{cmd:e(start_best)}}start with the largest log likelihood ({opt starts()} > 1){p_end}
{synopt:{cmd:e(n_best)}}starts within 0.1 of the largest log likelihood{p_end}
{synopt:{cmd:e(start_iter)}, {cmd:e(start_seed)}}EM iterations from each start and seed{p_end}
{synopt:{cmd:e(g_min)}, {cmd:e(g_avg)}, {cmd:e(g_max)}}estimation periods per individual{p_end}
{synopt:{cmd:e(lambda_lo)}, {cmd:e(lambda_hi)}}admissible interval of lambda at the fitted means{p_end}
{synopt:{cmd:e(n_bound)}}1 if lambda is on a bound of its interval{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and variance{p_end}
{synopt:{cmd:e(Pi)}}transition matrix of the states (with {opt xtrans()} or {opt tvtrans}, the mean over the transition rows){p_end}
{synopt:{cmd:e(Pi_t)}}with {opt xtrans()} or {opt tvtrans}: mean transition probabilities by period of arrival, one row per period{p_end}
{synopt:{cmd:e(omega)}}average posterior shares of the states{p_end}
{synopt:{cmd:e(omega_t)}}average posterior shares of the states by period of follow-up{p_end}
{synopt:{cmd:e(pi0)}}average initial probabilities of the states{p_end}
{synopt:{cmd:e(alpha)}}state intercepts, by equation{p_end}
{synopt:{cmd:e(delta)}}class shifts, by equation (zero for class 1){p_end}
{synopt:{cmd:e(pclass)}}average prior probabilities of the classes{p_end}
{synopt:{cmd:e(wclass)}}average posterior probabilities of the classes{p_end}
{synopt:{cmd:e(lambda_range)}}admissible interval of lambda{p_end}
{synopt:{cmd:e(ll_starts)}}log likelihood and EM iterations from each starting value ({opt starts()} > 1){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:hmdbivpois}{p_end}
{synopt:{cmd:e(depvar)}}the two counts{p_end}
{synopt:{cmd:e(dropout)}}the dropout indicator{p_end}
{synopt:{cmd:e(ivar)}, {cmd:e(tvar)}}identifier and time variable{p_end}
{synopt:{cmd:e(transitions)}}{cmd:homogeneous}, {cmd:heterogeneous} ({opt xtrans()}, {opt tvtrans}) or {cmd:none}{p_end}
{synopt:{cmd:e(xtrans)}, {cmd:e(tvtrans)}}options of the transitions{p_end}
{synopt:{cmd:e(lagtransform)}}transformation of the lagged counts{p_end}
{synopt:{cmd:e(dropcounts)}}counts in the dropout equation ({cmd:lag}, {cmd:current}, {cmd:both}){p_end}
{synopt:{cmd:e(curnames)}}names of the features of the current counts in the dropout equation{p_end}
{synopt:{cmd:e(classeq)}}equations shifted by the classes ({cmd:all}, {cmd:dropout}, {cmd:counts}){p_end}
{synopt:{cmd:e(statevary1)}, {cmd:e(statevary2)}, {cmd:e(statevaryd)}}terms with state-specific coefficients{p_end}
{synopt:{cmd:e(initcov)}, {cmd:e(classcov)}}covariates of the latent variables{p_end}
{synopt:{cmd:e(vce)}, {cmd:e(vcetype)}}variance estimator{p_end}

{pstd}
{cmd:estat lrtrans} stores {cmd:r(lr)}, {cmd:r(p)}, {cmd:r(p_chi2)}, {cmd:r(df)}, {cmd:r(ll_1)}, {cmd:r(ll_0)},
{cmd:r(reps)} and {cmd:r(weights)}. {cmd:estat stability} stores {cmd:r(G)}, {cmd:r(radius)}, {cmd:r(index)},
{cmd:r(det)}, {cmd:r(sim)} (share above the limit, paths, horizon, limit, probability of staying) and
{cmd:r(q)}.


{marker examples}{...}
{title:Examples}

{pstd}
The examples use {cmd:hmdbivpois_example.dta}, a simulated panel of 2,000 individuals observed in a
presample period and up to 6 estimation periods, generated from the model with three states, two
classes, informative dropout ({cmd:died} = 1 in the last observed period of those who leave) and
homogeneous transitions. The variables {cmd:state} and {cmd:class} hold the true latent state and
class. Download the data and the do-file with the examples into the current folder with{p_end}
{phang2}{cmd:. net get hmdbivpois, from("https://raw.githubusercontent.com/plidonni/hmdbivpois/main/")}{p_end}

{pstd}Setup{p_end}
{phang2}{cmd:. use hmdbivpois_example, clear}{p_end}
{phang2}{cmd:. xtset id time}{p_end}

{pstd}Three states, two classes, dropout equation{p_end}
{phang2}{cmd:. hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) k(3) classes(2) initcov(z1) classcov(zu)}{p_end}

{pstd}Postestimation{p_end}
{phang2}{cmd:. estat ic}{p_end}
{phang2}{cmd:. estat transition}{p_end}
{phang2}{cmd:. estat ame}{p_end}
{phang2}{cmd:. estat survival}{p_end}
{phang2}{cmd:. estat fit}{p_end}
{phang2}{cmd:. predict double pd, pdrop}{p_end}
{phang2}{cmd:. predict s_hat, viterbi}{p_end}
{phang2}{cmd:. predict c_hat, class}{p_end}
{phang2}{cmd:. tabulate state s_hat if time > 0}{p_end}
{phang2}{cmd:. tabulate class c_hat if time > 0}{p_end}
{phang2}{cmd:. estat lrtrans}{p_end}

{pstd}Choice of the number of states and classes with BIC{p_end}
{phang2}{cmd:. forvalues c = 1/2 {c -(}}{p_end}
{phang2}{cmd:.     forvalues k = 2/3 {c -(}}{p_end}
{phang2}{cmd:.         quietly hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) k(`k') classes(`c') initcov(z1) classcov(zu) nolog}{p_end}
{phang2}{cmd:.         display "classes `c'  k `k'  BIC " %10.2f e(bic)}{p_end}
{phang2}{cmd:.     {c )-}}{p_end}
{phang2}{cmd:. {c )-}}{p_end}

{pstd}Dropout ignored (missing at random): the model of {helpb hmbivpois}{p_end}
{phang2}{cmd:. hmdbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(3) initcov(z1)}{p_end}

{pstd}Multiple starting values{p_end}
{phang2}{cmd:. hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) k(3) classes(2) initcov(z1) classcov(zu) starts(5)}{p_end}
{phang2}{cmd:. matrix list e(ll_starts)}{p_end}

{pstd}Classes in the dropout equation only, against the default (2 restrictions){p_end}
{phang2}{cmd:. quietly hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) k(3) classes(2) initcov(z1) classcov(zu) nolog}{p_end}
{phang2}{cmd:. estimates store all}{p_end}
{phang2}{cmd:. quietly hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) k(3) classes(2) initcov(z1) classcov(zu) classeq(dropout) nolog}{p_end}
{phang2}{cmd:. lrtest all ., df(2)}{p_end}

{pstd}Counts of the period in the dropout equation, slope of x11 by state, log1p lags{p_end}
{phang2}{cmd:. hmdbivpois (y1 = x11 x12, statevary(x11)) (y2 = x21 x22) (died = xd), id(id) time(time) k(3) classes(2) initcov(z1) classcov(zu) dropcounts(both) lagtransform(log1p)}{p_end}
{phang2}{cmd:. estat stability}{p_end}


{marker references}{...}
{title:References}

{phang}
Bartolucci, F. 2006. Likelihood inference for a class of latent Markov models under linear hypotheses
on the transition probabilities. {it:Journal of the Royal Statistical Society, Series B} 68: 155-178.

{phang}
Bartolucci, F. and A. Farcomeni. 2015. A discrete time event-history approach to informative drop-out
in mixed latent Markov models with covariates. {it:Biometrics} 71: 80-89.

{phang}
Carrieri, V., D. Gyrd-Hansen and P. Li Donni. 2026. Dynamic healthcare demand with informative dropout:
evidence from register data. {it:Journal of the Royal Statistical Society Series A}.
doi:10.1093/jrsssa/qnag108.

{phang}
Famoye, F. 2010. On the bivariate negative binomial regression model. {it:Journal of Applied Statistics}
37: 969-981.

{phang}
Lakshminarayana, J., S. N. N. Pandit and K. Srinivasa Rao. 1999. On a bivariate Poisson distribution.
{it:Communications in Statistics - Theory and Methods} 28: 267-276.

{phang}
Laudicella, M. and P. Li Donni. 2022. The dynamic interdependence in the demand of primary and emergency
secondary care: a hidden Markov approach. {it:Journal of Applied Econometrics} 37: 521-536.


{title:Author}

{pstd}
Paolo Li Donni{break}
University of Palermo{break}
{browse "mailto:paolo.lidonni@unipa.it":paolo.lidonni@unipa.it}{break}
{browse "https://sites.google.com/site/plidonni/"}{break}
{browse "https://github.com/plidonni/hmdbivpois"}{break}
Please report bugs and suggestions by email.{p_end}
