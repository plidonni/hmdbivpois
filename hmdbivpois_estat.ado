*! version 1.0.0  18nov2026
*! estat after hmdbivpois. The program is in hmdbivpois.ado (subprogram Estat),
*! which holds the whole command: this file only hands the request over.
*!
*!   estat ame [, prior onestep filtered instate(#) inclass(#) level(#)]
*!       average marginal effects on E(y1), E(y2) and on the probability of
*!       dropout, with the probabilities of the latent variables held fixed
*!       (posterior by default): dy/dx for continuous regressors, discrete
*!       change from the base level for factor levels, one more count of each
*!       lagged outcome; delta-method standard errors
*!   estat survival [, over(varname) onestep prior filtered posterior]
*!       Kaplan-Meier survival by period of follow-up against the survival
*!       predicted by the model (one-step-ahead weights by default)
*!   estat fit [, max(#) classes posterior onestep prior filtered]
*!       observed distribution of each count (0, 1, ..., max or more) against
*!       the distributions predicted in each latent state (or class) and
*!       their mixture; Pearson chi2, dissimilarity index, largest cell
*!       difference
*!   estat transition [, period(#) level(#)]
*!       transition probabilities, expected durations, mean initial
*!       probabilities of the states and mean prior probabilities of the
*!       classes, with delta-method standard errors; with xtrans() or tvtrans
*!       the mean over the transition rows (into period # with period())
*!   estat lrtrans [, reps(#) seed(#)]
*!       LR test of no transitions between the states (notransition)
*!       against homogeneous transitions; chi-bar-squared p-values
*!   estat stability [, horizon(#) limit(#) nsim(#) seed(#)]
*!       explosive dynamics of the lagged counts: coefficients, spectral
*!       radius or index by latent type, deterministic and simulated paths
*!   estat ic [, n(#)]
*!       log likelihood, AIC, BIC with n = number of individuals (or n()) and
*!       with the number of person-periods, relative entropy of the
*!       classification into states and classes
*!   estat vce, estat summarize: as after any estimation command

program define hmdbivpois_estat
    version 16
    hmdbivpois __estat `0'
end
