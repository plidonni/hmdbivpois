*! version 1.0.0  18nov2026
*! predict after hmdbivpois. The program is in hmdbivpois.ado (subprogram
*! Predict), which holds the whole command: this file only hands the request
*! over.
*!
*!   predict [type] {newvar | stub* | newvarlist} [if] [in] [, statistic weights]
*!
*!   statistic   mu1        expected value of the first count (the default)
*!               mu2        expected value of the second count
*!               pdrop      probability of dropout after the period
*!               postprob   probabilities of the latent states (k variables)
*!               classprob  probabilities of the latent classes (classes() variables)
*!               state      most likely state of each period
*!               class      most likely class
*!               viterbi    state of each period in the most likely sequence of
*!                          (class, states) of the individual
*!               rho        correlation between the two counts given the latent
*!                          variables
*!               pr(a,b)    probability that the two counts equal (a, b)
*!   weights     (default)  posterior probabilities of the latent variables,
*!                          given all the observed outcomes of the individual
*!               filtered   given the outcomes up to the period
*!               onestep    given the outcomes up to the period before
*!               prior      given the covariates and the presample outcomes only
*!               instate(#) in state #, inclass(#) in class # (the other latent
*!                          variable keeps its weights)

program define hmdbivpois_p
    version 16
    hmdbivpois __predict `0'
end
