*******************************************************************************
* hmdbivpois_examples.do
* Runs the examples of help hmdbivpois on the simulated dataset
* hmdbivpois_example.dta (2,000 individuals, presample + up to 6 periods,
* 3 states, 2 classes, informative dropout).
* Install the command and get the files first:
*   net install hmdbivpois, from("https://raw.githubusercontent.com/plidonni/hmdbivpois/main/") replace
*   net get hmdbivpois, from("https://raw.githubusercontent.com/plidonni/hmdbivpois/main/") replace
*******************************************************************************

version 16
clear all
set more off

capture log close hmdbex
log using "hmdbivpois_examples.log", replace text name(hmdbex)

use hmdbivpois_example, clear
xtset id time
describe
summarize y1 y2 died x11 x12 x21 x22 xd z1 zu
bysort id (time): generate byte last = (_n == _N)
tabulate class died if last, row
drop last

* three states, two classes, dropout equation
hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) ///
    k(3) classes(2) initcov(z1) classcov(zu)

* postestimation
estat ic
estat transition
estat ame
estat survival
estat fit
predict double pd, pdrop
predict s_hat, viterbi
predict c_hat, class
tabulate state s_hat if time > 0
tabulate class c_hat if time > 0
estat lrtrans

* choice of the number of states and classes with BIC
forvalues c = 1/2 {
    forvalues k = 2/3 {
        quietly hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) ///
            k(`k') classes(`c') initcov(z1) classcov(zu) nolog
        display "classes `c'  k `k'  BIC " %10.2f e(bic)
    }
}

* dropout ignored (missing at random): the model of hmbivpois
hmdbivpois (y1 = x11 x12) (y2 = x21 x22), id(id) time(time) k(3) initcov(z1)

* multiple starting values
hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) ///
    k(3) classes(2) initcov(z1) classcov(zu) starts(5)
matrix list e(ll_starts)

* classes in the dropout equation only, against the default (2 restrictions)
quietly hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) ///
    k(3) classes(2) initcov(z1) classcov(zu) nolog
estimates store all
quietly hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) ///
    k(3) classes(2) initcov(z1) classcov(zu) classeq(dropout) nolog
lrtest all ., df(2)

* counts of the period in the dropout equation, slope of x11 by state, log1p lags
hmdbivpois (y1 = x11 x12, statevary(x11)) (y2 = x21 x22) (died = xd), id(id) time(time) ///
    k(3) classes(2) initcov(z1) classcov(zu) dropcounts(both) lagtransform(log1p)
estat stability

log close hmdbex
