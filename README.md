# hmdbivpois

**Hidden Markov bivariate Poisson regression with informative dropout, in Stata**

`hmdbivpois` fits the hidden Markov dropout bivariate Poisson model of Carrieri, Gyrd-Hansen and Li Donni (2026) to two counts observed on the same individuals over time, for example outpatient and emergency hospital care, when individuals may leave the panel (for example by death) and the dropout is informative. Two latent variables drive the counts and the dropout: a state that follows a first-order Markov chain, and a time-constant class (persistent frailty). Given them, the two counts have a bivariate Poisson distribution (Famoye form) whose means depend on regressors, on the lagged counts, on the state and on the class, and the probability of leaving after each period is a logit with the same latent variables. The dropout is missing not at random through the shared latent variables. Without the dropout equation and with one class the model is the hidden Markov bivariate Poisson of Laudicella and Li Donni (2022), fitted by [`hmbivpois`](https://github.com/plidonni/hmbivpois).

## Installation

```stata
net install hmdbivpois, from("https://raw.githubusercontent.com/plidonni/hmdbivpois/main/") replace
```

To download the example dataset and the do-file with the examples into the current folder:

```stata
net get hmdbivpois, from("https://raw.githubusercontent.com/plidonni/hmdbivpois/main/") replace
```

Requires Stata 16 or later. Then type `help hmdbivpois`. The technical guide is [`doc/hmdbivpois_guide.pdf`](doc/hmdbivpois_guide.pdf).

## Quick start

```stata
use hmdbivpois_example, clear
hmdbivpois (y1 = x11 x12) (y2 = x21 x22) (died = xd), id(id) time(time) ///
    k(3) classes(2) initcov(z1) classcov(zu)
estat ic
estat transition
estat survival
estat ame
predict double pd, pdrop
```

The data are in long form, one row per individual and period. The first observed period of each individual is the presample period: it provides the lagged counts of the first estimation period and the presample counts used in the prior probabilities of the states and of the classes. The dropout indicator is 1 only in the last observed period of those who leave the panel after it, and 0 elsewhere (also for those who are observed until the end of the follow-up).

## Features

- **Latent structure:** a Markov chain of `k()` states and `classes()` time-constant classes, with additive effects on the intercepts of the three equations; `classeq()` restricts the classes to the dropout equation or to the count equations.
- **Dropout:** logit of leaving after each period, on regressors, lagged counts and/or the counts of the period (`dropcounts()`), and the latent variables.
- **Lagged counts:** linear, `log1p`, `logc(#)`, `logzero`, indicator or top-coded (`lagtransform()`).
- **State-specific slopes** for chosen regressors and lags in each equation (`statevary()`, `lagvary()`).
- **Latent chain:** homogeneous transitions, transitions depending on covariates (`xtrans()`), transitions free by period (`tvtrans`), or no transitions (`notransition`).
- **Estimation:** EM algorithm with scaled forward-backward recursions on the expanded chain of states and classes, direct maximization with the analytic score when needed, multiple starting values (`starts()`); observed-information or cluster-robust standard errors.
- **Postestimation:** `predict` (expected counts, probability of dropout, posterior/filtered/one-step/prior probabilities of the states and of the classes, most likely state and class, Viterbi path, within-state correlation, joint probabilities); `estat ic`, `estat transition`, `estat ame` (average marginal effects on the counts and on the dropout), `estat survival` (Kaplan-Meier against the model), `estat fit` (observed and fitted distributions of the counts), `estat lrtrans` (chi-bar-squared test of no transitions), `estat stability` (explosive dynamics of the lags).

## Files

| File | Content |
|---|---|
| `hmdbivpois.ado` | the command (estimation, `predict` and `estat` routines) |
| `hmdbivpois_p.ado`, `hmdbivpois_estat.ado` | entry points for `predict` and `estat` |
| `hmdbivpois.sthlp` | help file |
| `hmdbivpois_example.dta` | simulated example data |
| `hmdbivpois_examples.do` | runs the examples of the help file |
| `doc/hmdbivpois_guide.pdf` | technical guide: model, estimation, postestimation, Monte Carlo evidence (LaTeX source in `doc/`) |
| `stata.toc`, `hmdbivpois.pkg` | installation files for `net install` / `net get` |

## Example data

`hmdbivpois_example.dta` is simulated from the model with 2,000 individuals, a presample period and up to 6 estimation periods, three states, two classes, the Famoye form and homogeneous transitions. The variables `state` and `class` hold the true latent variables. True values:

| Block | Parameters |
|---|---|
| Equation y1 | x11 = 0.4, x12 = −0.3, lag y1 = 0.03, lag y2 = −0.05, state intercepts (−0.3, 0.7, 1.5), class 2 shift 0.3 |
| Equation y2 | x21 = 0.3, x22 = −0.4, lag y1 = −0.10, lag y2 = 0.10, state intercepts (−2.5, −1.5, −0.6), class 2 shift 0.4 |
| Dependence | λ = −0.5 |
| Dropout (logit) | xd = 0.5, lag y1 = 0.05, lag y2 = 0.30, state intercepts (−4.2, −3.7, −3.2), class 2 shift 1.5 |
| Initial probabilities (states 2, 3 vs 1) | z1 = (0.5, 1.0), presample y1 = (0.15, 0.30), presample y2 = (0.30, 0.50), constant = (−0.5, −2.0) |
| Prior probability of class 2 (logit) | zu = 0.8, presample y1 = 0.1, presample y2 = 0.4, constant = −1.0 |
| Transition matrix | (0.85, 0.12, 0.03; 0.08, 0.84, 0.08; 0.03, 0.12, 0.85) |

## Citation

If you use `hmdbivpois`, please cite:

- Carrieri, V., D. Gyrd-Hansen, and P. Li Donni. 2026. Dynamic healthcare demand with informative dropout: evidence from register data. *Journal of the Royal Statistical Society Series A*. https://doi.org/10.1093/jrsssa/qnag108.
- Li Donni, P. 2026. hmdbivpois: Stata module to fit hidden Markov bivariate Poisson regression with informative dropout. https://github.com/plidonni/hmdbivpois.

## Development

The model and its extensions were specified by the author. Parts of the Stata and Mata code were improved with the assistance of an AI coding tool (Claude, Anthropic); the code was tested by the author with Monte Carlo simulations.

## License

MIT License, see `LICENSE`.

## Author

Paolo Li Donni, University of Palermo, [paolo.lidonni@unipa.it](mailto:paolo.lidonni@unipa.it), [personal webpage](https://sites.google.com/site/plidonni/). Please report bugs and suggestions as [issues](https://github.com/plidonni/hmdbivpois/issues).
