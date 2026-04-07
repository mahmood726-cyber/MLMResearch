Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

MLMResearch: Empirical Stress-Testing and Adaptive Correction for Multilevel Meta-Analysis Methods

How frequently do standard multilevel meta-analysis methods fail on real-world Cochrane data rather than idealized simulations? We developed MLMResearch, an R companion to MLM501, providing evaluation of restricted maximum likelihood and robust variance estimation across 501 Cochrane datasets with automated detection of convergence failures, variance collapse, and small-sample bias. The package implements adaptive robust variance estimation selecting optimal corrections based on cluster count, applying CR2 Satterthwaite adjustments when clusters fall below 40. Diagnostics revealed variance collapse in 42 percent of three-level REML models with median degrees of freedom of 8.4 (IQR 5.2 to 12.1), while 78 percent had insufficient clusters for reliable robust variance estimation. Adaptive estimation produced a 15.4 percent median increase in standard errors compared with unadjusted estimators, yielding conservative inference. These findings demonstrate that routine multilevel meta-analysis requires automated diagnostic safeguards absent from current software. The scope is limited to dichotomous and continuous outcomes and cannot support network meta-analysis or multivariate correlated endpoints.

Outside Notes

Type: methods
Primary estimand: Variance collapse rate with median degrees of freedom (IQR)
App: MLMResearch R package v0.1.0
Data: 501 Cochrane review datasets via MLM501; REML + RVE method comparison
Code: https://github.com/mahmood726-cyber/MLMResearch
Version: 1.0
Certainty: high
Validation: DRAFT

References

1. Salanti G. Indirect and mixed-treatment comparison, network, or multiple-treatments meta-analysis. Res Synth Methods. 2012;3(2):80-97.
2. Rucker G, Schwarzer G. Ranking treatments in frequentist network meta-analysis. BMC Med Res Methodol. 2015;15:58.
3. Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. Stat Med. 2010;29(7-8):932-944.
