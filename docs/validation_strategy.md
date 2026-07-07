# Validation Strategy

## Internal train/test split

The integrated dataset was split into training and held-out test patients. Elastic Net Cox, DeepSurv, and Random Survival Forest were evaluated on the same split for fair comparison.

## Leave-one-cohort-out validation

Leave-one-cohort-out validation was used to test cross-cohort generalization. In each fold, one complete cohort was held out as the test set and the model was trained on the remaining cohorts.

This avoids overestimating performance from random splits where samples from the same cohort can appear in both train and test sets.

## External validation

The frozen Elastic Net Cox model was externally validated on SCAN-B / GSE96058. The model was not retrained on SCAN-B. This tested both external cohort generalization and cross-platform transfer from microarray discovery data to RNA-seq validation data.
