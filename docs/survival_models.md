# Survival Models

## Elastic Net Cox

Elastic Net Cox was used as the main model because it is suitable for high-dimensional gene-expression data and performs regularized feature selection. It combines Ridge-like shrinkage and LASSO-like sparsity.

The final selected setting in the project was alpha = 0.1. The frozen full-data model selected 284 genes and was used for external validation.

## Random Survival Forest

Random Survival Forest was used as a non-linear survival benchmark. It can capture non-linear effects and does not require the proportional hazards assumption.

Because RSF is computationally expensive with thousands of genes, an in-fold univariate Cox filter was applied inside each training fold. Up to 500 genes were selected using only training data.

## DeepSurv

DeepSurv was used as a neural survival benchmark. It learns a non-linear risk score from expression features using a Cox partial-likelihood objective.

The model architecture used here was a feed-forward neural network with two hidden layers, batch normalization, dropout, and early stopping.

## Final interpretation

Elastic Net Cox remained the preferred model because it was sparse, interpretable, competitive with DeepSurv, stronger than RSF on the same split, and externally validated on SCAN-B.
