# Results Summary

## Discovery dataset

- 11 cohorts
- 3,818 patients
- 1,587 overall-survival events
- Top 5,000 MAD-filtered genes used for primary modelling

## Model comparison

Elastic Net Cox achieved an internal test C-index of 0.6513. DeepSurv achieved 0.6469 on the same split. Random Survival Forest achieved 0.6352.

## Leave-one-cohort-out validation

LOCO validation was used to test generalization to unseen cohorts. Performance varied by cohort, with the largest cohort METABRIC carrying substantial weight in event-weighted summaries.

## External SCAN-B validation

The frozen Elastic Net Cox model achieved an external C-index of 0.7163 on SCAN-B / GSE96058. Risk groups showed clear event-rate separation from low to high risk.

## Conclusion

The more complex models did not clearly outperform Elastic Net Cox. Elastic Net Cox was retained as the final model because it was interpretable, sparse, competitive internally, robust across cohort validation, and externally validated.
