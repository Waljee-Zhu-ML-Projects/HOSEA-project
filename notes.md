# Notes

Last updated: **September 23rd, 2021**

## Ideas

### Survival Analysis

- Is there a reason why we are not treating this as a survival analysis problem?
- Prediction horizon

### Feature engineering

- Why no mean difference? Instead of TV. TV shows just the amount of variability while mean difference would show the trend

### Imbalance

- Not sure if train/valid/test is done properly:
  - before or after downsampling?
- Different subsets every iteration?
- Stratified sampling / cluster controls into same amount of observations and take medioid

### Imputation

- I am not 100% convinced imputation is done correctly, i.e. within each dataset:
  - Shouldn't that be part of the model?
  - Then, you would resample from the training set
  - For example, what happens if you want to test on a single new observation (which is the use case here?)
  - The regression one probably don't work for this reason because the test set is pretty small
- Should it be done stratified? like estimating $p(x\mid 1)$ and $p(x\mid 0)$ or just $p(x)$, or even $p(x, y)$?
- Reduced model approach:
  - build models for each missingness pattern? 
  - need to check how many patterns there are

### `xgboost`

- Could use `xgb_model` argument to sequentially impute rather than parallel impute
  - I think this could prevent more overfitting
  - Would need a more efficient imputation implementation
- Could supply the AUC directly with `feval` ?

### Charlson scores imputation

- 