# eDITH 1.0.0.9000

- Functions `run_eDITH_optim`, `run_eDITH_optim_joint` have initial parameter set redefined.
- Function `eval_posterior_eDITH` exports quantiles and mean of the model parameters.

# eDITH 1.0.0

## Major changes

- Functions `run_eDITH_BT_joint`, `run_eDITH_optim_joint` added.
- Functions `sampling_strategy_eDNA`, `sampling_strategy_direct` added.
- Dataset `dataCD` added.
- Dependency to `rivnet` >= 0.4.2 updated.

# eDITH 0.3.0

## Major changes

- `run_eDITH_optim`: argument `n.restarts` is added.

## Minor changes

- `run_eDITH_optim`: `attempts.stats` is exported.
- `CITATION` added.

# eDITH 0.2.0

## Major changes

- `run_eDITH_optim`: log-posterior is now maximized. Default arguments for prior distributions
have been added.

## Minor changes

- `run_eDITH_BT`, `run_eDITH_optim`: updated error messages when unsuitable input is used.
- Vignette updated.
- BugReports link added in `DESCRIPTION`.

## Bugs fixed

- `run_eDITH_BT`, `run_eDITH_optim`: bug fixed in attribution of names to AEM covariates.

# eDITH 0.1.2

## Minor changes

- Vignette updated.

# eDITH 0.1.1

## Bugs fixed

- Bug fixed in the calculation of the detection probability when `ll.type = geom`.

