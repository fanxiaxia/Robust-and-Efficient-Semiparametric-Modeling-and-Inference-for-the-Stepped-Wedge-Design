# Simulation functions.

The following functions are used the conduct the simulation studies. Parameters shown are for Setting 2 (informative cluster sizes).

## Functions

### `libraries.R`

Loads all required R packages for the analysis pipeline.

**Required packages:**

-   `glmmTMB`: Generalized linear mixed models with Template Model Builder
-   `MASS`: General statistical functions
-   `Matrix`: Sparse and dense matrix operations
-   `stringr`: String manipulation
-   `dplyr`: Data manipulation
-   `tidyr`: Data tidying
-   `lme4`: Linear mixed-effects models
-   `geepack`: Generalized estimating equations
-   `doBy`: Groupwise operations
-   `fastDummies`: Dummy variable creation
-   `ggplot2`: Data visualization

### `effectmodels.R`

Provides four pre-configured effect curve models representing instant, lagged, curved (exponential), and piecewise-convex treatment effect patterns.

**Available models:**

-   `instant.model`: Immediate treatment effect
-   `lagged.model`: Delayed onset with 2-step lag
-   `curved.model`: Exponential growth pattern
-   `pconvex.model`: Piecewise-convex acceleration

### `clustersize.R`

**Function:** `generate_increasing_N_size()`

Generates a matrix of cluster sizes that vary across time and clusters, with some clusters experiencing growth over time while others remain stable.

**Parameters:**

-   `Ncluster`: Number of clusters
-   `Ntime`: Number of time periods
-   `basesize`: Starting cluster size (default: 10)

**Returns:** A matrix (Ncluster × Ntime) where each row represents a cluster and each column represents a time period. Half of the clusters show linear growth over time, and clusters are grouped into quartiles with different baseline size adjustments.

### `rdclusters.R`

**Function:** `randomize_clusters()`

Randomly assigns clusters to different crossover times in a stepped-wedge design, ensuring balanced allocation across available time points.

**Parameters:**

-   `Ncluster`: Number of clusters to randomize
-   `Time.vec`: Vector of available crossover time points

**Returns:**

A data frame with two columns:

-   `Cluster`: Cluster ID (1 to Ncluster)
-   `cross`: Assigned crossover time for each cluster

The function ensures approximately equal numbers of clusters cross over at each time point by dividing clusters evenly across the available times.

### `datagen_final.R`

**Function:** `generate_dataset()`

Generates simulated data for stepped-wedge cluster randomized trials with customizable treatment effect accumulation patterns, time trends, and cluster-level variation.

**Parameters:**

-   `mu`: Baseline mean outcome
-   `tau`: Standard deviation of random cluster effects
-   `theta`: treatment effect
-   `n_clusters`: Number of clusters in the trial
-   `n_time_points`: Number of time periods
-   `n_ind_per_cluster`: Matrix of individuals per cluster at each time point
-   `sigma`: Residual standard deviation
-   `delay_model`: Effect curve model (use pre-configured models from `effectmodels.R`)
-   `time_trend`: Temporal trend pattern ("none" or "incr" for increasing)
-   `rdtmsd`: Standard deviation of random time effects
-   `same.trend`: Logical; if TRUE, all clusters follow the same time trend
-   `time_size_int`: Logical; if TRUE, includes interaction between time and baseline cluster size

**Returns:**

A list containing:

-   `params`: Trial design parameters (n_clusters, n_time_points, crossover_times)
-   `beta_js`: Time trend coefficients
-   `theta_ls`: Treatment effect at each time step post-intervention
-   `data`: Data frame with columns for cluster (i), time (j), individual (k), time since intervention (l), treatment indicator (x_ij), outcome (y), and cluster sizes


### `est_gamma_final.R`

**Functions:** `est_gammaW_crc_adj()`, `est_gammaW_ctime_adj()`

#### `est_gammaW_crc_adj()`

Estimates time effects (gamma) and cluster-specific weight matrices using generalized linear mixed models with random cluster intercepts only.

**Parameters:**

- `deltahat`: Estimated treatment effect to be removed from outcomes
- `dataset`: Data frame containing cluster, time, individual, treatment (x), and outcome (y) variables
- `gamma.type`: Type of time effect specification ("cat" for categorical, "cont" for continuous linear, "zero" for no time trend)
- `adj.var`: Optional adjustment variables to include in the model (default: NULL)

**Returns:**
A list containing:

- `gammahat`: Data frame with predicted time effects for each observation (cluster, time, individual, gamma)
- `W`: List of inverse variance-covariance matrices, one per cluster
- `model_used`: Description of the fitted model structure
- `gammafit`: The fitted glmmTMB model object

#### `est_gammaW_ctime_adj()`

Estimates time effects (gamma) and cluster-specific weight matrices using generalized linear mixed models with random cluster intercepts and random time slopes.

**Parameters:**

- `deltahat`: Estimated treatment effect to be removed from outcomes
- `dataset`: Data frame containing cluster, time, individual, treatment (x), and outcome (y) variables
- `gamma.type`: Type of time effect specification ("cat" for categorical, "cont" for continuous linear, "zero" for no time trend)
- `adj.var`: Optional adjustment variables to include in the model (default: NULL)

**Returns:**
A list containing:

- `gammahat`: Data frame with predicted time effects for each observation (cluster, time, individual, gamma)
- `W`: List of inverse variance-covariance matrices, one per cluster
- `model_used`: Description of the fitted model structure
- `gammafit`: The fitted glmmTMB model object

### `est_prop_final.R`

**Functions:** `propx.cal()`, `West()`, `est_one_adj()`, `est_crc_newW()`, `est_ctime_newW()`

#### `propx.cal()`

Calculates treatment propensity scores and centered treatment indicators for $\tilde{L}_i$.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment (x), and outcome variables
- `ex.type`: Type of propensity score estimation ("cluster emp" for cluster-level, "individual emp" for individual-level)
- `adj.var.prop`: Formula string specifying adjustment variables for propensity model

**Returns:**
Dataset with added columns:

- `meanx`: Predicted propensity scores
- `L`: Centered treatment indicator (x - meanx)

#### `West()`

Constructs exchangeable correlation weight matrices for each cluster using a specified intracluster correlation coefficient.

**Parameters:**

- `dataset`: Data frame containing cluster identifiers
- `rhoest`: Estimated intracluster correlation coefficient

**Returns:**
A list of weight matrices (W), one per cluster, based on the exchangeable correlation structure.

#### `est_one_adj()`

Performs a single iteration of treatment effect estimation using the proposed estimator.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment, and outcome variables
- `gamma.df`: Data frame of estimated time effects
- `W.list`: List of cluster-specific weight matrices
- `ex.type`: Type of propensity score estimation (default: "individual emp")
- `adj.var.prop`: Formula string specifying adjustment variables for propensity model

**Returns:**
Single numeric value representing the estimated treatment effect.

#### `est_crc_newW()`

Iteratively estimates treatment effects (exchangeable working correlation), updating time effects and weight matrices at each iteration.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment, and outcome variables
- `delta.ini`: Initial value for treatment effect estimate
- `it.num`: Number of iterations to perform
- `ex.type`: Type of propensity score estimation (default: "individual emp")
- `gamma.type`: Type of time effect specification (default: "zero")
- `adj.var.out`: Adjustment variables for outcome model
- `adj.var.prop`: Adjustment variables for propensity model

**Returns:**
A list containing:

- `est`: Final treatment effect estimate
- `gamma`: Final estimated time effects
- `W`: Final cluster-specific weight matrices
- `gammafit`: Final fitted model object

#### `est_ctime_newW()`

Iteratively estimates treatment effects (working correlation incorporating random time effect), updating time effects and weight matrices at each iteration.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment, and outcome variables
- `delta.ini`: Initial value for treatment effect estimate
- `it.num`: Number of iterations to perform
- `ex.type`: Type of propensity score estimation (default: "individual emp")
- `gamma.type`: Type of time effect specification (default: "zero")
- `adj.var.out`: Adjustment variables for outcome model
- `adj.var.prop`: Adjustment variables for propensity model

**Returns:**
A list containing:

- `est`: Final treatment effect estimate
- `gamma`: Final estimated time effects
- `W`: Final cluster-specific weight matrices
- `gammafit`: Final fitted model object

### `seperm3_final.R`

**Function:** `se.perm3()`

Calculates permutation-based standard errors for treatment effect estimates in stepped-wedge cluster randomized trials.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment (x), and outcome (y) variables
- `deltahat`: Estimated treatment effect
- `gammahat`: Data frame of estimated time effects
- `W.list`: List of cluster-specific weight matrices
- `diffsize`: Indicator for differential cluster sizes
- `adj.var.prop`: Adjustment variables for propensity model
- `crosstime`: Vector of possible crossover times for permutation

**Returns:**
A numeric value representing the permutation-based standard error of the treatment effect estimate. The function performs `Nperm` permutations of cluster crossover times, recalculates L for each permutation, and computes variance components both within and between clusters.

### `seperm3loo_final.R`

**Function:** `se.perm3.loo()`

Calculates leave-one-out permutation-based standard errors for treatment effect estimates.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment (x), and outcome (y) variables
- `deltahat.list`: List of treatment effect estimates, one per cluster (from leave-one-out analysis)
- `gammahat.list`: List of time effect data frames, one per cluster (from leave-one-out analysis)
- `W.list`: List of cluster-specific weight matrices
- `diffsize`: Indicator for differential cluster sizes
- `adj.var.prop`: Adjustment variables for propensity model
- `crosstime`: Vector of possible crossover times for permutation

**Returns:**
A numeric value representing the leave-one-out permutation-based standard error. Similar to `se.perm3()` but uses cluster-specific estimates from a leave-one-out procedure, where each cluster's residuals are based on estimates derived from all other clusters.

### `sesand2_final.R`

**Function:** `se.sand2()`

Calculates sandwich-type robust standard errors for treatment effect estimates.

**Parameters:**

- `dataset`: Data frame with cluster, time, individual, treatment (x), and outcome (y) variables
- `deltahat`: Estimated treatment effect
- `gammahat`: Data frame of estimated time effects
- `W.list`: List of cluster-specific weight matrices
- `diffsize`: Indicator for differential cluster sizes
- `adj.var.prop`: Adjustment variables for propensity model

**Returns:**
A numeric value representing the sandwich-type robust standard error of the treatment effect estimate. Uses cluster-level summation of weighted residuals to construct the variance estimator with empirical sandwich covariance structure.

### `ci.R`

**Function:** `calculate_ci()`

Calculates confidence intervals for treatment effect estimates using normal approximation.

**Parameters:**

- `est`: Point estimate of the treatment effect
- `se`: Standard error of the estimate
- `conf_level`: Confidence level (default: 0.95 for 95% CI)

**Returns:**
A named vector with two elements:

- `lower`: Lower bound of the confidence interval
- `upper`: Upper bound of the confidence interval

### `doone_final.R`

**Function:** `do.one()`

Runs a single simulation trial comparing multiple estimation methods for stepped-wedge cluster randomized trials, including linear mixed models (LMM) and the proposed estimator with various variance estimators.

**Parameters:**

- `trials`: Simulation trial number (used as seed for reproducibility)

**Returns:**
A list containing:

- `cor.list`: Correlation between cluster size and crossover time
- `lmm_bias_s`: Bias for simple LMM models (4 variants: continuous/categorical time with random intercept/slope)
- `lmmcover_s`: Coverage indicators for simple LMM confidence intervals
- `lmm_bias_info`: Bias for LMM models adjusted for baseline cluster size
- `lmmcover_info`: Coverage indicators for adjusted LMM confidence intervals
- `prop_bias_simple`: Bias for proposed method without double baseline cluster size adjustment 
- `prop_bias_info`: Bias for proposed method with double baseline cluster size adjustment (4 variants)
- `prop_coverp_info`: Coverage using permutation-based standard errors
- `prop_covers_info`: Coverage using sandwich-type standard errors
- `coverp.loo`: Coverage using leave-one-out permutation standard errors (if `loo = TRUE`)

The function generates synthetic data with informative cluster sizes, fits multiple models, and evaluates their performance in terms of bias and confidence interval coverage.

### `main_info.R`

Main simulation script that orchestrates a Monte Carlo study comparing estimation methods for stepped-wedge cluster randomized trials with informative cluster sizes.

**Simulation Parameters:**

- `Ncluster`: Number of clusters (10)
- `Ntime`: Number of time periods (5)
- `mu_size`: Baseline mean outcome (3)
- `effectsize`: Maximum treatment effect (4)
- `rdcluster.var`: Variance of random cluster effects (0.25)
- `rderr.var`: Residual variance (4)
- `rdtimevar`: Variance of random time effects (0.25)

**Estimation Settings:**

- `it.num`: Number of iterations for the proposed method (10)
- `Nperm`: Number of permutations for variance estimation (20)
- `ex.type.set`: Propensity score type ("individual emp")
- `adjusted.var.out`: Adjustment variable for outcome model ("basesize")
- `adjusted.var.prop`: Adjustment variables for propensity model ("factor(basesize)*factor(time)")
- `est`: Methods to evaluate (LMM, proposed methods)
- `loo`: Whether to compute leave-one-out standard errors (TRUE)

**Execution:**

Runs 1,000 simulation trials in parallel using all available CPU cores and saves results to `reslist_final.RData`. Each trial calls `do.one()` to generate data and compare estimation methods.