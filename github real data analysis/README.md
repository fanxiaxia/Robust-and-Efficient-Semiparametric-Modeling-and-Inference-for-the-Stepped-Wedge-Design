# Data Analysis Functions

The following functions are used to analyze real data from Washington EPT trial. The analysis compares linear mixed models (LMM) with the proposed method.

## Functions

### `[1]preprocess.R`

Preprocesses raw trial data by creating time waves, handling missing values, calculating cluster sizes, and computing crossover times.

**Data transformations:**

- Aggregates continuous time (`timex`) into 5 discrete time waves (`timewave`)
- Removes observations with missing values
- Creates individual-level identifiers (`ind`)
- Calculates cluster size at each time point (`size`) and baseline cluster size (`basesize`)
- Identifies crossover times (first time point where treatment switches from 0 to 1)

**Output variables:**

- `dataset`: Preprocessed data frame with columns for cluster, time, treatment (x), outcome (y), individual ID, cluster size, and baseline size
- `diffsize`: Matrix of cluster sizes (clusters × time periods)
- `df_crosstime`: Data frame with crossover times and baseline sizes for each cluster
- Correlation between crossover time and baseline cluster size

### `[2]ept_main.R`

Main analysis script that fits both LMM and the proposed method without adjusting for baseline cluster size.

**Analysis Settings:**

- `it.num`: Number of iterations for the proposed method (10)
- `Nperm`: Number of permutations for variance estimation (50)
- `ex.type.set`: Propensity score type ("individual emp")
- `adjusted.var.out`: Adjustment for outcome model (NULL - no adjustment)
- `adjusted.var.prop`: Adjustment for propensity model ("factor(time)")
- `loo`: Whether to compute leave-one-out standard errors (TRUE)

**Models fitted:**

1. LMM with random cluster intercepts: `y ~ x + factor(time) + (1|cluster)`
2. LMM with random intercepts and time slopes: `y ~ x + factor(time) + (1 + factor(time)|cluster)`
3. Proposed method with permutation-based standard errors
4. Proposed method with leave-one-out permutation standard errors

**Output:**

- `res_simple.csv`: Results for proposed method (estimates, SE, confidence intervals)
- `res_lmm_simple.csv`: Results for LMM (estimates, SE, confidence intervals)

### `[3]ept_basesize.R`

Main analysis script that fits both LMM and the proposed method while adjusting for baseline cluster size (informative cluster size analysis).

**Analysis Settings:**

- `it.num`: Number of iterations for the proposed method (10)
- `Nperm`: Number of permutations for variance estimation (50)
- `ex.type.set`: Propensity score type ("individual emp")
- `adjusted.var.out`: Adjustment for outcome model ("log(basesize)")
- `adjusted.var.prop`: Adjustment for propensity model ("log(basesize)*factor(time)")
- `loo`: Whether to compute leave-one-out standard errors (TRUE)

**Models fitted:**

1. LMM with baseline size adjustment: `y ~ x + factor(time) + log(basesize) + (1|cluster)`
2. LMM with baseline size and random slopes: `y ~ x + factor(time) + log(basesize) + (1 + factor(time)|cluster)`
3. Proposed method with baseline size adjustment and permutation-based SE
4. Proposed method with baseline size adjustment and leave-one-out permutation SE

**Output:**

- `res_basesize.csv`: Results for proposed method with adjustment (estimates, SE, confidence intervals)
- `res_lmm_basesize.csv`: Results for adjusted LMM (estimates, SE, confidence intervals)


## Workflow

1. Run `[1]preprocess.R` to prepare the data
2. Run `[2]ept_main.R` for unadjusted analysis
3. Run `[3]ept_basesize.R` for analysis adjusting for informative cluster sizes
