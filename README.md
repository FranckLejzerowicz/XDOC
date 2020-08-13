# XDOC

Run and explore to results of a Dissimilarity Overlap Curve analysis, by implementing 
the R code from [Jakob Russel](https://github.com/Russel88/DOC).

## Installation

```
pip install -U git+https://github.com/FranckLejzerowicz/XDOC.git
```

## Input
    
* Data files:
  - `-i`: a feature table  with samples as columns and features as rows.
  - `-o`: a folder name (will be created)
  - `-m`: a sample metadata table for subsetting  

* Filtering and metadata-based subsetting:
  - `-c`: metadata columns to filter on.
  - `-v`: value in this column to filter on.
  - `-q`: one 0 - 100 value to use as quantile filter on the column (must be numerical).
  - `-fa`: feature abundance filter (i.e. features must have >= the given abundance in a sample to be kept 
  in this sample). If in range 0-1 used as min or equal percent of abundance in the sample, if >1 used as min 
  or equal number of reads in the sample.
  - `-fp`: feature prevalence filter at the `-fa` abundance (i.e. are kept features have >= `-fa` in `-fp` samples). 
  If in range 0-1 used as min or equal percent of samples, if >1 used as min or equal number of samples.   
  - `-f`: whether to first apply the metadata-based subsetting and the abundance/prevalence filters in the samples
  that remaing after potential dropouts, or first subset and then apply the filters.     
  
* bootstrapping :
  - `-r`: number of bootstrap draws.
  - `-subr`: 
  - `-pair`: 
  
* Loess regression and parameters:
  - `-mov_avg`: 
  - `-span`: 
  - `-degree`: 
  - `-family`: 
  - `-iterations`: 
  - `-surface`: 
  - `-ci`: 
  For more explanation, see [scikit-misc](https://has2k1.github.io/scikit-misc/generated/skmisc.loess.loess.html#skmisc.loess.loess)
  documentation, namely for [model](https://has2k1.github.io/scikit-misc/generated/skmisc.loess.loess_model.html#skmisc.loess.loess_model)
  and smoothing [control]([model](https://has2k1.github.io/scikit-misc/generated/skmisc.loess.loess_control.html#skmisc.loess.loess_control)). 

* Other parameters:
  `-cpus`: number of CPUs to use as both the pairwise DO calculations and bootstrapping are multipriocessed. 

* Other parameters:
  `--null`: whether to run null model(s).
  `-nulls`: the number of null models to run.
  `--non-zero`: only randomize the non-zero entries of each feature.

## Outputs



## Example


Running:

    ```
    ```
Would return file ``.

## Visualization

### Optional arguments

```
  -i, --i-otu TEXT                Input biom table.  [required]
  -o, --o-outdir TEXT             Output directory path.  [required]
  -m, --m-metadata TEXT           Metadata table.
  -c, --p-column TEXT             Column from metadata `-m` to use for
                                  filtering based on values of `-v`.

  -v, --p-column-value TEXT       Filtering value to select samples based on
                                  column passed to `-c`.

  -q, --p-column-quant INTEGER    Filtering quantile / percentile for samples
                                  based on column passed to `-c` (must be
                                  between 0 and 100).

  -fp, --p-filter-prevalence FLOAT
                                  Filter features based on their minimum
                                  sample prevalence (number >1 for sample
                                  counts: <1 for samples fraction).

  -fa, --p-filter-abundance FLOAT
                                  Filter features based on their minimum
                                  sample abundance (number >1 for abundance
                                  counts: <1 for abundance fraction).

  -f, --p-filter-order [meta-filter|filter-meta]
                                  Order to apply the filters: 'filter-meta'
                                  first the prevalence/abundance and then
                                  based on variable; 'meta-filter' first based
                                  on variable and then the
                                  prevalence/abundance on the remaining.
                                  [default: meta-filter]

  -r, --p-r INTEGER               Number of bootstraps.  [default: 100]
  -subr, --p-subr INTEGER         If ignored will do bootstrap, alternatively
                                  an integer denoting size of subsample.
                                  [default: 0]

  -pair, --p-pair TEXT            One of the two names pair matching samples
                                  names.

  -mov_avg, --p-mov-avg INTEGER   Moving average window to use for estimating
                                  where negative slope starts.  [default: 5]

  -ci, --p-ci FLOAT               Quantiles for confidence intervals.
                                  [default: 0.025, 0.5, 0.975]

  -span, --p-span FLOAT           Span of loess smoothing.  [default: 0.2]
  -degree, --p-degree FLOAT       Degree of loess smoothing (If 1 linear, >1
                                  polynomial).  [default: 1.0]

  -family, --p-family [gaussian|symmetric]
                                  'gaussian' for least-squares fitting,
                                  'symmetric' for robust fitting.  [default:
                                  symmetric]

  -iterations, --p-iterations INTEGER
                                  Number of iterations for robust fitting.
                                  [default: 4]

  -surface, --p-surface [direct|interpolate]
                                  'direct' estimation (slow exact) or
                                  'interpolate' estimation (fast approximate).
                                  [default: direct]

  -cpus, --p-cpus INTEGER         Number of cores to use.  [default: 1]
  -nulls, --p-nulls INTEGER       Number of null models.  [default: 1]
  --non-zero / --no-non-zero      Only shuffle non zero values for null model.
                                  [default: True]

  --null / --no-null              Perform DOC on null models.  [default:
                                  False]

  --verbose / --no-verbose
  --version                       Show the version and exit.
  --help                          Show this message and exit.

```

### Bug Reports

contact `flejzerowicz@health.ucsd.edu`