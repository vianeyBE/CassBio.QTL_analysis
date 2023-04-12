# Cassava Bioinformatic Platform: QTL analysis on bi-parental families

Pipeline for QTL analysis for Cassava using Rqtl R package

The pipeline has four main steps:

1. Single QTL mapping
2. MQM
3. Plots (Boxplots and Manhattan)
4. Annotation

## 1. Single QTL Mapping

This code performs QTL (Quantitative Trait Locus) analysis using Rqtl package. It computes single QTLs, LOD/Bayes intervals, scores, and effects using a mapqtl input format.

### Usage
```R

single_qtl(dir, dircross, dirfun, locfile, mapfile, phenofile, prefixResults, ncores, step = 0.5, off.end = 0, error.prob = 0.001, alpha = 0.1, n.perm = 1000, map.function = "kosambi", stepwidth = "fixed", model_scanone = "normal")

```

### Arguments
1. Single QTL function:
- `dir`: Directory (dir) to save results
- `dircross`: Directory where data is located. It can be empty if the cross files include the path and name
- `dirfun`: Directory were functions for plotting are located
- `locfile`: .loc file name. It can include the path
- `mapfile`: .map file name. It can include the path
- `phenofile`: .qua file name. It can include the path
- `prefixResults`: The prefix that the results will have
- `ncores`: Number of cores for permutation test

2. scanone function:
- `step`: Step size in cM used for interval mapping (default = 0.5)
- `off.end`: A value added to each end of the chromosome, to extend the genetic map. (default = 0)
- `error.prob`: The probability of a genotyping error. (default = 0.001)
- `alpha`: The significance level to use for peak detection (default = 0.1)
- `n.perm`: Number of permutations to use for significance testing. (default = 1000)
- `map.function`: The genetic map function to use. (default = "kosambi")
- `stepwidth`: Method used to compute step size. Can be "fixed" or "cov" (default = "fixed")
- `model_scanone`: Model used in the single-QTL scan. Can be "normal", "binary" or "poisson" (default = "normal")

### Details
The data is read using the ``read.cross()`` function. The cross data is in mapqtl format, which includes .map, .loc, and .qua files with the genetic map, markers, and phenotypes information, respectively.
The `calc.genoprob()` function is used to calculate the probabilities of each individual's genotype given their observed phenotype data.
The `scanone()` function is used to perform interval mapping using the Expectation-Maximization (EM) algorithm, which estimates the genetic effects of the markers in the cross.
A permutation test is conducted using `scanone()` to compute the genome-wide significance threshold.
Peaks are identified as locations in the genome where the genetic effect of the markers is above the genome-wide significance threshold.

### Output
*In progess*

### Dependencies
- `tidyverse`
- `snow`
- `qtl`

---

## Contact
For questions or feedback about this pipeline, please contact Vianey Barrera-Enriquez at v.barrera@cgiar.org.