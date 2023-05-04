# Cassava Bioinformatics Platform: QTL analysis Pipeline

Pipeline for QTL analysis on bi-parental families in Cassava

**Authors**: Vianey Barrera-Enriquez and Camilo E. Sanchez-Sarria

The pipeline has four main steps:

1. Single QTL mapping
2. Multiple QTL mapping (MQM) (*in progress*)
3. Customized plots (LOD profiles, Missing genotypes, etc)
4. Annotation of results

## 1. Single QTL Mapping

This code performs QTL (Quantitative Trait Locus) analysis using `Rqtl` package. It computes single QTLs, LOD/Bayes intervals, scores, and effects using a mapqtl input format.

### Usage

```R
single_qtl(dir, dircross, dirfun, locfile, mapfile, phenofile, prefixResults, ncores, step = 0.5, off.end = 0, error.prob = 0.001, alpha = 0.1, n.perm = 1000, map.function = "kosambi", stepwidth = "fixed", model_scanone = "normal")
```

### Arguments

1.1. For Single QTL function:
- `dir`: Directory to save results.
- `dircross`: Directory where data is located. It can be empty if the cross files include the path and name.
- `dirfun`: Directory were functions for plotting are located.
- `locfile`: .loc file name. It can include the path.
- `mapfile`: .map file name. It can include the path.
- `cphenofile`: .qua file name. It can include the path.
- `prefixResults`: The prefix that the results will have.
- `ncores`: Number of cores for permutation test.

1.2. For Single QTL function (`scanone` function part):
- `step`: Step size in cM used for interval mapping (default = 0.5).
- `off.end`: A value added to each end of the chromosome, to extend the genetic map (default = 0).
- `error.prob`: The probability of a genotyping error (default = 0.001).
- `alpha`: The significance level to use for peak detection (default = 0.1).
- `n.perm`: Number of permutations to use for significance testing (default = 1000).
- `map.function`: The genetic map function to use (default = "kosambi").
- `stepwidth`: Method used to compute step size. Can be "fixed" or "cov" (default = "fixed").
- `model_scanone`: Model used in the single-QTL scan (default = "normal". Options = "normal", "binary" or "poisson").

### Details

The data is read using the `read.cross()` function. The cross data is in mapqtl format, which includes .map, .loc, and .qua files with the genetic map, markers, and phenotypes information, respectively.
The `calc.genoprob()` function is used to calculate the probabilities of each individual's genotype given their observed phenotype data.
The `scanone()` function is used to perform interval mapping using the Expectation-Maximization (EM) algorithm, which estimates the genetic effects of the markers in the cross.
A permutation test is conducted using `scanone()` to compute the genome-wide significance threshold.
Peaks are identified as locations in the genome where the genetic effect of the markers is above the genome-wide significance threshold.

### Output

*In progess*

### Dependencies

- `qtl`
- `tidyverse`
- `snow`

## 2. Multiple QTL mapping (MQM)

*In progress*

## 3. Customized plots

*In progress*

## 4. Annotation of results

This code annotates the gen that overlaps and/or contain the significant QTLs obtained from the `rqtl` results.

### Usage

```R

QTL_Annotation(Wdir, Ddir, pat, wdyw)

```

### Arguments

- `Wdir`: Name of the directory that contains the GAPIT results. For example: home/user/folder.
- `Ddir`: Directory where is located the annotation files (annot, GFF files).
- `pat`: Enter the path of file names to look for. For example: QTL_LOD_Intervals. The path must finish with a point (.).
- `wdyw`: Enter what are you looking for to annotate. Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR.
- `annot`: Annotation details of the genes. txt file from the genome version used for alignment.
- `GFF`: gff3 file from the genome version used for alignment.
- `version`: (Optional) You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1. Default = 6.1).

## Details

*In progress*

### Examples

```R

# Set arguments
Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/"
Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
pat <- "QTL_LOD_Intervals."
wdyw <- "gene"

# Run function
QTL_Annotation(Wdir, Ddir, pat, wdyw)

```

### Output

A single CSV file containing relevant gene information plus QTLs' LOD, P-values, traits, and effects.

### Dependencies

- `tidyverse`

---

## Contact

For questions or feedback about this pipeline, please contact Vianey Barrera-Enriquez at v.barrera@cgiar.org.