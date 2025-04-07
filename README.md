# Cassava Bioinformatics Platform: QTL analysis pipeline 🧬

A modular, end-to-end pipeline to conduct **QTL (Quantitative Trait Loci)** mapping for cassava bi-parental populations. It performs SNP filtering, encoding, QTL single-marker analysis, high-resolution plotting, functional annotation, and validation through SNP-to-trait boxplots.

Designed with reproducibility and modularity in mind, this pipeline is implemented in `Python` and `R`, making it adaptable to different breeding projects.





## Authors 🙋

For questions or feedback about this pipeline, please contact:

- **Camilo E. Sánchez-Sarria** - [CGIAR]: c.e.sanchez@cgiar.org
- **Vianey Barrera-Enriquez** - [CGIAR]: v.barrera@cgiar.org





## Workflow overview 🚀

**| 1 | Quality control** 🧹 Filters the raw VCF file and converts it into a mapping-ready format

**| 2 | Encoding SNPs** 📉 Encodes bi-allelic genotypes for F1 QTL mapping

**| 3 | QTL analysis** ✂️ Performs single QTL scan using R/qtl

**| 4 | LOD profiles** 🧭 Plots LOD profiles from QTL scans

**| 5 | QTL regions** 📊 Visualizes converged QTL regions in the genome

**| 6 | QTL density** 🧮 Visualizes QTL density across chromosomes

**| 7 | Functional annotation** 🧾 Annotates QTL regions using upstream/downstream gene proximity

**| 8 | Marker validation boxplots** 📈 Validates candidate QTLs by SNP-to-trait boxplots





## Module descriptions 🧩

Each module contains a README.md describing:

- **Description**: A brief explanation of what the module does
- **Input arguments**: Parameters or flags used when running the module
- **Example usage**: Command-line or script usage example
- **Example output structure**: Example of the outputs produced by the module
- **Dependencies**: Required packages, tools, or environments





## Installation and dependencies 🛠️

To clone the repository and ensure reproducibility:

```bash

git clone https://github.com/vianeyBE/CassBio.QTL_analysis.git
cd CassBio.QTL_analysis

# OPTIONAL: Set up a conda environment
conda env create -f environment.yaml
conda activate cassbio-qtl

```

This pipeline relies on the following tools and packages. Make sure they are installed and accessible in your environment (or use the `environment.yaml` file to create a reproducible setup)

**Core tools**
- `R ≥ 4.0`: Required for QTL mapping, plotting, and annotation
- `Python ≥ 3.7`: Required for preprocessing and encoding

**R packages**
- `qtl`
- `tidyverse`
- `ggplot2`
- `ggsignif`
- `openxlsx`
- `gridExtra`
- `scales`

Install manually with

``` R

install.packages(c("tidyverse", "ggplot2", "ggsignif", "openxlsx", "gridExtra", "scales"))
if (!require("qtl")) install.packages("qtl")

```

**Python packages**
- `pandas`
- `numpy`
- `argparse`

Install via pip or Conda

``` bash

pip install pandas numpy argparse

```

**Optional tools**
- `PLINK`: VCF filtering or format conversion (optional)
- `TASSEL`: VCF to HapMap conversion (optional for downstream compatibility)





## Usage ▶️

Each script is designed to be run independently, allowing flexibility to modify or skip modules based on available data.

**Example workflow:**

``` bash

# Step 1: VCF filtering
python 01_FilterVCF_for_Map.py --input raw_data.vcf --output filtered.vcf

# Step 2: Encode markers for mapping
python 02_Encoding_SNP-F1-population.py --vcf filtered.vcf --labels sample_labels.csv --output encoded.csv

# Step 3: QTL single scan
Rscript 03_QTL_SingleMapping.R

# Step 4: Plot QTL scans
Rscript 04_QTL_SinglePlots.R

# Step 5: Create QTL convergence map
Rscript 05_QTL_MapPlot.R

# Step 6: Plot QTL density
Rscript 06_DensityMapPlot.R

# Step 7: Annotate genes near QTL
Rscript 07_QTL_Annotation.R

# Step 8: Validate QTLs with SNP x Trait boxplots
Rscript 08_QTL_Boxplots.R

```

See the header of each script for required arguments and usage examples.





## 1. Single QTL Mapping

This code performs QTL (Quantitative Trait Locus) analysis using `Rqtl` package. It computes single QTLs, LOD/Bayes intervals, scores, and effects using a mapqtl input format.

### Usage

``` R

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

QTL_Annotation(Wdir, Ddir, name, wdyw, annot, gff, version, recursive)

```

### Arguments

- `Wdir`: Name of the directory that contains the GAPIT results. For example: home/user/folder.
- `Ddir`: Directory where is located the annotation files (annot, GFF files).
- `name`: Enter the path or the name of file names to look for. For example: QTL_LOD_Intervals.
- `wdyw`: Enter what are you looking for to annotate (Options: CDS, five_prime_UTR, gene, mRNA, three_prime_UTR).
- `annot`: Annotation details of the genes. txt file from the genome version used for alignment.
- `gff`: gff3 file from the genome version used for alignment.
- `version`: You can choose between the genome of reference version 6.1 or 8.1 (Options: 6.1 or 8.1).
- `recursive`: A Boolean string that determines if the function should perform a recursive search or not (Default = F).


## Details

*In progress*

### Examples

``` R

# Set arguments
Ddir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/00_Data/"
wdyw <- "gene"
annot <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.annotation_info.txt"
gff <- "Mesculenta_305_v6.1/Mesculenta_305_v6.1.gene.gff3"
version <- "6.1"

# If recursive
Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/Everything/"
name <- "LodIntervals"
recursive <- "T"

# Non-recursive
Wdir <- "D:/OneDrive - CGIAR/Cassava_Bioinformatics_Team/01_ACWP_F1_Metabolomics/02_QTL_Analysis/CM8996/Significant_ones/"
name <- "QTL_results_heritability.csv"
recursive <- "F"


# Run function
QTL_Annotation(Wdir, Ddir, name, wdyw, annot, gff, version, recursive)

```

### Output

A single CSV file containing relevant gene information plus QTLs' LOD, P-values, traits, and effects.

### Dependencies

- `tidyverse`
- `openxlsx`

## 5. Boxplot of significant markers: genotypes vs. phenotype

This R function generates a boxplot for a given SNP. The function takes as inputs a CSV file with the phenotype values and a list of SNPs. The function can also receive an optional CSV file to add extra information about the samples to the plot (categories, family, ect.) The output of the function is a PDF file with the plot.

### Usage

```R

QTL_Boxplot(outputname, dir, phenofile, genofile, snp_list_file, labelfile = NULL)

``` 

### Arguments
- `outputname`: (required) character string with the base name for the output file.
- `dir`: (required) character string with the directory where the output file will be saved.
- `phenofile`: (required) character string with the name of the phenotype file in tabular format. The first column should contain the sample names, and the rest of the columns should contain the phenotypes.
- `genofile`: (required) character string with the name of the genotype file in hapmap format.
- `snp_list_file`: (required) character string with the name of the CSV file with three columns:
    - Column 01 - Name: SNPS. List of SNPS to plot. The name should be the same as in the `genofile` data.
    - Column 02 - Name: trait. Name of the trait as in the `phenofile` data.
    - Column 03 - Name: xlabel. Name of the trait to be included as a label.
- `labelfile`: (optional) character string with the name of the CSV file with two columns:
    - Column 01 - Name: Taxa. Sample names.
    - Column 02 - Name: label. Label or category to add to the plot.

### Details

The function loads all the necessary packages and data files. It then checks for the presence of an optional file with sample labels and prepares the data for the boxplot. The function generates a boxplot for each SNP specified in the input CSV file. It uses `ggplot2` package to generate the plot, with the genotypes on the x-axis and the phenotype on the y-axis. The function adds a label to the x-axis to indicate the trait being plotted. If `labelfile` is provided, it adds extra information about the samples, coloring the data by the levels in the provided file.

### Example

``` R

# Generate boxplot without extra labels
QTL_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv")

# Generate boxplot with extra labels
QTL_Boxplot("outputname", ".path/to/save/plots/", "phenotype.csv", "genotype.hmp", "snp_list.csv", "labelfile.csv")

```

### Output

A single PDF file containing the boxplot of the SNPs.

### Dependencies

- `tidyverse`
- `tibble`
- `dplyr`
- `janitor`
- `ggplot2`
- `Biostrings`
- `hrbrthemes`
- `forcats`
- `ggsignif`
- `RColorBrewer`





## Citing this pipeline 📌

If you use this pipeline in your research or publication, please cite it as:

> Authors: Sánchez-Sarria, C. E., and Barrera-Enríquez, V  
> Title: *QTL analysis modular pipeline for Cassava: A reproducible and scalable workflow for pre-processing steps, QTL analysis and annotation in Cassava (Manihot esculenta)*
> GitHub repository: https://github.com/vianeyBE/CassBio.QTL_analysis
> Version: v1.0  
> Year: 2025

Alternatively, cite this repository using the following BibTeX entry:

``` bibtex

@misc{cassava_gwas_pipeline,
  author       = {Sánchez-Sarria, C. E., and Barrera-Enríquez, V},
  title        = {QTL analysis modular pipeline for Cassava},
  year         = 2025,
  version      = {v1.0},
  url          = {ttps://github.com/vianeyBE/CassBio.QTL_analysis},
  note         = {QTL analysis modular pipeline for Cassava}
}

```





## License 📄

This pipeline is released under the MIT License (see LICENSE file)
