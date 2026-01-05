```markdown
# HrdDetective

## Overview

**HrdDetective** is an R package for detecting **Homologous Recombination Deficiency (HRD)** from whole-genome sequencing data. It implements an algorithm based on the PELT algorithm and Sequenza statistical model to identify HRD in human cancer cells. The package provides a complete workflow from BAM file processing to copy number segmentation, ploidy/cellularity estimation, and HRD scar score calculation.

The package performs:

- **BAM to seqz conversion** – extracts heterozygous SNPs and computes depth/BAF
- **Copy number segmentation** – identifies genomic segments with uniform copy number using PELT algorithm
- **Model optimization** – estimates tumor purity (cellularity) and ploidy based on Sequenza statistical model
- **Visualization** – generates publication-ready plots of segments and model fits

## Installation

```r
# Install from GitHub (if available)
# devtools::install_github("limeng12/HrdDetective")

# For development installation:
# install.packages("devtools")
devtools::install_github("limeng12/HrdDetective")
```

## Dependencies

### Required R Packages
- **R (>= 4.0.0)**
- **Rcpp (>= 1.0.0)**
- **plyr, dplyr, iotools, readr, tidyr**
- **slider, reshape2, callr, copynumber, squash**
- **ggplot2 (>= 3.4.0), stringr**

### Bioconductor Packages
- **Rsamtools, GenomicAlignments, Biostrings**
- **GenomicRanges, rtracklayer, BiocParallel**

### Suggested Packages
- **testthat (>= 3.0.0), knitr, rmarkdown, BiocStyle**
- **pheatmap, gridExtra, cowplot**

### System Requirements
- **C++11 compiler**
- **samtools** (optional, for BAM indexing)

## Quick Start

```r
# Install packages
# Install from CRAN
install.packages(c(
  "Rcpp", "plyr", "dplyr", "iotools", "readr", "tidyr",
  "slider", "reshape2", "callr", "copynumber", "squash",
  "ggplot2", "stringr", "testthat", "knitr", "rmarkdown",
  "pheatmap", "gridExtra", "cowplot"
))

# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "Rsamtools", "GenomicAlignments", "Biostrings",
  "GenomicRanges", "rtracklayer", "BiocParallel",
  "BiocStyle"
))

# Load the package
library(HrdDetective)
library(BSgenome.Hsapiens.UCSC.hg38)

# Run the complete workflow
t_output_dir_name <- "./test_results"

# Step 1: Convert BAM to seqz (using example data)
bamfile_normal <- system.file("extdata", "target.chr1.normal.dedup.bam", package = "HrdDetective")
bamfile_tumor <- system.file("extdata", "target.chr1.tumor.dedup.bam", package = "HrdDetective")

bam2seqz_r_snps(
  normal_bam = bamfile_normal,
  tumor_bam = bamfile_tumor,
  bsgenome = BSgenome.Hsapiens.UCSC.hg38,
  output_file = "sample.snps.seqz.txt"
)

# Step 2: Process seqz file
data.file <- "sample.snps.seqz.txt"
data <- read.table(data.file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
filtered_data <- data[!grepl("^(M|chrM)", data[, 1]), ]
write.table(filtered_data, data.file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Step 3: Run optimization and segmentation
seg_obj<-opti_and_segs(data.file=data.file, min.tumor.reads = 10);

# Step 4: Generate output
html_file <- generate_png_hrd_report(
  seg.tab = seg_obj$seg.tab,
  seqz_list = seg_obj$seqz_list,
  celluPloidyMat = seg_obj$celluPloidyMat,
  best_ploidy = seg_obj$best_ploidy,
  best_cellularity = seg_obj$best_cellularity,
  output_dir_name = t_output_dir_name
)

library(scarHRD);

scar_result<-scar_score(paste0(t_output_dir_name, "/scarHRD_input_test.tsv"),
                      reference = "grch38",
                      seqz=FALSE, outputdir=t_output_dir_name);
```

## Output Files

The analysis generates the following output files:

| File | Description |
|------|-------------|
| `raw_segments.txt` | Copy number segments with A/B allele counts |
| `scarHRD_input_test.tsv` | Formatted input for HRD scoring |
| `HRD_re.txt` | Estimated cellularity, ploidy, and sample/tumor ploidy |
| `hrd_plot_full.pdf` | Per-chromosome BAF, depth ratio, and copy number plots |
| `hrd_plot_seg.pdf` | Genome-wide horizontal segment plot |
| `model_fit_contour.pdf` | Likelihood landscape for cellularity/ploidy optimization |
| Various HRD score files | HRD, LST, and TAI scores |

## Key Functions

### Core Functions
- `bam2seqz_r_snps()` – Converts paired normal/tumor BAM files to seqz format using heterozygous SNPs
- `opti_and_segs()` – Main function for copy number segmentation and model optimization
- `generate_png_hrd_report()` – Generates all output files and visualization plots

### Supporting Functions
- `pileup_whole_bam()` – Fast pileup extraction from BAM files (C++ backend)
- `pileup_from_bam()` – Region-specific pileup extraction
- `draw_model_fit()` – Contour plot of cellularity/ploidy likelihood

## Parameters

### bam2seqz_r_snps()
| Parameter | Default | Description |
|-----------|---------|-------------|
| `normal_bam` | Required | Path to normal (germline) sample BAM file |
| `tumor_bam` | Required | Path to tumor sample BAM file |
| `genome_fasta` | NULL | Path to reference genome FASTA file |
| `bsgenome` | NULL | BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38) |
| `min_depth` | 10 | Minimum read depth per position |
| `min_af` | 0.25 | Minimum allele frequency for heterozygous SNPs |
| `max_af` | 0.75 | Maximum allele frequency for heterozygous SNPs |
| `gc_window` | 50 | Window size for GC content calculation (±250bp) |
| `output_file` | "output.seqz" | Path to output seqz file |

### opti_and_segs()
| Parameter | Default | Description |
|-----------|---------|-------------|
| `data.file` | Required | Input seqz file path |
| `number_of_cores` | 6 | CPU cores for parallel optimization |
| `CNt.max` | 8 | Maximum copy number considered |
| `min.tumor.reads` | 20 | Minimum tumor reads for SNP inclusion |
| `cellularity` | seq(0.1, 1, by=0.01) | Cellularity search range |
| `ploidy` | seq(1, 6, by=0.1) | Ploidy search range |

## Method Details

### Algorithm
1. **Copy number segmentation** using the PELT (Pruned Exact Linear Time) algorithm
2. **Model optimization** based on the Sequenza statistical model for estimating tumor cellularity and ploidy

## Example Data

The package includes example BAM files for chromosome 1:
- `target.chr1.normal.dedup.bam` – Normal sample BAM
- `target.chr1.tumor.dedup.bam` – Tumor sample BAM

These files are located in the `extdata` directory and can be accessed using `system.file()`.

## Performance

- **Parallel processing** for model optimization
- **Efficient C++ backend** for pileup extraction
- **Memory-efficient** chromosome-wise processing
- **Scalable** for whole-genome analysis

## Troubleshooting

### Common Issues
1. **BAM file errors**: Ensure BAM files are sorted and indexed (.bai files present)
2. **Memory issues**: Process chromosomes individually for large genomes
3. **No heterozygous SNPs**: Adjust `min_af` and `max_af` parameters
4. **Compilation errors**: Ensure C++11 compiler is available

## Citation

## License

MIT License - see LICENSE file for details.

## Maintainer

Meng Li <limeng49631@aliyun.com>

## Acknowledgments

- Built upon methods from Sequenza and PELT algorithm
- Uses efficient C++ implementation for pileup extraction
- Inspired by existing HRD detection tools

## Support

For bugs, feature requests, or questions:
- Open an issue on GitHub
- Contact: Meng Li <limeng49631@aliyun.com>
```
