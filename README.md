# ATACseq-preprocessing-ArchR
Pre-processing spatial ATAC-seq data using ArchR

## Install ArchR and required packages
<code> if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools") </code>

<code> if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") </code>

<code> devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories()) </code>

<code> library(ArchR) </code>

<code> ArchR::installExtraPackages() </code>

## Data
The data used in this project is from Dr. Rong Fan's lab at Yale University and can be found in the Gene Expression Omnibus (GEO) with accession code GSE171943 (reviewer token: itqreaggtxitnup).

**More information can be found here:**
<br>
https://www.biorxiv.org/content/10.1101/2021.06.06.447244v1
