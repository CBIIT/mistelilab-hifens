# mistelilab-hifens

This is a repo containing the R code used to generate some of the figures for a 2022 Nucleic Acid Research publication from the Misteli Lab at the National Cancer Institute/NIH.

-   [HiFENS: High-throughput FISH detection of endogenous pre-mRNA splicing isoforms](https://doi.org/10.1093/nar/gkac869)
-   Authors: Shilo, Asaf; Pegoraro, Gianluca; Misteli, Tom.
-   DOI: 10.1093/nar/gkac869

Once the repo is cloned, and to recreate locally the `renv` library used in the analysis, run the `renv::restore()` command at the R console *before* running any of the `.Rmd` scripts. This operation needs to be performed only once on each computer the analysis is run on.

The analysis for each figure in the manuscript is organized as a self-contained directory that contains:

-   An `input` folder containing the well and single cell results of the image analysis run in Columbus. The input folder contains a lot of data and is not under version control in Github. If the `input` folder is not present, the first time .Rmd analysis script is run the folder will be automatically downloaded from the [corresponding Figshare data repository for this project](<https://figshare.com/projects/HiFENS/138505>) and unzipped in the appropriate location.

-   An `.Rmd` script that contains the R code to run the analysis. The script can be "knitted " in RStudio to reproduce the output.

-   An `.md` output file that is rendered natively in Github as a webpage. This is the file that should be opened (Either in Github or locally using RStudio) to check the results of the analysis script.

-   A `metadata` folder containing the mapping of experimental treatments for each well in the 384-well imaging plates for that analysis.

-   An `output` folder that contains `.png` files containing single plots labelled according to the figure and panel number in the manuscript.

For information about this repo, please contact [Tom Misteli](mailto:mistelit@nih.gov) or [Gianluca Pegoraro](mailto:gianluca.pegoraro@nih.gov)
