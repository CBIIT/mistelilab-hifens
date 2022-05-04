# mistelilab-hifens

This is a repo containing the R code used to generate some of the figures for a preprint manuscript from the Misteli Lab at the National Cancer Institute/NIH.

-   Title: HiFENS: High-throughput FISH detection of endogenous pre-mRNA splicing isoforms
-   Authors: Shilo, Asaf; Pegoraro, Gianluca; Misteli, Tom
-   Link: <https://www.biorxiv.org/content/10.1101/2022.04.10.487792v1>
-   DOI: 10.1101/2022.04.10.487792

The analysis for each figure in the manuscript is organized as a self-contained directory that contains:

-   An `input` folder containing the well and single cell results of the image analysis run in Columbus. The input folder contains a lot of data and is not under version control in Github. If the `input` folder is not present, the first time .Rmd analysis script is run the folder will be downloaded from the Figshare data repository and unzipped in the appropriate location.

-   An `.Rmd` script that contains the R code to run the analysis. The script can be "knitted " in RStudio to reproduce the output.

-   An `.md` output file that is rendered natively in Github as a webpage. This is the file that should be opened to check the results of the analysis script.

-   An `output` folder that contains `.png` files containing single plots.

For information about this repo, please contact [Tom Misteli](mailto:mistelit@nih.gov) or [Gianluca Pegoraro](mailto:gianluca.pegoraro@nih.gov)
