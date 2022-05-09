Figures 4A - 4C
================
Asaf Shilo/Gianluca Pegoraro
May 9th 2022

### Analysis setup

Load required packages.

``` r
library(plyr)
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::arrange()   masks plyr::arrange()
    ## ✖ purrr::compact()   masks plyr::compact()
    ## ✖ dplyr::count()     masks plyr::count()
    ## ✖ dplyr::failwith()  masks plyr::failwith()
    ## ✖ dplyr::filter()    masks stats::filter()
    ## ✖ dplyr::id()        masks plyr::id()
    ## ✖ dplyr::lag()       masks stats::lag()
    ## ✖ dplyr::mutate()    masks plyr::mutate()
    ## ✖ dplyr::rename()    masks plyr::rename()
    ## ✖ dplyr::summarise() masks plyr::summarise()
    ## ✖ dplyr::summarize() masks plyr::summarize()

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:purrr':
    ## 
    ##     transpose

``` r
library(stringr)
library(ggthemes)
library(ggpubr)
```

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     mutate

``` r
library(curl)
```

    ## Using libcurl 7.79.1 with LibreSSL/3.3.5

    ## 
    ## Attaching package: 'curl'

    ## The following object is masked from 'package:readr':
    ## 
    ##     parse_date

Read and Change variable types in experimental metadata data.table.

``` r
dt_md <- fread("metadata/layout.txt")
dt_md <- dt_md[probe != "",]
dt_md[, `:=`(probe = factor(probe, levels = c("None", 
                                              "FGFR2",
                                              "FGFR2-IIIb",
                                              "PGK1",
                                              "FGFR2-D7-10",
                                              "All" )))]
```

### Download the data if needed

Download unzip the Columbus results of the experiments from Figshare if
they have not been already downloaded.

``` r
if(!dir.exists("ObjectLevelData")) {
  URL <- "https://figshare.com/ndownloader/files/35014492"
  curl_download(URL, "input.zip")
  unzip("input.zip")
}
```

### Read and process the cell level data

Set RegEx patterns for directory searches for cell level data.

``` r
pat_c <- ".*Nuclei Final\\[0\\]\\.txt$" # Pattern for Single Cell data files
```

Create a list of the RegEx patterns set in the previous chunk.
**Important:** the list names will be carried over all the next steps!!!

-   c = cell level data

``` r
pat_list <- list(c= pat_c)
```

Recursively search the `ObjectLevelData` directory and its
sub-directories for files whose name includes the RegEx patterns defined
two chunks above. The `path_list` function outputs absolute file names.
`path_list` is a list containing all the file names on a per cell-level.

``` r
list_files <- function(x) {
  dir(
  path = 'ObjectLevelData/',
  pattern = x,
  full.names = TRUE,
  recursive = TRUE,
  include.dirs = TRUE
  )
}

path_list <- llply(pat_list, list_files) 
```

Extract file names from absolut path and set them as list element names.

``` r
trim_names <- function(x) {
  names(x) <-
  basename(x) # This assigns the filename to the file that it is read
  y <- x ## This is necessary because of scoping issues
}

path_list <- llply(path_list, trim_names) 
```

Recursively read and merge object level data files as data.frames. Rows
are labeled with relative filenames (The `.id` variable). This and the
previous chunks are slightly modified tricks adopted from H. Wickam
[“Tidy Data” paper](http://vita.had.co.nz/papers/tidy-data.pdf).

``` r
read_merge <- function(x) {
  dt <- as.data.table(ldply(x, fread, sep = "\t", na.strings = NULL))
}

dt_list <- llply(path_list, read_merge)
```

Separate the cell level data.

``` r
dt_cell <- dt_list$c

rm(dt_list)
```

Change variables names.

``` r
setnames(
  dt_cell,
  c("Row",
    "Column",
    "ScreenName",
    "WellName",
    "Nuclei Final - Nucleus Roundness",
    "Nuclei Final - Number of Positive 488- per Cell",
    "Nuclei Final - Number of Positive 562- per Cell",
    "Nuclei Final - Number of Positive 647- per Cell"),
  c("row",
    "column",
    "plate",
    "well",
    "roundness",
    "spot_n_488",
    "spot_n_562",
    "spot_n_647")
  )
```

Join measurement data with experimental metadata.

``` r
setkey(dt_md, row, column)
setkey(dt_cell, row, column)

dt_full <- dt_cell[dt_md, nomatch = 0] 
```

Generate datasets for plotting.

``` r
dt_MP<- filter(dt_full, probe !='FGFR2')
```

### Plotting

Plotting total FGFR2.

![](output/Fig4_A-1.png)<!-- -->

Plotting FGFR2-IIIb.

![](output/Fig4_B-1.png)<!-- -->

Plotting PGK1 in B2-PGK1.

![](output/Fig4_C-1.png)<!-- -->

Create dataset with IIIb only for U2OS cells over-expressing IIIb.

``` r
dt_MP_IIIbOE <- filter(dt_MP, probe == 'All' & cell_line == 'U2OS FGFR2-IIIb')
```

Scatter plot of total FGFRII vs. IIIb only for U2OS cells
over-expressing IIIb.

![](output/Fig4_E-1.png)<!-- -->

Document the information about the analysis session.

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur/Monterey 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] curl_4.3.2        ggpubr_0.4.0      ggthemes_4.2.4    data.table_1.14.2
    ##  [5] forcats_0.5.1     stringr_1.4.0     dplyr_1.0.9       purrr_0.3.4      
    ##  [9] readr_2.1.2       tidyr_1.2.0       tibble_3.1.7      ggplot2_3.3.6    
    ## [13] tidyverse_1.3.1   plyr_1.8.7       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.8.3     lattice_0.20-45  lubridate_1.8.0  assertthat_0.2.1
    ##  [5] digest_0.6.29    utf8_1.2.2       R6_2.5.1         cellranger_1.1.0
    ##  [9] backports_1.4.1  reprex_2.0.1     evaluate_0.15    highr_0.9       
    ## [13] httr_1.4.2       pillar_1.7.0     rlang_1.0.2      readxl_1.4.0    
    ## [17] rstudioapi_0.13  car_3.0-13       Matrix_1.4-1     rmarkdown_2.14  
    ## [21] splines_4.2.0    labeling_0.4.2   munsell_0.5.0    broom_0.8.0     
    ## [25] compiler_4.2.0   modelr_0.1.8     xfun_0.30        pkgconfig_2.0.3 
    ## [29] mgcv_1.8-40      htmltools_0.5.2  tidyselect_1.1.2 fansi_1.0.3     
    ## [33] crayon_1.5.1     tzdb_0.3.0       dbplyr_2.1.1     withr_2.5.0     
    ## [37] grid_4.2.0       nlme_3.1-157     jsonlite_1.8.0   gtable_0.3.0    
    ## [41] lifecycle_1.0.1  DBI_1.1.2        magrittr_2.0.3   scales_1.2.0    
    ## [45] cli_3.3.0        stringi_1.7.6    carData_3.0-5    farver_2.1.0    
    ## [49] ggsignif_0.6.3   fs_1.5.2         xml2_1.3.3       ellipsis_0.3.2  
    ## [53] generics_0.1.2   vctrs_0.4.1      tools_4.2.0      glue_1.6.2      
    ## [57] hms_1.1.1        abind_1.4-5      fastmap_1.1.0    yaml_2.3.5      
    ## [61] colorspace_2.0-3 rstatix_0.7.0    rvest_1.0.2      knitr_1.39      
    ## [65] haven_2.5.0
