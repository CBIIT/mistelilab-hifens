Figures 3A, 3B, 3C, and 3E
================
Asaf Shilo/Gianluca Pegoraro
Aug 31st 2022

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

    ## Using libcurl 7.79.1 with LibreSSL/3.3.6

    ## 
    ## Attaching package: 'curl'

    ## The following object is masked from 'package:readr':
    ## 
    ##     parse_date

``` r
library(rstatix)
```

    ## 
    ## Attaching package: 'rstatix'

    ## The following objects are masked from 'package:plyr':
    ## 
    ##     desc, mutate

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(FSA)
```

    ## Registered S3 methods overwritten by 'FSA':
    ##   method       from
    ##   confint.boot car 
    ##   hist.boot    car

    ## ## FSA v0.9.3. See citation('FSA') if used in publication.
    ## ## Run fishR() for related website and fishR('IFAR') for related book.

    ## 
    ## Attaching package: 'FSA'

    ## The following object is masked from 'package:plyr':
    ## 
    ##     mapvalues

Set output folder.

Set the palette and the running theme for ggplot2.

``` r
theme_set(theme_minimal())
theme_update(axis.text.x = element_text(
angle = -90,
hjust = 0,
vjust = 0.5))

theme_update(plot.title = element_text(hjust = 0.5))
```

### Read the experimental metadata

Read and change variable types in experimental metadata data.table.

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

Extract file names from absolute path and set them as list element
names.

``` r
trim_names <- function(x) {
  names(x) <-
  basename(x) # This assigns the file name to the file that it is read
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
dt_IIIB <- filter(dt_full, probe == 'FGFR2-IIIb')
dt_Full <- filter(dt_full, probe == 'FGFR2')
dt_PGK1 <- filter(dt_full, probe == 'PGK1')
dt_D7_10 <-filter(dt_full, probe == 'FGFR2-D7-10')
dt_All <- filter (dt_full, probe == 'All')
dt_MP<- filter(dt_full, probe !='FGFR2')
```

### Statistcal tests

Kruskal-Wallis and Dunn test for FGFR2-Full

``` r
Full.kruskal <- dt_Full %>% kruskal_test(spot_n_488 ~ cell_line)
Full.kruskal
```

    ## # A tibble: 1 × 6
    ##   .y.            n statistic    df     p method        
    ## * <chr>      <int>     <dbl> <int> <dbl> <chr>         
    ## 1 spot_n_488  6848     2113.     3     0 Kruskal-Wallis

``` r
Full.kruskal_s <- dt_Full %>% kruskal_effsize(spot_n_488 ~ cell_line)

# Pairwise comparisons using Dunn’s test:
Full_pwc <- dt_Full %>% 
  dunn_test(spot_n_488 ~ cell_line) 
Full_pwc
```

    ## # A tibble: 6 × 9
    ##   .y.       group1 group2    n1    n2 statistic         p     p.adj p.adj.signif
    ## * <chr>     <chr>  <chr>  <int> <int>     <dbl>     <dbl>     <dbl> <chr>       
    ## 1 spot_n_4… MCF7   U2OS …  1913  1410    -17.4  4.59e- 68 1.38e- 67 ****        
    ## 2 spot_n_4… MCF7   U2OS …  1913  1945     23.5  6.02e-122 2.41e-121 ****        
    ## 3 spot_n_4… MCF7   U2OS …  1913  1580    -16.3  1.47e- 59 2.94e- 59 ****        
    ## 4 spot_n_4… U2OS … U2OS …  1410  1945     39.1  0         0         ****        
    ## 5 spot_n_4… U2OS … U2OS …  1410  1580      1.60 1.10e-  1 1.10e-  1 ns          
    ## 6 spot_n_4… U2OS … U2OS …  1945  1580    -38.7  0         0         ****

``` r
Full_pwc_p <- Full_pwc %>% add_xy_position(x = "cell_line")
```

Kruskal-Wallis and Dunn test for FGFR2-IIIb

``` r
IIIB.kruskal <- dt_IIIB %>% kruskal_test(spot_n_488 ~ cell_line)
IIIB.kruskal
```

    ## # A tibble: 1 × 6
    ##   .y.            n statistic    df     p method        
    ## * <chr>      <int>     <dbl> <int> <dbl> <chr>         
    ## 1 spot_n_488  7636     2737.     3     0 Kruskal-Wallis

``` r
IIIB.kruskal_s <- dt_IIIB %>% kruskal_effsize(spot_n_488 ~ cell_line)

# Pairwise comparisons using Dunn’s test:
IIIB_pwc <- dt_IIIB %>% 
  dunn_test(spot_n_488 ~ cell_line) 
IIIB_pwc
```

    ## # A tibble: 6 × 9
    ##   .y.       group1 group2    n1    n2 statistic         p     p.adj p.adj.signif
    ## * <chr>     <chr>  <chr>  <int> <int>     <dbl>     <dbl>     <dbl> <chr>       
    ## 1 spot_n_4… MCF7   U2OS …  1858  1678    -18.7  8.14e- 78 1.63e- 77 ****        
    ## 2 spot_n_4… MCF7   U2OS …  1858  2064     22.5  7.55e-112 2.27e-111 ****        
    ## 3 spot_n_4… MCF7   U2OS …  1858  2036    -23.6  1.68e-123 6.71e-123 ****        
    ## 4 spot_n_4… U2OS … U2OS …  1678  2064     41.0  0         0         ****        
    ## 5 spot_n_4… U2OS … U2OS …  1678  2036     -3.93 8.65e-  5 8.65e-  5 ****        
    ## 6 spot_n_4… U2OS … U2OS …  2064  2036    -47.3  0         0         ****

``` r
IIIB_pwc_p <- IIIB_pwc %>% add_xy_position(x = "cell_line")
```

Kruskal-Wallis and Dunn test for PGK1

``` r
PGK1.kruskal <- dt_PGK1 %>% kruskal_test(spot_n_562 ~ cell_line)
PGK1.kruskal
```

    ## # A tibble: 1 × 6
    ##   .y.            n statistic    df        p method        
    ## * <chr>      <int>     <dbl> <int>    <dbl> <chr>         
    ## 1 spot_n_562  7735      57.1     3 2.48e-12 Kruskal-Wallis

``` r
PGK1.kruskal_s <- dt_PGK1 %>% kruskal_effsize(spot_n_562 ~ cell_line)

# Pairwise comparisons using Dunn’s test:
PGK1_pwc <- dt_PGK1 %>% 
  dunn_test(spot_n_562 ~ cell_line) 
PGK1_pwc
```

    ## # A tibble: 6 × 9
    ##   .y.        group1  group2    n1    n2 statistic        p    p.adj p.adj.signif
    ## * <chr>      <chr>   <chr>  <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
    ## 1 spot_n_562 MCF7    U2OS …  1954  1745     6.65  2.91e-11 1.75e-10 ****        
    ## 2 spot_n_562 MCF7    U2OS …  1954  2101     6.39  1.65e-10 8.26e-10 ****        
    ## 3 spot_n_562 MCF7    U2OS …  1954  1935     4.06  4.83e- 5 1.93e- 4 ***         
    ## 4 spot_n_562 U2OS EV U2OS …  1745  2101    -0.563 5.74e- 1 5.74e- 1 ns          
    ## 5 spot_n_562 U2OS EV U2OS …  1745  1935    -2.69  7.19e- 3 2.16e- 2 *           
    ## 6 spot_n_562 U2OS F… U2OS …  2101  1935    -2.24  2.52e- 2 5.04e- 2 ns

``` r
PGK1_pwc_p <- PGK1_pwc %>% add_xy_position(x = "cell_line")
```

### Plotting

Plotting total FGFR2 in B1-FGFR2 System.

![](output/Fig3_A-1.png)<!-- -->

Plotting total FGFR2 in B1-FGFR2 System with p-values

![](output/Fig3_A_2-1.png)<!-- -->

Plotting FGFR2-IIIb in B1-IIIB System

![](output/Fig3_B-1.png)<!-- -->

Plotting FGFR2-IIIb in B1-IIIB System with p-values

``` r
spot_B1IIIb_violin_pv <- ggplot(dt_IIIB, aes(x = cell_line, y = spot_n_488 ))

spot_B1IIIb_violin_pv + geom_violin(fill = "green4") + geom_boxplot(outlier.shape = NA, width=0.03, color="black")+
  
            stat_pvalue_manual(IIIB_pwc_p, hide.ns = TRUE) +
            labs(y = "FGFR2-IIIb HCR Spots per cell") +
                theme_classic()+
                theme(
                legend.position ="none",
                plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1),
                axis.title = element_text(size = 20, color="black", face="bold"),
                axis.text = element_text(size=20,color="black", face="bold"))
```

![](output/Fig3_B_2-1.png)<!-- -->

Plotting PGK1 in B2-PGK1 System.

![](output/Fig3_C-1.png)<!-- -->

Plotting PGK1 in B2-PGK1 System with p-values

![](output/Fig3_C_2-1.png)<!-- -->

Document the information about the analysis session.

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
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
    ##  [1] FSA_0.9.3         rstatix_0.7.0     curl_4.3.2        ggpubr_0.4.0     
    ##  [5] ggthemes_4.2.4    data.table_1.14.2 forcats_0.5.1     stringr_1.4.0    
    ##  [9] dplyr_1.0.9       purrr_0.3.4       readr_2.1.2       tidyr_1.2.0      
    ## [13] tibble_3.1.7      ggplot2_3.3.6     tidyverse_1.3.1   plyr_1.8.7       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.8.3     lubridate_1.8.0  assertthat_0.2.1 digest_0.6.29   
    ##  [5] utf8_1.2.2       R6_2.5.1         cellranger_1.1.0 backports_1.4.1 
    ##  [9] reprex_2.0.1     evaluate_0.15    highr_0.9        httr_1.4.4      
    ## [13] pillar_1.7.0     rlang_1.0.4      readxl_1.4.0     rstudioapi_0.14 
    ## [17] car_3.0-13       rmarkdown_2.14   labeling_0.4.2   munsell_0.5.0   
    ## [21] broom_0.8.0      compiler_4.2.1   modelr_0.1.8     xfun_0.30       
    ## [25] pkgconfig_2.0.3  htmltools_0.5.2  tidyselect_1.1.2 fansi_1.0.3     
    ## [29] crayon_1.5.1     tzdb_0.3.0       dbplyr_2.1.1     withr_2.5.0     
    ## [33] grid_4.2.1       jsonlite_1.8.0   gtable_0.3.0     lifecycle_1.0.1 
    ## [37] DBI_1.1.2        magrittr_2.0.3   scales_1.2.0     cli_3.3.0       
    ## [41] stringi_1.7.6    carData_3.0-5    farver_2.1.0     ggsignif_0.6.3  
    ## [45] fs_1.5.2         xml2_1.3.3       ellipsis_0.3.2   generics_0.1.2  
    ## [49] vctrs_0.4.1      tools_4.2.1      glue_1.6.2       hms_1.1.1       
    ## [53] abind_1.4-5      fastmap_1.1.0    yaml_2.3.5       colorspace_2.0-3
    ## [57] rvest_1.0.2      knitr_1.39       haven_2.5.0
