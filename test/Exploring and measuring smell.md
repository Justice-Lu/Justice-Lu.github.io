---
layout: default
title: Exploring and measuring smell
nav_order: 3
---

## Exploring and measuring smell
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>


## Written Functions

Following chunk shows all the written function used in this vignette

## Differential gene expression data

The data table below is prepossessed by DEseq2 to calculate gene-level
differential expression.

``` r
head(get_OR_table())
```

    ##      ensembl_gene_id     logFC          FDR No1 No2 No3  ST1  ST2  ST3       id
    ## 1 ENSMUSG00000059069 2.5448640 6.676352e-65 144 298 249 1299 1764 1605  Olfr749
    ## 2 ENSMUSG00000053815 2.4720328 2.460729e-21  40 100  54  388  421  407  Olfr744
    ## 3 ENSMUSG00000075200 2.1649039 2.596337e-17  49  62  55  356  255  285 Olfr1044
    ## 4 ENSMUSG00000049528 1.7708260 2.680831e-09  29  57  46  183  189  155  Olfr429
    ## 5 ENSMUSG00000050028 1.0149519 1.436028e-06 187 265 264  500  547  672  Olfr745
    ## 6 ENSMUSG00000062782 0.9350916 5.008757e-04 291 304 348 1061  569  606  Olfr527
    ##                      name            odor concentration      odor_and_conc
    ## 1  olfactory receptor 749 2-Phenylethanol            1p 1% 2-Phenylethanol
    ## 2  olfactory receptor 744 2-Phenylethanol            1p 1% 2-Phenylethanol
    ## 3 olfactory receptor 1044 2-Phenylethanol            1p 1% 2-Phenylethanol
    ## 4  olfactory receptor 429 2-Phenylethanol            1p 1% 2-Phenylethanol
    ## 5  olfactory receptor 745 2-Phenylethanol            1p 1% 2-Phenylethanol
    ## 6  olfactory receptor 527 2-Phenylethanol            1p 1% 2-Phenylethanol
    ##   odor_category
    ## 1      Alcohols
    ## 2      Alcohols
    ## 3      Alcohols
    ## 4      Alcohols
    ## 5      Alcohols
    ## 6      Alcohols

From such data table it’s hard to understand what we’re looking for or
what’s actually different. Therefore it’s often a good idea to
restructure your data table into simpler form. Also keep in mind that
when trimming data table, different analysis may require different
measurements and structure.

## Volcano plot

A simple way to visualize and differential gene expression is through a
volcano plot. The volcano plot showcases the genes that are either
enriched in mouse given odor (right side) or control (left).

``` r
fdr_sigs_1p <- get_fdr(negativeOnly = TRUE, odor_file_name = "2m2t", concentration = "1p")
```

    ## [1] "filtered by FDR < 0.05 AND logFC < 0"

``` r
Vol_plot(Odor = "2-Methyl-2-thiazoline", conc = "1p", neg_sig_blue=TRUE)
```

    ## [1] "filtered concentration by  1p"

![](/resource/figure-markdown_github/unnamed-chunk-3-1.png)

However, there are limitations to such plots. As you can only observe
one category at a time.

## Simple Data Wrangling

From the DESeq2 output data table above is a classic example of how raw
data are usually not easily accessible for plotting. Here I use a very
simple example of data wrangling to reformat the data table into ORs by
logFC given its corresponding odor molecule using \< get_heat_table() \>
For more information on how the data table is transformed, refer to the
written functions section.

``` r
get_heat_table()
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "heat_table filtered by FDR < 1"

    ## # A tibble: 1,102 × 73
    ##    `10% (-)-2-Octanol` `1% (-)-Menthol` `10% (+)-2-Octanol` `1% (+)-Menthol`
    ##  *               <dbl>            <dbl>               <dbl>            <dbl>
    ##  1                3.42            0.266                2.84           0.425 
    ##  2                2.22            0.129                1.90          -0.193 
    ##  3                2.53            0.588                2.80          -0.176 
    ##  4                2.68           -0.246                3.03          -0.836 
    ##  5                1.90            0.312                2.12          -0.0861
    ##  6                1.80            0.156                1.62          -0.533 
    ##  7                1.81            0.192                1.41          -0.0169
    ##  8                2.30           -0.556                2.42          -0.638 
    ##  9                1.99            1.16                 2.03           0.369 
    ## 10                2.22            0.519                2.74           0.573 
    ## # … with 1,092 more rows, and 69 more variables: 1% 2-Phenylethanol <dbl>,
    ## #   1% Citronellol <dbl>, 1% Guaiacol <dbl>, 1% Linalool <dbl>,
    ## #   1% p-Cresol <dbl>, 1% 2-Methyl-2-pentenal <dbl>, 0.01% Anisaldehyde <dbl>,
    ## #   1% Benzaldehyde <dbl>, 0.01% Citral <dbl>, 1% Heptanal <dbl>,
    ## #   1% Octanal <dbl>, 1% Trans-cinnamaldehyde <dbl>, 1% Butyric acid <dbl>,
    ## #   1% Heptanoic acid <dbl>, 0.1% Isovaleric acid <dbl>,
    ## #   1% Ethyl butyrate <dbl>, 1% Isoamyl acetate <dbl>, …

## Ways to visualize data

Now that we have a more intuitive data frame that shows how each
specific odor molecules affect the individual activation of certain ORs,
we can utilize different data analysis and plots to visualize and more
importantly ask different questions

A classic way of visualizing how different certain samples are is to
utilize the principle component analysis (pca).

``` r
## Plotting pca plot for the pS6-IP experiments. 
pca_data <- get_pca_data(fdr_value = 1, scale=FALSE, center=TRUE, filter_conc = NULL)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "heat_table filtered by FDR < 1"

``` r
pca_percentage <- get_pca_data(fdr_value = 1, 
                               scale=FALSE,
                               center=TRUE,
                               filter_conc = NULL,
                               return_pca_percentage=TRUE)
```

    ## Warning: Setting row names on a tibble is deprecated.

    ## [1] "heat_table filtered by FDR < 1"

``` r
options(repr.plot.width=7, repr.plot.height=5)
ggplot(pca_data, aes(x=PC1, y=PC2, label=odor, colour=odor_category)) +
    geom_point(size = 7) + 
    #geom_text(size=7,position = position_nudge(y = +0.02)) +
    scale_color_manual(values = my_colors()) +
    theme_classic() + 
    xlab(pca_percentage[1]) + ylab(pca_percentage[2])
```

![](/resource/figure-markdown_github/unnamed-chunk-5-1.png) Here we
can broadly observe that most chemical groups does not have an global
impact in driving overall OR activation. However, certain groups of odor
are robust enough to deviate, such as the thiazole and Ketones category.

## Heatmap

Another interesting way of visualizing such data is to use heatmap.
Heatmap is extremely efficient in looking at large datasets and finding
certain trends or unique outliers.

``` r
## Plotting heatmap for pS6IP data
options(repr.plot.width=16, repr.plot.height=10)
p <- my_heatmap(fdr_value = 0.05,
                rm_0_rows=TRUE, 
                annotate=TRUE)
```

![](/resource/figure-markdown_github/unnamed-chunk-6-1.png)

## Dotplot search function

Say that we’re interested in a certain chemical, in this case
2-Methyl-2-Thiazoline (2mt). Since 2mt is one of the very few odors in
the natural environment that suppress or reduce activation of certain
ORs we’re interested in those corresponding ORs. One simple thing you
can do is to do a reverse search. Since we know the subset of ORs that
are inhibited by 2mt, we’re interested in seeing what odors do those
specific ORs also are responding to. We take those ORs are search our
data, and a good way of illustrating such idea is to not only look at
logFC but also how significant is such result depicted by FDR value.

``` r
## Extract list of Olfr with FDR < 0.05 for plotting 
fdr_2m2t_1p <- ps6_lookup(fdr_value = 0.05, 
           Odor = "2-Methyl-2-thiazoline", 
           conc = "1p", 
           negativeOnly = TRUE,
           select = TRUE)
```

    ## [1] "filtered odor by  2-Methyl-2-thiazoline"
    ## [1] "filtered fdr by  0.05"
    ## [1] "filtered by logFC < 0"
    ## [1] "filtered concentration by  1p"
    ## [1] "Selected columns "

``` r
fdr_2m2t_100p <- ps6_lookup(fdr_value = 0.05, 
           Odor = "2-Methyl-2-thiazoline", 
           conc = "100p", 
           negativeOnly = TRUE,
           select = TRUE)
```

    ## [1] "filtered odor by  2-Methyl-2-thiazoline"
    ## [1] "filtered fdr by  0.05"
    ## [1] "filtered by logFC < 0"
    ## [1] "filtered concentration by  100p"
    ## [1] "Selected columns "

``` r
options(repr.plot.width=13, repr.plot.height=15)
ps6_OR_dot(fdr_2m2t_100p$id[1:100], "100p 2MT pS6-IP top100")
```

![](/resource/figure-markdown_github/unnamed-chunk-7-1.png) A
narrower list of ORs

``` r
## Identify ligands for ORs that are inhibited by 1% 2mt
options(repr.plot.width=10, repr.plot.height=6)
ps6_OR_dot(fdr_2m2t_1p$id, "1p 2MT reduced OR's ligand pS6-IP")
```

![](/resource/figure-markdown_github/unnamed-chunk-8-1.png) This is
not the end of the analysis. But hopefully this vignette demonstrated
that data analysis is extremely fun and provides multiple angle in how
you can visualize and ask different questions in a complex biological
dataset.

## grouping of odor groups

``` r
## Write function to mutate and group odor groups manually 
## Updated and annotade pS6 master is saved 
Alcohols <- c("2-Phenylethanol",
               "Citronellol",
               "Guaiacol",
               "Linalool",
               "p-Cresol",
               "(+)-2-Octanol",
               "(-)-2-Octanol",
               "(+)-Menthol",
               "(-)-Menthol")
Aldehydes <- c("Anisaldehyde",
               "Citral",
               "2-Methyl-2-pentenal",
               "Benzaldehyde",
               "Heptanal",
               "Octanal",
               "Trans-cinnamaldehyde")
CarboxylicAcids <- c(
               "Isovaleric acid",
               "Butyric acid",
               "Heptanoic acid"
)
Esters <- c(
               "Ethyl butyrate",
               "Isoamyl acetate",
               "Methyl salicylate"
)
Ketones <- c(
               "Acetophenone",
               "2-Hydroxy acetophenone",
               "(+)-Carvone",
               "(-)-Carvone",
               "2-Heptanone",
               "2-Hexanone",
               "4-Methylacetophenone",
               "β-Damascone",
               "β-Ionone"
)
Pyradine <- c(
               "2-Ethyl-3-methylpyrazine",
               "2,5-Dimethylpyrazine",
               "Pyridine"
)
Sulfurous <- c(
               "2-Butene-1-thiol",
               "Cyclopentanethiol",
               "3-Methyl-1-butanethiol",
               "Dimethyl trisulfide"
)
Thiazole <- c(
               "SBT",
               "2-Methyl-2-thiazoline",
               "TMT",
               "nTMT"
)
Tiglates <- c(
               "Isopropyl tiglate",
               "Ethyl tiglate",
               "Hexyl tiglate"
)
Others <- c(
               "DHB",
               "(E)-β-farnesene",
               "α-Pinene",
               "Diacetyl",
               "(+)-Limonene"
)

odor_list <- c(
               "Isovaleric acid",
               "Pyridine",
               "3-Methyl-1-butanethiol",
               "2-Hexanone",
               "Cyclopentanethiol",
               "SBT",
               "(+)-Menthol",
               "(-)-Menthol",
               "(+)-Carvone",
               "2-Methyl-2-thiazoline",
               "TMT",
               "Dimethyl trisulfide",
               "Benzaldehyde",
               "Butyric acid",
               "(+)-2-Octanol",
               "2-Ethyl-3-methylpyrazine",
               "p-Cresol",
               "Anisaldehyde",
               "2,5-Dimethylpyrazine",
               "Isoamyl acetate",
               "β-Damascone",
               "Methyl salicylate",
               "(-)-Carvone",
               "(+)-Limonene",
               "Octanal",
               "Guaiacol",
               "Ethyl tiglate",
               "(E)-β-farnesene",
               "2-Methyl-2-pentenal",
               "DHB",
               "Isopropyl tiglate",
               "2-Phenylethanol",
               "nTMT",
               "Trans-cinnamaldehyde",
               "Hexyl tiglate",
               "Citral",
               "β-Ionone",
               "2-Butene-1-thiol",
               "Diacetyl",
               "Linalool",
               "α-Pinene",
               "Acetophenone",
               "Ethyl butyrate",
               "(-)-2-Octanol",
               "2-Heptanone",
               "Heptanoic acid",
               "Heptanal",
               "4-Methylacetophenone",
               "Citronellol"
)

Alcohols <- data.frame(Alcohols)  %>% mutate(odor_category="Alcohols") 
Aldehydes <- data.frame(Aldehydes)   %>% mutate(odor_category="Aldehydes")
CarboxylicAcids <- data.frame(CarboxylicAcids)   %>% mutate(odor_category="CarboxylicAcids")
Esters <- data.frame(Esters)   %>% mutate(odor_category="Esters")
Ketones <- data.frame(Ketones)   %>% mutate(odor_category="Ketones")
Others <- data.frame(Others)   %>% mutate(odor_category="Others")
Pyradine <- data.frame(Pyradine)   %>% mutate(odor_category="Pyradine")
Sulfurous <- data.frame(Sulfurous)   %>% mutate(odor_category="Sulfurous")
Thiazole <- data.frame(Thiazole)   %>% mutate(odor_category="Thiazole")
Tiglates <- data.frame(Tiglates)   %>% mutate(odor_category="Tiglates")

names(Alcohols)[1]<-"odor"
names(Aldehydes)[1]<-"odor"
names(CarboxylicAcids)[1]<-"odor"
names(Esters)[1]<-"odor"
names(Ketones)[1]<-"odor"
names(Others)[1]<-"odor"
names(Pyradine)[1]<-"odor"
names(Sulfurous)[1]<-"odor"
names(Thiazole)[1]<-"odor"
names(Tiglates)[1]<-"odor"

odor_list<-rbind(Alcohols,Aldehydes,CarboxylicAcids,Esters,Ketones,Pyradine,Sulfurous,Thiazole,Tiglates,Others)

#OR_table <- left_join(OR_table,odor_list, by="odor")
#write.csv(OR_table, file="./pS6IP_Annotadklafjasdkl;fjsdakl")

### For all the odors that did not get placed into a category manually above, they're sorted into others
#OR_table$odor_category[is.na(OR_table$odor_category)]<- "Others"
##Previously formatting odor_and_conc
#OR_table <- OR_table %>% mutate(odor_and_conc = paste(concentration, odor, sep = " "))
### Formatting percentage signs
#OR_table$odor_and_conc <- gsub("p1", "0.1%", OR_table$odor_and_conc)
#OR_table$odor_and_conc <- gsub("p01", "0.01%", OR_table$odor_and_conc)
#OR_table$odor_and_conc <- gsub("p ", "% ", OR_table$odor_and_conc)
```

``` r
sessionInfo()
```

    ## R version 4.1.1 (2021-08-10)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Big Sur 11.5.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] pheatmap_1.0.12   forcats_0.5.1     readr_2.0.1       tibble_3.1.4     
    ##  [5] tidyverse_1.3.1   purrr_0.3.4       data.table_1.14.0 tidyr_1.1.3      
    ##  [9] viridis_0.6.1     viridisLite_0.4.0 stringr_1.4.0     ggrepel_0.9.1    
    ## [13] plotly_4.9.4.1    ggplot2_3.3.5     dplyr_1.0.7      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.7         lubridate_1.7.10   assertthat_0.2.1   digest_0.6.27     
    ##  [5] utf8_1.2.2         R6_2.5.1           cellranger_1.1.0   backports_1.2.1   
    ##  [9] reprex_2.0.1       evaluate_0.14      highr_0.9          httr_1.4.2        
    ## [13] pillar_1.6.2       rlang_0.4.11       lazyeval_0.2.2     readxl_1.3.1      
    ## [17] rstudioapi_0.13    rmarkdown_2.10     labeling_0.4.2     htmlwidgets_1.5.4 
    ## [21] munsell_0.5.0      broom_0.7.9        compiler_4.1.1     modelr_0.1.8      
    ## [25] xfun_0.25          pkgconfig_2.0.3    htmltools_0.5.2    tidyselect_1.1.1  
    ## [29] gridExtra_2.3      fansi_0.5.0        crayon_1.4.1       tzdb_0.1.2        
    ## [33] dbplyr_2.1.1       withr_2.4.2        grid_4.1.1         jsonlite_1.7.2    
    ## [37] gtable_0.3.0       lifecycle_1.0.0    DBI_1.1.1          magrittr_2.0.1    
    ## [41] scales_1.1.1       cli_3.0.1          stringi_1.7.4      farver_2.1.0      
    ## [45] fs_1.5.0           xml2_1.3.2         ellipsis_0.3.2     generics_0.1.0    
    ## [49] vctrs_0.3.8        RColorBrewer_1.1-2 tools_4.1.1        glue_1.4.2        
    ## [53] hms_1.1.0          fastmap_1.1.0      yaml_2.2.1         colorspace_2.0-2  
    ## [57] rvest_1.0.1        knitr_1.33         haven_2.4.3
