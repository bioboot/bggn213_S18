---
layout: page
title: BGGN-213, Lecture 15
---

Pathway analysis from RNA-seq differential expression results
=============================================================

**BGGN-213 Lecture 15:**  
Barry Grant &lt; <http://thegrantlab.org> &gt;  
Date: 2018-05-22   (15:18:50 PDT on Tue, May 22)  
{:.message}


Overview
--------

Analysis of high-throughput biological data typically yields a list of genes or proteins requiring further interpretation - for example the ranked lists of differentially expressed genes we have been generating from our RNA-seq analysis to date.

Our intention is typically to use such lists to gain novel insights about genes and proteins that may have roles in a given phenomenon, phenotype or disease progression. However, in many cases these 'raw' gene lists are challenging to interpret due to their large size and lack of useful annotations. Hence, our expensively assembled gene lists often fail to convey the full degree of possible insight about the condition being studied.

Pathway analysis (also known as gene set analysis or over-representation analysis), aims to reduce the complexity of interpreting gene lists via mapping the listed genes to known (i.e. annotated) biological pathways, processes and functions.

> **Side-note**: Pathway analysis can actually mean many different things to different people. This includes analysis of Gene Ontology (GO) terms, protein–protein interaction networks, flux-balance analysis from kinetic simulations of pathways, etc. However, pathway analysis most commonly focuses on methods that exploit existing pathway knowledge (e.g. in public repositories such as GO or KEGG), rather than on methods that infer pathways from molecular measurements. These more general approaches are nicely reviewed in this paper:
>
> -   Khatri, et al. "*Ten years of pathway analysis: current approaches and outstanding challenges*." [PLoS Comput Biol 8.2 (2012): e1002375](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002375).

### Patway analysis with R and Bioconductor

There are many freely available tools for pathway or over-representation analysis. As of Nov 2017 Bioconductor alone has over [80 packages categorized under gene set enrichment](http://bioconductor.org/packages/release/BiocViews.html#___GeneSetEnrichment) and over [120 packages categorized under pathways](http://bioconductor.org/packages/release/BiocViews.html#___Pathways).

Here we play with just one, the [**GAGE** package](https://bioconductor.org/packages/release/bioc/html/gage.html) (which stands for **G**enerally **A**pplicable **G**ene set **E**nrichment), to do KEGG pathway enrichment analysis on RNA-seq based differential expression results.

The [KEGG pathway database](http://www.genome.jp/kegg/pathway.html), unlike GO for example, provides functional annotation as well as information about gene products that interact with each other in a given pathway, how they interact (e.g., activation, inhibition, etc.), and where they interact (e.g., cytoplasm, nucleus, etc.). Hence KEGG has the potential to provide extra insight beyond annotation lists of simple molecular function, process etc. from GO terms.

In this analysis, we check for coordinated differential expression over gene sets from KEGG pathways instead of changes of individual genes. The assumption here is that consistent perturbations over a given pathway (gene set) may suggest mechanistic changes.

### About our Input Data

The data for for hands-on session comes from GEO entry: [GSE37704](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37704), which is associated with the following publication:

-   Trapnell C, Hendrickson DG, Sauvageau M, Goff L et al. "*Differential analysis of gene regulation at transcript resolution with RNA-seq*". Nat Biotechnol 2013 Jan;31(1):46-53. [PMID: 23222703](https://www.ncbi.nlm.nih.gov/pubmed/23222703)

The authors report on differential analysis of lung fibroblasts in response to loss of the developmental transcription factor HOXA1. Their results and others indicate that HOXA1 is required for lung fibroblast and HeLa cell cycle progression. In particular their analysis show that "loss of HOXA1 results in significant expression level changes in thousands of individual transcripts, along with isoform switching events in key regulators of the cell cycle". For our session we have used their [Sailfish](https://www.nature.com/articles/nbt.2862) gene-level estimated counts and hence are restricted to protein-coding genes only.

Section 1. Differential Expression Analysis
-------------------------------------------

You can download the count data and associated metadata from here: [GSE37704\_featurecounts.csv](https://bioboot.github.io/bimm143_W18/class-material/GSE37704_featurecounts.csv) and [GSE37704\_metadata.csv](https://bioboot.github.io/bimm143_W18/class-material/GSE37704_metadata.csv). This is similar to our starting point for the last class where we used DESeq2 for the first time. We will use it again today!

``` r
library(DESeq2)
```

Load our data files

``` r
metaFile <- "data/GSE37704_metadata.csv"
countFile <- "data/GSE37704_featurecounts.csv"

# Import metadata and take a peak
colData = read.csv(metaFile, row.names=1)
head(colData)
```

    ##               condition
    ## SRR493366 control_sirna
    ## SRR493367 control_sirna
    ## SRR493368 control_sirna
    ## SRR493369      hoxa1_kd
    ## SRR493370      hoxa1_kd
    ## SRR493371      hoxa1_kd

``` r
# Import countdata
countData = read.csv(countFile, row.names=1)
head(countData)
```

    ##                 length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000186092    918         0         0         0         0         0
    ## ENSG00000279928    718         0         0         0         0         0
    ## ENSG00000279457   1982        23        28        29        29        28
    ## ENSG00000278566    939         0         0         0         0         0
    ## ENSG00000273547    939         0         0         0         0         0
    ## ENSG00000187634   3214       124       123       205       207       212
    ##                 SRR493371
    ## ENSG00000186092         0
    ## ENSG00000279928         0
    ## ENSG00000279457        46
    ## ENSG00000278566         0
    ## ENSG00000273547         0
    ## ENSG00000187634       258

Hmm... remember that we need the `countData` and `colData` files to match up so we will need to remove that odd first column in `countData` namely `contData$length`.

``` r
# Note we need to remove the odd first $length col
countData <- as.matrix(countData[,-1])
head(countData)
```

    ##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000186092         0         0         0         0         0
    ## ENSG00000279928         0         0         0         0         0
    ## ENSG00000279457        23        28        29        29        28
    ## ENSG00000278566         0         0         0         0         0
    ## ENSG00000273547         0         0         0         0         0
    ## ENSG00000187634       124       123       205       207       212
    ##                 SRR493371
    ## ENSG00000186092         0
    ## ENSG00000279928         0
    ## ENSG00000279457        46
    ## ENSG00000278566         0
    ## ENSG00000273547         0
    ## ENSG00000187634       258

This looks better but there are lots of zero entries in there so let's get rid of them as we have no data for these.

``` r
# Filter count data where you have 0 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)
```

    ##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000279457        23        28        29        29        28
    ## ENSG00000187634       124       123       205       207       212
    ## ENSG00000188976      1637      1831      2383      1226      1326
    ## ENSG00000187961       120       153       180       236       255
    ## ENSG00000187583        24        48        65        44        48
    ## ENSG00000187642         4         9        16        14        16
    ##                 SRR493371
    ## ENSG00000279457        46
    ## ENSG00000187634       258
    ## ENSG00000188976      1504
    ## ENSG00000187961       357
    ## ENSG00000187583        64
    ## ENSG00000187642        16

Nice now lets setup the `DESeqDataSet` object required for the **DESeq()** function and then run the DESeq pipeline. This is again similar to our last days hands-on session.

``` r
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds = DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing


A new `DESeqDataSet` is returned that contains all the results (and the input `DESeqDataSet` parameters) within it.  

Extracting out the tables of results we actually want from this object can be a bit tricky. The next section describes one common way to do this. 

### Extracting our results table
Calling the *DESeq* packages **results()** function on our `DESeqDataSet` without any arguments will extract the estimated log2 fold changes and p values like so:


```r
res <- results(dds)
res
```

```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 15280 rows and 6 columns
##                   baseMean log2FoldChange      lfcSE        stat
##                  <numeric>      <numeric>  <numeric>   <numeric>
## ENSG00000279457   29.91358     0.17927483 0.32459294   0.5523066
## ENSG00000187634  183.22965     0.42644724 0.14017817   3.0421802
## ENSG00000188976 1651.18808    -0.69272061 0.05484412 -12.6307173
## ENSG00000187961  209.63794     0.72975918 0.13178350   5.5375609
## ENSG00000187583   47.25512     0.04055411 0.27169055   0.1492658
## ...                    ...            ...        ...         ...
## ENSG00000273748  35.302652      0.6743994  0.3034582   2.2223801
## ENSG00000278817   2.423024     -0.3889516  1.1295943  -0.3443286
## ENSG00000278384   1.101796      0.3328870  1.6590966   0.2006435
## ENSG00000276345  73.644956     -0.3561673  0.2075751  -1.7158482
## ENSG00000271254 181.595903     -0.6096640  0.1412340  -4.3166951
##                       pvalue         padj
##                    <numeric>    <numeric>
## ENSG00000279457 5.807383e-01 6.846746e-01
## ENSG00000187634 2.348712e-03 5.109223e-03
## ENSG00000188976 1.429690e-36 1.745815e-35
## ENSG00000187961 3.067131e-08 1.109758e-07
## ENSG00000187583 8.813439e-01 9.191354e-01
## ...                      ...          ...
## ENSG00000273748 2.625763e-02 4.756160e-02
## ENSG00000278817 7.305992e-01 8.086868e-01
## ENSG00000278384 8.409773e-01 8.927559e-01
## ENSG00000276345 8.618983e-02 1.389975e-01
## ENSG00000271254 1.583827e-05 4.470014e-05
```

The returned `res` object is not a standard R data.frame but one that carries extra meatadata on the meaning of the columns:


```r
mcols(res, use.names = TRUE)
```

```
## DataFrame with 6 rows and 2 columns
##                        type
##                 <character>
## baseMean       intermediate
## log2FoldChange      results
## lfcSE               results
## stat                results
## pvalue              results
## padj                results
##                                                                description
##                                                                <character>
## baseMean                         mean of normalized counts for all samples
## log2FoldChange log2 fold change (MLE): condition hoxa1 kd vs control sirna
## lfcSE                  standard error: condition hoxa1 kd vs control sirna
## stat                   Wald statistic: condition hoxa1 kd vs control sirna
## pvalue              Wald test p-value: condition hoxa1 kd vs control sirna
## padj                                                  BH adjusted p-values
```

The column `log2FoldChange` is the effect size estimate. It tells us how much the gene's expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples.  This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene's expression is increased by a multiplicative factor of \(2^{1.5} \approx 2.82\).

*DESeq2* performs for each gene a *hypothesis test* to see whether evidence is sufficient to decide against the *null hypothesis* that there is zero effect of the treatment on the gene and that the observed difference between treatment and
control was merely caused by experimental variability (i.e., the type of variability that you can expect between different
samples in the same treatment group). As usual in statistics, the result of this test is reported as a *p* value, and it is found in the column `pvalue`. Remember that a *p* value indicates the probability that a fold change as strong as the observed one, or even stronger, would be seen under the situation described by the null hypothesis.  

We can also summarize the results with the *DESeq2* specific version of the **summary()** function. This will report some additional useful information: 


```r
summary(res)
```

```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.1
## LFC > 0 (up)     : 4352, 28% 
## LFC < 0 (down)   : 4400, 29% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 590, 3.9% 
## (mean count < 1)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

Note that there are many many genes with differential expression reported above. Let's therefore be more strict about which set of genes are considered 'significant'. There are two main ways we can do this:

* Lower the false discovery rate threshold (i.e. the threshold on the adjusted p-value (`padj`) in the results table)
* Raise the log2 fold change threshold from 0 to a higher value.



> **Q.** In the summary of our results printed above (and by default) the FDR level is set to 10% (i.e. adjusted p-value < 0.1) and the log2 fold change threshold is set to 0. Use the `alpha` and `lfcThreshold` input arguments to the **results()** function to change these to an FDR of 5% and a log2 fold change threshold of 2. Then use the **summary()** function to find out how many genes are up and down at these thresholds.  




```
## 
## out of 15280 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)     : 99, 0.65% 
## LFC < 0 (down)   : 134, 0.88% 
## outliers [1]     : 0, 0% 
## low counts [2]   : 1482, 9.7% 
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

You could also use the ever useful **table()** function on your output of the **results()** function like so:


```r
table(res$padj < 0.05)
```

```
## 
## FALSE  TRUE 
## 13565   233
```


```r
table(res$log2FoldChange > 2)
```

```
## 
## FALSE  TRUE 
## 14723   557
```

Then combining to determine the number of genes that meet both the *p* value and log2 fold change thresholds (*UP* genes: 99; and *DOWN* genes: 134): 


```r
table( res$padj < 0.05, res$log2FoldChange > 2)
```

```
##        
##         FALSE  TRUE
##   FALSE 13292   273
##   TRUE    134    99
```


> **Side-Note:** In high-throughput biology, we are careful to not use the *p* values directly as evidence against the null, but to correct for **multiple testing**.  
>
> What would happen if we were to simply threshold the *p* values at a low value, say 0.05? There are 326 genes with a *p* value below 0.05 among the 15280 genes for which the test succeeded in reporting a *p* value:


```r
table(res$pvalue < 0.05)
```

```
## 
## FALSE  TRUE 
## 14954   326
```

> *DESeq2* uses the Benjamini-Hochberg (BH) adjustment as implemented in the base R **p.adjust()** function; in brief, this method calculates for each gene an adjusted *p* value that answers the following question: if one called significant all genes with an adjusted *p* value less than or equal to this gene's adjusted *p* value threshold, what would be the fraction of false positives (the *false discovery rate*, FDR) among them, in the sense of the calculation outlined above? These values, called the BH-adjusted *p* values, are given in the column `padj` of the `res` object.  
>
> The FDR is a useful statistic for many high-throughput experiments, as we are often interested in reporting or focusing on a set of interesting genes, and we would like to put an upper bound on the percent of false positives in this set.  
>
> Hence, if we consider a fraction of 5% false positives acceptable, we can consider all genes with an adjusted *p* value below 5% = 0.05 as significant. How many such genes are there?  


```r
table(res$padj < 0.05)
```

```
## 
## FALSE  TRUE 
## 13565   233
```


We can now subset the results table to extract the genes with adjusted *p* value less than 0.05 and then sort them by their log2 fold change estimate to get the significant genes with the strongest down-regulation:


```r
# Make a new results object 'resSig' with only significant genes
resSig <- subset(res, padj < 0.05)

# Print the first 10 strongest DOWN genes
ord.down <- order(resSig$log2FoldChange)
head(resSig[ ord.down, ], 10)
```

```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 10 rows and 6 columns
##                  baseMean log2FoldChange     lfcSE       stat       pvalue
##                 <numeric>      <numeric> <numeric>  <numeric>    <numeric>
## ENSG00000152779  20.47585      -4.159774 0.6372338  -3.389295 7.007247e-04
## ENSG00000139269  34.67082      -3.707145 0.4489612  -3.802434 1.432814e-04
## ENSG00000100526 104.74410      -3.670758 0.2538387  -6.581965 4.642702e-11
## ENSG00000162520 248.71515      -3.522675 0.1645404  -9.254111 2.160200e-20
## ENSG00000117650 126.83830      -3.436330 0.2201916  -6.523092 6.887258e-11
## ENSG00000183287 603.48721      -3.436112 0.1129688 -12.712465 5.041941e-37
## ENSG00000013293 186.24913      -3.419521 0.1869698  -7.592247 3.144050e-14
## ENSG00000183763  46.40721      -3.398924 0.3586012  -3.901059 9.577266e-05
## ENSG00000157456 291.43640      -3.380308 0.1501209  -9.194645 3.762361e-20
## ENSG00000171320  49.47850      -3.368015 0.3462612  -3.950819 7.788429e-05
##                         padj
##                    <numeric>
## ENSG00000152779 4.259295e-02
## ENSG00000139269 9.643890e-03
## ENSG00000100526 6.536735e-09
## ENSG00000162520 6.931729e-18
## ENSG00000117650 9.599029e-09
## ENSG00000183287 3.864928e-34
## ENSG00000013293 5.633974e-12
## ENSG00000183763 6.955111e-03
## ENSG00000157456 1.128545e-17
## ENSG00000171320 5.746778e-03
```

> **Q.** Do the same as above but print out the top 10 strongest up-regulated genes. HINT: see the help for the **order()** function to see how to return the decreasing ordered indices you will want for accesing your `resSig` result.  



```
## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
## Wald test p-value: condition hoxa1 kd vs control sirna 
## DataFrame with 10 rows and 6 columns
##                   baseMean log2FoldChange     lfcSE      stat       pvalue
##                  <numeric>      <numeric> <numeric> <numeric>    <numeric>
## ENSG00000128052  158.24955       8.822086 1.0347214  6.593162 4.305562e-11
## ENSG00000141668   43.94156       8.820638 1.2040490  5.664751 1.472384e-08
## ENSG00000004799   97.26023       8.117454 1.0373919  5.896956 3.702688e-09
## ENSG00000109321 3907.19054       7.479166 0.2880675 19.020423 1.155483e-80
## ENSG00000162892   25.27893       7.052106 1.2210031  4.137669 3.508527e-05
## ENSG00000073756 1357.62707       6.337747 0.3399843 12.758671 2.789195e-37
## ENSG00000125845 1172.23868       5.525337 0.2077712 16.967397 1.431378e-64
## ENSG00000163739   89.96993       5.356228 0.8223964  4.081034 4.483576e-05
## ENSG00000163734   33.68060       5.218578 0.6984361  4.608263 4.060459e-06
## ENSG00000072041 2369.28278       5.097353 0.1868148 16.579804 9.754639e-62
##                         padj
##                    <numeric>
## ENSG00000128052 6.124551e-09
## ENSG00000141668 1.751376e-06
## ENSG00000004799 4.561579e-07
## ENSG00000109321 3.188670e-77
## ENSG00000162892 2.719700e-03
## ENSG00000073756 2.263842e-34
## ENSG00000125845 2.821451e-61
## ENSG00000163739 3.399142e-03
## ENSG00000163734 3.591424e-04
## ENSG00000072041 1.495495e-58
```




## Annotating our genes and mapping to Entrez IDs

Since we mapped and counted against the Ensembl annotation, our results only have information about Ensembl gene IDs. However, our pathway analysis downstream will use KEGG pathways, and genes in KEGG pathways are annotated with Entrez gene IDs. So lets add them as we did the last day.

``` r
library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
    ##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
    ## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
    ## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [25] "UNIGENE"      "UNIPROT"

``` r
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")

res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")

res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

head(res, 10)
```

    ## log2 fold change (MLE): condition hoxa1_kd vs control_sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 10 rows and 9 columns
    ##                  baseMean log2FoldChange      lfcSE      stat    pvalue
    ##                 <numeric>      <numeric>  <numeric> <numeric> <numeric>
    ## ENSG00000117519  4483.627      -2.422719 0.06001850 -40.36620         0
    ## ENSG00000183508  2053.881       3.201955 0.07241968  44.21388         0
    ## ENSG00000159176  5692.463      -2.313737 0.05757255 -40.18820         0
    ## ENSG00000150938  7442.986      -2.059631 0.05386627 -38.23601         0
    ## ENSG00000116016  4423.947      -1.888019 0.04318301 -43.72134         0
    ## ENSG00000136068  3796.127      -1.649792 0.04394825 -37.53942         0
    ## ENSG00000164251  2348.770       3.344508 0.06907610  48.41773         0
    ## ENSG00000124766  2576.653       2.392288 0.06171493  38.76352         0
    ## ENSG00000124762 28106.119       1.832258 0.03892405  47.07264         0
    ## ENSG00000106366 43719.126      -1.844046 0.04194432 -43.96415         0
    ##                      padj      symbol      entrez
    ##                 <numeric> <character> <character>
    ## ENSG00000117519         0        CNN3        1266
    ## ENSG00000183508         0      FAM46C       54855
    ## ENSG00000159176         0       CSRP1        1465
    ## ENSG00000150938         0       CRIM1       51232
    ## ENSG00000116016         0       EPAS1        2034
    ## ENSG00000136068         0        FLNB        2317
    ## ENSG00000164251         0       F2RL1        2150
    ## ENSG00000124766         0        SOX4        6659
    ## ENSG00000124762         0      CDKN1A        1026
    ## ENSG00000106366         0    SERPINE1        5054
    ##                                                        name
    ##                                                 <character>
    ## ENSG00000117519                                  calponin 3
    ## ENSG00000183508 family with sequence similarity 46 member C
    ## ENSG00000159176         cysteine and glycine rich protein 1
    ## ENSG00000150938 cysteine rich transmembrane BMP regulator 1
    ## ENSG00000116016            endothelial PAS domain protein 1
    ## ENSG00000136068                                   filamin B
    ## ENSG00000164251                 F2R like trypsin receptor 1
    ## ENSG00000124766                                   SRY-box 4
    ## ENSG00000124762        cyclin dependent kinase inhibitor 1A
    ## ENSG00000106366                    serpin family E member 1

Great, this is looking good so far. Now lets see how pathway analysis can help us make further sense out of this ranked list of differentially expressed genes.

Section 2. Pathway Analysis
---------------------------

Here we are going to use the [**gage**]() package for pathway analysis. Once we have a list of enriched pathways, we're going to use the [**pathview**]() package to draw pathway diagrams, shading the molecules in the pathway by their degree of up/down-regulation.

### KEGG pathways

The **gageData** package has pre-compiled databases mapping genes to KEGG pathways and GO terms for common organisms. `kegg.sets.hs` is a named list of 229 elements. Each element is a character vector of member gene Entrez IDs for a single KEGG pathway. (See also `go.sets.hs`). The `sigmet.idx.hs` is an index of numbers of signaling and metabolic pathways in `kegg.set.gs`. In other words, KEGG pathway include other types of pathway definitions, like "Global Map" and "Human Diseases", which may be undesirable in pathway analysis. Therefore, `kegg.sets.hs[sigmet.idx.hs]` gives you the "cleaner" gene sets of signaling and metabolic pathways only.

> **Side-Note**: While there are many freely available tools to do pathway analysis, and some like gage are truly fantastic, many of them are poorly maintained or rarely updated. The DAVID tool that a lot of folks use for simple gene set enrichment analysis was not updated at all between Jan 2010 and Oct 2016.

First we need to do our one time install of these required bioconductor packages:

``` r
source("http://bioconductor.org/biocLite.R")
biocLite( c("pathview", "gage", "gageData") )
```

Now we can load the packages and setup the KEGG data-sets we need.

``` r
library(pathview)
```

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

``` r
library(gage)
library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)

kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
```

    ## $`hsa00232 Caffeine metabolism`
    ## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
    ## 
    ## $`hsa00983 Drug metabolism - other enzymes`
    ##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
    ##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
    ## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
    ## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
    ## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
    ## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
    ## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
    ## [50] "8833"   "9"      "978"   
    ## 
    ## $`hsa00230 Purine metabolism`
    ##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
    ##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
    ##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
    ##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
    ##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
    ##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
    ##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
    ##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
    ##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
    ##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
    ##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
    ##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
    ##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
    ##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
    ##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
    ## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
    ## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
    ## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
    ## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
    ## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
    ## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
    ## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
    ## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
    ## [162] "9615"

The main **gage()** function requires a named vector of fold changes, where the names of the values are the Entrez gene IDs.

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

    ##      1266     54855      1465     51232      2034      2317 
    ## -2.422719  3.201955 -2.313737 -2.059631 -1.888019 -1.649792

Now, let’s run the pathway analysis. See help on the gage function with `?gage`. Specifically, you might want to try changing the value of `same.dir`. This value determines whether to test for changes in a gene set toward a single direction (all genes up or down regulated) or changes towards both directions simultaneously (i.e. any genes in the pathway dysregulated).

Here, we're using `same.dir=TRUE`, which will give us separate lists for pathways that are upregulated versus pathways that are down-regulated. Let’s look at the first few results from each.

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
```

Lets look at the result object. It is a list with three elements ("greater", "less" and "stats").

``` r
attributes(keggres)
```

    ## $names
    ## [1] "greater" "less"    "stats"

So it is a list object (you can check it with `str(keggres)`) and we can use the dollar syntax to access a named element, e.g.

``` r
head(keggres$greater)
```

    ##                                         p.geomean stat.mean       p.val
    ## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
    ## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
    ## hsa04740 Olfactory transduction       0.014450242  2.239717 0.014450242
    ## hsa04916 Melanogenesis                0.022339115  2.023074 0.022339115
    ##                                           q.val set.size        exp1
    ## hsa04640 Hematopoietic cell lineage   0.3847887       49 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.3847887      103 0.005655916
    ## hsa04142 Lysosome                     0.3847887      117 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.3847887       26 0.009619717
    ## hsa04740 Olfactory transduction       0.4624078       39 0.014450242
    ## hsa04916 Melanogenesis                0.5297970       85 0.022339115

``` r
head(keggres$less)
```

    ##                                      p.geomean stat.mean        p.val
    ## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
    ## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
    ## hsa03013 RNA transport            1.471026e-03 -3.007785 1.471026e-03
    ## hsa04114 Oocyte meiosis           1.987557e-03 -2.915377 1.987557e-03
    ## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
    ## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
    ##                                         q.val set.size         exp1
    ## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
    ## hsa03030 DNA replication          0.007127774       36 8.909718e-05
    ## hsa03013 RNA transport            0.078454709      143 1.471026e-03
    ## hsa04114 Oocyte meiosis           0.079502292       98 1.987557e-03
    ## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
    ## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03

Each `keggres$greater` and `keggres$less` object is data matrix with gene sets as rows sorted by p-value. Lets look at both up (greater), down (less), and statistics by calling **head()** with the **lapply()** function. As always if you want to find out more about a particular function or its return values use the R help system (e.g. `?gage` or `?lapply`).

``` r
lapply(keggres, head)
```

    ## $greater
    ##                                         p.geomean stat.mean       p.val
    ## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
    ## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
    ## hsa04740 Olfactory transduction       0.014450242  2.239717 0.014450242
    ## hsa04916 Melanogenesis                0.022339115  2.023074 0.022339115
    ##                                           q.val set.size        exp1
    ## hsa04640 Hematopoietic cell lineage   0.3847887       49 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.3847887      103 0.005655916
    ## hsa04142 Lysosome                     0.3847887      117 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.3847887       26 0.009619717
    ## hsa04740 Olfactory transduction       0.4624078       39 0.014450242
    ## hsa04916 Melanogenesis                0.5297970       85 0.022339115
    ## 
    ## $less
    ##                                      p.geomean stat.mean        p.val
    ## hsa04110 Cell cycle               1.004024e-05 -4.353447 1.004024e-05
    ## hsa03030 DNA replication          8.909718e-05 -3.968605 8.909718e-05
    ## hsa03013 RNA transport            1.471026e-03 -3.007785 1.471026e-03
    ## hsa04114 Oocyte meiosis           1.987557e-03 -2.915377 1.987557e-03
    ## hsa03440 Homologous recombination 2.942017e-03 -2.868137 2.942017e-03
    ## hsa00240 Pyrimidine metabolism    5.800212e-03 -2.549616 5.800212e-03
    ##                                         q.val set.size         exp1
    ## hsa04110 Cell cycle               0.001606438      120 1.004024e-05
    ## hsa03030 DNA replication          0.007127774       36 8.909718e-05
    ## hsa03013 RNA transport            0.078454709      143 1.471026e-03
    ## hsa04114 Oocyte meiosis           0.079502292       98 1.987557e-03
    ## hsa03440 Homologous recombination 0.094144560       28 2.942017e-03
    ## hsa00240 Pyrimidine metabolism    0.138500584       95 5.800212e-03
    ## 
    ## $stats
    ##                                       stat.mean     exp1
    ## hsa04640 Hematopoietic cell lineage    2.857393 2.857393
    ## hsa04630 Jak-STAT signaling pathway    2.557207 2.557207
    ## hsa04142 Lysosome                      2.384783 2.384783
    ## hsa00140 Steroid hormone biosynthesis  2.432105 2.432105
    ## hsa04740 Olfactory transduction        2.239717 2.239717
    ## hsa04916 Melanogenesis                 2.023074 2.023074

Now, let's process the results to pull out the top 5 upregulated pathways, then further process that just to get the IDs. We’ll use these KEGG pathway IDs downstream for plotting.

``` r
## Sanity check displaying all pathways data
pathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(pathways)
```

    ##                                                                          id
    ## hsa04640 Hematopoietic cell lineage     hsa04640 Hematopoietic cell lineage
    ## hsa04630 Jak-STAT signaling pathway     hsa04630 Jak-STAT signaling pathway
    ## hsa04142 Lysosome                                         hsa04142 Lysosome
    ## hsa00140 Steroid hormone biosynthesis hsa00140 Steroid hormone biosynthesis
    ## hsa04740 Olfactory transduction             hsa04740 Olfactory transduction
    ## hsa04916 Melanogenesis                               hsa04916 Melanogenesis
    ##                                         p.geomean stat.mean       p.val
    ## hsa04640 Hematopoietic cell lineage   0.002709366  2.857393 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.005655916  2.557207 0.005655916
    ## hsa04142 Lysosome                     0.008948808  2.384783 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.009619717  2.432105 0.009619717
    ## hsa04740 Olfactory transduction       0.014450242  2.239717 0.014450242
    ## hsa04916 Melanogenesis                0.022339115  2.023074 0.022339115
    ##                                           q.val set.size        exp1
    ## hsa04640 Hematopoietic cell lineage   0.3847887       49 0.002709366
    ## hsa04630 Jak-STAT signaling pathway   0.3847887      103 0.005655916
    ## hsa04142 Lysosome                     0.3847887      117 0.008948808
    ## hsa00140 Steroid hormone biosynthesis 0.3847887       26 0.009619717
    ## hsa04740 Olfactory transduction       0.4624078       39 0.014450242
    ## hsa04916 Melanogenesis                0.5297970       85 0.022339115


Now, let's try out the **pathview()** function from the [pathview package](https://bioconductor.org/packages/release/bioc/html/pathview.html) to make a pathway plot with our result shown in color. To begin with lets manually supply a `pathway.id` (namely the first part of the `"hsa04110 Cell cycle"`) that we could see from the print out above.

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

    ## Info: Downloading xml files for hsa04110, 1/1 pathways..
    ## Info: Downloading png files for hsa04110, 1/1 pathways..
    ## 'select()' returned 1:1 mapping between keys and columns
    ## Info: Writing image file hsa04110.pathview.png

This downloads the patway figure data from KEGG and adds our results to it. You can play with the other input arguments to **pathview()** to change the dispay in various ways including generating a PDF graph. For example:

``` r
# A different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa04110", kegg.native=FALSE)
```

Here is the default low resolution raster PNG output from the first pathview() call above:

![]({{ site.baseurl }}/class-material/hsa04110.pathview.png)

Note how many of the genes in this pathway are pertubed (i.e. colored) in our results.

Now, let's process our results a bit more to automagicaly pull out the top 5 upregulated pathways, then further process that just to get the IDs needed by the **pathview()** function. We'll use these KEGG pathway IDs for plotting below.

``` r
## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

    ## [1] "hsa04640" "hsa04630" "hsa04142" "hsa00140" "hsa04740"

Finally, lets pass these IDs in `keggresids` to the **pathview()** function to draw plots for all the top 5 pathways.

``` r
pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")
```



Here are the plots:  

![]({{ site.baseurl }}/class-material/hsa00140.pathview.png)

![]({{ site.baseurl }}/class-material/hsa04142.pathview.png)

![]({{ site.baseurl }}/class-material/hsa04630.pathview.png)

![]({{ site.baseurl }}/class-material/hsa04640.pathview.png)

![]({{ site.baseurl }}/class-material/hsa04740.pathview.png)



Section 3. Gene Ontology (GO)
=============================

We can also do a similar procedure with gene ontology. Similar to above, go.sets.hs has all GO terms. go.subs.hs is a named list containing indexes for the BP, CC, and MF ontologies. Let’s only do Biological Process.

``` r
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)
```

    ## $greater
    ##                                              p.geomean stat.mean
    ## GO:0007156 homophilic cell adhesion       4.893044e-05  3.971869
    ## GO:0060429 epithelium development         6.727999e-05  3.834578
    ## GO:0007610 behavior                       2.171759e-04  3.534089
    ## GO:0048729 tissue morphogenesis           2.471263e-04  3.498950
    ## GO:0002009 morphogenesis of an epithelium 3.227727e-04  3.429293
    ## GO:0016337 cell-cell adhesion             8.194676e-04  3.163087
    ##                                                  p.val     q.val set.size
    ## GO:0007156 homophilic cell adhesion       4.893044e-05 0.1337863      107
    ## GO:0060429 epithelium development         6.727999e-05 0.1337863      478
    ## GO:0007610 behavior                       2.171759e-04 0.2457053      404
    ## GO:0048729 tissue morphogenesis           2.471263e-04 0.2457053      403
    ## GO:0002009 morphogenesis of an epithelium 3.227727e-04 0.2567334      326
    ## GO:0016337 cell-cell adhesion             8.194676e-04 0.3753986      318
    ##                                                   exp1
    ## GO:0007156 homophilic cell adhesion       4.893044e-05
    ## GO:0060429 epithelium development         6.727999e-05
    ## GO:0007610 behavior                       2.171759e-04
    ## GO:0048729 tissue morphogenesis           2.471263e-04
    ## GO:0002009 morphogenesis of an epithelium 3.227727e-04
    ## GO:0016337 cell-cell adhesion             8.194676e-04
    ## 
    ## $less
    ##                                             p.geomean stat.mean
    ## GO:0000279 M phase                       1.582159e-16 -8.314874
    ## GO:0048285 organelle fission             8.120979e-16 -8.149796
    ## GO:0000280 nuclear division              2.314155e-15 -8.024006
    ## GO:0007067 mitosis                       2.314155e-15 -8.024006
    ## GO:0000087 M phase of mitotic cell cycle 6.404776e-15 -7.881237
    ## GO:0007059 chromosome segregation        1.055849e-11 -6.988384
    ##                                                 p.val        q.val
    ## GO:0000279 M phase                       1.582159e-16 6.292245e-13
    ## GO:0048285 organelle fission             8.120979e-16 1.614857e-12
    ## GO:0000280 nuclear division              2.314155e-15 2.300848e-12
    ## GO:0007067 mitosis                       2.314155e-15 2.300848e-12
    ## GO:0000087 M phase of mitotic cell cycle 6.404776e-15 5.094359e-12
    ## GO:0007059 chromosome segregation        1.055849e-11 6.998521e-09
    ##                                          set.size         exp1
    ## GO:0000279 M phase                            492 1.582159e-16
    ## GO:0048285 organelle fission                  373 8.120979e-16
    ## GO:0000280 nuclear division                   349 2.314155e-15
    ## GO:0007067 mitosis                            349 2.314155e-15
    ## GO:0000087 M phase of mitotic cell cycle      359 6.404776e-15
    ## GO:0007059 chromosome segregation             141 1.055849e-11
    ## 
    ## $stats
    ##                                           stat.mean     exp1
    ## GO:0007156 homophilic cell adhesion        3.971869 3.971869
    ## GO:0060429 epithelium development          3.834578 3.834578
    ## GO:0007610 behavior                        3.534089 3.534089
    ## GO:0048729 tissue morphogenesis            3.498950 3.498950
    ## GO:0002009 morphogenesis of an epithelium  3.429293 3.429293
    ## GO:0016337 cell-cell adhesion              3.163087 3.163087


Section 4. Reactome Pathway Analysis
------------------------------------

Reactome is database consisting of biological molecules and their relation to pathways and processes. Reactome, such as many other tools, has an online software available (<https://reactome.org/>) and R package available (<https://bioconductor.org/packages/release/bioc/html/ReactomePA.html>).  

If you would like more information, the documentation is available here: <https://reactome.org/user/guide>

Let's now conduct over-representation enrichment analysis and pathway-topology analysis with Reactome using the previous list of significant genes generated from our differential expression results above.

First, Using R, output the list of significant genes at the 0.05 level as a plain text file:

``` r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:", length(sig_genes)))
```

    ## [1] "Total number of significant genes: 8151"

``` r
write.table(sig_genes, file="significant_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Then, to perform pathway analysis online go to the Reactome website (<https://reactome.org/PathwayBrowser/#TOOL=AT>). Select “choose file” to upload your significant gene list. Then, select the parameters “Project to Humans”, then click “Analyze”.  

> **Question**: What pathway has the most significant “Entities p-value”? Do the most significant pathways listed match your previous KEGG results? What factors could cause differences between the two methods?  


Section 5. GO Analysis
----------------------

Gene Set Gene Ontology (GO) Enrichment is a method to determine over-represented or under-represented GO terms for a given set of genes. GO terms are formal structured controlled vocabularies (ontologies) for gene products in terms of their biological function. The goal of this analysis is to determine the biological process the given set of genes are associated with.

To perform Gene Set GO Enrichment online go to the website <http://www.geneontology.org/page/go-enrichment-analysis>. Paste your significant gene list from section 4. Then, select “biological process” and “homo sapiens”, and click submit.

> **Question**: What pathway has the most significant “Entities p-value”? Do the most significant pathways listed match your previous KEGG results? What factors could cause differences between the two methods?  



Bounus: Gene clustering, heatmaps and PCA
----------------------

Many statistical methods for analysis of multidimensional data, for example *clustering* and *principal components analysis* (PCA), work best for data that generally has the same range of variance at different ranges of the mean values.  

However, for counts from RNA-seq the expected variance grows with the mean. For example, if one performs PCA or clustering directly on a matrix of counts then the results will be heavely influenced by the genes with the highest counts (because they show the largest absolute differences between samples).  To address this problem the *DESeq2* package offers the **vst()** function (that stands for Variance Stabilizing Transformation) that stabilizes the variance of count data across different mean values. We will use **vst()** here as input for our clustering.


```r
vsd <- vst(dds, blind = FALSE)
```

Since gene clustering is only really relevant for genes that actually carry a signal, one usually would only cluster a subset of the most highly variable genes. Here, for demonstration purposes we select the 20 genes with the highest variance across samples.



```r
library("genefilter")

#row.variance <- apply(assay(vsd), 1, var)
row.variance <- rowVars(assay(vsd))
ord.variance <- order( row.variance, decreasing = TRUE) 

# Focus on top 20 most variable genes for demo purposes
mat  <- assay(vsd)[ ord.variance[1:20], ]
```

The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s average across all samples. To do this we center each genes’ values across samples by subtracting their mean values, and then plot the heatmap (figure below). 


```r
library(pheatmap)
mat.center  <- mat - rowMeans(mat)
pheatmap(mat.center)
```

![]({{ site.baseurl }}/class-material//unnamed-chunk-27-1.png)<!-- -->


> **Side-note:** We can also do PCA with our `vsd` object.


```r
pcaData <- plotPCA(vsd, intgroup="condition", returnData = TRUE)
pcaData
```

```
##                 PC1        PC2         group     condition      name
## SRR493366 -12.74667  0.3046597 control_sirna control_sirna SRR493366
## SRR493367 -13.04690 -0.0278821 control_sirna control_sirna SRR493367
## SRR493368 -13.28760 -0.3782467 control_sirna control_sirna SRR493368
## SRR493369  13.38297 -0.1131647      hoxa1_kd      hoxa1_kd SRR493369
## SRR493370  11.48330  1.1139575      hoxa1_kd      hoxa1_kd SRR493370
## SRR493371  14.21490 -0.8993238      hoxa1_kd      hoxa1_kd SRR493371
```


```r
library(ggplot2)

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition) ) +
  geom_point(size =3) 
```

![]({{ site.baseurl }}/class-material//unnamed-chunk-28-1.png)



Session Information
-------------------

The `sessionInfo()` prints version information about R and any attached packages. It's a good practice to always run this command at the end of your R session and record it for the sake of reproducibility in the future.

``` r
sessionInfo()
```

    ## R version 3.4.1 (2017-06-30)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] gageData_2.16.0            gage_2.28.2               
    ##  [3] pathview_1.18.2            org.Hs.eg.db_3.5.0        
    ##  [5] AnnotationDbi_1.40.0       DESeq2_1.18.1             
    ##  [7] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
    ##  [9] matrixStats_0.52.2         Biobase_2.38.0            
    ## [11] GenomicRanges_1.30.1       GenomeInfoDb_1.14.0       
    ## [13] IRanges_2.12.0             S4Vectors_0.16.0          
    ## [15] BiocGenerics_0.24.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.3.1             bit64_0.9-7            splines_3.4.1         
    ##  [4] Formula_1.2-2          latticeExtra_0.6-28    blob_1.1.0            
    ##  [7] GenomeInfoDbData_1.0.0 yaml_2.1.16            RSQLite_2.0           
    ## [10] backports_1.1.2        lattice_0.20-35        digest_0.6.14         
    ## [13] RColorBrewer_1.1-2     XVector_0.18.0         checkmate_1.8.5       
    ## [16] colorspace_1.3-2       htmltools_0.3.6        Matrix_1.2-12         
    ## [19] plyr_1.8.4             XML_3.98-1.9           pkgconfig_2.0.1       
    ## [22] genefilter_1.60.0      zlibbioc_1.24.0        xtable_1.8-2          
    ## [25] scales_0.5.0           BiocParallel_1.12.0    htmlTable_1.11.2      
    ## [28] tibble_1.3.4           annotate_1.56.1        KEGGREST_1.18.0       
    ## [31] ggplot2_2.2.1          nnet_7.3-12            lazyeval_0.2.1        
    ## [34] survival_2.41-3        magrittr_1.5           memoise_1.1.0         
    ## [37] evaluate_0.10.1        KEGGgraph_1.38.0       foreign_0.8-69        
    ## [40] graph_1.56.0           tools_3.4.1            data.table_1.10.4-3   
    ## [43] stringr_1.2.0          munsell_0.4.3          locfit_1.5-9.1        
    ## [46] cluster_2.0.6          Biostrings_2.46.0      compiler_3.4.1        
    ## [49] rlang_0.1.6            grid_3.4.1             RCurl_1.95-4.10       
    ## [52] rstudioapi_0.7         htmlwidgets_1.0        bitops_1.0-6          
    ## [55] base64enc_0.1-3        rmarkdown_1.8          gtable_0.2.0          
    ## [58] DBI_0.7                R6_2.2.2               gridExtra_2.3         
    ## [61] knitr_1.18             bit_1.1-12             Hmisc_4.1-1           
    ## [64] rprojroot_1.3-2        Rgraphviz_2.22.0       stringi_1.1.6         
    ## [67] Rcpp_0.12.15           png_0.1-7              geneplotter_1.56.0    
    ## [70] rpart_4.1-12           acepack_1.4.1
