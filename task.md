Immunology project for computation biologists
================

The goal of this project is to design a pooled CRISPR screen to study an
immunological process of your
choice.

## Introduction to pooled CRISPR screens combined with single-cell RNA sequencing

In a pooled CRISPR screen, cells are transfected by a library of vectors
that all target different genes in different ways. Each of these vectors
is barcoded so that we can later derive which cell was targeted by which
vector. A couple of hours after transfection, cells are profiled using
single-cell (RNA) sequencing techniques (Wagner, Regev, and Yosef 2016).
After some computational analyses, we derive for each cell it’s
expression profile and the vectors by which the cells were transfected.

Most often, a gene is knocked-out (loss-of-function), for example
Perturb-Seq (Dixit et al. 2016) or CROP-Seq (Datlinger et al. 2017)).
But genes can also be targeted in other ways, for example by
overexpression (gain-of-function) using CRISPRa (Wang et al. 2019). The
possibilities are basically endless; one might imagine that libraries
will be designed to target individual point mutations (Anzalone et al.
2019).

Most screens until now were performed *in vitro*, but exploring *in
vivo* screens is probably not far off. The main limiting factor here is
making sure that the correct cells are transfected; one solution for
immunology might be to do an adoptive transfer of bone marrow cells
which were transfected.

## The project

Doing a pooled screen can be very expensive. To get enough signal over
noise, you need to profile enough cells (thousands) for each vector, so
checking the effect of 100 vectors requires you to profie hundreds of
thousands of cells, and with about $0.5-$1 per cell, …. you can do the
calculation\!

To lower the cost, we typically want to do a targeted screen, where we
narrow down to a subset of genes, between 10 and 50. This may be based
on function (e.g. transcription factors) and expression pattern
(e.g. genes specific to alveolar macrophages).

In this example, we’re going to design a screen for transcription
factors that may be important for alveolar macrophages. For the project,
you should design a screen for a different topic. Feel free to explore
and be creative\!

## Selecting potential genes based on function

``` r
library(org.Mm.eg.db)
```

<GO:0003700> is the Gene Ontology identifier for “DNA-binding
transcription factor activity”

You can explore the Gene Ontology at
<https://www.ebi.ac.uk/QuickGO/term/GO:0003700>

Our expression dataset uses gene symbols, so we have to convert the
entrez ids

``` r
tfs_entrezid <- as.list(org.Mm.egGO2ALLEGS)$`GO:0003700`
tfs <- unique(AnnotationDbi::select(org.Mm.eg.db, keys=tfs_entrezid, columns=c("SYMBOL"), keytype="ENTREZID")$SYMBOL)
```

    ## 'select()' returned many:1 mapping between keys and columns

``` r
length(tfs)
```

    ## [1] 1001

## Selecting potential genes based on expression

``` r
library(limma)
```

### Using ImmGen RNA-seq

The ImmGen consortium has generated bulk gene expression dataset of a
lot of immune cells. These datasets are generated using a standardized
protocol to limit the batch effects, and the exact FACS gating
strategies are (in most cases) provided so that the sorting can be
reproduced
(e.g. <http://rstats.immgen.org/Skyline/resources/Sorting_PDFs_OSMNP/Mo.6Chi11b+.APAP.36h.Lv.pdf>).

The most extensive dataset was generated using microarrays:

  - To browse the data: <http://www.immgen.org/databrowser/index.html>
  - Data can be downloaded at
    <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15907> and
    <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37448>

Nowadays, these is also an extensive RNA-seq dataset available:

  - Browsing: <http://rstats.immgen.org/Skyline/skyline.html>
  - ImmGen ULI RNASeq, containing an overview of all immune populations:
    <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127267>
  - ImmGen MNP OpenSource, a community-based dataset because of the
    extensive heterogeneity and the expertise required for correct
    sorting (“Open-Source ImmGen: Mononuclear Phagocytes” 2016):
    <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122108>

We download the MNP dataset here:

``` r
# geo <- getGEO("GSE122108")
if (!file.exists("data/GSE122108_Gene_count_table.csv.gz")) {
  curl::curl_download("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122108&format=file&file=GSE122108_Gene_count_table.csv.gz", "data/GSE122108_Gene_count_table.csv.gz")
}
```

``` r
# read in the counts matrix
counts_df <- read.table("data/GSE122108_Gene_count_table.csv.gz", header = TRUE, stringsAsFactors = FALSE)
counts <- as.matrix(counts_df[, 2:ncol(counts_df)])
rownames(counts) <- counts_df$gene_symbol
# colnames(counts) <- gsub('-', 'min', colnames(counts)) # R doesnt like - or + signs in names
# colnames(counts) <- gsub('\\+', 'plus', colnames(counts))

# remove genes that are very lowly expressed
counts <- counts[apply(counts, 1, max) >= 10, ]

# create a dataframe with information on each sample
sample_info <- data.frame(
  sample_id =  colnames(counts),
  population_id = gsub("(.*)\\.[0-9]*", "\\1", colnames(counts)),
  stringsAsFactors = FALSE
)
```

``` r
head(sample_info, 20)
```

    ##              sample_id     population_id
    ## 1  MF.64pLYVEpIIn.Ao.1 MF.64pLYVEpIIn.Ao
    ## 2  MF.64pLYVEpIIn.Ao.2 MF.64pLYVEpIIn.Ao
    ## 3  MF.64pLYVEpIIn.Ao.3 MF.64pLYVEpIIn.Ao
    ## 4  MF.64pLYVEpIIp.Ao.1 MF.64pLYVEpIIp.Ao
    ## 5  MF.64pLYVEpIIp.Ao.2 MF.64pLYVEpIIp.Ao
    ## 6  MF.64pLYVEpIIp.Ao.3 MF.64pLYVEpIIp.Ao
    ## 7  MF.64pLYVEnIIp.Ao.1 MF.64pLYVEnIIp.Ao
    ## 8  MF.64pLYVEnIIp.Ao.2 MF.64pLYVEnIIp.Ao
    ## 9  MF.64pLYVEnIIp.Ao.3 MF.64pLYVEnIIp.Ao
    ## 10            MF.PC.44             MF.PC
    ## 11            MF.PC.45             MF.PC
    ## 12            MF.PC.46             MF.PC
    ## 13           MF.F.PC.1           MF.F.PC
    ## 14           MF.F.PC.2           MF.F.PC
    ## 15           MF.F.PC.3           MF.F.PC
    ## 16            MF.PC.47             MF.PC
    ## 17            MF.PC.48             MF.PC
    ## 18            MF.PC.49             MF.PC
    ## 19            MF.PC.50             MF.PC
    ## 20            MF.PC.51             MF.PC

We can use limma-voom to find genes that are differentially expressed
between two populations. Let’s first explore how to get differentially
expressed genes between two populations of interest, for example Kupffer
cells and alveolar
macrophages:

``` r
population_ids_oi <- c("MF.alv.11cp64pSiglecFp.Lu","MF.KC.Clec4FpTim4p64p.Lv")
```

Plot the genes in a heatmap. We will scale the expression per gene.
![](task_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

We can do a differential expression for many different populations. The
genes of interest are then for example the genes that are upregulated in
our population of interest compared to all other “reference”
populations. In our case, the reference are all types of
macrophages.

``` r
populations_reference <- c("MF.KC.Clec4FpTim4p64p.Lv", "MF.F.PC", "MF.480p.SP")
population_oi <- "MF.alv.11cp64pSiglecFp.Lu"
differential_expression <- get_marker_genes(counts, sample_info, populations_reference, population_oi)
```

![](task_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

## Designing the CRISPR library

``` r
genes_oi <- intersect(markers, tfs)
cat(paste0(genes_oi, collapse = "\n"))
```

    ## Bhlhe41
    ## Cebpa
    ## Cebpb
    ## Cebpd
    ## Ctnnb1
    ## Foxf1
    ## Gtf3a
    ## Hes6
    ## Hes7
    ## Maff
    ## Mkx
    ## Ovol2
    ## Plscr1
    ## Rara
    ## Rfx2
    ## Runx1
    ## Runx2
    ## Snai3
    ## Spi1
    ## Tfeb
    ## Thap11
    ## Trerf1
    ## Zbtb7a
    ## Zfp296
    ## Zfp358

We’re lucky that we here get about 20 genes. What would you do if you
get 500 genes? You might for example be a bit stricter and select the 20
genes that have the highest expression in our population of interest.
What would you do if you get 2 genes? You can loosen the restrictions a
bit and allow that the genes are also expressed in one or two of the
reference populations.

We can submit this list to …
<https://portals.broadinstitute.org/gpp/public/analysis-tools/sgrna-design>
(or any related tool). This will search for appropriate gRNA sequences.
It will try to avoid off-targets by doing a BLAST search against the
genome of interest, and selecting those gRNAs that are only found once,
i.e. at the location of interest.

``` r
sgrna_design <- read_tsv("results/sgrna_designs.txt")
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_character(),
    ##   Quota = col_double(),
    ##   `Target Taxon` = col_double(),
    ##   `Target Gene ID` = col_double(),
    ##   `Target Alias` = col_logical(),
    ##   `Initial Spacing Requirement` = col_double(),
    ##   `Off-Target Match Ruleset Version` = col_double(),
    ##   `Off-Target Tier Policy` = col_double(),
    ##   `sgRNA Cut Position (1-based)` = col_double(),
    ##   `Exon Number` = col_double(),
    ##   `Target Cut Length` = col_double(),
    ##   `Target Total Length` = col_double(),
    ##   `Target Cut %` = col_double(),
    ##   `On-Target Efficacy Score` = col_double(),
    ##   `On-Target Rank` = col_double(),
    ##   `Off-Target Rank` = col_double(),
    ##   `On-Target Rank Weight` = col_double(),
    ##   `Off-Target Rank Weight` = col_double(),
    ##   `Combined Rank` = col_double(),
    ##   `Pick Order` = col_double(),
    ##   `Picking Round` = col_double()
    ## )

    ## See spec(...) for full column specifications.

Let’s choose for each target gene one guide RNA:

``` r
sgrna_design %>% 
  group_by(Input) %>% 
  top_n(1, `Combined Rank`) %>% 
  dplyr::select(`Target Gene Symbol`, `sgRNA Sequence`)
```

    ## Adding missing grouping variables: `Input`

    ## # A tibble: 19 x 3
    ## # Groups:   Input [19]
    ##    Input  `Target Gene Symbol` `sgRNA Sequence`    
    ##    <chr>  <chr>                <chr>               
    ##  1 Bcl6b  Bcl6b                GCCTTGGGGTCTGGGCTGGC
    ##  2 Erg    Erg                  GAGGGTGGGGCTGCAGGGCC
    ##  3 Gata4  Gata4                CTGCCGCCGCTGCCGCAGCC
    ##  4 Hey1   Hey1                 CTATCGGAGTTTGGGGTTTC
    ##  5 Hic1   Hic1                 GGCGGCGGCGGTGGCCCGGC
    ##  6 Hnf4a  Hnf4a                GGTGAGGGTGCAGGGGGTGG
    ##  7 Id3    Id3                  GTCCTGGCAGAGCCGGCGCC
    ##  8 Ifi205 Ifi205               GGGCTGTGGAAGTCTCTTCC
    ##  9 Ifi208 Ifi208               ACACTGCTGGGCTCTGTTTT
    ## 10 Meis2  Meis2                CCATGGCTGGGTGGTGGGGA
    ## 11 Nfib   Nfib                 TGCAGGAAGGATGGGTCTCT
    ## 12 Nfxl1  Nfxl1                CCCGGGCCCCCGGGGGCTGC
    ## 13 Nr1h4  Nr1h4                GTCTGTGGAGACAGGGCCTC
    ## 14 Nr2f2  Nr2f2                CAGGGCGGCCCTGGCGGCCC
    ## 15 Ppara  Ppara                CTTCAGATAAGGGACTTTCC
    ## 16 Rorc   Rorc                 TCTCTGTGGGGCCCTGTCCA
    ## 17 Smad6  Smad6                GGAGTCGGGGGCCGGGGCTG
    ## 18 Sox18  Sox18                GCTGGACGGGGAGGCGGGCG
    ## 19 Tead4  Tead4                CGGTGCGGAGGGTGAGGGGG

Time to start cloning\!

## Analysis

To analyze the data, you need to know for each read the cell to which it
belongs (cell barcode) and how this cell was perturbed (gRNA barcode).

Analyzing the data happens in two steps: preprocessing and modelling.

In the preprocessing step, you count the number of molecules that are
part of each cell. You also assign to each cell with which vector it was
transfected, so that you know which genes were affected in your cell.

In the modelling step, you have a look at the effect of the perturbation
on your cell of interest. These analyses can be:

  - Finding modules of genes that are affected by different
    perturbations. These modules can have functions, e.g. a module that
    is related to lipid metabolism.
  - If you combine perturbations of different genes, you can also
    analyze interactions between genes, e.g. whether they work
    synergistically, antagonistically, …

The current state-of-the-art of CRISPR screen modelling is probably:
<https://science.sciencemag.org/content/365/6455/786/tab-figures-data>
together.

## References

<div id="refs" class="references">

<div id="ref-anzaloneSearchandreplaceGenomeEditing2019">

Anzalone, Andrew V., Peyton B. Randolph, Jessie R. Davis, Alexander A.
Sousa, Luke W. Koblan, Jonathan M. Levy, Peter J. Chen, et al. 2019.
“Search-and-Replace Genome Editing Without Double-Strand Breaks or
Donor DNA.” *Nature*, October, 1–1.
<https://doi.org/10.1038/s41586-019-1711-4>.

</div>

<div id="ref-datlingerPooledCRISPRScreening2017">

Datlinger, Paul, André F. Rendeiro, Christian Schmidl, Thomas
Krausgruber, Peter Traxler, Johanna Klughammer, Linda C. Schuster,
Amelie Kuchler, Donat Alpar, and Christoph Bock. 2017. “Pooled CRISPR
Screening with Single-Cell Transcriptome Readout.” *Nature Methods* 14
(3): 297–301. <https://doi.org/10.1038/nmeth.4177>.

</div>

<div id="ref-dixitPerturbseqDissectingMolecular2016">

Dixit, Atray, Oren Parnas, Biyu Li, Jenny Chen, Charles P. Fulco, Livnat
Jerby-Arnon, Nemanja D. Marjanovic, et al. 2016. “Perturb-Seq:
Dissecting Molecular Circuits with Scalable Single Cell RNA Profiling of
Pooled Genetic Screens.” *Cell* 167 (7): 1853–1866.e17.
<https://doi.org/10.1016/j.cell.2016.11.038>.

</div>

<div id="ref-OpensourceImmGenMononuclear2016">

“Open-Source ImmGen: Mononuclear Phagocytes.” 2016. *Nature Immunology*
17 (7): 741–41. <https://doi.org/10.1038/ni.3478>.

</div>

<div id="ref-wagner_revealing_2016">

Wagner, Allon, Aviv Regev, and Nir Yosef. 2016. “Revealing the Vectors
of Cellular Identity with Single-Cell Genomics.” *Nature Biotechnology*
34 (11): 1145–60. <https://doi.org/10.1038/nbt.3711>.

</div>

<div id="ref-wangMultiplexedActivationEndogenous2019">

Wang, Guangchuan, Ryan D. Chow, Zhigang Bai, Lvyun Zhu, Youssef Errami,
Xiaoyun Dai, Matthew B. Dong, et al. 2019. “Multiplexed Activation of
Endogenous Genes by CRISPRa Elicits Potent Antitumor Immunity.” *Nature
Immunology* 20 (11): 1494–1505.
<https://doi.org/10.1038/s41590-019-0500-4>.

</div>

</div>
