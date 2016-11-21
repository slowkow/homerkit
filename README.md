# homerkit

`homerkit` is an R package that implements functions to read [HOMER] output
files.

[HOMER]: http://homer.salk.edu/homer/

## Installation

How to install HOMER: <http://homer.salk.edu/homer/download.html>

```r
install.packages("devtools")
devtools::install_github("slowkow/homerkit")
```

## Usage

### 1. Run HOMER findMotifs.pl on your target genes

```bash
head -n3 target_genes.txt
ENSG00000003989
ENSG00000017427
ENSG00000028277
```

```bash
gene_file="target_genes.txt"
bg_file="background_genes.txt"
out_dir="output"

mkdir -p $out_dir
# Find motifs that are enriched in the promoters of your target genes.
findMotifs.pl $gene_file human $out_dir \
  -bg $bg_file &> ${out_dir}/run_homer.log
```

### 2. Run HOMER annotatePeaks.pl on every motif

```bash
# Find the target genes for each motif.
for motif in $out_dir/*/*.motif; do
  if [[ ! -f ${motif}.tsv ]]
  then
    annotatePeaks.pl tss hg38 \
      -size -500,250 -m $motif -list $gene_file \
      1> ${motif}.tsv 2> ${motif}.tsv.log
  fi
done
```

### 3. Read all of the HOMER output files with homerkit

```r
# install.packages("devtools")
# devtools::install_github("slowkow/homerkit")

library(homerkit)
h <- read_homer_output("output")
```

Novel motif target genes:

```r
head(split(h$novel_motif_peaks$gene_name, h$novel_motif_peaks$motif), 3)
```

```
$motif1
[1] "RERG"  "CSF3"  "CXCL6" "CXCL1" "CXCL5" "CXCL3" "CXCL2" "CSF2"  "ELF3" 

$motif10
[1] "IER3"  "MT1X"  "MMP3"  "CCL20"

$motif11
 [1] "IL6"    "CCL7"   "CXCL6"  "CXCL1"  "CXCL5"  "CXCL3"  "CXCL2"  "GPR183"
 [9] "NR4A2"  "PLD1" 
```

Possible transcription factors that match `motif1`:

```r
subset(h$novel_motif_tfs, motif == "motif1")
```

```
# A tibble: 10 × 8
                                                 match_name match_rank offset
                                                      <chr>      <dbl>  <dbl>
1  NFkB-p65-Rel(RHD)/ThioMac-LPS-Expression(GSE23622)/Homer          1      2
2                                      RELA/MA0107.1/Jaspar          2      2
3                                 MF0003.1_REL_class/Jaspar          3      2
4        NFkB-p65(RHD)/GM12787-p65-ChIP-Seq(GSE19485)/Homer          4      1
5                                       REL/MA0101.1/Jaspar          5      2
6                                     NFKB2/MA0778.1/Jaspar          6      1
7                                    PB0012.1_Elf3_1/Jaspar          7      4
8                                    NFATC1/MA0624.1/Jaspar          8      5
9                                    NFATC3/MA0625.1/Jaspar          9      5
10                                    NFKB1/MA0105.4/Jaspar         10      1
# ... with 5 more variables: orientation <chr>, score <dbl>, motif <chr>,
#   alignment1 <chr>, alignment2 <chr>
```

Known motif target genes:

```r
head(split(h$known_motif_peaks$gene_name, h$known_motif_peaks$motif), 3)
```

```
$known1
 [1] "MAP3K8" "CFB"    "CSF3"   "CXCL8"  "CXCL6"  "CXCL1"  "CXCL5"  "CXCL3" 
 [9] "CXCL2"  "NR4A2"  "ELF3"   "PID1"  

$known10
 [1] "RERG"            "SPECC1L-ADORA2A" "IER3"            "CFB"            
 [5] "SLC11A2"         "NR4A1"           "IL23A"           "MT1L"           
 [9] "CXCL8"           "CXCL1"           "CXCL3"           "CXCL2"          
[13] "FLVCR2"          "STEAP1"          "SERPINA9"        "AVPI1"          
[17] "GPR183"          "MMP3"            "PTGS2"           "ELF3"           
[21] "HSD11B1"         "CCL20"          

$known11
 [1] "CSF3"     "PIM2"     "MT1X"     "GAB2"     "SERPINA9" "IGF1"    
 [7] "IL1B"     "TNFAIP6"  "STAT4"    "ELF3"     "ACKR3"
```

Known transcription factors:

```r
head(unique(h$known_motif_peaks[,c("motif", "best_guess")]), 3)
```

```
# A tibble: 3 × 2
   motif                                               best_guess
   <chr>                                                    <chr>
1 known1 NFkB-p65-Rel(RHD)/ThioMac-LPS-Expression(GSE23622)/Homer
2 known2       NFkB-p65(RHD)/GM12787-p65-ChIP-Seq(GSE19485)/Homer
3 known3                             TATA-Box(TBP)/Promoter/Homer
```

## Contributing

Please [submit an issue][issues] to report bugs or ask questions.

Please contribute bug fixes or new features with a [pull request][pull] to this
repository.

[issues]: https://github.com/slowkow/homerkit/issues
[pull]: https://help.github.com/articles/using-pull-requests/

## Related work

- https://github.com/MalteThodberg/homeR
