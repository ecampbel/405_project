# Meeting Minutes & Agenda items, weekly

[Oct 20](#oct-20) ; [Oct 27](#oct-27) ; [Oct 31](#oct-31)

## Oct 31

To-Do:

- **Eleanor**: Learn [edgeR](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) for differential expression
- **Eleanor**: Recount data in HTSeq
- **Imogen**: Run STAR to align the unaligned sequences to chloroplast genome
- **Fares**: Create visualisations for pathway analysis – within a pathway, which proteins are up/downregulated (can look at KEGG outputs, last column in Fares' code includes geneID **KO terms** for everything in pathway)

### Thus far

Our study focuses on the juvenile aposymbiotic (APO) and stable symbiosis (10D) RNA-seq data from Chan et al (2018). The reads from of the three biological replicates were trimmed using trimmomatic to remove Illumina adaptors and low-quality reads. They were then aligned using STAR to the E. chlorotica genome (assembled by Cai et al., 2019).

The unmapped sequences resulting from the original STAR alignment were also mapped to the chloroplast genome (source?) to investigate whether kleptoplasty stabilisation resulted in expression of plastid genes.


## Oct 27

Imogen/Eleanor: 

- find gene name (EGW) –> accession number data in order to parse Fares' gene table.
- draft intro?
- R fiddling with gene counts 

Joe:

- Try mapping reads to chloroplast reference genome next week some time
- In RNAseq file there is an ummapped reads file

## Oct 20

Assigned tasks to be completed by Oct 27:

- Eleanor: find comparison slug data (non photosynthetic?), make github repository
- Fares & Joe: trim reads
- Imogen: look through gene count data
- ??? : figure out data normalisation

**Overview of timeline pushing forwards:**

Ot 30th: Have draft pipeline in R ready for visualization and analysis

Nov 15th: Have data good enough to begin writing draft

Nov 27th: First draft of report, slides due

Dec 7th: Final report due

