# BacBuster

Foodborne bacterial pathogens are a serious public health concern [[1](https://doi.org/10.3934%2Fmicrobiol.2017.3.529)].  Particularly concerning are bacterial strains which have acquired resistance to antibiotics, antiseptics, and disinfectants.  Many public health surveillance programs, such as [PulseNet](https://www.cdc.gov/pulsenet/index.html) at the CDC, exist to "help identify unsafe foods and production processes" [[2](https://www.cdc.gov/pulsenet/next-gen-wgs.html)].  These programs take samples from humans (usually feces) and the environment (such as food processing machinery).  Typically [[3](https://www.cdc.gov/foodnet/reports/cidt-questions-and-answers.html)], isolates from these samples are analyzed via whole genome sequencing (WGS).  

WGS allows scientists to identify genetic similarities among new and previously characterized isolates.  WGS provides the resolution needed to cluster samples and to infer bacterial phylogenies. [[4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5608882/)] Large-scale genomic information allows outbreaks to be linked to their source, even if that source is hundreds of miles away. [[5](https://www.cdc.gov/ncezid/dfwed/keyprograms/tracking-foodborne-illness-wgs.html)] (The assumption is that if two organisms share very similar genomes, then they must be close in evolutionary history.  In practice, [bacterial recombination](https://en.wikipedia.org/wiki/Bacterial_recombination) must also be considered.)  Additionally, WGS allows for the characterization of the problematic pathogens because it allows genes for antibacterial resistance and virulence factors to be identified.

The challenge of the bioinformatician is to analyze WGS data.  This is no easy task considering the many steps required.  **BacBuster** is presented as an end-to-end genomic epidemiology pipeline. [[6](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5572866/)]  

The goals of **BacBuster** are two-fold:
1. Clustering of isolates by genomic features and source
2. Characterization of those clusters

**BacBuster** gives epidemiologists the knowledge they need to control bacterial pathogens.

The workflow of **BacBuster** is diagrammed below: 

# Running BacBuster
Once Nextflow has been installed on your system, BacBuster can be run using the command

```
nextflow run BacBuster.nf
```
Nextflow channels will use skesa for genome assembly, prodigal for gene prediction, AMR finder plus for gene annotation, and snippy for comparative genomics functions,

**Skesa**
 [[7](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1540-z)] is a DeBruijn graph-based de-novo assembler designed for assembling microbial genomes. They have a higher sequence quality assemblies compared to other genome assembly tools. NCBI has used SKESA to assemble 272,000+ reads from in SRA.  

We can download skesa in out local directory using the following command 
```
$ git clone https://github.com/ncbi/SKESA
```
For more instructions, refer the  [documentation](https://github.com/ncbi/SKESA) 

**Prodigal**
 [[8](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)] is a gene finding algorithm that uses a combination of Critica and Glimmer, BLAST to locate missing genes and correct errors, followed by manual expert curation. At the start of the program, it scans through the entire input sequence and checks the prevalence of G's and C's in each of the three codon positions in every open reading frame (ORF). The codon position with the highest GC content for each ORF is labeled the "winner", and a running total for that position is added up. After all ORFs have been evaluated, the sum provides an estimated indication of the tendency of each codon position to favor G's and C's. 

**Snippy**
 [[9](https://github.com/tseemann/snippy)] finds SNPS between a haploid reference genome and NGS sequence reads. The Snippy-multi script can be used to simplify running a set of isolate sequences against the same reference. The script requires a tab-separated input file containing a sequence IDs and paths for the isolate sequences.

Overall, our nextflow pipeline will accomodate all these above-mentioned tools in various channels and the framework will be connected to the various functionalities in our predictive web server, created using [Streamlit](https://streamlit.io/). The user will be presented with a GUI that enables them to run these steps separately or in order and will have the option to mail the results to their email.
