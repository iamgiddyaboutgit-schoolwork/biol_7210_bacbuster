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
