I downloaded mambaforge from [here](https://github.com/conda-forge/miniforge#mambaforge).

I added a default channel to mamba: 
```
conda config --add channels bioconda
```

The bioconda channel hosts a lot of biology-relevant software. 

Conda/mamba environments can be exported using a command similar to
```
mamba env export | head -n -1 > ./Team3-WebServer_falco_env.yml
```

Note that the head command is used to remove the last line which contains prefix information unique to your local setup which may not be relevant for other users.

Useful information on [parallelization](https://www.nextflow.io/docs/latest/faq.html?highlight=parallel#how-do-i-process-multiple-input-files-in-parallel).  

# Installing Java using [SDKMAN!](https://sdkman.io/install)
This will install sdkman to a hidden file in the ```${HOME}``` directory.
```
curl -s "https://get.sdkman.io" | bash
```
and then 
```
source "$HOME/.sdkman/bin/sdkman-init.sh"
```
Verify installation with
```
sdk version
```
Finally, install Amazon Corretto (a distribution of the Open Java Development Kit).
```
sdk install java 17.0.6-amzn
```

# Nextflow Notes
Upon invocation within a directory, nextflow creates a project specific .nextflow.log file, .nextflow cache directory as well as a work directory.  The output of the pipeline will be stored in the work directory. 

# Prodigal Notes
Prodigal pairs nicely with skesa (which tends to produce
a lot of contigs) because "By default, Prodigal's parameters 
are ideal for scaffolds and/or multiple FASTA with many contigs."

"TIP: You should be careful using the -c option with draft 
genomes in many contigs, as this will prevent Prodigal from 
predicting partial genes. Similarly, if you have a single 
scaffold with many gaps in it, you should be careful using 
the -e option, as you may also lose many partial genes."

How does Prodigal handle multiple sequences in the draft genome?
"If we encounter multiple sequences, we insert TTAATTAATTAA between 
each one to force stops in all six frames."
"If the genome consists of multiple chromosomes, you can analyze 
them together or separately. Chromosomes should only be separated 
if (1) each chromosome is at least 500kb, and (2) you have reason 
to believe the chromosomes are quite different in terms of GC content,
RBS motif usage, and other parameters."
By default, "partial genes are allowed to run into gaps of N's, 
which means you should get the same results analyzing 1000 contigs 
in one file, or analyzing one scaffold with the 1000 contigs joined 
together by runs of N's."

"Prodigal contains no special routines to deal with viruses. 
As such, it cannot handle certain phenomena that occur sometimes 
in viruses, such as translational frame shifts. Viruses should 
generally be analyzed as above, with short genomes analyzed in 
anonymous mode and longer ones in normal mode."

# AMRFinderPlus Notes
After installing into an environment, the database can be updated using ```amrfinder --update```.  The database will be installed to the same directory as your environment.

# General, Uncategorized Notes
The raw read files can be checked for corruption using the md5 check sums using the following command within the appropriate directory.
```
md5sum -c *.md5
``` 
