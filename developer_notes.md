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

# Nextflow Notes
Upon invocation within a directory, nextflow creates a project specific .nextflow.log file, .nextflow cache directory as well as a work directory.  The output of the pipeline will be stored in the work directory. 

# General, Uncategorized Notes
The raw read files can be checked for corruption using the md5 check sums using the following command within the appropriate directory.
```
md5sum -c *.md5
``` 
