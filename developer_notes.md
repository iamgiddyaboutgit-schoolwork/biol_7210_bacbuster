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