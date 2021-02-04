# TARDiS for Phylogenetics
Temporal And diveRsity Distribution Sampler (TARDiS) for Phylogenetics

# Quick start
Download TARDis and make sure [dependencies](#dependencies-and-os) are installed. For the quickest start, just run our example:

`[path/to/tardis]/tardis -s`

and the TARDiS explorer GUI will open in your default browser. Retrieve example data in `data/example` and click on `Run Tardis`.

# Yes but what does TARDiS do exactly?
TARDiS subsamples genetic data sets optimizing genetic diversity and temporal sampling according to parameters set by the users. The optimization is driven by a genetic algorithm.

# Citation
A paper describing TARDiS principles and application is

* [Marini S, et al. Optimizing viral genome subsampling by genetic diversity and temporal distribution (TARDiS) for Phylogenetics](https://www.biorxiv.org/content/10.1101/2021.01.15.426832v1)

In an early, simplified version, TARDiS principles were also applied in:

* [Lednicky J, et al. Earliest Detection to Date of SARS-CoV-2 in Florida: Identification Together With Influenza Virus on the Main Entry Door of a University Building, February 2020](https://www.researchsquare.com/article/rs-87486/v1)

# Running TARDiS
To run TARDiS, you will need the following inputs:

* A genetic sequence alignment in fasta format (example: `data/example/aln.fa`)
* A distance matrix, i.e., a square matrix where you stored the genetic distances for you genetic sequence pairs. It can be a csv or rds file. Rows and columns should be named as the fasta headers (example: `data/example/jc.distance.precalc.csv`). Note that this function works for *aligned* fasta files
* A metadata file in csv format. This file should include two columns, `Accession.ID`, with the fasta headers, and `Collection.date`, with the sampling date in the dd/mm/YYYY format (example: `data/example/metadata.csv`)

## GUI
For experimenting with TARDiS, you can run the GUI as explained [above](#quick-start). You can retrieve example data in `data/example`.

![GUI](/shiny_local/gui.png)

Important: **this GUI is intended for experimenting with small sets**, and small GA populations. All the GUI outputs are stored as `shiny_local/output/jc.distance.precalc.rds`. If you don't have a distance file, the `Jukes-Cantor distance` calculates the genetic distance. This distance files is to used as part of the TARDiS input. For larger data sets, please use the command line version instead. 

## Command line (requires [Nextflow](https://www.nextflow.io/))
You can run TARDiS from the command line as follows:

1. Run `[path/to/tardis]/tardis [path/to/tardis]/example/example.config`
2. Profit! You results are in `[current/directory]/output/example`

Run `tardis -h` to display all available options. Note that to ensure reproducibility of the results, the user can specify the randomization seeds to be used in `data/seeds.txt`.

### Config file
You must specify the parameters for your run in a config file. Use the format of the following example:
```
params.data_set = "example"
params.nsamples = 4
params.gensize = 100
params.nbatches = 1
params.ncores = 2
params.ngenerations = 10
params.fracnew = 0.14
params.fracevolved = 0.85
params.fracelite = 0.01
params.wdiv = 1
params.wtem = 1
params.distopt = "max"
```
Settable parameters in the config file are: params.data_set (name of the data set), params.nsamples (number of genomes in the subsample), params.gensize (size of the generation per batch, see [below](#batches), default 1 batch), params.nbatches (number of batches, see below), params.ncores (number of cores for parallel computing), params.ngenerations (number of generations), params.fracnew (fraction of newly generated individuals per generation), params.fracevolved (fraction of evolved individuals per generation), params.fracelite (fraction of elite individuals per generation; elite individuals are the ones with the highest fitness, to be copied as they are in the new generation), params.wdiv (weight of the genetic diversity), params.wtem (weight of the time distribution), params.distopt (target of genetic diversity optimization, "max", "median", or "mean" of the initial population).

### HPC and resource allocation
The default NextFlow execution profile (option -p) is "local", which uses the local machine directly. In an HPC environment, you can use the "small" (5 GB), "medium" (30 GB), or "large" (128 GB) profiles, which assume the presence of the SLURM scheduler.

### Group mode
TARDiS can be run in group mode (option -g) from the command line. This is useful when the user has a large genomes file, with genomes pertaining to different groups, to be subsampled independently and combined into a single output at the end. For example, groups could correspond to different geographical regions. 
In this case, the config file will be in comma-delimited format, with one row for each group. Note that for group mode config files:
* Column names are [parameter names](#config-file)
* A special group column, called `group`, needs to be present, with a different value in each line
* The execution profile can be specified in the group mode file config in the "profile" column (this is not possible in a single run config file
* `NAs` are accepted (default values will be used)
* Groups that are NOT listed in the group-mode config file will be included in the final output as a whole, without being subsampled

Also note that the metadata file needs to include a `group` column to identify the group of each genome.
A working toy example is provided in `data/example_group`. There are three groups, a, b, and c. While groups a and b will be subsampled, the whole group c (absent from the config file) will be included in the output without being subsampled.

### Batches
To ease the calculation burden for large populations, data can be split into batches. Remember that params.gensize defines the number of individuals per batch, while params.nbatches
defines the number of batches. So to have a 500K population split into 50 batches of 10K individuals each, you can set `params.gensize = 10000`, and `params.nbatches = 50`. Note that this will submit 50 jobs (1 per batch) for each generation to your workload manager. When all jobs in the first generation are  complete, the 50 jobs for the next generation will be submitted, and so on.

# Dependencies and OS
To run TARDis, please install
* [R >= 3.6.1](https://www.r-project.org/)
* [Pyhton >= 3.7](https://www.python.org/) (works with Pyhton 2.7 as well)
* [optparse >= 1.6.6](https://cran.r-project.org/web/packages/optparse/index.html)
* [doRNG >= 1.8.2](https://cran.r-project.org/web/packages/doRNG/index.html)
* [dplyr >= 1.0.0](https://cran.r-project.org/web/packages/dplyr/index.html)
* [ggplot2 > 3.3.1](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [gridExtra >= 2.3](https://cran.r-project.org/web/packages/gridExtra/index.html)

For the GUI/explorer version, please install
* [Shiny >= 1.4.0.2](https://www.r-project.org/nosvn/pandoc/shiny.html)
* [directoryInput >= 137dc69](https://github.com/wleepang/shiny-directory-input)

For the local/hpc command line version, please install
* [Nextflow >= 20.01.0](https://www.nextflow.io/docs/latest/getstarted.html)

To calculate Jukes-Cantor distances
* [DECIPHER >= 2.14](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)

TARDiS has been successfully used on Linux Ubuntu (local, commandline), Chrome (GUI) and SLURM (hpc). Please let use know if you are using it on other platforms.

# Contacts
* [Simone Marini](https://github.com/smarini)
* [Alberto Riva](https://github.com/albertoriva)

# License
TARDiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See <https://www.gnu.org/licenses/>.

