# tardis-phylogenetics
Temporal And diveRsity Distribution Sampler (TARDiS) for Phylogenetics

## Quick start
Download TARDiS and make sure [dependencies](#Dependendencies-and-OS) are installed. For the quickest start, just run our example:

`[path/to/tardis]/tardis -s`

and the TARDiS explorer GUI will open in your default browser. Retrieve example data in `data/example` and click on `Run Tardis`.

### Input
To run TARDiS, you will need the following inputs:

* A genomic data set in fasta format (example: `data/example/aln.fa`)
* A distance matrix, i.e., a square matrix where you stored the genomic distances for you genomic sequence pairs. It can be a csv or rds file. Rows and columns should be named as the fasta headers (example: `data/example/jc.distance.precalc.csv`). Note that this function works for *aligned* fasta files.
* A metadata file in csv format. This file should include two columns, `Accession.ID`, with the fasta headers, and `Collection.date`, with the sampling date in the dd/mm/YYYY format (example: `data/example/metadata.csv`)

### GUI
For experimenting with TARDiS, you can run the GUI as explained [above](#Quick-start). You can retrieve example data in `data/example`.

![GUI](/shiny_local/gui.png)

Important: **this GUI is intended for experimenting with small sets**, and small GA populations. All the GUI outputs are stored as `shiny_local/output/jc.distance.precalc.rds`. If you don't have a distance file, the `Jukes-Cantor distance` calculates the genomic distance. This distance files is to used as part of the TARDiS input. For large data sets, please use the command line instead. 

### Local command line ([Nextflow](https://www.nextflow.io/))
On a single machine, you can run TARDiS as

1. Run `[path/to/tardis]/tardis example`
2. Profit! You results are in `[current/directory]/output/example`

Run `tardis --help` to check out all the running options. Note that for reproducibility of the results, the user can specify the seeds to be used in `data/seeds.txt`.

Depending on your resources, if you are considering a large solution population, you can use more than one batch, i.e, each GA generation will be processed in batches. Parameters can be specified by command line (see `tardis --help`), or by a configuration file, named after your data set (example `example.config`).

### hpc Command line ([Nextflow](https://www.nextflow.io/))
Particularly large data sets can be run on hpc. To do so, you need generate an a nextflow configuration file for the specific workload manager in use. We provide an example file for SLURM in `hpc.nextflow.config`.

#### Batches
To ease the calculation burden for large data sets, data can be split into batches. In the configuration file, `params.gensize` defines the individual per batch, while `params.nbatches`
defines the number of batches. So to have 500K population split into 50 batches of 10K individuals each, you can set `params.gensize = 10000`, and `params.nbatches = 50`. Note that this will submit to your workload manager
50 jobs (1 per batch) for each generation. When all the jobs in the first generation are  complete, the 50 batches for the new generation will be submitted, and so on.

## Yes but what TARDiS does exactly?
TARDiS subsamples genomic data sets optimizing genomic diversity and temporal sampling according to parameters set by the users.
For a detailed discussion please be patient, a preprint and/or a user manual will be out soon.

## Dependendencies and OS
To run TARDiS, please install
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

## Contacts
* [Simone Marini](https://github.com/smarini)
* [Alberto Riva](https://github.com/albertoriva)

## License
TARDiS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See <https://www.gnu.org/licenses/>.
