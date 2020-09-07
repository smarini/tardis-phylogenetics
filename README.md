# tardis-phylogenetics
Temporal And diveRsity Distribution Sampler (TARDiS) for Phylogenetics

## Quick start
Make sure [dependencies](#Dependendencies-and-OS) are installed. For the quickest start, just run our example

2. Download TARDiS
3. Run `[path to tardis-phylogenetics]/tardis example`
4. Profit! You results are in `[path to tardis-phylogenetics]/output/example`


### Input
To run TARDiS, you will need the following inputs:

* A genomic data set in fasta format (example: `data/example/aln.fa`)
* A distance matrix, i.e., a square matrix where you stored the genomic distances for you genomic sequence pairs. It can be a csv or rds file. Rows and columns should be named as the fasta headers (example `data/example/jc.distance.precalc.csv`). Note that this function work for *aligned* fasta files
* A metadata file in csv format. This file should include two columns, `Accession.ID`, with the fasta headers, and `Collection.date`, with the sampling date in the dd/mm/YYYY format (example in `data/example/metadata.csv`)

### Shiny GUI
For experimenting with TARDiS, you can run the Shiny app in `shiny_local`, your results will be in `shiny_local/output`.

![Shiny GUI](/shiny_local/gui.png)

Important: **this GUI is intended for experimenting with small sets**, and small GA populations. All the GUI outputs are stored as `shiny_local/output/jc.distance.precalc.rds`. If you don't have a distance file, the `Jukes-Cantor distance` calculates the genomic distance. This distance files is to used as part of the TARDiS input. For large data sets, please use the command line instead. 

### Command line ([Nextflow](https://www.nextflow.io/), local)
On a single machine, you can run TARDiS as explained above. Depending on your resources, if you are considering a large solution population, you can use more than one batch, i.e, each GA generation will be processed in batches. Parameters can be specified by command line (see `tarids --help`), or by a configuration file, named after your data set (example `example.config`).

### Command line ([Nextflow](https://www.nextflow.io/), hpc)
Particularly large data sets can be run on hpc. To do so, you need generate an a nexflow configuration file for the specific workload manager in use. We provide an example file for SLURM in `hpc.nextflow.config`.

### Options
Run `tarids --help` to check out all the running options. Note that for reproducibility of the results, the user can specify the seeds to be used in `data/seeds.txt`.

## Yes but what TARDiS do exactly?
TARDiS subsamples genomic data sets optimizing genomic diversity and temporal sampling according to parameters set by the users. For a detailed discussion please be patient, a preprint will be out soon.

## Dependendencies and OS
To run local/hpc command line TARDiS, please install
* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
* [DECIPHER](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html)
* [optparse](https://cran.r-project.org/web/packages/optparse/index.html)
* [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)

* [Shiny](https://www.r-project.org/nosvn/pandoc/shiny.html) is needed for the GUI.

TARDiS has been successfully used on Ubuntu (local) and SLURM (hpc). Please let use know if you are using on other platforms.

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
