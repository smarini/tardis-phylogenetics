#!/usr/bin/env nextflow

params.data_set = "USA"
params.nsamples = 10
params.gensize = 100
params.nbatches = 5
params.ncores = 2
params.ngenerations = 10

params.fracnew = 0.14
params.fracevolved = 0.85
params.fracelite = 0.01

params.distances = "jc.distance.precalc.csv"
params.metadata = "metadata.csv"

init = Channel.from(0)
incr = Channel.create()

batches = Channel.from(1..params.nbatches)

process Init {
	executor "local"

	input:
	val zero from init

	output:
	val zero into start

	"""
	#!/bin/bash
	if [ ! -d "${workflow.launchDir}/data/${params.data_set}" ];
	then
	  echo "Error: this program expects a subdirectory called data/${params.data_set} in the current directory."
	  exit 1
	fi
	mkdir -p ${workflow.launchDir}/output/${params.data_set}
	"""
}

Channel 
	.from incr
	.map { it + 1 }
	.until { it >= params.ngenerations }
	.set { igens }
	
generations = start.concat(igens)

process Generation {
	executor "local"

	input: 
	val gen from generations

	output:
	val gen into combo

	"""
	#!/bin/bash
	echo Starting generation: $gen
	"""
}

Channel
	.from combo
	.combine(batches)
	.set { genround }

process DoIt {
	memory "25G"

	input:
	tuple gen, batch from genround

	output:
	tuple gen, batch into genend

	"""
	#!/bin/bash
	echo "Generation $gen, Batch $batch"
	if [ "$gen" == "0" ];
	then
	  fnew=1; fevo=0; feli=0
        else
          fnew=${params.fracnew}; fevo=${params.fracevolved}; feli=${params.fracelite}
        fi
	Rscript ${workflow.projectDir}/bin/make.gen.R \
          --frac.new \$fnew \
	  --frac.evo \$fevo \
	  --frac.eli \$feli \
	  --nsamples ${params.nsamples} \
	  --gen.size ${params.gensize} \
	  --tot.batches ${params.nbatches} \
	  --n_cores ${params.ncores} \
	  --data_set ${params.data_set} \
          --w_div 1 \
          --w_tem 1 \
          --generation $gen \
          --batch $batch \
	  --basedir ${workflow.launchDir} \
	  --distance ${workflow.launchDir}/${params.distances} \
	  --metadata ${workflow.launchDir}/${params.metadata}
	"""
}

Channel
	.from genend
	.groupTuple(size: params.nbatches)
	.set { gensummary }

process Collect {
	executor "local"

	input:
	tuple gen, batches from gensummary

	output:
	val gen into incr

	"""
	#!/bin/bash
	echo "Collecting generation $gen, Batches $batches"
	"""
}


