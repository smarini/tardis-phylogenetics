#!/usr/bin/env nextflow

params.data_set = "USA"
params.n_samples = 10
params.gensize = 100
params.nbatches = 5
params.n_cores = 2
params.ngenerations = 10

params.frac_new = 0.14
params.frac_evolved = 0.85
params.frac_elite = 0.01

params.w_div = 1
params.w_tem = 1
params.dist_opt = "max"

params.seeds = "data/seeds.txt"
params.distances = "jc.distance.precalc.csv"
params.metadata = "metadata.csv"
params.outdir = "output/${params.data_set}"
params.afile = ""

init = Channel.from(0)
incr = Channel.create()
endg = Channel.from(-1)

batches = Channel.from(1..params.nbatches)

process Init {
	executor "local"

	input:
	val zero from init

	output:
	val zero into start

	"""
	#!/bin/bash
	#if [ ! -d "${workflow.launchDir}/data/${params.data_set}" ];
	#then
	#  echo "Error: this program expects a subdirectory called data/${params.data_set} in the current directory."
	#  exit 1
	#fi

	checkParams.py ${params.frac_new} ${params.frac_evolved} ${params.frac_elite} ${params.gensize} ${params.dist_opt}
	"""
}

Channel 
	.from incr
	.map { it + 1 }
	.until { it >= params.ngenerations }
	.set { igens }
	
generations = start.concat(igens).concat(endg)

process Generation {
	executor "local"

	input: 
	val gen from generations

	output:
	val gen into combo
	val gen into terminate

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

	when:
	gen >= 0

	"""
	#!/bin/bash
	echo "Generation $gen, Batch $batch"
	if [ "$gen" == "0" ];
	then
	  fnew=1; fevo=0; feli=0
        else
          fnew=${params.frac_new}; fevo=${params.frac_evolved}; feli=${params.frac_elite}
        fi
	Rscript ${workflow.projectDir}/bin/make.gen.R \
          --frac.new \$fnew \
	  --frac.evo \$fevo \
	  --frac.eli \$feli \
	  --n.samples ${params.n_samples} \
	  --gen.size ${params.gensize} \
	  --tot.batches ${params.nbatches} \
	  --n.cores ${params.n_cores} \
	  --data.set ${params.data_set} \
          --w.div ${params.w_div} \
          --w.tem ${params.w_tem} \
          --generation $gen \
          --batch $batch \
	  --dist.opt ${params.dist_opt} \
	  --basedir ${workflow.launchDir} \
          --seeds ${workflow.launchDir}/${params.seeds} \
	  --distance ${workflow.launchDir}/${params.distances} \
	  --metadata ${workflow.launchDir}/${params.metadata} \
	  --out.dir ${workflow.launchDir}/${params.outdir}
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

process Conclusion {
	executor "local"

	input:
	val gen from terminate

	when:
	gen == -1

	script:
	"""
	#!/bin/bash
	extractSeqs.py ${workflow.launchDir}/${params.outdir} ${workflow.launchDir}/${params.afile} ${params.ngenerations} ${params.nbatches} ${params.n_samples} GA.${params.data_set}
	"""
}
