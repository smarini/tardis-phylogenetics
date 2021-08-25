#!/usr/bin/env nextflow

params.data_set = "USA"
params.nsamples = 10
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
	memory "1G"

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

	checkParams.py ${params.fracnew} ${params.fracevolved} ${params.fracelite} ${params.gensize} ${params.distopt}
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
	memory "1G"

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
          fnew=${params.fracnew}; fevo=${params.fracevolved}; feli=${params.fracelite}
        fi
	# echo Rscript ${workflow.projectDir}/bin/make.gen.R \
        #   --frac.new \$fnew \
	#   --frac.evo \$fevo \
	#   --frac.eli \$feli \
	#   --n.samples ${params.nsamples} \
	#   --gen.size ${params.gensize} \
	#   --tot.batches ${params.nbatches} \
	#   --n.cores ${params.ncores} \
	#   --data.set ${params.data_set} \
        #   --w.div ${params.wdiv} \
        #   --w.tem ${params.wtem} \
        #   --generation $gen \
        #   --batch $batch \
	#   --dist.opt ${params.distopt} \
        #   --seeds ${workflow.launchDir}/${params.seeds} \
	#   --distance ${params.distances} \
	#   --metadata ${params.metadata} \
	#   --out.dir ${workflow.launchDir}/${params.outdir} > ${workflow.launchDir}/make.gen.${gen}.cmd
	Rscript ${workflow.projectDir}/bin/make.gen.R \
          --frac.new \$fnew \
	  --frac.evo \$fevo \
	  --frac.eli \$feli \
	  --n.samples ${params.nsamples} \
	  --gen.size ${params.gensize} \
	  --tot.batches ${params.nbatches} \
	  --n.cores ${params.ncores} \
	  --data.set ${params.data_set} \
          --w.div ${params.wdiv} \
          --w.tem ${params.wtem} \
          --generation $gen \
          --batch $batch \
	  --dist.opt ${params.distopt} \
          --seeds ${workflow.launchDir}/${params.seeds} \
	  --distance ${params.distances} \
	  --metadata ${params.metadata} \
	  --out.dir ${workflow.launchDir}/${params.outdir}
	"""
}

Channel
	.from genend
	.groupTuple(size: params.nbatches)
	.set { gensummary }

process Collect {
	executor "local"
	memory "1G"

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
	memory "1G"

	input:
	val gen from terminate

	when:
	gen == -1

	script:
	"""
	#!/bin/bash
	extractSeqs.py ${workflow.launchDir}/${params.outdir} ${params.afile} ${params.ngenerations} ${params.nbatches} ${params.nsamples} GA.${params.data_set}
	Rscript ${workflow.projectDir}/bin/print.results.R -s ${params.data_set} -d ${params.distances} --metadata ${params.metadata} \
         --outdir ${workflow.launchDir}/${params.outdir} --ngenerations ${params.ngenerations} --n.batches ${params.nbatches}
	"""
}
