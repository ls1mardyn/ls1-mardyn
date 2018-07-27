#!/bin/bash
# pp
# postprocessing template for Hermit (http://www.hlrs.de/systems/platforms/cray-xe6-hermit/)
# Martin Bernreuther <bernreuther@hlrs.de>
#
#PBS -l mppwidth=1
#PBS -l walltime=0:05:00
#PBS -N pp
##PBS -k eo
##PBS -M <emailadress>
##PBS -m abe

if [ $# -ge 1 ]; then
	# prepare to submit job
	jobids=`echo "$GENCMDOUTPUT" | tr '\f' ':'`
	qsub -W depend=afterany:${jobids} pp
else
	# submitted job
	#cd ${PBS_O_WORKDIR}
	cd $DSTROOTPATH
	pwd
	# get runtime for all jobs
	echo -e '# job\thost\treturncode'
	echo -e '# job\tnumproc\tjob runtime\truntime' >runtimes.dat
	for j in `printf "$CREATEDJOBS" | tr '\f' ' '`; do
		output=`ls $j/*.o*`
		if [ -r "${output}" ]; then
			errors=$((`grep -ci 'ERROR' ${output}`))
			jobruntime=$((`grep 'running for' ${output} | sed -e 's/^\s*running for\s*//' -e 's/\s.*$//'`))
			host=$((`grep 'Execution host:' ${output} | sed -e 's/.*Execution host:\s*//' -e 's/\s*$//'`))
			numproc=$((`grep 'Running with .* processes' ${output} | sed -e 's/.*Running with\s*//' -e 's/\sprocesses.*//'`))
			computationtime="$(grep 'Computation in main loop took:' ${output} | sed -e 's/.*Computation in main loop took:\s*//' -e 's/\s.*$//')"
			echo -e "$j\t${host}\t${errors}"
			echo -e "$j\t${numproc}\t${jobruntime}\t${computationtime}" >>runtimes.dat
		fi
	done
	gnuplot runtimes.gplt
fi

