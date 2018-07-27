bins=" mar-gcc-default-no-flto mar-gcc-default-wi-flto mar-gcc-marchhaswell-no-flto mar-gcc-marchhaswell-wi-flto mar-icc-15-0-4 mar-icc-16-0-3 mar-icc-17-0-2 mar-icc-18-0-0"
for rep in {0..4} ; do
	echo $rep
	for bin in $bins ; do
		echo $bin
		outfile=run-$rep/$bin-out.txt
		echo $outfile
		./$bin comptest.xml --steps=11 --final-checkpoint=0 >$outfile
	done
done

