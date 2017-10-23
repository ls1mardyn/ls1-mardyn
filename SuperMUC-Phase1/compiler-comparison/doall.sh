bins=" mar-gcc-7 mar-icc-16 mar-icc-17 mar-icc-18"
for rep in {0..4} ; do
	echo $rep
	for bin in $bins ; do
		echo $bin
		outfile=run-$rep/$bin-out.txt
		echo $outfile
		./$bin comptest.xml --steps=11 --final-checkpoint=0 >$outfile
	done
done

