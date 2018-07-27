for d in ./*/*.cfg ; do
	sed -i.bak '/output/d' $d
done

