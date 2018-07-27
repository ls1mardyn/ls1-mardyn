#!/bin/bash
numNodes=$1
# Step 2: generate smaller .xml files for reducing memory utilisation
echo "script $0 generating input file for $1 nodes."

input_file_base="$numNodes-nodes.xml"
echo "using $input_file_base as base"

# generate up to 4096 times smaller files
for ((iSize=1; iSize <= 32768; iSize *= 2)) ; 
do
	onenodesize=`cat $input_file_base | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`
	size=`echo "" | awk "END {print $onenodesize / ($iSize ^ (1/3)) }"`
	
	#for loop over schemes
	schemes="slice c08"
	for iScheme in $schemes ;
	do
		configname="$numNodes-nodes-$iSize-$iScheme.xml"
		sed -e "s/slice/$iScheme/g; s/$onenodesize/$size/g" $input_file_base > $configname
	done

done
