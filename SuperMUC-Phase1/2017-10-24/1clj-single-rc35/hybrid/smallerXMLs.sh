#!/bin/bash
numNodes=$1
# Step 2: generate smaller .xml files for reducing memory utilisation
echo "script $0 generating input file for $1 nodes."

input_file_base="$numNodes-nodes.xml"
echo "using $input_file_base as base"

# generate up to 4096 times smaller files
for ((i=1; i <= 4096; i *= 2)) ; 
do
	onenodesize=`cat $input_file_base | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`
	size=`echo "" | awk "END {print $onenodesize / ($i ^ (1/3)) }"`
	
	#change size of new input:
	configname="$numNodes-nodes-$i.xml"
	sed -e "s/$onenodesize/$size/g" $input_file_base > $configname

done
