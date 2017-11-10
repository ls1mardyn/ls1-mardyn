#!/bin/bash
numNodes=$1
# Step 1: generate base .xml file from weakscaling start, which fills whole memory
echo "script $0 generating input file for $1 nodes."

input_file_name1node="ljfluid1node.xml"
onenodesize=`cat $input_file_name1node | grep "<lx>" | sed "s/>/</g" | cut -d '<' -f 3`
size=`echo "" | awk "END {print $onenodesize * $numNodes ^ (1/3) }"`

#change size of new input:
configname="$numNodes-nodes.xml"
sed -e "s/$onenodesize/$size/g" $input_file_name1node > $configname
