#!/bin/bash
numNodes=$1
echo "generating input for $1 number of nodes"

# Step 1: generate base .xml file from weakscaling start, which fills whole memory
bash baseXML.sh $1
# Step 2: generate smaller .xml files for reducing memory utilisation
bash smallerXMLs.sh $1
# Step 3: generate jobscripts hybrid MPI x OpenMP variation of number of MPI ranks
bash varyMPI.sh $1
