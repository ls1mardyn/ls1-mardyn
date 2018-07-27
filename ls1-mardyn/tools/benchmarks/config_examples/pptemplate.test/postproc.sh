#!/bin/sh
# postproc
# postprocessing template
# M. Bernreuther <bernreuther@hlrs.de>

echo -e "GENCMDOUTPUT:\t$GENCMDOUTPUT"
echo -e "CREATEDJOBS:\t$CREATEDJOBS"
echo
echo -e "DSTROOTNAME:\t$DSTROOTNAME"
echo -e "DSTROOTPATH:\t$DSTROOTPATH"
echo -e "DELIMITER:\t$DELIMITER"
echo -e "GENTEMPLATENAME:\t$GENTEMPLATENAME"
echo -e "GENTEMPLATEPATH:\t$GENTEMPLATEPATH"
echo -e "PPTEMPLATENAME:\t$PPTEMPLATENAME"
echo -e "PPTEMPLATEPATH:\t$PPTEMPLATEPATH"
echo -e "LOGFILENAME:\t$LOGFILENAME"
echo -e "LOGFILEPATH:\t$LOGFILEPATH"
echo

cd $DSTROOTPATH
pwd
for j in `printf "$CREATEDJOBS" | tr '$DELIMITER' ' '`; do
	#echo $j
	cd $j
	pwd
	cd ..
done
