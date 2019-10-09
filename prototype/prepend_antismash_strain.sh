#!/bin/bash

AS_INPUT=$1

# test for input empty
if [ -z "$AS_INPUT" ]; then
	echo "Please provide path to antismash folder";
	exit
fi

for STRAIN_FOLDER in $(ls $AS_INPUT);
do 
	echo Process $STRAIN_FOLDER
	STRAIN=$(echo $STRAIN_FOLDER | cut -d. -f1)
	for BGC in $(ls $AS_INPUT/$STRAIN_FOLDER/*.gbk); 
	do
		FILENAME=$(basename $BGC)
		if ! [[ $FILENAME =~ ^"$STRAIN"* ]]; then
			mv $AS_INPUT/$STRAIN_FOLDER/$FILENAME $AS_INPUT/$STRAIN_FOLDER/$STRAIN.$FILENAME
		fi
	done
	for BLAST in $(ls "$AS_INPUT"/$STRAIN_FOLDER/knownclusterblast/*.txt); 
	do
		FILENAME=$(basename $BLAST)
		if ! [[ $FILENAME =~ ^"$STRAIN"* ]]; then
			mv $AS_INPUT/$STRAIN_FOLDER/knownclusterblast/$FILENAME $AS_INPUT/$STRAIN_FOLDER/knownclusterblast/$STRAIN.$FILENAME
		fi
	done
done

