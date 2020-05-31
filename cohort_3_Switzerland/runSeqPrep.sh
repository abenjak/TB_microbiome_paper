#!/bin/bash
###############################################################################
#
#	Filename:	runSeqPrep
#
#	Description:	runs SeqPrep on given fastq.gz files.
#	
#	Usage:		runSeqPrep 'path to R1.fastq.gz' 'path to R2.fastq.gz' 'min overlapping length' 'min final read length' 'output base name'
#
#	Note:		Will run with these parameters: -A AGATCGGAAGAGCACACGTCT -B AGATCGGAAGAGCGTCGTGTA -y I (=max Phred33 score 41, default is 61, which can cause downstream problems in some programs)
###############################################################################

usage="
Usage: runSeqPrep 'path to R1.fastq.gz' 'path to R2.fastq.gz' 'min overlapping length' 'min final read length' 'output base name (to which '_SeqPrep-merged.fastq.gz', '_SeqPrep_adapt-trimmed.fastq.gz' and '_SeqPrep-NOTmerged.fastq.gz'  will be appended)'

Note: Will run with these parameters: -A AGATCGGAAGAGCACACGTCT -B AGATCGGAAGAGCGTCGTGTA -y I (=max Phred33 score 41 instead of the default 61)
You can provide additional SeqPrep parameters as suffix to the 'min final read length' argument, just make sure all is enclosed with quotes"

while getopts ':h:' option; do
	case "$option" in
  	  h) echo "$usage"
    	   exit
    	   ;;
esac
done

#if no arguments are given show usage:
if [ $# -eq 0 ]
then
 echo "$usage"
 exit
fi

# show usage if "-h" argument given:
if [ $# = "-h" ]
 then 
  echo "$usage"
  exit
fi

#rename input parameters for clarity and define variables

inR1=$1
inR2=$2
OVERLAP=$3
MINLEN=$4
BASENAME=$5

#in case some variables are missing report it and exit:
if [ -z "$5" ]
then
	echo "
Some arguments are missing. It should be: <R1> <R2> <min overlap> <min read length> <basename>

Exiting now...
"
	exit
fi

fnameR1=$(basename $inR1 .fastq.gz)
fnameR2=$(basename $inR2 .fastq.gz)

# first write the command to be executed to the log file:
echo 'SeqPrep -S -f '$inR1' -r '$inR2' -1 '$BASENAME'_SeqPrep_adapt-trimmed.R1.fastq.gz -2 '$BASENAME'_SeqPrep_adapt-trimmed_R2.fastq.gz -3 '$BASENAME'_SeqPrep-NOTmerged_R1.fastq.gz -4 '$BASENAME'_SeqPrep-NOTmerged_R2.fastq.gz -A AGATCGGAAGAGCACACGTCT -B AGATCGGAAGAGCGTCGTGTA -y I -L '$MINLEN' -o '$OVERLAP' -s '$BASENAME'_SeqPrep-merged.fastq.gz' -E $BASENAME'_SeqPrep-merged_alignments.txt.gz' > $BASENAME'_SeqPrep.log'

#now run SeqPrep and save the log (use -S to enable the spinner progress, else the Pairs Processed is always 0. I filed the bug at Github).
echo "
Running SeqPrep on $BASENAME...
"
SeqPrep -S -f $inR1 -r $inR2 -1 $BASENAME'_SeqPrep_adapt-trimmed.R1.fastq.gz' -2 $BASENAME'_SeqPrep_adapt-trimmed_R2.fastq.gz' -3 $BASENAME'_SeqPrep-NOTmerged_R1.fastq.gz' -4 $BASENAME'_SeqPrep-NOTmerged_R2.fastq.gz' -A AGATCGGAAGAGCACACGTCT -B AGATCGGAAGAGCGTCGTGTA -y I -L $MINLEN -o $OVERLAP -s $BASENAME'_SeqPrep-merged.fastq.gz' -E $BASENAME'_SeqPrep-merged_alignments.txt.gz' 2>> $BASENAME'_SeqPrep.log'

#delete the unnecessary progress line in the log file:
sed -i '/Processing reads/d' $BASENAME'_SeqPrep.log'

#parse the log file and calculate the merging rate:
MERGED=$(grep 'Pairs Merged' $BASENAME'_SeqPrep.log' | sed -r 's/Pairs Merged..([0-9])/\1/')
TOTAL=$(grep 'Pairs Processed' $BASENAME'_SeqPrep.log' | sed -r 's/Pairs Processed..([0-9])/\1/')
perl -e 'print (($ARGV[0] / $ARGV[1])*100);' -- $MERGED $TOTAL | sed -r 's/(^.+\..{1,2}).+/merging is \1%/' >> $BASENAME'_SeqPrep.log'
#print the merging ratio to screen:
tail -n1 $BASENAME'_SeqPrep.log'

exit 0
