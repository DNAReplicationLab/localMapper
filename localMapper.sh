#!/bin/bash -l

#~ VERSION
scriptVersion="1.0"

#~ DESCRIPTION
#~ shell script to map Illuminia paired-end or single reads to reference, sort resulting bam file, mark duplicates 
#~ (for paired-end), index the bam file, extract number of reads mapped in bins (5', first in pair for paired end).
#~ The script creates a bash file and executes it (unless used with the -t | --test argument)

#~ AUTHORS
#~ Dzmitry Batrakou, Conrad Nieduszynski

#~ REQUIREMENTS
#~ You need to have all the necessary programs installed and available for the script. There are several options:
#~	1) - easy: all programs are installed into your search path (you still may want to define $BOWTIE2_INDEXES);
#~			if running Ubuntu (or other Debian-based Linux, try running the following command:
#~				$ sudo apt-get install bowtie2 samtools bedtools picard-tools
#~			if successful, you don't need to worry about the variables below.
#~	2) - painful: the programs are downloaded manually and cannot be run by typing $ bowtie2 (or similar) into the command line. You must use one of the options below:
#~		A) installation folders are added manually to your search path (again, mind $BOWTIE2_INDEXES);
#~			see https://askubuntu.com/questions/60218/how-to-add-a-directory-to-the-path
#~		B) the variables are specified in this script (see below).
#~		C) variables $BOWTIE2, $BOWTIE2_INDEXES, $SAMTOOLS, $BEDTOOLS, $PICARDTOOLS are specified elsewhere (e.g your ~/.profile file)
#~			see https://help.ubuntu.com/community/EnvironmentVariables

#~ VARIABLES
#~  If the required programs are not found in your search path, use lines below to specify and uncomment as needed
#~ (make sure that the files are executable):
# BOWTIE2="/path/to/bowtie2"
# BOWTIE2_INDEXES="/path/to/bowtie2/indexFolder"  #~ either specify this variable or use absolute path with the -g parameter.
# SAMTOOLS="/path/to/samtools"
# BEDTOOLS="/path/to/bedtools"
# PICARDTOOLS="/path/to/picard.jar" #~ only required for paired-end reads

#~ NOTE
#~ surround multiple read files with quotes and separate by commas without spaces, like so:
#~			"firstFile.fastq.gz,secondFile.fastq.gz"

# collect command line options
while [ "$1" != "" ]; do
	case $1 in
		-1 | --first )		shift				# file(s) with first mate pair sequences in fastq format
							firstMATEfiles=$1
									;;
		-2 | --second)		shift				# file(s) with second mate pair sequences in fastq format
							secondMATEfiles=$1	
									;;
		-c | --cores )		shift				# maximum number of cores to use [1]
							theNUMBERofCORES=$1
									;;
		-g | --genome )		shift				# name of reference genome (must be present in bowtie2 index directory
							genomePath=$1
							genomeName=${1##*/}
									;;
		-h | --help )
						echo "$(basename $0 .sh) $scriptVersion"
						echo "Usage: $(basename $0) [OPTIONS]..."
						echo "Map and count unique reads using bowtie2, samtools and bedtools."
						echo ""
						echo "PARAMETERS"
						echo "	-1, --first	first mate pair sequences in fastq format"
						echo "	-2, --second	second mate pair sequences in fastq format"
						echo "	-c, --cores	number of cores to use, defaults to 1"
						echo "	-g, --genome	path to the indexed genome"
						echo "	-h, --help	display this help and exit"
						echo "	-m, --memory	amount of RAM to use per core"
						echo "	-s, --sample	name base for the resulting files"
						echo "	-t, --test	dry run (only creates the bash file)"
						echo "	-U, --reads	single end sequence reads in fastq format"
						echo "	-w, --window	the window size in basepairs, defaults to 1000"
						echo ""
						echo "EXAMPLE USAGE"
						echo "	$(basename $0) -g sacCer3 -U 'singleReadFile1.fastq.gz,singleReadFile2.fastq.gz' -s sampleName -w 1000 -c 8 -m 500M"
						echo "		(single end reads, multiple files, BOWTIE2_INDEXES specified)"
						echo ""
						echo "	$(basename $0) -g /path/to/bowtie2_index/sacCer3 -1 firstMateFile.fastq.gz -2 secondMateFile.fastq.gz -s sampleName -w 1000 -c 8 -m 500M"
						echo "		(paired-end reads, BOWTIE2_INDEXES not specified):"
						exit 1
									;;

		-m | --memory)		shift				# amount of memory to use per CPU core [100M]
							memory=$1
									;;

		-s | --sample)		shift				# basename for resulting files
							sampleName=$1
									;;
		-t | --test )		testMode=TRUE			# testing mode: don't execute the bash script (file modifying commands are ommitted)
									;;
		-U | --reads)		shift				# file(s) with single end sequence reads in fastq format
							readFILES=$1
									;;
		-w | --window )		shift				# the window size in basepairs [1000]
							theWINDOWsize=$1
									;;
        esac
        shift
done

if ! hash bowtie2 2>/dev/null; then
	if [ -z "$BOWTIE2" ] ; then
		echo "Could not find bowtie2. You must define BOWTIE2 variable (in this script or your ~/.profile file)"
		exit 3
	fi
else
	BOWTIE2="bowtie2"
fi
if ! hash samtools 2>/dev/null; then
	if [ -z "$SAMTOOLS" ] ; then
		echo "Could not find samtools. You must define SAMTOOLS variable (in this script or your ~/.profile file)"
		exit 3
	fi
else
	SAMTOOLS="samtools"
fi
if ! hash bedtools 2>/dev/null; then
	if [ -z "$BEDTOOLS" ] ; then
		echo "Could not find bedtools. You must define BEDTOOLS variable (in this script or your ~/.profile file)"
		exit 3
	fi
else
	BEDTOOLS="bedtools"
fi
if [ -z "$secondMATEfiles" ]; then
	if [ -z "$readFILES" ]; then
		echo "Missing read file(s)"
		exit 4
	else
		echo "Single end sequence file(s) detected: $readFILES"
	fi
else
	if ! hash PicardCommandLine 2>/dev/null; then
		if [ -z "$PICARDTOOLS" ] ; then
			echo "Could not find picard-tools. Define path to the jar file in PICARDTOOLS variable (in this script or your ~/.profile file)"
			exit 3
		fi
	fi
	matePAIR=yes
	echo "Mate pair detected: $firstMATEfiles and $secondMATEfiles"
fi

# set default options
WORKING_DIRECTORY=`pwd`
if [ ! -w $WORKING_DIRECTORY ] ; then WORKING_DIRECTORY="/home/"${USER} ; echo "Not a writable directory. The files will be written to /home/${USER} instead." ; fi
theNUMBERofCORESdefault=1
theWINDOWsize=1000
memoryDefault=100M

# create subdirectories
mkdir $WORKING_DIRECTORY/${sampleName}
mkdir $WORKING_DIRECTORY/${sampleName}/processed
mkdir $WORKING_DIRECTORY/${sampleName}/raw
mkdir $WORKING_DIRECTORY/${sampleName}/tmp
# set up variables
sampleDir=$WORKING_DIRECTORY/${sampleName}
processed=$WORKING_DIRECTORY/${sampleName}/processed
raw=$WORKING_DIRECTORY/${sampleName}/raw
tmp=$WORKING_DIRECTORY/${sampleName}/tmp
theBASHfile=${sampleDir}/${sampleName}_$(basename $0 .sh).bash
echo "	The number of cores to use is ${theNUMBERofCORES:=$theNUMBERofCORESdefault}"
echo "	The amount of memory to use per core is ${memory:=$memoryDefault}"
if [ "$testMode" = TRUE ] ; then
	echo "	Running in testing mode"
fi

# define the comands
if [ -n "$matePAIR" ] ; then
	command1="$BOWTIE2 -p $theNUMBERofCORES -X 750 -q --phred33 -x $genomePath -1 \"$firstMATEfiles\" -2 \"$secondMATEfiles\" -S ${tmp}/${sampleName}.sam"
	command2="mv {$firstMATEfiles,$secondMATEfiles} ${raw} && chmod 0444 ${raw}/*"
	command3="$SAMTOOLS sort -l 9 -m $memory -O bam -o ${tmp}/${sampleName}.bam -T ${tmp}/sort.tmp -@ $theNUMBERofCORES ${tmp}/${sampleName}.sam"
	command4="rm ${tmp}/${sampleName}.sam"
	if ! hash PicardCommandLine 2>/dev/null; then
		command5="java -XX:ParallelGCThreads=$((theNUMBERofCORES-1)) -jar $PICARDTOOLS MarkDuplicates I=${tmp}/${sampleName}.bam O=${processed}/${sampleName}.bam M=${sampleDir}/${sampleName}_picard_metrics.txt TMP_DIR=${tmp}"
	else
		command5="PicardCommandLine MarkDuplicates I=${tmp}/${sampleName}.bam O=${processed}/${sampleName}.bam M=${sampleDir}/${sampleName}_picard_metrics.txt TMP_DIR=${tmp}"
	fi
else
	command1="$BOWTIE2 -a -p $theNUMBERofCORES -q --phred33 -x $genomePath -U \"$readFILES\" -S ${tmp}/${sampleName}.sam"
	command2="mv {$readFILES} ${raw} && chmod 0444 ${raw}/*"
	command3="$SAMTOOLS sort -l 9 -m $memory -O bam -o ${processed}/${sampleName}.bam -T ${tmp}/sort.tmp -@ $theNUMBERofCORES ${tmp}/${sampleName}.sam"
	command4="rm ${tmp}/${sampleName}.sam"
	command5="sleep 1"
fi
command6="$SAMTOOLS index ${processed}/${sampleName}.bam"

# check if the genome divided into windows already exists, and if not - make one:
command7="$SAMTOOLS idxstats ${processed}/${sampleName}.bam | \
awk 'BEGIN {OFS=\"\\t\"} {if (\$2>0) print (\$1,\$2)}' > ${tmp}/${genomeName}.txt  && \
$BEDTOOLS makewindows -g ${tmp}/${genomeName}.txt -w $theWINDOWsize > ${tmp}/${genomeName}.${theWINDOWsize}bp.bed"

# determine the number of 5’ ends of reads from the sample that map to each position in the genome (you should check the -F and -f flags in the samtools comment and the grep -v XS:i: which we use to remove reads mapping to multiple genomic locations (how these are recorded depends on the mapping software you used)
if [ -n "$matePAIR" ]; then
	command8="$SAMTOOLS view -h -@ $theNUMBERofCORES -q 30 -F 3840 -f 64 -L ${tmp}/${genomeName}.${theWINDOWsize}bp.bed ${processed}/${sampleName}.bam | grep -v XS:i: | $SAMTOOLS view -@ $theNUMBERofCORES -b - | $BEDTOOLS genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS=\"\\t\"} {if (\$3>0) print \$1,\$2,\$2,\"name\",\$3}' > ${tmp}/${sampleName}.bed"
	else
	command8="$SAMTOOLS view -h -@ $theNUMBERofCORES -q 30 -F 3840 -L ${tmp}/${genomeName}.${theWINDOWsize}bp.bed ${processed}/${sampleName}.bam | grep -v XS:i: | $SAMTOOLS view -@ $theNUMBERofCORES -b - | $BEDTOOLS genomecov -5 -d -ibam stdin | awk 'BEGIN {OFS=\"\\t\"} {if (\$3>0) print \$1,\$2,\$2,\"name\",\$3}' > ${tmp}/${sampleName}.bed"
fi
# this sums reads for the control sample within the specified windows (e.g. 1000 bp) and uses awk to convert to bed format
command9="$BEDTOOLS map -a ${tmp}/${genomeName}.${theWINDOWsize}bp.bed -b ${tmp}/${sampleName}.bed -null 0 -o sum | awk 'BEGIN {OFS=\"\\t\"} {if (\$4>0) print \$1,\$2,\$3,\"name\",\$4}' > ${processed}/${sampleName}.bed"
# clean up
command10="rm -r ${tmp}"

# write the bash file
echo "#!/bin/bash" > $theBASHfile
echo "#$(basename $0 .sh) version: $scriptVersion" >> $theBASHfile
echo $command1 >> $theBASHfile
if [ "$testMode" != TRUE ] ; then
	echo $command2 >> $theBASHfile
fi
echo $command3 >> $theBASHfile
if [ "$testMode" != TRUE ] ; then
	echo $command4 >> $theBASHfile
fi
echo $command5 >> $theBASHfile
echo $command6 >> $theBASHfile
echo $command7 >> $theBASHfile
echo $command8 >> $theBASHfile
echo $command9 >> $theBASHfile
if [ "$testMode" != TRUE ] ; then
	echo $command10 >> $theBASHfile
fi
# execute the bash file
if [ "$testMode" != TRUE ] ; then
	/bin/bash $theBASHfile 2> ${sampleDir}/${sampleName}_$(basename $0 .sh)_errors.txt
fi
