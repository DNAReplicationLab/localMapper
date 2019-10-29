read -p "This script creates replication profile for S.cerevisiae.
It needs ~20 Gb of free hard drive space,
good internet connection and a few hours to run.
Are you sure you want to run it? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then

read -p "This script creates replication profile for S.cerevisiae.
It needs ~20 Gb of free hard drive space,
good internet connection and a few hours to run.
Are you sure you want to run it?  (y/n)?" choice
case "$choice" in 
  y|Y ) echo "yes";;
  n|N ) echo "no";;
  * ) echo "invalid";;
esac



	## Get the latest sra-toolkit
	echo "Downloading the latest SRA Toolkit..."
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
	tar xvzf sratoolkit.current-ubuntu64.tar.gz
	rm sratoolkit.current-ubuntu64.tar.gz
	cd sratoolkit.*/bin/
	## Fetch G2 sample from SRA and convert it to fastq files (slow)
	echo "Fetching G2 phase sample from SRA. This may take a while..."
	./fastq-dump -B --gzip --split-files SRR926345
	mv SRR926345_*.fastq.gz ../..
	cd ../..
	## Get sacCer3 genome
	echo "Downloading and indexing sacCer3 genome"
	mkdir sacCer3
	cd sacCer3
	wget ftp://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
	tar xvzf chromFa.tar.gz
	cat chrI.fa chrII.fa chrIII.fa chrIV.fa chrV.fa chrVI.fa chrVII.fa chrVIII.fa chrIX.fa chrX.fa chrXI.fa chrXII.fa chrXIII.fa chrXIV.fa chrXV.fa chrXVI.fa chrM.fa > sacCer3.fasta
	rm chromFa.tar.gz *.fa
	bowtie2-build -f sacCer3.fasta sacCer3
	echo "Processing G2 sample"
	/home/dzmitry/Dropbox/myGitHub/localMapper.sh -g sacCer3/sacCer3 -1 SRR926345_1.fastq.gz -2 SRR926345_2.fastq.gz -s T9475_G2



	/home/dzmitry/Dropbox/myGitHub/localMapper 

	cd sratoolkit.*/bin/
	echo "Fetching S phase sample from SRA. This will take a wee bit..."
	./prefetch -o SRR926346.sra SRR926346
	./fastq-dump -B --gzip --split-files SRR926346.sra


fi
