#! usr/bin/bash
input = "/home/sanjay/sidrah/md1"
while read -r line
do
	sampleName="$(cut -d'_' -f1 <<<$line)"
	mv ${sampleName}_1.fastq.gz "${sampleName}_S${i}_L001_R1_001.fastq.gz"
	mv ${sampleName}_2.fastq.gz "${sampleName}_S${i}_L001_R2_001.fastq.gz"
	# echo "${sampleName}_S${i}_L001_R1_001.fastq.gz" "${sampleName}_S${i}_L001_R2_001.fastq.gz"
	((i+=1))
done<"SRR.txt"
