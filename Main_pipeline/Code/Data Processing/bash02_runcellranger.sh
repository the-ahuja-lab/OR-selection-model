#!/usr/bin/bash

# transcriptome="/storage/debarka/princey18246/debarkaWS/refdata-cellranger-mm10-3.0.0"
# fastqs="/storage/debarka/princey18246/debarkaWS/skin/fq"
#transcriptome="/home/sidrah19220/princy/crispr/reference_hg19/hg19"
#fastqs="/home/sidrah19220/princy/crispr/fq"

transcriptome = "/home/sidrah19220/mouse_a/ref/mm10"
fastqs = "/home/sidrah19220/nmd/mouse/md1"
while read line
do
	sample="$(cut -d'_' -f1 <<<$line)"
	cellranger count --id=$sample --transcriptome /home/sidrah19220/mouse_a/ref/mm10 --fastqs /home/sidrah19220/nmd/mouse/md1 --sample=$sample --expect-cells=8000 --localcores=12
done<"SRR.txt"
