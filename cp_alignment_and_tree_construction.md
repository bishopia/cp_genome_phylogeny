### goal: convert set of complete genomes into concatenated gene multifasta for mafft alignment and raxml tree building

The analytical pipeline i'm trying to follow here comes from Yu's 2017 dissertation where many diatom chloroplast genomes are aligned and analyzed. The following paragraph is most relevant:

"Phylogenetic analysis: Sequences of 20 plastid genes (psaA, psbC, petD, petG, atpA, atpG, rbcL, rbcS,
rpoA, rpoB, rps14, rpl33,rnl, rns, ycf89, sufB, sufC, dnaK, dnaB, clpC) from 22 diatom
taxa were aligned with MAFFT (Katoh et al., 2005). This included 15 published diatom
plastid genomes and the seven genomes sequenced in this study. All sequences were
included, and protein-coding genes were partitioned by gene and codon position. A
maximum likelihood tree was constructed with RAxML7.2.8 (Stamatakis, 2006a), using
the substitution model GTR+G+I and “-f a” option, and 1000 bootstrap replicates were
performed to evaluate support for clades." - From Yu (2017) dissertation: Tempo and mode of diatom plastid genome evolutionM Yu - 2017 - repositories.lib.utexas.edu

Here's my script:

```
#!/bin/bash

#create list of genes to extract for each file
for FILE in `ls S*complete*fasta`
do
	#create variable for sample name, e.g. S31, S32, etc.
	y=`echo "$FILE" | cut -c1-3`

	#remove unwanted duplicate or fragmented genes, save to list
	grep -v 'fragment' $FILE | grep -v "'_" - | grep -v "'-I" - | sort - | grep -v "II" - | grep -f gene.list - > ${y}_genes

	#report gene count per new file
	echo $y
	wc -l ${y}_genes
done

#create multigene fasta for each file, then create another that concatenates and renames 
for FILE in `ls S*complete*fasta`
do
	#create variable for sample name, e.g. S31, S32, etc.
	y=`echo "$FILE" | cut -c1-3`

	#create variable for sample specific gene list
	GENE_LIST=${y}_genes

	#remove multigene.fasta if already present in directory; if you don't you contenate to get get copies
	rm ${y}_multigene.fasta

	#loop to construct new multigene fastas
	for GENE in `cat ${GENE_LIST}`
	do
		grep -A1 "${GENE}" $FILE >> ${y}_multigene.fasta
		#echo $GENE
	done

	#transform multigene fastas to concatenate and rename
	grep -v '>' ${y}_multigene.fasta | tr -d '\n\r' | sed "1 i\>${y}_yu_geneset" - > ${y}_yu_geneset.fa
	echo "" >> ${y}_yu_geneset.fa
done
```

From here, concatenate files for each sample and run through mafft:
```
#!/bin/bash

#SBATCH --time=XX
#SBATCH --mem=XXG #default is 1 core with 2.8GB of memory
#SBATCH -n 12
##SBATCH -p bigmem
#SBATCH --account=epscor-condo
#SBATCH -J mafft_alignment

#activate alignment env, includes muscle, mafft, etc.
source activate /gpfs/home/ibishop/data/ibishop/condas/mafft

MULTIFASTA=all_yu.fa

mafft --localpair --maxiterate 1000 --thread -1 $MULTIFASTA > output_with_stel_outgroup

#and run RAxML; note that the first time I ran the above alignment it told me to toss identical sequence, and gave mem a reduced alignment, which i used for tree building instead.
~$CONDA_DIR/bin/raxmlHPC -m GTRGAMMAI -f a -p 2 -x 2 --bootstop-perms=1000 -n output_raxml -N 100 -s output_with_stel_outgroup.reduced
```

From here, you can copy over to laptop and import into Geneious for visualization



