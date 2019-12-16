# DTPcourse2019
RNA-seq tutorial

# Start of the session:

Download the files for this tutorial from 
https://www.dropbox.com/sh/q07p0xb3qqcahny/AADrZHwb9UoHHmB6YjKNdLZBa?dl=0
(chicken+NDV) and
https://www.dropbox.com/sh/2aff7egibvev79f/AABXNTGfJ3w5P5Q0SM0DVP2ca?dl=0
(mosquitoes) in the directory where you choose to work.
You can also download only a few of the files initially, then start working while you download the others.

Then open a terminal in the same directory. 

Decompress all gzipped FASTA files (extension .fa.gz or .fasta.gz) and annotation files (extension .gtf.gz) with the command

> gunzip FILENAME.gz

where FILENAME.gz is the name of the compressed file. 
There is no need to decompress FASTQ files (extension .fq.gz or .fastq.gz).

If you are interested in visualising the content of compressed files without decompressing them, you can use the command

> zcat FILENAME.gz | less

----------------------
----------------------

# Short reads data from RNA-seq of chicken embryo fibroblast (CEF) cells:

In this part of the tutorial, we will look at the impact of viral infections on the cell state. 
The aim is to compare the patterns of gene expression in infected and uninfected chicken cells.

Three RNA-seq replicates are provided both for naive CEF cells and for CEF cells infected with Newcastle Disease Virus (NDV). 
The corresponding filenames have the form CEFnaive_"replicate number"_R"read number".fq.gz and are in FASTQ format.
These RNA-seq samples are all paired-end: "_R1" and "_R2" denote files containing the first or the second read of each pair.

Since the reads contain fragments of viral transcripts as well as chicken transcripts, both the annotation chicken_NDV_annotation.gtf and the FASTA sequence chicken_NDV_genome.fasta contain also the virus. 
All chicken genes begin with the standard Ensemble code for chicken "ENSG".
The six NDV genes are called "NP", "P", "M", "F", "HN", "L" in order of distance from the transcription start site of the virus. 

----------------------

# Initialise alignments by indexing reference genomes:

> kallisto index -i Gallus_gallus.Gallus_gallus-5.0.cdna.all.idx Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz 

(The next steps could take ~15 minutes)

> gemtools index -i chicken_NDV_genome.fasta

> gemtools t-index -a chicken_NDV_annotation.gtf -i chicken_NDV_genome.gem

----------------------

# Pseudoalignment of short reads with kallisto:

We take as an example the sample CEFNDV_1, but the same procedure can be repeated for all samples. 

> kallisto quant -i Gallus_gallus.Gallus_gallus-5.0.cdna.all.idx -o kallistoOutput_CEFNDV_1 -b 100 CEFNDV_1_R1.fq.gz CEFNDV_1_R2.fq.gz

The results can be found in the directory kallistoOutput_CEFNDV_1 chosen as in the above command. 
The file abundance.tsv contains all read counts and TPMs for the sample.

To extract just the counts for differential expression analysis, you can use commands like

> cut -f 1,4 kallistoOutput_CEFNDV_1/abundance.tsv | sed '1d' > CEFNDV_1.kallisto.counts.txt

----------------------

# Alignment of short reads with GEM:

First align the reads

> gemtools rna-pipeline -i chicken_NDV_genome.gem -a chicken_NDV_annotation.gtf -f CEFNDV_1_R1.fq.gz CEFNDV_1_R2.fq.gz -q 33 -t 2 --no-bam

then convert the resulting file to SAM

> zcat CEFNDV_1_R1.map.gz | gem-2-sam -q offset-33 -o CEFNDV_1.sam -I chicken_NDV_genome.gem -l

then convert to BAM

> samtools view -bS CEFNDV_1.sam > CEFNDV_1.sam

then sort by position along the genome

> samtools sort CEFNDV_1.bam -T tempname -O 'bam' > CEFNDV_1.sorted.bam

and the final BAM file should be indexed for later purposes:

> samtools index CEFNDV_1.sorted.bam
 
----------------------

The files *.gtf.counts.txt contain the read counts for each gene therein.

You can use Excel/Openoffice or R to compare and analyse these data.

In case you would have problems with the ordering, you can use the "sort" command in Ubuntu before opening the files in Excel/Openoffice:

> sort -k 1 FILE1.gtf.counts.txt > FILE1.gtf.counts.sorted.txt

replacing FILE1 with the correct name.

You could also join multiple files in a single table before analysing them, using e.g. the join command:

> join -o 1.1,2.1,1.2,2.2 -a 1 -a 2 FILE1.gtf.counts.sorted.txt FILE2.gtf.counts.sorted.txt > jointcount.sorted.csv

replacing again FILE1 and FILE2 with the correct names.

To join multiple files FILE1 FILE2 ... FILEN at once, use the script

> echo FILE1 FILE2 ... FILEN | awk '{for(i=1;i<=NF;i++){print $i}}' | awk '{n=n+1; while ((getline line < $1)>0){ nf[n]=split(line,v,"\t"); key=v[0+$2]; listkeys[key]=listkeys[key]+1; s=""; for(i=1;i<=nf[n];i++){if(i!=$2){s=s"\t"v[i]}}; iskeyin[n" "key]=1; res[n" "key]=s} }END{ for(key in listkeys){ if(listkeys[key]==n){ s=key; for(i=1;i<=n;i++){s=s""res[i" "key]}; print s} else { s=key; for(i=1;i<=n;i++){if(iskeyin[i" "key]+0==1){s=s""res[i" "key]} else {for(j=1;j<=(nf[i]-1);j++){s=s"\tNO_RECORD"}}}; print s} } }' | awk '{gsub("NO_RECORD","0"); print}' | sort -k 1b,1 > jointcount.sorted.csv

----------------------

# Differential Expression Analysis:

We will use DESeq2, one of the standard tools for differential expression analysis.

First, prepare a TAB-separated text file (let's call it "data_table.txt") in the following format:

                condition	 type
               
     CEFnaive_1	naive	paired-end

     CEFnaive_2	naive	paired-end

     CEFnaive_3	naive	paired-end

     CEFNDV_1	infected	paired-end

     CEFNDV_2	infected	paired-end

     CEFNDV_3	infected	paired-end

Then, open R:

> R

And run these commands

>> library(DESeq2)

>> counts<-read.table("jointcount.sorted.csv",header=T,stringsAsFactors=F)

>> cts<-round(counts)

>> keep <- rowSums(cts) >= 10 ; cts <- cts[keep,]

>> coldata <- read.csv("data_table.txt", row.names=1, sep="\t")

>>  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

>> dds$condition <- relevel(dds$condition, ref = "naive")

>> dds <- DESeq(dds)

>> results_DEA <- results(dds, contrast=c("condition","infected","naive"))

>> results_DEA_shrinklFC <- lfcShrink(dds, coef="condition_infected_vs_naive", type="apeglm")

(the latter is better in practice, since log-fold changes are shrunk when the significance is low).

The results can be then plotted/analysed within R or exported as tab-separated file as

> write.table(results_DEA_shrinklFC, file="results_DEA.tsv", ,sep="\t", quote=F)

The two more useful columns are log2FoldChange (which is the log-fold change in base 2) and padj (the p-value adjusted for multiple testing).

A good explanation and introduction to DESeq2 can be found here: http://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

----------------------

# SNP calling from viral deep sequencing:

The content and position of the aligned read against the genome can be seen graphically from

> samtools tview -p NDV_LaSotaChina:1 CEFNDV_1.sorted.bam chicken_NDV_genome.fasta

where each "." corresponds to a base that is identical to the reference. Press q or CTRL-c to exit.

Note how noisy the data are!

To call variants from these data, first extract the information base-by-base (pileup format) with samtools:

> samtools mpileup -f chicken_NDV_genome.fasta -p NDV_LaSotaChina:1-16000 CEFNDV_1.sorted.bam > ViralReads.pileup

then run Varscan2 on the pileup file to generate a tabular file with all the variants and related information

> java -jar VarScan.jar pileup2snp ViralReads.pileup --min-coverage 100 --min-reads2 50 --min-avg-qual20 --min-var-freq 0.001 --p-value 0.05 > listVariants.tabular

You can also try to change the above parameters (see manual at http://dkoboldt.github.io/varscan/using-varscan.html#v2.3_pileup2snp) according to reasonable criteria, andd have have a look at the difference in the variants in the output file.

Then examine the resulting output file containing only variants that passed the filter, by doing

> less listVariants.tabular

The final question of this tutorial is:

these sequences belong to a single viral 'population', or viral swarm/quasispecies. In this context, we expect all sequences to form a 'cloud' of genotypes differing by a few random mutations. Is it the case? Can you find reliable nucleotide variants?

And can you spot something that looks odd, like a series of insertions of +G in multiple copies at a given position in the middle of the P sequence (location ~2300)? 
That is an example of RNA editing occurring during transcription and generating three different proteins (P, V and W) from the same gene through insertions of +G and +GG and the corresponding frameshifts.

----------------------
----------------------

# Testis-specific differences in gene expression from RNA-seq data of Anopheles

----------------------

# Alignment of short reads from RNA-sequencing of mosquitoes

List of files/datasets (Anopheles gambiae and merus):

- Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa : reference genome of A.gambiae

- Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.gtf : annotation of A.gambiae genome

- Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.genes.bed : annotation of genes in BED format, generated by the command

> cat Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.gtf | sed '1,5d' | grep -v transcript_id | cut -f 1,4,5,9 | cut -f 1 -d ";" | sed 's/gene_id //' | sed 's/"//g' > Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.genes.bed

- Anopheles"species"Testis_"replicate number"_R"read number".fastq.gz : reads from transcriptome of adult male testis of Anopheles species gambiae or merus (two replicates each, paired-end)

- Anopheles"species"Carcass_"replicate number"_R"read number".fastq.gz : reads from transcriptome of adult male whole body/carcass of Anopheles species gambiae or merus (two replicates each, paired-end)

----------------------

Index reference transcriptome for GEM (takes about 15 minutes)

> gemtools index -i Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

> gemtools t-index -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.gtf -i Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.gem

Align with GEM the first replicate of A.gambiae carcass

> gemtools rna-pipeline -i Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.gem -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.gtf -f AnophelesgambiaeCarcass_1_R1.fastq.gz AnophelesgambiaeCarcass_1_R2.fastq.gz -q 33 -t 2 --no-bam

then convert the resulting file to SAM

> zcat AnophelesgambiaeCarcass_1_R1.map.gz | gem-2-sam -q offset-33 -o AnophelesgambiaeCarcass_1.sam -I Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.gem -l

and follow the steps above to generate a sorted, indexed BAM file. 

Repeat the above steps for the other mosquito datasets as well.

And finally check the mapping statistics:

> samtools flagstat AnophelesgambiaeCarcass_1.sorted.bam

> samtools flagstat AnophelesmerusCarcass_1.sorted.bam

How do they compare? How many reads would you lose by mapping to the "wrong" reference genome? And how this depends on choosing to consider all mapped reads or properly mapped pairs only?

You can use

> samtools mpileup AnophelesmerusCarcass_1.sorted.bam | less

or

> samtools tview AnophelesmerusCarcass_1.sorted.bam Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa

to explore the data by hand. You can move around using arrows on the keyboard; press q to exit.

----------------------

# Expression track visualization

If you are interested in visualizing the data as expression profiles on the genome:

Bam files can be converted to traces in bedGraph format via commands like

> bedtools genomecov -trackline -trackopts autoscale=on -trackopts graphType=bar -bg -ibam AnophelesgambiaeCarcass_1.sorted.bam > mytrack.bg

but for practical reasons the size of the file should be less than 20 Megabytes, which can be done by 

> head -n 500000 mytrack.bg > mytrack_partial.bg

then they can be visualised at https://www.vectorbase.org/Anopheles_gambiae/Location/View?r=2L%3A2686000-2746000 by clicking on Custom Tracks on the left and loading your bedGraph file mytrack_partial.bg there.

----------------------

The files *.gtf.counts.txt contain the read counts for each gene therein.

See above for how to put them together.

----------------------

# Differential expression analysis

Follow the steps above to perform the differential expression analysis with DESeq2 in R.

This way you can compare the difference in expression in male testes vs body with the absolute expression of the gene or with the log-fold change.

Which genes are more expressed? Which genes have the largest differences in expression levels in testis, compared with their average expression? Which genes are almost only expressed in testis?

Finally, you can run a list of the most interesting genes that you found through easy-to-use webservers for functional annotation and gene enrichment analysis such as http://www.pantherdb.org or http://www.flymine.org in order to understand which biological functions are overrepresented among these genes.

----------------------

SNP calling from the previous alignments using SAMtools:

> samtools mpileup -uf Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa AnophelesmerusCarcass_1.sorted.bam | bcftools call -mO v > AnophelesmerusCarcass_1.sorted.samtools.vcf

then filter for low quality polymorphisms and select either polymorphic variants

> bcftools filter -s LowQual -e '%QUAL<20 || DP<10 || AC>1 || AC<1' AnophelesmerusCarcass_1.sorted.samtools.vcf | grep -v LowQual > AnophelesmerusCarcass_1.sorted.samtools.filtered.variants.vcf

or select fixed differences with the reference

> bcftools filter -s LowQual -e '%QUAL<20 || DP<10 || AC<2' AnophelesmerusCarcass_1.sorted.samtools.vcf | grep -v LowQual > AnophelesmerusCarcass_1.sorted.samtools.filtered.diffs.vcf

To find the number of fixed differences per gene (including introns), you can use

> bedtools coverage -counts -a Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.8.genes.bed -b AnophelesmerusCarcass_1.sorted.samtools.filtered.diffs.vcf > genes.AnophelesmerusCarcass_1.sorted.samtools.filtered.diffs.vcf.counts.txt

and similarly for the number of variants.

----------------------

You could also join one of the previous count files with the file genes.AnophelesmerusCarcass_1.sorted.samtools.filtered.diffs.vcf.counts.txt as well.

If you do this, you could analyse in Excel/Openoffice or R the resulting table to answer an interesting question:

is there a correlation between the absolute expression level of a gene and its evolutionary properties?

More specifically, do highly expressed genes in testis tend to diverge more between different species, or less?
