#!/usr/bin/make -f
## mapping dSMF-gw sequences with bwa-meth
## required software: fastqc, cutadapt, trimmomatic, bwa-meth, samtools, picard, qualimap

###############################
########### VARIABLES #########
###############################

genomefile:= ${HOME}/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa
trimmomaticDIR := ${HOME}/Trimmomatic-0.36
trimAdapterFile := ${trimmomaticDIR}/adapters/TruSeq_2-3_PE.fa
methIndGenomeFiles := $(addsuffix ${genomefile}.bwameth.ct2, .sa .amb .ann .pac .bwt)
bname :=  $(addprefix 180126_SNK268_A_L001_JIB-, 1 2 3 4)
longbname := $(addsuffix _R1, $(bname)) $(addsuffix _R2, $(bname))
PRESEQ:=/Users/semple/mySoftware/preseq


#fileList=( $(cat fileList.txt) )

#list of the final output files
objects := 	${methIndGenomeFiles} \
	$(addsuffix .filt.bam, $(addprefix aln/, $(bname))) 

#list of the various reports and stats produced during analysis
statsObjects := $(addsuffix _fastqc.html, $(addprefix rawData/fastQC/, $(longbname))) \
	$(addsuffix _fastqc.html, $(addprefix cutadapt/fastQC/, $(longbname))) \
	$(addsuffix _forward_paired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addsuffix _forward_unpaired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addsuffix _reverse_paired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addsuffix _reverse_unpaired_fastqc.html, $(addprefix trim/fastQC/, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _stats.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _picard_insert_size_metrics.txt, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _picard_insert_size_histogram.pdf, $(bname))) \
	$(addprefix aln/prefilt/report_, $(addsuffix _qualimap.pdf, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _picard.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _stats.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _flagstats.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _picard_insert_size_metrics.txt, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _picard_insert_size_histogram.pdf, $(bname))) \
	$(addprefix aln/postfilt/report_, $(addsuffix _qualimap.pdf, $(bname))) \
	$(addprefix aln/postfilt/, $(addsuffix _depthStats.txt, $(bname))) \
	$(addprefix plots/moreSeqProj_, $(addsuffix .pdf, $(bname))) 
	
	
#	$(addprefix bed/preseq/readCounts_, $(addsuffix .txt, $(bname))) \
#	$(addprefix bed/preseq/c_curve_output_, $(addsuffix .txt, $(bname))) \
#	$(addprefix bed/preseq/lc_extrap_output_, $(addsuffix .txt, $(bname))) \	
#	$(addprefix bed/preseq/bound_pop_output_, $(addsuffix .txt, $(bname))) \

	
#list of files to delete at the end of the analysis
intermediateFiles := $(addsuffix .fastq.gz, $(addprefix cutadapt/, $(longbname))) \
	$(addsuffix _forward_paired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix _forward_unpaired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix _reverse_paired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix _reverse_unpaired.fq.gz, $(addprefix trim/, $(bname))) \
	$(addsuffix .sam, $(addprefix aln/, $(bname))) \
	$(addsuffix .dup.bam, $(addprefix aln/, $(bname))) \
	$(addprefix aln/postfilt/, $(addsuffix _depthCol.txt, $(bname))) \
	$(addprefix bed/, $(addsuffix _sort.bed, $(bname))) 
	

#list of secondary files to keep
secondaryFiles :=   $(addsuffix .sorted.bam, $(addprefix aln/, $(bname)))


###############################
########### RULES  ############
###############################

all:  $(objects) $(statsObjects) $(secondaryFiles)

#use cleanall when you want to force rerunning all the analysis
cleanall:
	rm -f $(intermediateFiles)
	rm -f $(secondaryFiles)
	
#use clean if the intermediate files are not properly removed (should not be required)
clean:
	rm -f $(intermediateFiles)

cleanall4rerun:
	rm -f $(objects)
	rm -f $(statsObjects)
	rm -f $(secondaryFiles)

.PHONY: all clean cleanall cleannall4rerun
.INTERMEDIATE: $(intermediateFiles)
.SECONDARY: $(secondaryFiles)




#######################################################
## get initial read stats                            ##
#######################################################

#run fastqc on sequences
rawData/fastQC/%_fastqc.html: rawData/%.fastq.gz
	mkdir -p rawData/fastQC
	fastqc $^ -o rawData/fastQC 


#######################################################
## trim adaptors with cutadapt                       ##
#######################################################

# use cutadapt to trim
cutadapt/%_R1.fastq.gz cutadapt/%_R2.fastq.gz: rawData/%_R1.fastq.gz rawData/%_R2.fastq.gz
	mkdir -p cutadapt
	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
		-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
		-o cutadapt/$*_R1.fastq.gz -p cutadapt/$*_R2.fastq.gz \
		rawData/$*_R1.fastq.gz rawData/$*_R2.fastq.gz

#redo fastQC on trimmed reads
cutadapt/fastQC/%_fastqc.html: cutadapt/%.fastq.gz
	mkdir -p cutadapt/fastQC
	fastqc $^ -o cutadapt/fastQC 


#######################################################
## quality trim reads with Trimmomatic               ##
#######################################################

# use trimmomatic to trim
trim/%_forward_paired.fq.gz trim/%_forward_unpaired.fq.gz trim/%_reverse_paired.fq.gz trim/%_reverse_unpaired.fq.gz: cutadapt/%_R1.fastq.gz cutadapt/%_R2.fastq.gz
	mkdir -p trim
	java -jar ${trimmomaticDIR}/trimmomatic-0.36.jar PE cutadapt/$*_R1.fastq.gz cutadapt/$*_R2.fastq.gz trim/$*_forward_paired.fq.gz trim/$*_forward_unpaired.fq.gz trim/$*_reverse_paired.fq.gz trim/$*_reverse_unpaired.fq.gz ILLUMINACLIP:${trimAdapterFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2> trim/report_$*_trimmomatic.txt

# redo fastQC on trimmed reads	
trim/fastQC/%_fastqc.html: trim/%.fq.gz
	mkdir -p trim/fastQC
	fastqc $^ -o trim/fastQC


#######################################################
## align to genome with BWA-meth and convert to bam  ##
#######################################################
	
# convert and index genome file for bwameth alignment
${methIndGenomeFiles}: ${genomefile}	
	bwameth.py index ${genomefile}

# align sequences to meth converted genome with bwameth
aln/%.sam: trim/%_forward_paired.fq.gz trim/%_reverse_paired.fq.gz ${methIndGenomeFiles}
	mkdir -p aln
	bwameth.py --threads 3 --reference ${genomefile} trim/$*_forward_paired.fq.gz trim/$*_reverse_paired.fq.gz > aln/$*.sam


# use samtools to convert to bam and sort
aln/%.sorted.bam: aln/%.sam
	samtools view -u $^ | samtools sort -o $@
	rm $^


#######################################################
## Get alignment stats                               ##
#######################################################

# 	get alignment stats
aln/prefilt/report_%_flagstats.txt: aln/%.sorted.bam
	mkdir -p aln/prefilt
	samtools flagstat $^ > $@

aln/prefilt/report_%_stats.txt: aln/%.sorted.bam
	samtools stats $^ > $@
	
#samtools view -cF 0x100 accepted_hits.bam

# Get insert size statistics and plots with picard and qualimap (set path to $QUALIMAP in .bash_profile)
aln/prefilt/report_%_picard_insert_size_metrics.txt aln/prefilt/report_%_picard_insert_size_histogram.pdf \
	aln/prefilt/report_%_qualimap.pdf: aln/%.sorted.bam ${PICARD} ${QUALIMAP}
	java -Xms1g -Xmx5g -jar ${PICARD} CollectInsertSizeMetrics I=aln/$*.sorted.bam \
	O=aln/prefilt/$*_picard_insert_size_metrics.txt \
    H=aln/prefilt/$*_picard_insert_size_histogram.pdf
	${QUALIMAP} bamqc -bam aln/$*.sorted.bam -c -outdir aln/prefilt -outfile $*_report_qualimap.pdf -outformat PDF


#######################################################
## Mark duplicates and filter reads. Then redo stats ##
#######################################################

# mark duplicates with picard (path to picard should be set in $PICARD variable in .bash_profile or in session)
aln/%.dup.bam aln/postfilt/report_%_picard.txt: aln/%.sorted.bam ${PICARD}
	mkdir -p aln/postfilt
	java -Xmx5g -jar ${PICARD} MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.dup.bam M=aln/postfilt/report_$*_picard.txt

# 	remove mitochondrial reads
aln/%.filt.bam: aln/%.dup.bam
	samtools view -q 30 -F 1804 -b $^ > $@
	rm $^ 	
# NOTE: sam flag 1804 means the following:
# read unmapped
# mate unmapped
# not primary alignment
# read fails platform/vendor quality checks
# read is PCR or optical duplicate

# 	get alignment stats again post-filtering
aln/postfilt/report_%_flagstats.txt: aln/%.filt.bam
	samtools flagstat $^ > $@

aln/postfilt/report_%_stats.txt: aln/%.filt.bam
	samtools stats $^ > $@

# Get insert size statistics and plots with picard and qualimap post-filtering
aln/postfilt/report_%_picard_insert_size_metrics.txt aln/postfilt/report_%_picard_insert_size_histogram.pdf \
	aln/postfilt/report_%_qualimap.pdf: aln/%.filt.bam ${PICARD} ${QUALIMAP}
	java -Xms1g -Xmx5g -jar ${PICARD} CollectInsertSizeMetrics I=aln/$*.filt.bam \
	O=aln/postfilt/$*_picard_insert_size_metrics.txt \
    H=aln/postfilt/$*_picard_insert_size_histogram.pdf
	${QUALIMAP} bamqc -bam aln/$*.filt.bam -c -outdir aln/postfilt -outfile $*_report_qualimap.pdf -outformat PDF



#######################################################
## get median coverage 						         ##
#######################################################

aln/postfilt/%_depthCol.txt: aln/%.filt.bam
	samtools depth -a $^ | cut -f3  > $@

aln/postfilt/%_depthStats.txt: aln/postfilt/%_depthCol.txt
	echo "min\tmax\tmedian\tmean" > $@
	./R/mmmm.r < $^ >> $@
	#depthStats=`./R/mmmm.r < $^`
	#echo ${depthStats}
	#echo "${depthStats}" >> $@



#######################################################
## predict unique seqs with greater seq depth 	     ##
####################################################### 

#convert to bed
bed/%_sort.bed: aln/%.filt.bam
	mkdir -p bed
	bedtools bamtobed -i $^ | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 > $@


#preform predictions
plots/moreSeqProj_%.pdf: bed/%_sort.bed
	mkdir -p bed/preseq
	mkdir -p plots
	${PRESEQ}/preseq c_curve -o bed/preseq/c_curve_output_$*.txt bed/$*_sort.bed
	#${PRESEQ}/preseq lc_extrap -o bed/preseq/lc_extrap_output_$*.txt bed/$*_sort.bed
	#${PRESEQ}/preseq bound_pop -o bed/preseq/bound_pop_output_$*.txt bed/$*_sort.bed
	wc -l bed/$*_sort.bed >> bed/preseq/readCounts_$*.txt
	#Rscript R/sequenceMore.R $* $PWD/bed/preseq
	touch plots/moreSeqProj_$*.pdf


