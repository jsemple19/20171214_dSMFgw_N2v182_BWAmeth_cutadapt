#! /bin/bash
bnames=( N2_DE_gwV006_forward_paired_1cfa438c62d9 N2_DE_gwV007_forward_paired_1cfa5510d503 F2_DE_gwV008_forward_paired_1cfa54f2b6cd F2_DE_gwV009_forward_paired_1cfa73e0160a	)
PRESEQ=/Users/semple/mySoftware/preseq



for bname in ${bnames[@]}
do
	echo $bname
	myFile=${bname}.bam
	depthFile=${bname}_depthStats.txt
	
	samtools depth -a $myFile > $depthFile
		
	depthCol=${bname}_depthCol.txt
	cut -f3 ${depthFile} > ${depthCol}
	
	depthStats=`./mmmm.r < ${depthCol}`
	
	if [ ! -e "depthStats.txt" ]
	then
		echo "sample min max median mean" > depthStats.txt
	fi
	
	echo $bname $depthStats >> depthStats.txt

	# now project the duplication rate for higher sequencing depth
	bedFile=${bname}.bed
	bedtools bamtobed -i $myFile > $bedFile

	bedSorted=${bname}_sort.bed
	sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 $bedFile > $bedSorted

	${PRESEQ}/preseq c_curve -o c_curve_output_${bname}.txt $bedSorted

	${PRESEQ}/preseq lc_extrap -o lc_extrap_output_${bname}.txt $bedSorted

	${PRESEQ}/preseq bound_pop -o bound_pop_output_${bname}.txt $bedSorted

	wc -l $bedSorted >> readCounts_${bname}.txt
	
	Rscript ../../sequenceMore.R ${bname} $PWD

	rm $bedFile $bedSorted c_curve_output_${bname}.txt lc_extrap_output_${bname}.txt \
	bound_pop_output_${bname}.txt readCounts_${bname}.txt $depthFile $depthCol
done
