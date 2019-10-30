for i in *.bam
do
	file=$(echo $i | cut -d . -f1)
	temp_filename=$file.filtered.sam
	sam_header=$(samtools view -H $i)
	echo $temp_filename
	echo "$sam_header$(cat $temp_filename)" > $temp_filename
	#echo $sam_header
	# echo "task goes here
	# $(cat $tepm_filename)" > $temp_filename
done

mkdir verified_sams 
mv *.sam verified_sams
cd verified_sams
mkdir bams 
for i in *.sam
do
	samtools view -Sb $i > $i.bam
mv *.bam bams
