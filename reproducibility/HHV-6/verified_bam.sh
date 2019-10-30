suffix=$2
echo $suffix
iter=0
while read line 
do
	read_id=$(echo $line | cut -d , -f1 | cut -d - -f2)
	hit_count=$(echo $line | cut -d , -f2)
	sample_name=$(echo $line | cut -d . -f1)
	hit_count_trimmed=${hit_count%$'\r'}
	if (( "$hit_count_trimmed"> "10")); then
		temp_filename=$sample_name.$suffix
		#echo $temp_filename
		#if (("$iter" == "0" )); then
		#	samtools view -H $temp_filename >> $sample_name.filtered.sam
		#fi
	#iter=$(expr $iter + 1)
	#USCOUNTER=$(expr $USCOUNTER + 1)
	samtools view $temp_filename | grep -m 1 $read_id | head -n1 >> $sample_name.filtered.sam
	fi
done < $1 

for i in *.bam
do
	file=$(echo $i | cut -d . -f1)
	temp_filename=$file.filtered.sam
	sam_header=$(samtools view -H $i)
	echo $temp_filename
	echo -e "$sam_header$(cat $temp_filename)" > $temp_filename
	#echo $sam_header
	# echo "task goes here
	# $(cat $tepm_filename)" > $temp_filename
done

mkdir sams 
mv *.sam sams/
mv sams/ .. 