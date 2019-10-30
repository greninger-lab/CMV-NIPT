cd $1 
for i in *.bam; do echo $i; samtools view -H  $i >> $i.filtered.sam ; samtools view $i | grep '37M' >> $i.filtered.sam 
samtools view $i | grep '36M' >> $i.filtered.sam
samtools view $i | grep '35M' >> $i.filtered.sam
samtools view $i | grep '34M' >> $i.filtered.sam
 done

mkdir ../filtered_sams
mv *.sam ../filtered_sams
cd ../filtered_sams
mkdir ../filtered_bams
for i in *.sam 
do
	samtools view -Sb -@8 $i > $i.bam
done

mv *.bam ../filtered_bams
cd ..

rscript --vanilla makefastasfromsams.r filtered_bams/
cd filtered_bams

find . -name "*.fasta" | xargs -n 100 -P 8 cat >> resequenced_combined.txt

mv resequenced_combined.txt .. 

echo '34 finished'