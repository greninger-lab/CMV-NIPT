#Workflow for generating CMV figures 

#Bowtie2 alignment to CMV Merlin genome (NC_006273.2)
bowtie2 -x <cmv_merlin_reference> -1 <r1.fastq> -2 <r2.fastq> --local --no-unal -p <cores> -S <aligned.sam> 

#Convert SAM to BAM 
samtools view -@ <cores> -Sb <aligned.sam> -o <aligned.bam> 

#Deduplicate aligned bam file 
picard MarkDuplicates I=<aligned.bam> O=<deduplicated.bam> M=<metrics.txt> REMOVE_DUPLICATES= TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY= SILENT

#Only keep reads with 34 or more matches and extract reads for repeatmasking and BLAST later.
./34m.sh <deduplicated_bam_folder> 

#Repeatmask output of previous command
Repeatmasker -int -pa <cores> cmv_combined.txt
mv cmv_combined.txt.masked cmv_combined_masked.fasta 

#remove all repeat masked lines with NNNNN in it 
for file in cmv_combined_masked.fasta
do
sed '$!N;/NNNNN/!P;D' "$file"
done > cmv_combined_masked_n_removed.fasta

#Local BLAST against NT (requires blast NT database)
blastn -query cmv_combined_masked_n_removed.fasta -db /db/blast_db/nt -num_threads 42 -perc_identity 95 -evalue 1e-5 -out cmv_vs_full_nt.txt

#Creates count table of BLAST hits to "Human hepesvirus 5" from the BLAST results.
python blast_hits.py cmv_masked_blastn_out.txt 'Human herpesvirus 5'

#creates all_sample_data.csv- A table of sample name, FPM, FPKM, read counts, and classification (strong vs intermediate positive based on the FPM cutoff of .3)
#file paths inside the script require editing based on how the above commands were run
rscript --vanilla fpm_calculate.R

#figure 1B
rscript --vanilla figure_1b.r

#figure 3
#run the first script to create a frequency table of insert sizes 
samtools view -@ 40 <deduplicated_121r04_resequenced.bam> | awk '{print $9}' | sort | uniq -c > duplicates_removed_read_lengths.txt
rscript --vanilla figure_3.r