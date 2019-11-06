#Workflow for generating CMV figures 

#Bowtie2 alignment to CMV Merlin genome (NC_006273.2)
bowtie2 -x <hhv-6a/b_reference> -1 <r1.fastq> -2 <r2.fastq> --local --no-unal -p <cores> -S <aligned.sam> 

#Convert SAM to BAM 
samtools view -@ <cores> -Sb <aligned.sam> -o <aligned.bam> 

#sort BAM for picard
samtools sort <aligned.bam> <sorted.bam> 

#Deduplicate aligned bam file 
picard MarkDuplicates I=<sorted.bam> O=<deduplicated.bam> M=<metrics.txt> REMOVE_DUPLICATES= TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY= SILENT

#Only keep reads with 34 or more matches and extract reads for repeatmasking and BLAST later.
./34m.sh <deduplicated_bam_folder> 

#Repeatmask output of previous command
Repeatmasker -int -pa <cores> hhv6_combined.txt
mv hhv6_combined.txt.masked hhv6_combined_masked.fasta 

#remove all repeat masked lines with NNNNN in it 
for file in hhv6_combined_masked.fasta
do
sed '$!N;/NNNNN/!P;D' "$file"
done > hhv6_combined_masked_n_removed.fasta

#Local BLAST against NT (requires blast NT database)
blastn -query hhv6_combined_masked_n_removed.fasta -db /db/blast_db/nt -num_threads 42 -perc_identity 97 -evalue 1e-5 -out hhv6_vs_full_nt.txt

#Creates count table of BLAST hits to "Human hepesvirus 5" from the BLAST results.
python blast_hits.py hhv6_masked_blastn_out.txt '<Human herpesvirus 6>'

#Creates BAMs with just the verified reads 
#run from within folder of filtered bams for both HHV-6A and HHV-6B (./34m.sh step)
uw-virongs-vikas:verified_sams gerbix$ python verified_bam.py blast_hits<HHV6a/b>.csv <suffix for filtered bams> 


#creates all_sample_data.csv- A table of sample name, FPM, FPKM and read counts
#file paths inside the script require editing based on how the above commands were run
rscript --vanilla fpm_calculate.R

#for the scripts below file paths will have to be changed from the ones used on our local machine
#figure 5 inidividual figures
#path to any bams should be from the verified bams folder
rscript --vanilla figure_5A.r
rscript --vanilla figure_5B.r
rscript --vanilla figure_5C.r

#figure 5 panel
rscript --vanilla figure_5.r

