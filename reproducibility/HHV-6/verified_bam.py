import sys 
import glob
import csv
import subprocess 
blast_hits_file = sys.argv[1]
suffix = sys.argv[2]

#suffix = "*" + suffix

file_list = glob.glob( "*" + suffix)

duplicate_list=[]

with open(blast_hits_file, 'r') as hitsfile:
	reader = csv.reader(hitsfile)
	for row in reader:
		#print(row)
		read_id = row[0].split('-')[1]
		sample = row[0].split('-')[0].split('.')[0]
		hit_count = int(row[1])
		check= sample + read_id
		#print(duplicate_list)
		if hit_count > 5 and check not in duplicate_list: 
			temp_filename = sample + suffix 
			#print(temp_filename)
			grep_cmd = "samtools view " + temp_filename + " | grep " + read_id + " >> " + sample + ".filtered.sam"
			duplicate_list.append(check)
			subprocess.run(grep_cmd, shell = True)

subprocess.call('/Users/gerbix/Documents/vikas/NIPT/nipt_git_repo/reproducibility/HHV-6/append_sam_header.sh', shell=True)