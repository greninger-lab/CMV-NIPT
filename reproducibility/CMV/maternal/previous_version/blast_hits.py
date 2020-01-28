#!/usr/bin/env python
####MODIFIED FOR LOCAL BLAST#### 
#args[1] = blast file
#args[2] = virus identifier

import sys 
import csv
import re

blastfile = open(sys.argv[1], "r").readlines()
#print blastfile

#deletelines=[]
print(sys.argv[2])
count = 0 
for i in range(1, len(blastfile)): 
	#print(type(blastfile[i]))
	if(re.search('-$', blastfile[i])): 
		count +=1
		#print(blastfile[i].rstrip("\n"))
		blastfile[i] = blastfile[i].rstrip("\n") + blastfile[i+1].rstrip()
		#print(blastfile[i])
		#i +=1
		#deletelines.append((i+1))
		# print(blastfile[i])
		#print(blastfile[i + 1])
		#print(blastfile[i+1])
#print(count)


querylist= []
count_start=14
matches_start=0
matches_end=0
alignments_start=[]
alignments_end=[]
names=[]
for i in range(1,len(blastfile)):
	numlines=0
	if "Query=" in blastfile[i]:
		#print(blastfile[i])
		# should use a single if not statement here 
		# if blastfile[i+6][0] == '*':
		# 	print ('fuck')
		if blastfile[i+5][0]!='*' and blastfile[i+6][0]!='*':
			#print(blastfile[i])
			#print(blastfile[i+6])

			# default argument of split is whitespace 
			temp_name=blastfile[i].split(" ")[1]
			#print (temp_name.rstrip())
			names.append(temp_name.rstrip())
			#print(temp_name)
			continue
			#print(temp_name.rstrip())
		#print(temp_name)

	if "alignments:" in blastfile[i]:
		matches_start=(i+3)
		alignments_start.append(matches_start)
		continue
		#print(matches_start)
		#print('ok')
		#print(matches_start)
	if blastfile[i][0] == ">" and "\n" in blastfile[i-1] and (blastfile[i-3][0]!='S'):
		#print(blastfile[0])
		#print (blastfile[i])
		matches_end=(i-2)
		alignments_end.append(matches_end)
		continue
	if(len(alignments_end) > len(alignments_start)): 
		print(i)
		print('fuck')
		#print(matches_end)
		#print (temp_name)
		#print ('\n')
		# print ('new')
		# print(blastfile[i-3])
		# print(blastfile[i-2])
		# print(blastfile[i-1])
		# print(blastfile[i])
		#print ('ok')
	#print ("\n")
		#print(matches_end)
	#print(range(matches_end,matches_start))
	#hhv6_count=0
	# for j in range(matches_end,matches_start): 
	# 	if ('Herpesvirus 6') in blastfile[j]:
	# 		hhv6_count+=1
	# 	if(hhv6_count>1):
	# 		print(hhv6_count)
		# count_end=i
		# #print(blastfile[count_end])
		# numlines=count_end-count_start
		# print(numlines)

		#count_start=count_end
#print(blastfile[11])
#print(len(names))
data=[]
perc_total=[]


for i in range(len(alignments_start)):
	#print (len(alignments_start))
	#print (len(alignments_end))
	#print('\n')
	count=0
	temp=[]
	# if (alignments_start[i]==alignments_end[i]):
	# 	print ("ok")
	#print(alignments_start[i])
	#print(alignments_end[i])
	for j in range(alignments_start[i],alignments_end[i]):
		if sys.argv[2] in blastfile[j]:
			count+=1
			# print(blastfile[j])
			# print(names[i])
			# print(j+1)
			# print(alignments_start[i])
			# print(alignments_end[i])
			# print('\n')
			# print(blastfile[j])
	#temp_perc = 100* (count/(alignments_end[i] - alignments_start[i]))
	#perc_total.append(temp_perc)
	temp=[names[i].strip(),count]
	#print (names[i].strip())
	#print(temp)
	data.append(temp)
#print(data)
myFile = open('blast_hits.csv', 'w')
with myFile:
    writer = csv.writer(myFile)
    writer.writerows(data)
print("file written as blast_hits.csv")
	#print(i)		
	# 		hhv6_count+=1
	# 	if(hhv6_count>1):
	# 		print(hhv6_count)
#print(len(names))
print(len(alignments_start))
print(len(alignments_end))
print(len(names))

#print(names)

###make a temp list of lines between alignment_end and alignment_start 
###take that list and grep for herpesvirus 6 
###count that and put in a line with the name in a .csv












