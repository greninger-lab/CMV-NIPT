import random
output = open("MultiSampleRead3Probe350_alex.txt", "w")

probelen = 350
probestart = random.randint(1, 235149 - probelen)
probestop = probestart + probelen - 1

for a in range (0, 100):
        print(str(a))
        input1 = open("Sample2Read3_alex.txt", "r")  #1st column read length, 2nd colum number of reads

        readstart = []
        readstop = []
        readcount = 0

        for line in input1:
                line = line.strip()
                line = line.split("\t")
                count = int(line[1])
                for x in range (0, count):
                        z = 235149 - int(line[0])
                        y = random.randint(1, z)
                        readstart.append(y)
                        y = y + int(line[0]) - 1
                        readstop.append(y)
                        readcount = readcount + 1

        input1.close()

        load = 0        

        for x in range (0, readcount):
                if probestart >= readstart[x] and probestop <= readstop[x]:
                        load = load + 1
        output.write(str(load))
        output.write("\n")

output.close()          
