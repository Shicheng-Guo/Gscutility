#!/usr/bin/env python
# from Devon Ryan in 2016
f = open("GSE17972_HUMtg5lib.qmap.chrY.txt") # Change me!
line = f.next().strip()
chrom = line[1:]
of = open("{}.bedGraph".format(chrom), "w")
pos = 0
for line in f:
    cols = line.split()
    M = int(cols[1])
    UM = int(cols[2])
    if M + UM >=5:
        of.write("{}\t{}\t{}\t{}\n".format(chrom, pos, pos + 1, 100. * float(M) / (M + UM)))
    pos += 1
of.close()
f.close()
