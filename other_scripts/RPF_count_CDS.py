import HTSeq
import sys
from collections import Counter
"""
usage: python RPF_counts_CDS.py bamFile gtfFile
change the following parameters for specific instance: 
	exclude_start_distance = 45 # 15 codons
	exclude_stop_distance = 15 # 5 codons
	minLen = 26
	maxLen = 34
only used for strand specific mapping.
"""
# if there are multiple TIS, use the most 5' end start codon and the most 3' end stop codon
# read the gtf file
def readGTF(gtfFile):
	gtf = HTSeq.GFF_Reader(gtfFile)
	start_codon_sites = {}  
	stop_codon_sites = {} 
	CDS_features = HTSeq.GenomicArrayOfSets("auto", stranded="yes") 
	i = 0
	for f in gtf:
		i += 1
		if i % 10000 == 0:
			sys.stderr.write("%d GFF lines processed.\r" % i)
		gname = f.attr['gene_id']
		if f.type == "CDS":
			CDS_features[f.iv] += gname
		if f.type == "start_codon":
			if gname not in start_codon_sites:
				start_codon_sites[gname] = f.iv.start_d 
			else:
				if f.iv.strand == "+":
					start_codon_sites[gname] = min(f.iv.start, start_codon_sites[gname])
				else:
					start_codon_sites[gname] = max(f.iv.start_d, start_codon_sites[gname])
		if f.type == "stop_codon":
			if gname not in stop_codon_sites:
				stop_codon_sites[gname] = f.iv.end_d
			else:
				if f.iv.strand == "+":
					stop_codon_sites[gname] = max(f.iv.end, stop_codon_sites[gname])
				else:
					stop_codon_sites[gname] = min(f.iv.end_d, stop_codon_sites[gname])
	return start_codon_sites, stop_codon_sites, CDS_features

#main function
##define the paramaters:
bamFile = sys.argv[1]
gtfFile = sys.argv[2]

exclude_start_distance = 45 # 15 codons
exclude_stop_distance = 15 # 5 codons
minLen = 26
maxLen = 34
#read the gtfFile
start_codon_sites, stop_codon_sites, CDS_features = readGTF(gtfFile)

#read the bamFile
#only unique mapped reads are used
#strand specific, intersection_straict mode.
counts = Counter()
empty = 0
ambiguous = 0
lowqual = 0
notaligned = 0
nonunique = 0
bam = HTSeq.BAM_Reader(bamFile)
for r in bam:
	if not r.aligned:
		notaligned += 1 
		continue
	if r.optional_field("NH") > 1:
		nonunique += 1
		continue
	if r.iv.chrom in ["MT","chrM","chrMT"]:   #skip the mitochondria chrosome, change this if your name is different.
		continue
	if minLen<= len(r.read.seq) <=maxLen:
		pass
	else:
		continue 

	iv_seq = (co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0)

	fs = None
	for iv in iv_seq:
		for iv2, fs2 in CDS_features[iv].steps():
			if fs is None:
				fs = fs2.copy()
			else:
				fs = fs.intersection(fs2)


	if fs is None or len(fs) == 0:
		empty += 1
	elif len(fs) >1:
		ambiguous += 1
	else:
		gname = list(fs)[0]
		try:
			if abs(start_codon_sites[gname] - r.iv.start_d) < exclude_start_distance:
				continue
			elif abs(r.iv.end_d - stop_codon_sites[gname]) < exclude_stop_distance:
				continue
			else:
				counts[gname] += 1
		except:
			counts[gname] += 1

for g in sorted(counts):
	print "%s\t%d" % (g, counts[g])

print "__no_feature\t%d" % empty
print "__ambiguous\t%d" % ambiguous
print "__too_low_aQual\t%d" % lowqual
print "__not_aligned\t%d" % notaligned
print "__alignment_not_unique\t%d" % nonunique
