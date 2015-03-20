#!/usr/bin/env python

# this program makes dictionaries for UCE probe and contig files
# and matches the sequences

# function to read in the uce file
def read_uce_file(fname):
    lines = open(fname).read().splitlines()
    uce = {}
    for i in range(0, len(lines) - 1, 2):
        tok = lines[i].split()[0].replace(">", "")
        uce[tok] = lines[i + 1]
    return uce

# function to read contigs file
def read_contigs_file(fname):
    lines = open(fname).read().splitlines()
    contigs = {}
    for i in range(0, len(lines) - 1, 2):
        toks = lines[i].split()
        tok = toks[11] + " " + toks[12] + " " + toks[16]
        contigs[tok] = lines[i + 1]
    return contigs

# function to write a dictonary. Each line contains the key and the value
def write_dict(d, fname):
    fhand = open(fname, "w")
    for key in d:
        fhand.write(key + " " + d[key] + "\n")
    fhand.close()

# read in the probes file
uce = read_uce_file("uce-2.5k-probes.fasta")

# read in the contigs file
contigs = read_contigs_file("birds-contigs-assembled-from-uce-loci.fasta")

# write the uce dictionary
write_dict(uce, "uce-file.txt")

# write the contigs dictionary
write_dict(contigs, "contigs-file.txt")

# match UCE probes to contigs
out = open("script-output.txt", "w")
# loop over each probe seq
for probe_seq in uce:
    # now loop over each seq from the contigs file
    for bird_seq in contigs.keys():
        # if we get a match print the uce and contigs info plus
        # where the probe fits into the contigs seq
        if uce[probe_seq] in contigs[bird_seq]:
            contigs_split = contigs[bird_seq].split(uce[probe_seq])
            out.write(probe_seq + " " +
                      bird_seq + " " +
                      contigs_split[0] + " " +
                      uce[probe_seq] + " " +
                      contigs_split[1] + "\n")
out.close()
