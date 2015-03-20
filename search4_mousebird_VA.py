#!/usr/bin/env python

# this program takes a pre-existing alignment file and a UCE probe sequence and
# splits up the contig match sequences into three pieces:
# pre-sequence, probe/core, and post-sequence 
# and makes comparisons across species

from __future__ import division

# create a dictionary for a UCE alignment file

def read_nex(fname):
    nex = {} # creating a dictionary
    lines = open(fname).read().splitlines()
    for line in lines:
        # skip comments
        if line[0] == "#":
            continue
        # skip the discriptions
        if ";" in line:
            continue
        # skip the data start
        if line == "matrix":
            continue
        toks = line.split() # splits the line at any white space into species and sequence
        nex[toks[0]] = toks[1] # for the key / species I'm assigning a value / sequence

    return nex



# begin main program
print("UCE probe-contig matches: flanks and cores")

# using a probe sequence from the original uce-2.5k-probes.fasta file and a nexus alignment file created during McCormack et al's analysis
### potential issue: what if not all species are a perfect match with the probe sequence? Will they be skipped over?

probe = "CTCTGACTCCTGGTGTTCTTTATTCCAATAATTTAAATTGGAATTGATTTTTCAATTTAGCCGTTGACTTGCTATGGACACACAGCTGAAAAGGCCGAATGGAATTGCATGGCTGAAATA"

nex = read_nex("chr1_32397.nex")

out = open("pre&post-output.txt", "w")
for key in nex:
    if probe in nex[key]:
        seq = nex[key].split(probe) # split command finds all locations of the specified split string (default is white space), the probe sequence in this case
        print "*********************" # prints stars
        print "name =", key
        print "pre-seq =", seq[0] 
        print "probe/core =", probe
        print "post-seq =", seq[1]
        print "pre-seq length =", len(seq[0])
        print "post-seq length =", len(seq[1])
        out.write("*********************" + "\n"
                  "name = " + key + " \n" +
                  "pre-seq = " + seq[0] + " \n" +
                  "probe/core = " + probe + " \n" +
                  "post-seq = " + seq[1] + " \n")
        
        

# splitting the DNA sequence by the probe will return an array with 2 elements.  [0] is the sequence before the probe and [1] is the sequence after the probe.
# since we lost the probe sequence as part of the splitting we need to add it back when we print the result.
# meaning that seq[0] + probe + seq[1] is just the original DNA sequence
# lengths should be consistent




# cross-species comparisons:

# how do I quantify differences pre-seq's and post-seq's?
# want to use sequences from gallus_gallus as the reference (basis of comparison)
# came across levenshtein distance on stack overflow (http://stackoverflow.com/questions/24204715/matching-2-strings-and-allowing-a-5-mismatch-rate)
# using code from stack overflow and wikibooks (http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python)

print("\n")
print("Cross-species comparisons of UCE flanks")
print("\n")

def levenshtein(s, t):
        ''' From Wikipedia article; Iterative with two matrix rows. '''
        if s == t: return 0
        elif "----------" in t: return "NA" # potential issue: what if species have a lot of no scores (-)? Should their levenshtein value be NA?
        elif len(s) == 0: return len(t)
        elif len(t) == 0: return len(s)
        v0 = [None] * (len(t) + 1)
        v1 = [None] * (len(t) + 1)
        for i in range(len(v0)):
            v0[i] = i
        for i in range(len(s)):
            v1[0] = i + 1
            for j in range(len(t)):
                cost = 0 if s[i] == t[j] else 1
                v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
            for j in range(len(v0)):
                v0[j] = v1[j]
 
        return v1[len(t)]



print("5' pre-seq flanks:")
print("\n")

# create a dictonary for the pre- probe/core sequences
pre_seq = {}
for key in nex:
    if probe in nex[key]:
        seq = nex[key].split(probe)
    pre_seq[key] = seq[0]

# save the urocolius_indicus pre sequence
urocolius_indicus = pre_seq["urocolius_indicus"]

# levenshtein comparison loop
for key in pre_seq:
    print key, levenshtein(urocolius_indicus, pre_seq[key])
    
# levenshtein comparison loop again to
# print the distance and create a dictionary
# and save the data to a file for plotting in R
out = open("species-predist-urocolius.txt", "w") # create the file
out.write("species distance\n") # write column headers
preseq_dist = {}
for key in pre_seq:
    lev_dist = levenshtein(urocolius_indicus, pre_seq[key])
    preseq_dist[key] = lev_dist
    out.write(key + " " + str(lev_dist) + "\n")
    
print("\n")
print("3' post-seq flanks:")
print("\n")

# create a dictonary of the post- probe/core sequences
post_seq = {}
for key in nex:
    if probe in nex[key]:
        seq = nex[key].split(probe)
    post_seq[key] = seq[1]

# save the urocolius_indicus post sequence
urocolius_indicus = post_seq["urocolius_indicus"]

# levenshtein comparison loop
for key in post_seq:
    print key, levenshtein(urocolius_indicus, post_seq[key])
    
# levenshtein comparison loop again to
# print the distance and create a dictionary
# and save the data to a file for plotting in R
out = open("species-postdist-urocolius.txt", "w") # create the file
out.write("species distance\n") # write column headers
postseq_dist = {}
for key in post_seq:
    lev_dist = levenshtein(urocolius_indicus, post_seq[key])
    postseq_dist[key] = lev_dist
    out.write(key + " " + str(lev_dist) + "\n")

out.close()
    
#
# Make a simple plot of the distances
#