#! /usr/bin/env python3
import matplotlib.pyplot as plt

def readfq(fp):  # this is a generator function
    # From https://github.com/lh3/readfq/blob/master/readfq.py
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last: break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

def second_parse(s):
    inseq = False
    new_s = []
    t = ""
    for i in range(len(s)):
        for x in s[i]:
            if (x=='A' or x == 'T' or x == 'C' or x == 'G'):
                if (not inseq): 
                    inseq = True
                t = t + str(x)
            elif inseq and (not (x=='A' or x == 'T' or x == 'C' or x == 'G')):
                if(t != ""):
                    new_s.append(t)
                t = ""
                inseq = False
        if(t != ""):
            new_s.append(t)
    return new_s

def cumul_seq_length(s):
    n = 0
    for t in s:
        n = n+len(t)
    return n

s = []

file = open('file1.fa', 'r')
read = readfq(file)
read_save = list(read)
s.append(second_parse([read_save[i][1] for i in range(len(read_save))]))

file = open('file2.fa', 'r')
read = readfq(file)
read_save = list(read)
s.append(second_parse([read_save[i][1] for i in range(len(read_save))]))

file = open('file3.fa', 'r')
read = readfq(file)
read_save = list(read)
s.append(second_parse([read_save[i][1] for i in range(len(read_save))]))

file = open('file4.fa', 'r')
read = readfq(file)
read_save = list(read)
s.append(second_parse([read_save[i][1] for i in range(len(read_save))]))

file = open('file5.fa', 'r')
read = readfq(file)
read_save = list(read)
s.append(second_parse([read_save[i][1] for i in range(len(read_save))]))

file = open('file6.fa', 'r')
read = readfq(file)
read_save = list(read)
s.append(second_parse([read_save[i][1] for i in range(len(read_save))]))

print("file 1 : nb of seq :" + str(len(s[0])) + ", cumulative seq length : " + str(cumul_seq_length(s[0])) + '\n')
print("file 2 : nb of seq :" + str(len(s[1])) + ", cumulative seq length : " + str(cumul_seq_length(s[1])) + '\n')
print("file 3 : nb of seq :" + str(len(s[2])) + ", cumulative seq length : " + str(cumul_seq_length(s[2])) + '\n')
print("file 4 : nb of seq :" + str(len(s[3])) + ", cumulative seq length : " + str(cumul_seq_length(s[3])) + '\n')
print("file 5 : nb of seq :" + str(len(s[4])) + ", cumulative seq length : " + str(cumul_seq_length(s[4])) + '\n')
print("file 6 : nb of seq :" + str(len(s[5])) + ", cumulative seq length : " + str(cumul_seq_length(s[5])) + '\n')

def reverse_dna(dna):
    new_dna = dna
    n = len(dna)
    for i in range(n):
        if(dna[i] == 'A'):
            new_dna[n-1-i] = 'T'
        elif(dna[i] == 'T'):
            new_dna[n-1-i] = 'A'
        elif(dna[i] == 'C'):
            new_dna[n-1-i] = 'G'
        elif(dna[i] == 'G'):
            new_dna[n-1-i] = 'C'
        else:
            return list("The given string does not correspond to DNA.")
    return new_dna

def canonical_kmers(k, t):
    can_kmer_set = set()
    for i in range(len(t)-k):
        can_kmer_set.add(min("".join(t[i:i+k]), "".join(reverse_dna(list(t[i:i+k])))))
    return can_kmer_set

def canonical_kmers_list(k, t):
    total = set()
    for i in range(len(t)):
        total = total.union(canonical_kmers(k,t[i]))
    return total

def jaccard_index(st1, st2):
    i = len(st1.intersection(st2))
    u = len(st1.union(st2))
    return i/u

#for i in range(6):
#    for j in range(i):
#        print("jaccard index between file " + str(i+1) + " and " + str(j+1) + " : " + str(jaccard_index(canonical_kmers_list(20, s[i]), canonical_kmers_list(20, s[j]))))

# Looking at the dirrerent Jaccard indexes, I can infer that the genomes are not similar at all, except for 4 and 1 which are very similar, and (3 and 0) and (5 and 3)
# which are not as different from each other as the other pairs. 


def num_kmer(w):
    num_w = 0
    four_power = 1
    for i in range(len(w)):
        if(w[i] == 'A'):
            num_w = num_w + four_power * 0
        elif(w[i] == 'C'):
            num_w = num_w + four_power * 1
        elif(w[i] == 'G'):
            num_w = num_w + four_power * 2
        elif(w[i] == 'T'):
            num_w = num_w + four_power * 3
        four_power = four_power * 4
    return num_w

def kmer_counting(w,k):
    allwords = dict()
    F = 0
    for i in range(len(w)-k+1):
        if(not (num_kmer(w[i:i+k]) in allwords)):
            F = F+1
            allwords[num_kmer(w[i:i+k])] = True
    return F


x1 = []
y1 = []
for i in range(1,101):
    x1.append(i)
    y1.append(kmer_counting(s[1][0], i))

k = 100
plt.hist(x1,bins=k, weights=y1)
plt.show()
