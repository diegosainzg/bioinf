from Bio import SeqIO
fh=open("alineated.fa")
ids=[]
sequences=[]
for sequence in SeqIO.parse(fh,"fasta"):
   ids.append(sequence.id)
   sequences.append(sequence.seq)


def startstop_codon(dna ,frame):
    dna=dna.upper()
    for i in range(frame,len(dna),3):
        codon1=dna[i:i+3]
        if codon1=="ATG":
            pos1=i
            for j in range(pos1,len(dna),3):
                codon2=dna[j:j+3]
                if codon2 in ["TAG","TGA","TAA"]:
                    pos2=j
                    yield(pos2-pos1+3,dna[pos1:pos2+3]) #returns a generatior with the length and the orf
                    break
def resume_data(iterable):
    dist_lists=[]
    seq_list=[]
    for orflen, orf in iterable:
        dist_lists.append(orflen)
        seq_list.append(orf)
    return(dist_lists,seq_list)


def dnaToProtein(dna):
    map = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
           "UCU": "S", "UCC": "s", "UCA": "S", "UCG": "S",
           "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
           "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
           "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
           "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
           "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
           "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
           "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
           "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
           "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
           "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
           "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
           "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
           "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
           "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
    seq=dna.replace('T','U')
    protein_seq=""
    for i in range(0,len(seq),3):
        codon=seq[i:i+3]
        aa=map[codon]
        if aa =="STOP": # Stop codon
            break
        protein_seq+=aa
    return(protein_seq)

#distlist, seqlist= resume_data(startstop_codon(sequences[0],0))
#max_seq=seqlist[distlist.index(max(distlist))]
#print(dnaToProtein(str(max_seq)))
def consensus(alignment,threshold=0.5):
    n = len(alignment[0]) # Number of aa each seq
    nSeq = float(len(alignment))
    consensus = ""
    for i in range(n): # each column
        counts={}
        for seq in alignment:
            letter = seq[i]
            if letter == "-":
                continue # Space detected
            counts[letter] = counts.get(letter,0) +1 # 0 by default
        fractions = []
        for letter in counts:
            frac = counts[letter]/nSeq
            fractions.append([frac,letter])
        fractions.sort() #sort the list of each by the fraction
        bestFraction, bestLetter = fractions[-1]
        if bestFraction<threshold:
            consensus += "X"
        else:
            consensus += bestLetter
    return consensus
print ( consensus(sequences),"s")








