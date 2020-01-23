# read fasta file
def open_fasta(fastaname):
    fh = open(fastaname)
    sequences = []  # gobal list contatining all sequences
    seq = ""  # partial seq
    for line in fh:
        if line.startswith(">"):
            if seq:  # Not the first sequence to read, so seq was not empty
                sequences.append(seq)  # append last read seq
            seq = ""  # First sequence to read
        else:  # it is a line with bases
            seq += line.rstrip()  # remove newline characters
    if seq:  # for the last seq read
        sequences.append(seq)
    fh.close()
    return sequences

sequences =open_fasta("chromo.fasta")
stop_seq=["TAG","TAA","A","T"]
my_seq= sequences[0]
orfs=[]
orf=""
for i in range(len(my_seq)-2):
    codon=my_seq[i:i+3]
    if orf: #an ATG already started
        if codon in stop_seq: #there is a stop
            orfs.append(orf)
            orf = ""
            continue
        else: #Not a stop, so add the codon
            orf+=codon
            i = i + 3
        continue
    else: #not any ATG
        if codon=="ATG": #Initial codon
            if orf:#we append the last orf read
                orfs.append(orf)
                orf=""
                orf += codon
                i +=3

            else:
                orf=""
                orf+=codon
                i+=3
                continue

        else:
            continue
print(orfs)






