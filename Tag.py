#1/bin/python3

from Bio import SeqIO
import sys

input_files = []

def read_fasta_file(file_path):
    sequence = ""
    # Parse the FASTA file
    records = SeqIO.parse(file_path, "fasta")

    # Iterate through each record in the file
    for record in records:
        # Access the sequence object
        sequence += str(record.seq).lower()

        # Print or manipulate the sequence as needed
        #print(f"Header: {record.id}")
        #print(f"Sequence: {sequence}")
        #print()
    return sequence
    

def read_Zero_Reference_File(file_path):
    try:
        with open(file_path, 'r') as f:
            ZeroDNAref = f.read()
        return ZeroDNAref
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None


for i, arg in enumerate(sys.argv):
    if i >= 1:
        input_files.append(arg)
        print(arg)

for file in input_files[:]:

    #fasta_file_path = "GCF_000001405.40_GRCh38.p14_genomic.fna" #Human Genome - February 2022
    #fasta_file_path = "GCF_000001635.27_GRCm39_genomic.fna" #House Mouse Genome (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)
    #fasta_file_path = "Physarum_polycephalum_AB027295.1.fasta"
    #fasta_file_path = "Monkey_GCF_003339765.1_Mmul_10_genomic.fna"

    fasta_file_path = file #"YOUR FASTA FILENAME HERE!"

    DNAZeroRefFile = "Qubits_Split.DNA"
    
    OutFile = file.split(".")[0] + "_Freq_Matched_DNA.Zero"

    # Call the function with the file path
    DNAsequence = read_fasta_file(fasta_file_path)
    ZeroDNAref = read_Zero_Reference_File(DNAZeroRefFile)

    ZeroDNAref = ZeroDNAref[22000:len(ZeroDNAref)-500000] # stripping leading and trailing leaving DNA Zero

    for i in range(len(ZeroDNAref)): #Adjusting the loop.
        if ZeroDNAref[i] == ZeroDNAref[-1]:
            ZeroDNAref = ZeroDNAref[i:len(ZeroDNAref)-1]
            break

    Aligned = False
    AlignedIndex = 0
    Count = 0
    Speed_of_Light = 299792458
    NucleotideSize = 0.33 / 1e9 #nm to meters
    Frequencies = []

    nucleotideIndex = 0
    f = open(OutFile, "w")
    for nucleotide in DNAsequence:
        if not Aligned:
            if nucleotide != "n":
                for i in range(len(ZeroDNAref)):
                    if nucleotide == ZeroDNAref[i]: #aligned to DNA Zero
                        Aligned = True
                        AlignedIndex = i
        elif Aligned:
            if ((nucleotideIndex + 1) < (len(DNAsequence) -1)) and (AlignedIndex < len(ZeroDNAref) -1):
                if nucleotide == ZeroDNAref[AlignedIndex] and DNAsequence[nucleotideIndex+1] == ZeroDNAref[AlignedIndex+1]:
                    WaveLength = Count * NucleotideSize
                    Frequencies.append(Speed_of_Light/WaveLength)
                    #f.write("%d=|%s-%s|" % (Count, nucleotide, DNAsequence[nucleotideIndex+1]))
                    if AlignedIndex + 1 == len(ZeroDNAref) -1: #loop it back
                        AlignedIndex = 0
                    else:
                        AlignedIndex = AlignedIndex + 1
                    Count = 0
            elif AlignedIndex == len(ZeroDNAref)-1:
                AlignedIndex = 0

            #else:
                #f.write("%d=|%s|--------|%s|" % (Count, nucleotide, DNAsequence[nucleotideIndex+1]))
        Count = Count + 1
        nucleotideIndex = nucleotideIndex + 1

    Spectrum = {}
    for Frequency in Frequencies:
        if Frequency not in Spectrum:
            Spectrum.update({Frequency:0})
        else:
            Spectrum.update({Frequency:Spectrum[Frequency]+1})

    sortedSpectrum = dict(sorted(Spectrum.items()))

    for Frequency in sortedSpectrum:
        f.write("Frequency=%f Count=%d\n" % (Frequency, sortedSpectrum[Frequency]))

    f.close()
    print("Saved...%s" % (OutFile))