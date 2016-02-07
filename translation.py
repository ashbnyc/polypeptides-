__author__ = 'anb2126'
# Translation of genetic code to protein
# Step 1) Split up sequence in to codons
# Step 2) Translate each codon to amino acid
# Step 3) Assemble chain of amino acids
# Step 4) Convert polypeptide to list of janin hydrophobicity values (integers)
# Step 5) Plot janin integers on plot of aa position (x) vs janin value
# Step 6) Find electrodensity size scale (Pymol?)
# Step 7) Repeat steps 4-6 with size scale
# Step 8) Find method to overlay/compare other picorna Vp1s using the above

# Code for method which does not print non-codons at end

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

janin_code = {
    'I': 0.7, 'F': 0.5, 'V': 0.6, 'L': 0.5, 'W': 0.3, 'M': 0.4,
    'A': 0.3, 'G': 0.3, 'C': 0.9, 'Y': -0.4, 'P': -0.3, 'T': -0.2,
    'S': -0.1, 'H': -0.1, 'E': -0.7, 'N': -0.5, 'Q': -0.7, 'D': -0.6, 'K': -1.8, 'R': -1.4}

kd_code = {
    'I': 4.5, 'F': 2.8, 'V': 4.2, 'L': 3.8, 'W': -0.9, 'M': 1.9, 'A': 1.8, 'G': -0.4, 'C': 2.5,
    'Y': -1.3, 'P': -1.6, 'T': -0.7, 'S': -0.8, 'H': -3.2, 'E': -3.5, 'N': -3.5, 'Q': -3.5,
    'D': -3.5, 'K': -3.9, 'R': -4.5}

ew_code = {
    'I' : 0.73, "F" : 0.61, 'V': 0.54, 'L': 0.53, 'W': 0.37, 'M': 0.26, 'A': 0.25,
     'G': 0.16, 'C': 0.04, 'Y': 0.02, 'P': -0.07, 'T': -0.18, 'S': -0.26, 'H': -0.40,
     'E': -0.62, 'N': -0.64, 'Q': -0.69, 'D': -0.72, 'K': -1.10, 'R': -1.8}



def translate_dna(dna):
    last=len(dna)-2
    print(last)
    polypep=""
    for base in range(0,last,3):
        codon=dna[base:base+3]
        print("Codon:",codon)
        amino=gencode.get(codon.upper(), 'X')
        print(amino)
        polypep+=amino
    print(polypep)
    return(polypep)

def janin_scale(polypep):
    janin_points = []
    for aa in polypep:
        janin = janin_code.get(aa)
        print(janin)
        janin_points.append(janin)
    print("\n Janin_points:", janin_points )
    return [janin_points]

def kd_scale(polypep):
    kd_points = []
    for aa in polypep:
        kd = kd_code.get(aa)
        print(kd)
        kd_points.append(kd)
    print("\n Kyte and Doolittle points:", kd_points)
    return [kd_points]

def ew_scale(polypep):
    ew_points = []
    for aa in polypep:
        ew = ew_code.get(aa)
        print(ew)
        ew_points.append(ew)
    print("\n Eisenberg Weis points:", ew_points)
    return [ew_points]

dna=input("Paste Sequence for Translation Here:")

dna=dna.upper()

translate_dna(dna)

pp = translate_dna(dna)

janin_scale(pp)

kd_scale(pp)

ew_scale(pp)
