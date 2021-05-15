#Origionally modified from code by Jonathan Pekar
import json

#dict from https://pythonforbiologists.com/dictionaries
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
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', "---": "-"}


input_seq_path = "/mnt/c/Users/hjenn/Documents/1_UCSD/Research/COVID/Muller_Lineage/15_updatedData/GISAID_NYCPHL_04_29/gisaid_2021_04_30_00_rename_aln.fasta"
outPrefixPath = "/mnt/c/Users/hjenn/Documents/1_UCSD/Research/COVID/Muller_Lineage/15_updatedData/GISAID_NYCPHL_04_29/renameAln_codingRegions"
annotation_path = '/mnt/c/Users/hjenn/Documents/1_UCSD/Research/COVID/epiMuller/data/genemap_SARSCoV2.gff'


fasta = {}
input_seq = open(input_seq_path, "r")
idLine = 'na'
for line in input_seq:
    if ">" in line:
        if idLine != "na":
            fasta[idLine] = seq
        idLine = line.strip()
        seq = ""
    else: 
        seq = seq + line.strip()
fasta[idLine] = seq

fasta_coding_d = {}
fasta_aa_d = {}
stop_codon_dict = {}
gene_boundry_d = {}
aa_index = 0

for l in open(annotation_path).readlines():
    if '#' in l:
        continue
    l = l.strip().split('\t')
    start = int(l[3])
    end = int(l[4])
    gene = l[-1].split('"')[1]
    if gene == 'ORF9b' or gene == 'ORF14': 
        continue
    

    #fasta_path = 'PATHTODFASTA_%s.fasta' % gene
    seq_set = set()
    seq_keys = []
    stop_codon_dict[gene] = {}
    
    #with open(fasta_path, 'w') as f:
    for k in fasta.keys():
        seq = fasta[k][start-1:end][:-3] # exclude stop codon

        seq_noStop = ''
        aaSeq_noStop = ""
        hadStop = False
        for i in range(int(len(seq)/3)):
            codon = seq[i*3:i*3+3].upper()
            if codon == 'TGA' or codon == 'TAA' or codon == 'TAG' or hadStop:
                hadStop = True
                seq_noStop += '---'
                aaSeq_noStop += "-"
            else:

                seq_noStop += codon
                if "N" in codon:
                    aa = "-"
                elif "-" in codon:
                    aa = "-"
                else:
                    if codon not in gencode:
                        aa = "-"
                    else:
                        aa = gencode[codon]
                aaSeq_noStop += aa

        #f.write('%s\n%s\n' % (k, seq_noStop)) 
        if k not in fasta_coding_d:
            fasta_coding_d[k] = ""
        fasta_coding_d[k] += seq_noStop

        if k not in fasta_aa_d:
            fasta_aa_d[k] = ""
        fasta_aa_d[k] += aaSeq_noStop

    gene_boundry_d[gene] = [aa_index, aa_index+int(len(seq)/3)]
    aa_index += int(len(seq)/3)


fasta_path = outPrefixPath + '_allCDS.fasta'
f = open(fasta_path, "w")
for k in fasta_coding_d:
    f.write('%s\n%s\n' % (k, fasta_coding_d[k])) 


fasta_aa_path = outPrefixPath + '_allAA.fasta'
f = open(fasta_aa_path, "w")
for k in fasta_aa_d:
    f.write('%s\n%s\n' % (k, fasta_aa_d[k])) 


json_path = outPrefixPath + '_geneAAboundries.json'
with open(json_path, 'w') as fp:
    json.dump(gene_boundry_d, fp)
    print("dumped json")