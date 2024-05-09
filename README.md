## TL-Seq Analysis for transcripts




### Description 
This application estimates new RNA transcripts using metabolic labelling which converts T to C 


#### **Installation:**


#### **Usage:**
Usage:

python fragment_preprocessing.py [Options]


Options:
- -h, --help: usage
- -f, --folder_location: Directory to save the output for the analysis
- -p, --pear_merged_loc: Directory for merged output files
- -b, --base_quality: Directory of the base quality
- -m, --metadata: Path to the sample metadata CSV file
- -v, --vcf_file: Path to the VCF file for SNPs
- -r, --ref_loc: Path to the reference sequence file for the sample
- -chr, --chromosome: Chromosome in which the gene of interest is located (choices: 1-23, X, Y, x, y)
- -gs, --gene_pos_start: Start position of the gene of interes
- -ge, --gene_pos_end: End position of the gene of interes
- -s, --strandedness: Strandedness of the gene, either '+' or '-' 



**Sample metadata:**


|Sample name | Region | Time | Treatment | F | R | Product sequence|
|------------|--------|------|-----------|---|---|-----------------|
|c_con_exon2 | exon2  |con   |c|TCACTGTCCTTCTGCCATGG|CCCGCACACTAGGTAGAGAG|TCACTGTCCTTCTGCCATGGccctgtggatgcgcctcctgcccctgctggcgctgctggccctctggggacctgacccagccgcagcctttgtgaaccaacacctgtgcggctcacacctggtggaagCTCTCTACCTAGTGTGCGGG|
|woc_con_exon2 | exon2  |con   |woc          |TCACTGTCCTTCTGCCATGG|CCCGCACACTAGGTAGAGAG|TCACTGTCCTTCTGCCATGGccctgtggatgcgcctcctgcccctgctggcgctgctggccctctggggacctgacccagccgcagcctttgtgaaccaacacctgtgcggctcacacctggtggaagCTCTCTACCTAGTGTGCGGG |
