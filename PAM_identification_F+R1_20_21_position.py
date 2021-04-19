#!/usr/bin/python
from Bio import SeqIO
from Bio import motifs
from Bio.Seq import reverse_complement
import re
import os

Seq_fastq_R1 = "2_S2_L001_R1_001.fastq"                         # - *.fastq files to analyze
Seq_fastq_R2 = "2_S2_L001_R2_001.fastq"                         #
                                                                #
Flank5 = "CTTCCGATCT"                                           # - 10 nt sequence from adapter
Flank3 = "TGTCCTCTTC"                                           # - 10 nt sequence from plasmid backbone
                                                                #
TAM = "[ATGC]{7}"                                               # - number of N in the randomized sequence (7)
TAM_len = 7                                                     #
                                                                #
Target = "[ATGC]{20,21}"                                        # 
start = 20                                                      # - region in the target to extract TAM sequences (20-21)
stop = 21                                                       #
                                                                #
data_folder = "data"                                            # - create "data" folder in the current directory for the output file
current_directory = os.getcwd()                                 #
data = os.path.join(current_directory,data_folder)              # 
if not os.path.exists(data):                                    #
   os.makedirs(data)                                            #
os.chdir(data)                                                  #
TAMs = open("TAMs_for_WebLogo.txt", "w")                        # - output file with extracted TAM sequences (for visualization using WebLogo)
os.chdir("../")                                                 #


for record in SeqIO.parse(Seq_fastq_R1, "fastq"):                                
    m = str(record.seq)                                                                                                                    
    if re.search(Flank5 + Target + TAM + Flank3, m):                        
        n = re.search(Flank5 + Target + TAM + Flank3, m)                    
        o = n.group(0)                                                                                                                   
        tam = o[-(TAM_len + len(Flank3)):-(len(Flank3))]                                               
        tam_rc = reverse_complement(tam)                                                                                   
        TAMs.write(str(tam_rc) + '\n')                                      

for record in SeqIO.parse(Seq_fastq_R2, "fastq"):                                
    m = str(record.seq)                                                                                                                   
    if re.search(Flank5 + Target + TAM + Flank3, m):                        
        n = re.search(Flank5 + Target + TAM + Flank3, m)                    
        o = n.group(0)                                                                                                
        tam = o[-(TAM_len + len(Flank3)):-(len(Flank3))]                                            
        tam_rc = reverse_complement(tam)                                    
        TAMs.write(str(tam_rc) + '\n')                                      

TAMs.close()
                                                                             





