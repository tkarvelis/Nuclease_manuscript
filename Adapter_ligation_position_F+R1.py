#!/usr/bin/python
from Bio import SeqIO
import re
import os

Seq_fastq_R1 = "XXX_R1.fastq"                         # - *.fastq files to analyze
Seq_fastq_R2 = "XXX_R2.fastq"                         #
                                                                #
Flank5 = "CTTCCGATCT"                                           # - 10 nt sequence from adapter
Flank3 = "TGTCCTCTTC"                                           # - 10 nt sequence from randomized region
                                                                #
TAM = "[ATGC]{7}"                                               # - number of N in the randomized sequence (7)
TAM_len = 7                                                     #
                                                                #
Target = "[ATGC]{0,30}"                                         # 
start = 0                                                       # - region in the target to test for adapter ligation (0-30)
stop = 30                                                       #
                                                                #
data_folder = "data"                                            # - create "data" folder in the current directory for the output files
current_directory = os.getcwd()                                 #
data = os.path.join(current_directory,data_folder)              # 
if not os.path.exists(data):                                    #
   os.makedirs(data)                                            #
os.chdir(data)                                                  #
Positions = open("Positions_F+R1.txt", "w")                     # - output file with counted cleavage positions
Statistics = open("Statistics_F+R1.txt", "w")                   # - output file with statistics
os.chdir("../")                                                 #

r = 0
f = 0
positions = []

for record in SeqIO.parse(Seq_fastq_R1, "fastq"):                                
    m = str(record.seq)                                                     
    r = r + 1                                                               
    if re.search(Flank5 + Target + TAM + Flank3, m):                        
        n = re.search(Flank5 + Target + TAM + Flank3, m)                    
        o = n.group(0)                                                       
        f = f + 1                                                           
        spacer = o[len(Flank5):-(TAM_len + len(Flank3))]                                           
        spacer_len = len(spacer)                                                                     
        positions.append(spacer_len)                                                                   

for record in SeqIO.parse(Seq_fastq_R2, "fastq"):                                
    m = str(record.seq)                                                     
    r = r + 1                                                               
    if re.search(Flank5 + Target + TAM + Flank3, m):                        
        n = re.search(Flank5 + Target + TAM + Flank3, m)                    
        o = n.group(0)                                                       
        f = f + 1                                                           
        spacer = o[len(Flank5):-(TAM_len + len(Flank3))]                                           
        spacer_len = len(spacer)                                                                     
        positions.append(spacer_len)                                                                   

while start <= stop:                                                        
    cpc = positions.count(start)                                             
    Positions.write(str(start) + '\t' + str(cpc / f * 100) + '\n')          
    start = start + 1                                                       

Statistics.write("Total reads: " + str(r) + '\n')                                                       
Statistics.write("Reads with adapter ligated at the target: " + str(f) + " (" + (str(f / r * 100)) + ("%)") + '\n')   #  

Positions.close()
Statistics.close()





