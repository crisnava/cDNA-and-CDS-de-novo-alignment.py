#!/usr/bin/env python

inp_f = '/home/.../de_novo_RZE037.algs.tab' 
out_f = '/home/.../RZE037_outp.txt'
out_f2 = '/home/.../RZE037_cov_cDNA_outp.txt'
out_f3 = '/home/.../RZE037_cov_CDS_outp.txt'

import matplotlib.pyplot as plt

cdna_covs, cds_covs = [], []
fo = open(out_f, "w") 
fcov_cDNA = open(out_f2, "w")
fcov_CDS = open(out_f3, "w")
for line in open(inp_f):
    if line.startswith('#'):
        continue   
    align = line.split() 
   
    overlap = float(align[3])/float(align[5])
    if overlap < 0.3:
        continue
    
    target_header = align[1].split('|')
    ref_transcript = target_header[0]
    CDS_start, CDS_end = map(int, target_header[1:3])
    CDS_length = CDS_end-CDS_start
    
    # Aligned cDNA and CDS start
    align_cDNA_start = int(align[2])
    mis_CDS_start = align_cDNA_start - CDS_start
    #alg starts before CDS start
    if mis_CDS_start<0:
        mis_CDS_start=0
    
    # Aligned cDNA and CDS end
    align_cDNA_end = int(align[2]) + int(align[3])
    
    # Bases missed in 3'
    mis_cDNA_end = int(align[5]) - align_cDNA_end
    
    if align_cDNA_end > CDS_end:
        mis_CDS_end = 0
    else:
        mis_CDS_end = CDS_end - align_cDNA_end
                
    # Total bases missed in the alignment and coverage
    mis_cDNA_total = int(align_cDNA_start) + int(mis_cDNA_end)
    cov_cDNA = 1 - float(mis_cDNA_total)/int(align[5])
    
    mis_CDS_total = mis_CDS_start + mis_CDS_end
   
    if align_cDNA_start>CDS_end or align_cDNA_end<CDS_start: 
        cov_CDS = 0.0
    else:
        cov_CDS = 1 - float(mis_CDS_total)/int(CDS_length)
        
    cds_covs.append(cov_CDS)
    cdna_covs.append(cov_cDNA)
    
cumulative=-1
plt.hist(cds_covs, bins=100, label="CDS", cumulative=cumulative)
plt.hist(cdna_covs, bins=70, label="cDNA", cumulative=cumulative)
plt.title('Reversed cumulative histogram of de novo transcriptome', fontsize=15)
plt.xlabel('Fraction', fontsize=12,)
plt.ylabel('No. of transcripts', fontsize=12)
plt.legend()
plt.show()
