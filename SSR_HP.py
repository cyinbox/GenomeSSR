# -*- coding: utf-8 -*-
"""
Changchuan Yin
Department of Mathematics, Statistics, and Computer Science
University of Illinois at Chicago
USA
Email cyin1@uic.edu

Citation

Yin., C. (2022). Evolutionary trend of SARS-CoV-2 inferred by the homopolymeric nucleotide repeats. Computational and Mathematical Biophysics.

Last updated 07/09/2022
"""

from fuzzysearch import find_near_matches
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------------
def seqs_fasta(seqFasta):  
    seqs = []
    for record in SeqIO.parse(seqFasta, 'fasta'):
        idx = record.description
       
        seq = {}
        seq['seq'] = record.seq
        seq['seqId'] = idx
        
        seqs.append(seq)
 
    return seqs

#-----------------------------------------------------------------------------
def getATCG(seq): 
   num_A = seq.count("A") 
   num_T = seq.count("T") 
   num_C = seq.count ("C") 
   num_G = seq.count ("G") 
   
   return [num_A, num_T, num_G, num_C, len(seq)]

#-----------------------------------------------------------------------------
def removeConsecutives(data):
    """
    remove consecutives (these consectives are the subset of the longer mers)
    # Example
    # data = [802, 941, 1066, 13123, 13124,14965, 15954, 16242, 18314,18315, 18778, 21421, 23515, 28372, 28880, 29131]
    # For example GGGG: 13123, and 13124 are the same as GGGGG at 13123  
    """
    cons = []
    for i in range(0, len(data)):
       if (data[i] - data[i-1]) == 1:
           cons.append(data[i-1])
           cons.append(data[i])
    
    dataX = []
    for item in data:
      if item not in cons:
         dataX.append(item)
          
    return dataX   

#-----------------------------------------------------------------------------
def SSRFreqs(genomeSeq):
    """
    Parameters
    ----------
    genomeSeq : a genome sequence
    Example: virusSeq='ATCGGC'

    Returns
    -------
    ssrFreqs : dictionary of SSR
    
    at2cg : SSR metric     

    """
    ntMers = ['AAA','AAAA','AAAAA','AAAAAA','AAAAAAA','TTT','TTTT','TTTTT','TTTTTT','TTTTTTT','CCC','CCCC','CCCCC','CCCCCC','CCCCCCC','GGG','GGGG','GGGGG','GGGGGG','GGGGGGG']
    maxDist = 0
    mersDict = {}
    
    L = len(genomeSeq)
    mersList = []
    
    for seqQ in ntMers:
       posStart = []  
       matches = find_near_matches(seqQ, genomeSeq, max_l_dist=maxDist)
       
       for match in matches:
           start = match.start
           posStart.append(start)
           
       mersDict[seqQ] = posStart
       mer = len(seqQ)
       
       R_pos = L*[0]
       for i in posStart:
          R_pos[i] = mer
       
       mersList.append(R_pos)
       
    positionData = []
    freqs = []
    ssrFreqs = {}
    
    for mers,posStart in mersDict.items():
        posStart = removeConsecutives(posStart)
        positionData.append(posStart)
        
        freq = len(posStart)
        freqs.append(freq)
        ssrFreqs[mers] = freq

    num = 0
    denum = 0
    for key,value in ssrFreqs.items():
        if 'AAA' in key or 'TTT' in key:
            num = num+value
        else:
            denum = denum+value
    hp = round(num/denum,4)
    
    return ssrFreqs, hp

#-----------------------------------------------------------------------------
def getPolyNTFreqs(virusSeq):
    ntMers = ['AAA','AAAA','AAAAA','AAAAAA','TTT','TTTT','TTTTT','TTTTTT',\
              'CCC','CCCC','CCCCC','CCCCCC','GGG','GGGG','GGGGG','GGGGGG']
    maxDist = 0 
    mersDict = {}
    
    L = len(virusSeq)
    
    mersList = []
    for seqQ in ntMers:
       posStart = []  
       matches = find_near_matches(seqQ, virusSeq, max_l_dist=maxDist)
       
       for match in matches:
           start = match.start
           posStart.append(start)
           
       mersDict[seqQ] = posStart
       mer = len(seqQ)
       
       R_pos = L*[0]
       for i in posStart:
          R_pos[i] = mer
        
       mersList.append(R_pos)
       
    positionData = []
    freqs = []
    for mers,posStart in mersDict.items():
        posStart = removeConsecutives(posStart)
        positionData.append(posStart)
        
        freq = len(posStart)*1000/L
        freqs.append(freq)
 
    s = sum(freqs[8:16])
    if s > 0:
       R_HP = round(sum(freqs[0:8])/sum(freqs[8:16]), 4) 
       
    else:
       R_HP = 0
       print('outlier')

    return [freqs,R_HP,mersDict]

#------------------------------------------------------------------------------
def eventPlot(idVirus, xLabel, yLabels, positionData, rowColors=None, lineWidth=2):
    m = len(positionData)   

    positionData = np.array(positionData)
    
    if rowColors is None:
       rowColors = ['C{}'.format(i) for i in range(m)]
    
    lineSize = m*[0.55]                                 
    linewidths1 = m*[lineWidth]
    
    plt.figure(figsize=(10,8))
    plt.eventplot(positionData, color=rowColors, linelengths=lineSize, linewidths=linewidths1)    
    plt.xticks(fontsize=14, fontweight='bold')
    plt.yticks(fontsize=12, fontweight='bold')
    plt.yticks(np.arange(m), tuple(yLabels), rotation=45)
    plt.xlabel(xLabel, fontsize=14, fontweight='bold')
    plt.show()

#------------------------------------------------------------------------------
def plotPolyNT(seqName, mersDict, freqs):
    ntMers = ['AAA','AAAA','AAAAA','AAAAAA','TTT','TTTT','TTTTT','TTTTTT','CCC','CCCC','CCCCC','CCCCCC','GGG','GGGG','GGGGG','GGGGGG']
 
    positionData = []
    for mers,posStart in mersDict.items():

        posStart = removeConsecutives(posStart)
        positionData.append(posStart)
    
    # Plotting the distribution of the poly NTs  
    xLabel = 'Nucleotide position (' + seqName +')'
    rowColors = ['C0','C0','C0','C0','C2','C2','C2','C2','C1','C1','C1','C1','C4','C4','C4','C4']
    eventPlot(seqName, xLabel, ntMers, positionData, rowColors=rowColors, lineWidth=2)
    
    NTs = ['A', 'T', 'C', 'G']
    colors = ['C0', 'C2', 'C1', 'C4']
    freuencies = []
    
    freuencies.append(freqs[0:4])
    freuencies.append(freqs[4:8])
    freuencies.append(freqs[8:12])
    freuencies.append(freqs[12:16])
    
    number_groups = len(NTs) 
    bin_width = 0.50/(number_groups+1)
    
    fig, ax = plt.subplots(figsize=(10,8))
    
    for i in range(number_groups):
        ax.bar(x=np.arange(len(NTs)) + i*bin_width, 
               height=freuencies[i],
               width=bin_width,
               color=colors[i],alpha=0.90,
               align='center')
    group = ['3A 3T 3C 3G', '4A 4T 4C 4G', '5A 5T 5C 5G', '6A 6T 6C 6G']
    
    ax.set_xticks(np.arange(len(NTs)) + 1/(2*(number_groups+1)))
    ax.set_xticklabels(group)
    ax.legend(NTs, facecolor='w')
    
    plt.xlabel(seqName, fontsize=12, fontweight='bold')
    plt.ylabel('Frequencies in 1 kb', fontsize=12, fontweight='bold')
    plt.ylim(0, 20)
    plt.show()

#------------------------------------------------------------------------------
if __name__ == "__main__":
    
    seqName = 'SARS-CoV-2'
    genomeFasta = 'NC_045512.fasta'
    seqs = seqs_fasta(genomeFasta)
    seq = seqs[0]['seq']
  
    ssrFreqs, hp = SSRFreqs(seq)
    print('Frequency counts', ssrFreqs, hp)
    
    [nA, nT, nG, nC, L] = getATCG(seq)
    at2cg = round((nA+nT)/(nC+nG),4)
    print(at2cg)
    
    [nfreqs, R_HP, mersDict] = getPolyNTFreqs(seq)
    plotPolyNT(seqName, mersDict, nfreqs)
    
    print('The normalized frequencies:', nfreqs)