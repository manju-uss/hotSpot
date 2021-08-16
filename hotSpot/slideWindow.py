import pandas as pd
from Bio import SeqIO
import os
from termcolor import colored
import progressbar

#---------------------------------------------------------------
# Running sliding window of 9 amino-acid and finding overlapping epitopes present in above CSV file 
# Format of file will be: proteinName Start End WindowSequence Epitopes HowManyEpitopes MergedEpitopeSequence HLAs
#---------------------------------------------------------------
def slidingwindow (protein_name, proteome, inputfile, windowSize):
    
    print(colored('-> Performing sliding window calculation:' + protein_name,
                  "green",
                  attrs = ["bold"]
                 )
         )
    
    proteinName=[] 
    Start=[]
    End=[]
    WindowSequence=[]
    Epitopes=[]
    NumberOfEpitopes=[]
    MergedEpitopeSequence=[]
    HLAs=[]
    csvfile=inputfile.replace('.out','.csv')
    mycsv=pd.read_csv(csvfile,header=0)
    mycsvrange=mycsv.shape[0]
    for record in SeqIO.parse(open(proteome,'r'),'fasta'):
        if record.id == protein_name:
            mylast = len(record.seq) - (windowSize - 1)
            bar = progressbar.ProgressBar(maxval=mylast, 
                                          widgets=[progressbar.Bar('=', 
                                                                   '[', ']'), 
                                                   ' ',
                                                   progressbar.Percentage()]
                                         )
            bar.start()
            for window in range(0,mylast):
                bar.update(window + 1)
                bar.finish()
                mhepitope=[]
                mhstart=[]
                mhend=[]
                mhalleles=[]
                for lines in range(0,mycsvrange):
                    if (
                        (mycsv.iloc[lines,3] < window 
                        and mycsv.iloc[lines,4] < window+(windowSize-1)
                        and mycsv.iloc[lines,4] > window)
                    or (mycsv.iloc[lines,3] == window
                        and mycsv.iloc[lines,4] == window+(windowSize-1))
                    or (mycsv.iloc[lines,3] >= window
                        and mycsv.iloc[lines,4] <= window+(windowSize-1))
                    or (mycsv.iloc[lines,3] > window 
                        and mycsv.iloc[lines,4] > window+(windowSize-1) 
                        and mycsv.iloc[lines,3] < window+(windowSize-1))
                    ):
                        
                        mhepitope.append(mycsv.iloc[lines,0])
                        mhstart.append(mycsv.iloc[lines,3] - 1)
                        mhend.append(mycsv.iloc[lines,4])
                        mhalleles=mhalleles + mycsv.iloc[lines,5].replace("\'","" ).split(';')
                        
                if len(mhstart) > 0:
                    myminimum=min(mhstart)
                    mymaximum=max(mhend)
                    MergedEpitopeSequence.append(str(record.seq[myminimum:mymaximum]))
                    HLAs.append(','.join( repr(e) for e in mhalleles))
                    Epitopes.append(','.join( repr(e) for e in mhepitope))
                    proteinName.append(protein_name)
                    Start.append(window)
                    End.append(window+(windowSize-1))
                    WindowSequence.append(str(record.seq[window:window+(windowSize-1)]))
                    NumberOfEpitopes.append(len(mhepitope))
                    
        frame = pd.DataFrame({'Protein': proteinName,
                              'Start':Start,
                              'End':End,
                              'Sequence': WindowSequence,
                              'Epitopes': Epitopes,
                              'NumberOfEpitopes': NumberOfEpitopes,
                              'MergedEpitopeSequence':MergedEpitopeSequence, 
                              'HLAs': HLAs})
    print('')
    
    outfile=inputfile.replace('.out','.window')
    frame.to_csv(outfile,sep='\t',index=False)


