import pandas as pd
from Bio import SeqIO
import os
from termcolor import colored

def out2csv (filename,proteome,protein_name,window):
    
    print(colored('-> Finding location of Epitopes in the proteome:' + protein_name,
                  "green",
                  attrs = ["bold"]
                 )
         )
    
    hla_number=[]
    data=pd.read_csv(filename)
    epitope=data.iloc[:,0].values.tolist()
    #	print(epitope)
    coverage=data.iloc[:,1].values.tolist()
    hla=data.iloc[:,2].values.tolist()
    for i in range(0,len(hla)):
        hla_number.append(len(hla[i].split( ',')))
        hla[i]=hla[i].replace(',',';')
        
    #--- Identifying start and end of epitopes 
    start=[]
    end=[]
    record_found = 0
    for i in range(0,len(epitope)):
        for record in SeqIO.parse(open(proteome,"r"),"fasta"):
            if record.id == protein_name:
                record_found = record_found + 1
                values= str(record.seq).lower().find(epitope[i].lower())
                start_value = values + 1
                end_value = start_value + (window - 1)
                start.append(start_value)
                end.append(end_value)
    if record_found == 0:
        print(colored("The name of the proteins is not found in the given Fasta file",
                      "red",
                      attrs = ["bold"]
                     )
             )
                
    frame = pd.DataFrame({'epitope':epitope,
                          'coverage':coverage,
                          'hla_number':hla_number,
                          'start': start, 
                          'end':end,
                          'hlas': hla})
    
    outputname=filename.replace('.out','.csv')
    frame.to_csv(outputname,index=False)
