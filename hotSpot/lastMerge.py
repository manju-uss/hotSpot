import pandas as pd
import os
from termcolor import colored
from predivac.PopCoverage import PopulationCoverage
import predivac
from predivac import Resources

#-------------- remaining

#--- Finaly merging and calculating overall cover: Because overlapping hotspots were identified and I need to merge them together which are nearby
#--------------
def lastmerge (myfile,mhc_class,population):
    
    print(colored('-> merging hotspot',
                  "green",
                  attrs = ["bold"]
                 )
         )
    
    infile = myfile.replace('.txt','.Hotspot.txt')
    df=pd.read_csv(infile,sep='\t')
    df.sort_values(['Protein','Start'],inplace=True,ascending=True)
    proteins=df.loc[:,'Protein'].unique()
    allprotein=[]
    allstart = []
    allend = []
    allepitopes=[]
    allhlas = []
    mhc_class= [mhc_class]
    population = [population]
    pcal = PopulationCoverage(population=population,mhc_class=mhc_class)
    for prot in proteins:
        mydf=df[df['Protein']==prot].reset_index(drop=True)
        start = mydf.loc[0,'Start']
        end   = mydf.loc[0,'End']
        epitopes   = mydf.loc[0,'Epitopes'].split(',')
        hlas   =mydf.loc[0,'HLAs'].split(',')
        for i in range(1,mydf.shape[0]):
            inepitope=mydf.loc[i,'Epitopes'].split(',')
            inhlas=mydf.loc[i,'HLAs'].split(',')
            if mydf.loc[i,'Start'] <= mydf.loc[i-1,'End']:
                end = mydf.loc[i,  'End' ]
                epitopes=epitopes+inepitope
                hlas=hlas+inhlas
            elif mydf.loc[i,'Start'] > mydf.loc[i-1,'End']:
                epitopes=list(set(epitopes))
                hlas=list(set(hlas))
                allprotein.append(prot)
                allstart.append(start)
                allend.append(end)
                allepitopes.append(epitopes)
                allhlas.append(hlas)
                start = mydf.loc[i,'Start']
                end = mydf.loc[i,'End']
                epitopes   = mydf.loc[i,'Epitopes'].split(',')
                hlas   =mydf.loc[i,'HLAs'].split(',')
        epitopes=list(set(epitopes))
        hlas=list(set(hlas))
        allprotein.append(prot)
        allstart.append(start)
        allend.append(end)
        allepitopes.append(epitopes)
        allhlas.append(hlas)
    total = {'protein': allprotein, 'start': allstart, 'end': allend, 'epitopes': allepitopes, 'hlas': allhlas}
    p_total=pd.DataFrame(total)
    p_total['hlas']=p_total['hlas'].astype(str).str[1:-1].str.replace(r"[\']", r"").str.replace(" ", "")
    p_total['epitopes']=p_total['epitopes'].astype(str).str[1:-1].str.replace(r"[\']", r"").str.replace(" ", "")
    p_total.loc[:,['hlas']].to_csv('temp.txt',sep='\t',header=None)
    #---- 
    mymy=[]
    mycoverage=pcal.calculate_coverage(filename="temp.txt")
    mylist=mycoverage[1].values.tolist()
    for j in mylist:
        mymy.append(float(str(j).replace('[','').replace(']','')))
    p_total['coverage']=mymy
    outfile = myfile.replace('.txt','.Merged.Hotspot.txt')
    p_total.to_csv(outfile,sep='\t',index=False)
