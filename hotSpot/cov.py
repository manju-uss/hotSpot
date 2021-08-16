import pandas as pd
import predivac
from predivac import Resources
from predivac.PopCoverage import PopulationCoverage
import progressbar
from termcolor import colored
#---------------------------------------------------------------
# coverage calculation
#-------------------------------------------------------------
def mycoverage (inputfile,hla,population):

    print(colored('-> Calculating coverage for sliding windows:',
                  "green",
                  attrs = ["bold"]
                 )
         )
    
    outfileSlidingWindow = inputfile.replace('.out','.window')
    outfile = inputfile.replace('.out','.window2Coverage.txt')
    formycoveragefile=pd.read_csv(outfileSlidingWindow,sep='\t')
    formycoveragefile.iloc[:,7] = formycoveragefile.iloc[:,7].str.replace(r"[\']", r"")
    mhc_class = hla.split('=')
    population = population.split(',')
    pcal = PopulationCoverage(population=population,mhc_class=mhc_class)
    mymy=[]
    bar = progressbar.ProgressBar(maxval=formycoveragefile.shape[0], 
                                  widgets=[progressbar.Bar('*', 
                                                           '[', ']'), ' ', 
                                           progressbar.Percentage()]
                                 )
    bar.start()
    if formycoveragefile.shape[0] >= 10:
        for myrows in range(10,formycoveragefile.shape[0],10):
            bar.update(myrows + 1)
            bar.finish()
            s = myrows - 10
            e = myrows
            formycoveragefile.iloc[list(range(s,e)),[3,7]].to_csv('temp.txt',sep='\t',index=None,header=None)
            mycoverage=pcal.calculate_coverage(filename="temp.txt")
            mylist=mycoverage[1].values.tolist()
            for j in mylist:
                mymy.append(float(str(j).replace('[','').replace(']','')))
        formycoveragefile.iloc[list(range(e,formycoveragefile.shape[0])),[3,7]].to_csv('temp.txt',sep='\t',index=None,header=None)
        mycoverage=pcal.calculate_coverage(filename="temp.txt")
        mylist=mycoverage[1].values.tolist()
        for j in mylist:
            mymy.append(float(str(j).replace('[','').replace(']','')))
    if formycoveragefile.shape[0] < 10:
        formycoveragefile.iloc[:,[3,7]].to_csv("temp.txt",sep="\t",index=None,header=None)
        mycoverage=pcal.calculate_coverage(filename="temp.txt")
        mylist=mycoverage[1].values.tolist()
        for j in mylist:
            mymy.append(float(str(j).replace('[','').replace(']','')))
    formycoveragefile['coverage']=mymy
    formycoveragefile.to_csv(outfile,sep='\t',index=False)
