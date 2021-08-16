import pandas as pd
import os
from termcolor import colored
import progressbar
from scipy import stats
import numpy
import matplotlib.pyplot as plt
import random
#--------------
#--- simulation
#--------------
def mysimulation (myfile, coverage, quantile_Value):
    
    print(colored('-> Performing Simulation to find cut-off for hotspot',
                  "green",
                  attrs = ["bold"]
                 )
         )
    
    df = pd.read_csv(myfile,sep="\t")
    Q = df.loc[:,'NumberOfEpitopes'].quantile(quantile_Value)
#    print(str(Q )+ ":" + str(coverage))
    filtr = df.loc[(df['NumberOfEpitopes'] >=Q ) & (df['coverage'] >= coverage)] 
    myrows=filtr.shape[0] - 1
    print(myrows)
    
    if myrows < 10:
        print(colored("After filtering < 10 regios were selected for simulation. Decrease the coverage or quantile value",
                      "red",
                      attrs = ["bold"]
                     )
             )
        exit()
        
    else:
        myselectionOfrows=int(myrows*0.1)
        mean_value=[]
        Coverage_mean_value=[]
        bar = progressbar.ProgressBar(maxval=10000, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
        bar.start()
        for i in range (0,10000):
            bar.update(i + 1)
            bar.finish()
            num1 = []
            for i in range(0,myselectionOfrows):
                num1.append(random.randint(0,myrows))
            value = filtr.iloc[num1,:]
            mean_value.append(value.loc[:,'NumberOfEpitopes'].mean())
            Coverage_mean_value.append(value.loc[:,'coverage'].mean())
        print ('')
        data = mean_value
        density = stats.kde.gaussian_kde(data)
        x = numpy.arange(3, 5, .1)
        plt.plot(x, density(x))
        plt.show()
        #plt.save('out.png')
        mynumpy=numpy.array(mean_value)
        mycutoff_simulation=numpy.percentile(mynumpy,0.975)
        mynumpy_coverage=numpy.array(Coverage_mean_value)
        mycutoff_simulation_coverage=numpy.percentile(mynumpy_coverage,0.975)
        print ('My cut-off values are: ')
        print('Coverage: ' + str(mycutoff_simulation_coverage))
        print('Number of the Epitopes: ' + str(mycutoff_simulation))
        final_filtr = filtr.loc[(filtr['NumberOfEpitopes'] > mycutoff_simulation ) & (filtr['coverage'] > mycutoff_simulation_coverage )]
        outfile = myfile.replace('.txt','.Hotspot.txt')
        final_filtr.to_csv(outfile,sep='\t',index=False)
    
