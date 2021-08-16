import pandas as pd
import argparse
import os
from hotSpot import out2csv
from hotSpot import cov
from hotSpot import lastMerge
from hotSpot import sim
from hotSpot import slideWindow

parser = argparse.ArgumentParser('hotSpot')

parser.add_argument("--window",  
                    help="Size of window (9/30/... anything but should be integer)",
                    default = 9,
                    type = int
                   )

parser.add_argument("--hla",
                    help="Class I or II",
                    choices = ["I","II"],
                    type = str,
                    required = True
                   )

parser.add_argument("--proteome",  
                    help="Proteome fasta file",
                    type = str,
                    required = True
                   )

#parser.add_argument("--protein",  help="Name of Protein (as give in proteome fasta file)")
parser.add_argument("--protein_file",  
                    help="Tab seperated files for: Predivacl_outfile predivac_outfile_inititals proteinNAmeAsGivenInProteome",
                    type = str,
                    required = True
                   )

parser.add_argument("--population",  
                    help="Name of population (India)",
                    type=str,
                    required = True,
                    choices = ["Algeria","American Samoa","Amerindian","Arab","Argentina","Asian",
                               "Australia","Australian Aborigines","Austria",
                               "Austronesian","Belarus","Belgium","Berber","Black","Bolivia",
                               "Borneo","Brazil","Bulgaria","Burkina Faso",
                               "Cameroon","Canada","Cape Verde","Caucasoid","Central Africa",
                               "Central African Republic","Central America",
                               "Chile","China","Colombia","Congo","Cook Islands","Costa Rica",
                               "Croatia","Cuba","Czech Republic","Denmark",
                               "East Africa","East Asia","Ecuador","England","Equatorial Guinea",
                               "Ethiopia","Europe","Fiji","Finland",
                               "France","Gabon","Gambia","Georgia","Germany","Ghana","Greece",
                               "Guatemala","Guinea-Bissau","Hispanic",
                               "Hong Kong","India","Indonesia","Iran","Ireland Northern",
                               "Ireland South","Israel","Italy","Ivory Coast",
                               "Jamaica","Japan","Jew","Jordan","Kenya","Kiribati","Kurd",
                               "Lebanon","Liberia","Macedonia","Malaysia","Mali",
                               "Martinique","Melanesian","Mestizo","Mexico","Micronesian",
                               "Mixed","Mongolia","Morocco","Mulatto","Nauru",
                               "Netherlands","New Caledonia","New Zealand","Nigeria","Niue",
                               "North Africa","North America","Northeast Asia",
                               "Norway","Oceania","Oman","Oriental","Other","Pakistan",
                               "Papua New Guinea","Paraguay","Persian","Peru","Philippines",
                               "Poland","Polynesian","Population Area","Population Country",
                               "Population Ethnicity","Portugal","Romania","Russia",
                               "Rwanda","Samoa","Sao Tome and Principe","Saudi Arabia",
                               "Scotland","Senegal","Serbia","Siberian","Singapore","Slovakia",
                               "Slovenia","South Africa","South America","South Asia",
                               "Southeast Asia","Southwest Asia","Spain","Sri Lanka","Sudan",
                               "Sweden","Switzerland","Taiwan","Thailand","Tokelau","Tonga",
                               "Trinidad and Tobago","Tunisia","Turkey","Uganda",
                               "Ukraine","United Arab Emirates","United Kingdom","United States",
                               "Venezuela","Vietnam","Wales","West Africa","West Indies","Zambia","Zimbabwe"]
                   )

parser.add_argument("--coverage",  
                    help="Coverage cut-off for simulation. Default is 0.6",
                    type=float,
                    default = 0.6
                   )

parser.add_argument("--quantile",  
                    help="quantile of HLA for simulation. Default is 0.9",
                    type=float,
                    default = 0.9
                   )

#parser.add_argument("--input",  help="Predivac output file name (generated in discovery mode) (xyz.txt)")

args = parser.parse_args()

window=args.window
hla_class=args.hla
proteome=args.proteome
#protein_name=args.protein
myfilename=args.protein_file
population_name=args.population
cov_simulation=args.coverage
quantile_hla=args.quantile

pd.set_option('display.max_rows', 
              None)     #The line is enough to display all rows from dataframe
pd.set_option('display.max_columns', 
              None)
pd.set_option('display.max_columns', 
              None) # The line is enough to display all columns from dataframe.  


def main ():
    print('')
    #_____________________
    if not os.path.exists(myfilename):
        print(colored("--protein_file: file does not exists",
                      "red",
                      attrs = ["bold"]
                     )
             )
        exit()

    if not os.path.exists(proteome):
        print(colored("--proteome: file does not exists",
                      "red",
                      attrs = ["bold"]
                     )
             )
        exit()
    # Running all functions
    #_____________________        
    myfile=pd.read_csv(myfilename,
                       sep='\t',
                       header=None)
    
    o_exist = 0
    for i in range(0,myfile.shape[0]):
        predivacOutFile=myfile.iloc[i,0]
        if not os.path.exists (predivacOutFile):
            print(colored("predivac out file given in the --protein_file"+
                          " does not exist",
                          "red",
                          attrs = ["bold"]
                         )
                 )
            o_exist = o_exist + 1
    if o_exist > 0:
        exit()
        
    for i in range(0,myfile.shape[0]):
        
        predivacOutFile=myfile.iloc[i,0]
        
        protein_name=myfile.iloc[i,1]
        
        print('Hotspot calculation for Predivac Output files (discovery mode)')
        print('')
        
        out2csv.out2csv(predivacOutFile,
                        proteome,
                        protein_name,
                        window
                       )
        
        slideWindow.slidingwindow(protein_name,
                                  proteome,
                                  predivacOutFile,
                                  window)
        
        cov.mycoverage(predivacOutFile,
                       hla_class,
                       population_name)
        
        
    total_proteome = pd.DataFrame({'Protein' : [],
                                   'Start' : [],
                                   'End' : [],
                                   'Sequence' : [],
                                   'Epitopes' : [],
                                   'NumberOfEpitopes' : [],
                                   'MergedEpitopeSequence' : [],
                                   'HLAs' : [],
                                   'coverage' : []})
    
    for i in range(0,myfile.shape[0]):
         mycoverageFiles=myfile.iloc[i,0].replace('.out','.window2Coverage.txt')
         my = pd.read_csv(mycoverageFiles,sep= '\t')
         total_proteome=pd.concat([total_proteome,my])
    #_____________________
    
    totalWindowCoverage='Total_window_coverage.'+ hla_class + "-" + population_name + '.txt'
    
    total_proteome.to_csv(totalWindowCoverage,
                          sep='\t',
                          index=None)
    
    sim.mysimulation(totalWindowCoverage,
                     cov_simulation,
                     quantile_hla) 
    
    lastMerge.lastmerge(totalWindowCoverage,
                        hla_class,
                        population_name)


if __name__ == "__main__":
    main()
