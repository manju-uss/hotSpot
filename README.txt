Installation:
____________
cd hotSpot
pip install ./

Usage:
_____
hotSpot --window 9 --hla II --proteome check.fa --protein_file proteinFile.txt --population Chile --coverage 0.6 --quantile 0.9

Options:
_______
  --window WINDOW       Size of window (9/30/... anything but should be integer)
  --hla {I,II}          Class I or II
  --proteome PROTEOME   Proteome fasta file
  --protein_file PROTEIN_FILE
                        Tab seperated files for: Predivacl_outfile predivac_outfile_inititals proteinNAmeAsGivenInProteome
  --population {Algeria,American Samoa,Amerindian,Arab,Argentina,Asian,Australia,Australian Aborigines,Austria,Austronesian,Belarus,Belgium,Berber,Black,Bolivia,Borneo,Brazil,Bulgaria,Burkina Faso,Cameroon,Canada,Cape Verde,Caucasoid,Central Africa,Central African Republic,Central America,Chile,China,Colombia,Congo,Cook Islands,Costa Rica,Croatia,Cuba,Czech Republic,Denmark,East Africa,East Asia,Ecuador,England,Equatorial Guinea,Ethiopia,Europe,Fiji,Finland,France,Gabon,Gambia,Georgia,Germany,Ghana,Greece,Guatemala,Guinea-Bissau,Hispanic,Hong Kong,India,Indonesia,Iran,Ireland Northern,Ireland South,Israel,Italy,Ivory Coast,Jamaica,Japan,Jew,Jordan,Kenya,Kiribati,Kurd,Lebanon,Liberia,Macedonia,Malaysia,Mali,Martinique,Melanesian,Mestizo,Mexico,Micronesian,Mixed,Mongolia,Morocco,Mulatto,Nauru,Netherlands,New Caledonia,New Zealand,Nigeria,Niue,North Africa,North America,Northeast Asia,Norway,Oceania,Oman,Oriental,Other,Pakistan,Papua New Guinea,Paraguay,Persian,Peru,Philippines,Poland,Polynesian,Population Area,Population Country,Population Ethnicity,Portugal,Romania,Russia,Rwanda,Samoa,Sao Tome and Principe,Saudi Arabia,Scotland,Senegal,Serbia,Siberian,Singapore,Slovakia,Slovenia,South Africa,South America,South Asia,Southeast Asia,Southwest Asia,Spain,Sri Lanka,Sudan,Sweden,Switzerland,Taiwan,Thailand,Tokelau,Tonga,Trinidad and Tobago,Tunisia,Turkey,Uganda,Ukraine,United Arab Emirates,United Kingdom,United States,Venezuela,Vietnam,Wales,West Africa,West Indies,Zambia,Zimbabwe}
                        Name of population (India)
  --coverage COVERAGE   Coverage cut-off for simulation. Default is 0.6
  --quantile QUANTILE   quantile of HLA for simulation. Default is 0.9

Format of the input files:
_________________________

Please see the predivac.out file given in the test folder

Contact:
_______

k2007.manju@gmail.com

