import pandas as pd
import json
import os
import shutil

#get all files named 'UnipeptResponse.json' in the directory and subdirectories
def GetUnipeptJsonFiles(directory):
    files_list = []
    for subdir, dirs, files in os.walk(directory):
        for file in files:
            if file == 'UnipeptPeptides.json':
                files_list.append(os.path.join(subdir, file))
    return files_list



def ExportPeptidesForUnipeptCli(input,out):
#load dictionary from json
    with open(input) as json_file:
        print(input)
        peptide_to_psms = json.load(json_file)
        #print one peptide per line in new txt file as many times as it is in the dictionary sub entry 'key'
        with open(out, 'w') as f:
            for key in peptide_to_psms.keys():
                for i in range(peptide_to_psms[key]['psms']):
                    f.write("%s\n" % key)



UnipeptResults = GetUnipeptJsonFiles('/home/tanja/Peptonizer2000/Peptonizer2000/results')

for path in UnipeptResults:
    ExportPeptidesForUnipeptCli(path,'/home/tanja/Peptonizer2000/Peptonizer2000/results/UnipeptPeptides/'+path.split('/')[-3]+'_UnipeptPeptides.txt')

