
import argparse


parser = argparse.ArgumentParser(description = 'Filter host spectra from mgf or mzML file')


parser.add_argument('--Pout', type = str, required = True, help = 'path to PeptideShaker simple psm report file')
parser.add_argument('--out',type = str, required = True, help= 'output filepath')


args = parser.parse_args()




def filterPout(Poutfile,output):

    '''
    filter peptides with 'Cont_'  in the pout file
    input: list of spectrum titles, .mzML file of spectra
    output: filtered .mzML file
    '''
    LinesToWrite = []
    Spectra = open(Poutfile,'r')
    cleanedSpectra = open(output,'w')


    
    
    for line in Spectra:

        writeLine = True
            
        if line.find('Cont_') != -1:
            writeLine = False
         

        else:
            LinesToWrite.append(line)      
                              

    cleanedSpectra.writelines(LinesToWrite)

    Spectra.close()
    cleanedSpectra.close()




filterPout(args.Pout,args.out)



