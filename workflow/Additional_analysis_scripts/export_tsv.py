import re

outfilee = 'S03_sample.tsv'
input = '/home/tanja/Peptonizer2000/Peptonizer2000/results/s03/CAMPI_SIHUMIx/S03/ms2rescore/rescored/rescored.filtered.psms.tsv'
fdr_threshold = 0.01
def MS2RescoreOutParser(pout_file, fdr_threshold,outfile):
    '''
    Parses the ms2rescore pout file for peptides, psm numbers and peptide scores
    :param pout_file: str, path to pout file(s)
    :param fdr_threshold: float, fdr threshold below which psms are kept
    :param decoy_flag: str, can be emtpy string, decoy flag in pout file
    :return: dict, peptides:[score,#psms]
    '''

    with open(outfile, "w") as outfile:
        outfile.write("Peptide\t1-q\n")  # Write header
    

        with open(pout_file, "r") as f:
            next(f)  # skip header
            for line in f:
                # skip empty lines
                if line.rstrip() == "":
                    continue
                splitted_line = line.rstrip().split("\t")[0:8]
                    #assert len(splitted_line) >= 6, "Input file is wrongly formatted. Make sure that the input is a valid .pout file."
                peptide, psm_id, run, colelction, collection , score, q, pep = splitted_line
                if float(q) < fdr_threshold:
                    peptide = re.sub(r"\[.*?\]", "", peptide).split("/")[0]
                    one_minus_q = 1 - float(q)
                            
                        # Write peptide and (1-q) to the output file
                    outfile.write(f"{peptide}\t{one_minus_q:.6f}\n")

MS2RescoreOutParser(input, fdr_threshold,outfilee)