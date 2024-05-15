###rules that run ms2rescore

def InputMod(name, unimod_accession, mass_shift, amino_acid, n_term, c_term):
    return '{"name":"'+name+'", "unimod_accession":'+str(unimod_accession)+', "mass_shift":'+str(mass_shift)+', "amino_acid":'+str(amino_acid)+', "n_term":'+str(n_term)+', "c_term":'+str(c_term)+'},'

rule createMS2RescoreConfig:
    input: ExperimentDir+'{spectrum_name}/Xtandem/tandem.output.xml'
    output:
        ExperimentDir+'{spectrum_name}/ms2rescore/config.json'
    run:
        f_out = open(output[0], "w")

        lines = [' {\"$schema\":\"./config_schema.json\",']

        lines.append('"ms2rescore":{')
        lines.append('"feature_generators":{')
        lines.append('"basic":{},')
        lines.append('"ms2pip":{')
        lines.append('"model":'+FragModel+',')
        lines.append('"ms2tolerance":'+str(XTANDEM_FMME))
        lines.append('}')
        lines.append('},')
        lines.append('"psm_id_pattern":"'+SpectrumTitlePattern+'",')
        lines.append('"spectrum_id_pattern": null ,')
        lines.append('"modification_mapping":{')
        lines.append('"+57":"Carbamidomethyl",')
        lines.append('"+15.99":"Oxidation",')
        lines.append('"+39.9945":"Pyro-Carbamidomethyl",')
        lines.append('"-17.0266":"Gln->pyro-Glu",')
        lines.append('"-18.0106":"Glu->pyro-Glu",')
        lines.append('"+42.0106":"Acetyl"')
        lines.append('},')
        lines.append('"fixed_modifications":{},')
        lines.append('"processes":16,')
        lines.append('"rename_to_usi": false,')
        lines.append('"fasta_file":null,')
        lines.append('"write_report": false,')
        lines.append('"rescoring_engine":{')
        lines.append('"percolator":{')
        lines.append('}},')
        lines.append('"psm_file_type" : "xtandem" ,')
        lines.append('"id_decoy_pattern": "DECOY"')
        lines.append('}}')

        f_out.writelines([line + "\n" for line in lines])
        f_out.close()


rule RunMS2Rescore:
    input:
        ExperimentDir+'{spectrum_name}/Xtandem/tandem.output.xml',
        ExperimentDir+'{spectrum_name}/ms2rescore/config.json',
        SpectraDir+'{spectrum_name}.mgf',
        InputDB
    conda: 'envs/graphenv.yml'
    log: ExperimentDir+'{spectrum_name}/ms2rescore/ms2rescore.log'
    params: OutputName = ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored'
    output: ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.psms.tsv'
    shell: 'cp config/config.yaml '+ResultsDir +' && ms2rescore -c {input[1]} -s {input[2]} -o {params.OutputName} -p {input[0]} -f {input[3]} '

rule FilterPouT:
    input:
        ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.psms.tsv'
    output:
        ExperimentDir+'{spectrum_name}/ms2rescore/rescored/rescored.psms.filtered.tsv'
    shell: 'python3 filterHostSpectra.py --PoutFile {input} --out {output}'



