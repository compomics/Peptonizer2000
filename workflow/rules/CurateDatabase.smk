rule AddCrap:
    input:
       DatabaseDirectory+'crap.fasta',
       DatabaseDirectory+ReferenceDBName+ '.fasta'

    output: DatabaseDirectory+ ReferenceDBName+'+crap.fasta'

    shell:'cat {input} > {output}'  

rule RemoveDuplicatesCrap:
    input: DatabaseDirectory+ ReferenceDBName+'+crap.fasta'
    output: DatabaseDirectory+ ReferenceDBName+'+crap_UNI.fasta'
    conda: 'envs/graphenv.yml'
    shell:  'seqkit rmdup -s {input} > {output}'

rule AddDecoysCrap:
     input: DatabaseDirectory+ ReferenceDBName+'+crap_UNI.fasta'
     output: 
          DatabaseDirectory+ ReferenceDBName+'+crap_UNI_concatenated_target_decoy.fasta',
          DatabaseDirectory+ReferenceDBName+'+crap_UNI_decoy.fasta'
     conda: 'envs/decoy.yml'
     shell: 'decoypyrat {input} -o {output[1]} -d DECOY && cat {input} {output[1]} > {output[0]}'

#rules to produce files necessary for searching after filtering host spectra or to search all spectra but whithout host or crap DB added
rule RemoveDuplicates:
     input: DatabaseDirectory+ ReferenceDBName+'.fasta'
     output: DatabaseDirectory+ReferenceDBName+'_UNI.fasta'
     conda: 'envs/graphenv.yml'
     shell:  'seqkit rmdup -s {input} > {output}'

rule AddDecoys:
     input: DatabaseDirectory+ReferenceDBName+'_UNI.fasta'
     output: 
          DatabaseDirectory+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta',
          DatabaseDirectory+ReferenceDBName+'_UNI_decoy.fasta'
     conda: 'envs/decoy.yml'
     shell: 'decoypyrat {input} -o {output[1]} -d DECOY && cat {input} {output[1]} > {output[0]}'
          
#if the spectra aren't beeing filtered, check whether host and crap should be added to the search DB
def DBToUse(condition):
     if condition:
        return DatabaseDirectory+ReferenceDBName+'+crap_UNI_concatenated_target_decoy.fasta'
     else:
        return DatabaseDirectory+ReferenceDBName+'_UNI_concatenated_target_decoy.fasta'

InputDB = DBToUse(AddCrapToDB)

