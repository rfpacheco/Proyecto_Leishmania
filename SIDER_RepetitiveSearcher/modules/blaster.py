import os


# 1) Creacion del diccionario Blast

# Creation of the Blast Dic of our whole genome
def BlastN_Dic (Path_Input):
    try:
        os.system("makeblastdb -in " + Path_Input + " -dbtype nucl -parse_seqids")
        print("\nBlast Dictionary created in", Path_Input)
    except:
        print("\nError: Blast Dictionary couldn't be created")

# BlastN_Dic(Path_Input)

    #Arg 0: STRING. Directorio del archivo FASTA del cual queremos realizar el diccionario BLAST. Los archivos generados se colocaran en ese misma directorio, por ello se recomienda que este dentro de una carpeta unicamente dedicada a estos archivos.