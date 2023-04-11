#Script made by R. F. Pacehco Hinojosa

import os
import csv
import subprocess
import shutil
from sys import exit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#________________________________________________________________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#---------------------------PARTE 1: CREACION Y FILTRACION--------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#________________________________________________________________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#__________________________________MODULO: CREATORS______________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#Aqui estan las funciones dedicadas a la creaciones de archivos, documentos, etc.

#1) Creacion de la carpeta donde se guardara todo
def Folder_Creator(Folder_Name):#Es una string
    if not os.path.exists(Folder_Name):
        os.mkdir(Folder_Name)
        print("\nDirectory", Folder_Name, "created")
    else:
        print("\nDirectory", Folder_Name, "already exists. Be carefule before continuing")

# Folder_Creator(Folder_Name)

    #Arg 0: STRING. Nombre de la carpeta que se creara en el directorio donde sea ejecutado el código

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#2.1) Creamos el script creador de archivos CSV, que usaran varios modulos
def CSV_Creator(Writing_Path_Input, Writing_Input):
    with open(Writing_Path_Input, "w") as OutCSV:
        writer = csv.writer(OutCSV)
        writer.writerows(Writing_Input)
        print("\nCSV:", Writing_Path_Input, "has been created.")

#CSV_Creator(Writing_Path_Input, Matrix_Length_Filter)

    #Arg 0: STRING. Directorio en la que quiero que el archivo CSV se cree, y recordar que termine con la extension .csv
    #Arg 1: Lo que queremos meterle a el archivo CSV, generalmente es una Matrix/Array 3D, no se escribe en STRING
    
#2.2) Mezcla dos archivos CSV y utiliza el CSV_Creator para formar uno que abarque los dos

def CSV_Mixer(Path_Input1, Path_Input2, Writing_Path_Input):
    CSV_Mixer_Matrix = []
    with open(Path_Input1, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            CSV_Mixer_Matrix.append(row)
    with open(Path_Input2, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            CSV_Mixer_Matrix.append(row)
            
    CSV_Creator(Writing_Path_Input, CSV_Mixer_Matrix)
    
#CSV_Mixer(Path_Input1, Path_Input2, Writing_Path_Input)

    #Arg 0: STRING. Directorio del primer archivo a juntar, ocuparia las primeras filas del archivo final
    #Arg 1: STRING. Directorio del segundo archivo a juntar, se escribiria destras de las filas del Arg 0 en el archivo final
    #Arg 2: STRING. Directorio del archivo de salida que contendra los 2 archivos anteriores

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#3) Creacion Fasta para el siguiente BLASTN
def Fasta_Creator(Path_Input, Fasta_Output_Path):
    Matrix_Fasta_Creator = []
    Numbering = 0
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            Numbering += +1
            rec = SeqRecord(
                    Seq(row[15]),
                    id = "Seq_" + str(Numbering) + "_" + row[1] + "_" + row[14],#Que tenga aqui el sentido es esencial para luego filtrarlos
                    description = "Leishmania infantum " + row[14]
            )
            Matrix_Fasta_Creator.append(rec)
            
    SeqIO.write(Matrix_Fasta_Creator, Fasta_Output_Path, "fasta")
    print("\nFasta created at:", Fasta_Output_Path)

# Fasta_Creator(Path_Input, Fasta_Output_Path)

    #Arg 0: STRING. Directorio del archivo CSV a leer de donde queremos extraer las secuencias FASTA
    #Arg 1: STRING. Directorio del archivo FASTA que contiene las secuencias del archivo CSV. Recordar terminar en la extension .fasta

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#__________________________________MODULO : BLASTER______________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#1) Creacion del diccionario Blast
def BlastN_Dic (Path_Input):##Es una string
    try:
        os.system("makeblastdb -in " + Path_Input + " -dbtype nucl -parse_seqids")
        print("\nBlast Dictionary created in", Path_Input)
    except:
        print("\nError: Blast Dictionary couldn't be created")

# BlastN_Dic(Path_Input)

    #Arg 0: STRING. Directorio del archivo FASTA del cual queremos realizar el diccionario BLAST. Los archivos generados se colocaran en ese misma directorio, por ello se recomienda que este dentro de una carpeta unicamente dedicada a estos archivos.

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#2) Ejecucion del programa Blaster
def BLASTN_Blaster(Query_Path, Dict_Path, OutFile_Path, Perc_Identity): #Todos los ponemos en string ya que el programa os.system los incorpora ya sin las comillas. Puede que se conveniente sustituir en un futuro el modulo de os por el  subprocess
    try:
        os.system("blastn -word_size 11 -query " 
                  + Query_Path + " -db " 
                  + Dict_Path + " -out " 
                  + OutFile_Path + " -perc_identity " 
                  + Perc_Identity + " -outfmt '10 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore sstrand sseq'") #Uso el outfmt 10 para que las separaciones sean por comas y no por tabulaciones o espacios
        print("\nBlaster succesful", OutFile_Path, "created.")
    except:
        print("\nError: Blaster couldn't be loaded, somthing happened")
        
# BLASTN_Blaster(Query_Path,
#                Dict_Path,
#                OutFile_Path,
#                Perc_Identity)

    #Arg 0: STRING. Directorio del archivo FASTA query que queremos lanzar en BLASTN a nuestro diccionario.
    #Arg 1: STRING. Directorio del Diccionario en formato FASTA, i.e., es el mismo directorio que para la creacion del diccionario. El programa ejecutara el Fasta y se apoyara en otros archivos con el mismo nombre pero con extensiones .nhr, .nin, .nog, . nsd, .nsi y .nsq generados al ejecutar el creador de diccionario.
    #Arg 2: STRING. Directorio en formato CSV (recordar escribir la extension .csv) en donde se guardara los resultados del BLASTN.
    #Arg 3: STRING. Numero en el que le indicamos el porcentaje de homologia que queremos para realziar el BLASTN.
    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#__________________________________MODULO: FILTERS_______________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
  
#1) Filtro pora la obtencion de las etiquetas/titulos para cada cromosoma del archivo, en el caso de Leishmania obtendriamos asi las etiquetas LinJ.01, LinJ.02, etc. El objetivo de esta funcion es poder nombrar correctamente los archivos ya que los lenguajes de programacion no suelen admitir en un formano de no string, numero que empiecen en 0, de esta forma podemos automatizar su correcto etiquetado sobre todo par los numeros del 01 al 09.
def Chromosome_IDs(Path_Input, Name):#Todo STRING
    Max_Chr = len(list(SeqIO.parse(Path_Input, "fasta"))) #Leemos el archivo fasta y obtenemos su numero de cromosomas.
    Chromosome_Number = []
    List = (list(range(1, Max_Chr + 1)))
    for Number in List:
        Number = str(Number)
        if len(Number) == 1:
            Chromosome_Number.append(Name + ".0" + Number)#Para que coincida con los CSV
        else:
            Chromosome_Number.append(Name + "." + Number)
    return(Chromosome_Number)

#Chromosome_IDs(Path_Input, Name)

    #Arg 0: STRING. Directorio del archivo FASTA a leer.
    #Arg 1: STRING. Nombre para dar a los resultados, en el caso de Leihsmania infantum, seria "LinJ"

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
     
#2)Version mejorada para utilizarlo segun la columna que se desee editar, sea por longitud de las secuencias como por porcentaje de homologia
def Filter_by_Column(Path_Input, Column, Size_Filter, Writing_Path_Input):#Todo STRING menos Size_Filter
    if Column == "length":
        Column = 3
    elif Column == "percent":
        Column = 2
    
    Matrix_Filter_by_Column = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",") #Recordar que antes al poner outfmt 10, ahora estan separados por comas.
        for row in reader:
            if Column == 3:
                if 1000 >= int(row[Column]) >= Size_Filter:###1000 por el problema de duplicaciones que me dijo Requena
                    Matrix_Filter_by_Column.append(row)
            elif Column == 2:
                if float(row[Column]) >= Size_Filter:#Necesario para pasar de STRING a FLOAT
                    Matrix_Filter_by_Column.append(row)
    CSV_Creator(Writing_Path_Input, Matrix_Filter_by_Column)

#Filter_by_Column(Path_Input, Column, Size_Filter, Writing_Path_Input)

    #Arg 0: STRING. Directorio del archivo en formato CSV al que queremos filtrar los datos.
    #Arg 1: STRING. Puede ser "length" o "percent", dependiendo de lo que quereamos filtrar.
    #Arg 2: INT. Numero que nos indica el minimo para filtrar dependiendo del Arg 1.
    #Arg 3: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, recordar poner la extension .csv
    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#3) Filtrado de los "-" en las secuencias
def Dash_Filter(Path_Input, Writing_Path_Input):# Todo son strings
    Matrix_Dash_Filter = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            row[15] = row[15].replace("-","")
            Matrix_Dash_Filter.append(row)
    CSV_Creator(Writing_Path_Input, Matrix_Dash_Filter)

# Dash_Filter(Path_Input, Writing_Path_Input)

    #Arg 0: STRING. Directorio del archivo en formato CSV al que queremos filtrar los datos. Tiene que estar en STRING.
    #Arg 1: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Tiene que estar en formato STRING.

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------   

#3) Eliminacion de duplicado. En este primer paso se define la funcion para crear la matrix de los resultados filtrados.

#3.1) El primer paso elimina los duplicados segun la hebra que se le haya impuesto en el segundo argumento
def Pre_Duplicate_Filter(Path_Input, DNA_sense, Max_Diff):#Todos STRING menos Max_Diff
    Location_Start = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if DNA_sense in row[14]:
                Location_Start.append(int(row[10]))
                
    Matrix_Filter = []#Aqui guardara las ROWS que sacaremos la final     
    Position_Global = []##!!!!!MUY IMPORTANTE¡¡¡¡Con esto le decimos que no repita localizaciones, es vital
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if DNA_sense in row[14] and int(row[10]) not in Position_Global:
            #Arriba le hemos dicho que revise las localizaciones en el list GLOBAL, para que no se repita. Se hace aquí porque Position_RecPlus va cambiando constantemente.
                Position_Rec = []
                for Position in Location_Start:
                    if abs(int(row[10]) - Position) < Max_Diff:
                        if Position not in Position_Rec:#Esta parte es para que en lo que dura una organizacion de Position_RecPlus, no se repitan ningun valor dentro de ella y por consiguiente, dentro del Global
                            Position_Rec.append(Position)
                            Position_Global.append(Position)
                            
                DNASeq_Filter = []
                with open(Path_Input, "r") as Main_File:#Tengo que abrirlo de nuevo para empezar en la primera row siempre
                    reader = csv.reader(Main_File, delimiter = ",")
                    for row in reader:
                        if int(row[10]) in Position_Rec:#Comprobamos si la localizacion de START esta dentro de Position_RecPlus de este momento -recordar que va cambiando-. Si esta dentro de estas localizaciones, entonces nos centramos en mirar duplicaciones solo dentro de estas localizaciones
                            if row[15] in DNASeq_Filter:
                                continue #Si la secuencia esta dentro de nuestra base de datos, el CONTINUE salta el resto del codigo y vuelve al anterior loop FOR, evitando asi añadir duplicados
                            else:#Al no estar la secuencia dentro de nuestra base de datos, se le añade a nuestra base de datos Y a la matrix final global.
                                DNASeq_Filter.append(row[15])
                                Matrix_Filter.append(row)                        
    return(Matrix_Filter)

#Pre_Duplicate_Filter(Path_Input, DNA_sense, Max_Diff)

    #Arg 0: STRING. Directorio del archivo en formato CSV al que queremos filtrar los datos.
    #Arg 1: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Recordar la extension .csv    
    #Arg 2: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE


#3.2)En esta segunda funcion, se llama la anterior para que lo junte en una matrix y luego, posteriormente se fusionan la matrix de minus con la de plus. Llamaremos al CSV_Creator para que nos de el siguiente CSV a editar. La funcion trabaja a nivel de un elemento dado, no funciona a nivel del genoma entero, para eso, mirar los filtros a nivel genomico.

def Duplicate_Filter(Path_Input, Max_Diff, Writing_Path_Input):#Todo STRING menos Max_Diff
    DNA_sense = ["plus", "minus"]
    Matrix_Main1 = Pre_Duplicate_Filter(Path_Input, DNA_sense[0], Max_Diff)
    Matrix_Main2 = Pre_Duplicate_Filter(Path_Input, DNA_sense[1], Max_Diff)
    Matrix_Main = Matrix_Main1 + Matrix_Main2
        
    CSV_Creator(Writing_Path_Input, Matrix_Main)
    
# Duplicate_Filter(Path_Input, Max_Diff, Writing_Path_Input)

    #Arg 0: STRING. Directorio del archivo en formato CSV al que queremos filtrar los datos.
    #Arg 1: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    #Arg 2: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Recordar la extension .csv
    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#5)Ahora nos toca analizar las posiciones solapantes. La funcion trabaja a nivel de un elemento dado, no funciona a nivel del genoma entero, para eso, mirar los filtros a nivel genomico.

#5.1)Aqui metemos todos los START y END para tenerlos localizados.
def Solap_Location_Filter(Path_Input):
    Pos_Plus_Start = []
    Pos_Plus_End = []
    Pos_Minus_Start = []
    Pos_Minus_End = []
    with open (Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if "plus" in row[14]:
                Pos_Plus_Start.append(int(row[10]))
                Pos_Plus_End.append(int(row[11]))
            else:
                Pos_Minus_Start.append(int(row[10]))
                Pos_Minus_End.append(int(row[11]))
            
    return(Pos_Plus_Start, Pos_Plus_End, Pos_Minus_Start, Pos_Minus_End)

#Solap_Location_Filter(Path_Input)
    #Arg 0: STRING. Directorio del documento .csv a filtrar

#5.2)Aqui Metemos los START y END dentro de grupos por proximidad de coordenadas. Usamos la funcion anterior.
def Solap_Location_Grouping(Path_Input, DNA_sense, Max_Diff):#Todo STRING menos Max_Diff
    #Number es: 0 para Plus:Start, 1 para Plus:End, 2 para Minus:Start, 3 para Minus:End
    Position_List = Solap_Location_Filter(Path_Input)
    if DNA_sense == "plus":
        Position_List_Main = [Position_List[0], Position_List[1]]
    elif DNA_sense == "minus":
        Position_List_Main = [Position_List[2], Position_List[3]]
        
    Matrix1 = []
    Matrix2 = []
    
    for List in Position_List_Main:
        if List == Position_List_Main[0]:
            Matrix = Matrix1
        elif List == Position_List_Main[1]:
            Matrix = Matrix2
        #De esta forma me aseguro que cada List vaya a una Matrix diferente y al final las junto
        
        for Position in List:
            Main_Statement = False
            for Group in Matrix:
                for Member in Group:
                    if abs(Member - Position) <= Max_Diff:
                        Group.append(Position)
                        Main_Statement = True
                        break#Romperia codigo y saldriamos del loop "for Member" para coontinuar a "if not Found_Matrix"
            if not Main_Statement:#Es lo mismo que si "If Found_Matrix == False". Si Group esta vacio entonces vendra aqui directamente a poner el numero en ([]), es decir, en un array 3D
                Matrix.append([Position])
                
    Matrix_Main = [Matrix1, Matrix2]
    return(Matrix_Main)

#Solap_Location_Grouping(Path_Input, DNA_sense, Max_Diff)
    #Arg 0: STRING. Directorio del documento .csv a filtrar
    #Arg 1: STRING. Sentido de la hebra, puede ser "plus" o "minus"
    #Arg 2: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE

#5.3)Aqui obtenemos los minimos y maximos dependiendo de la hebra que tengamos. Se usa la funcion anterior
def Solap_MinMax(Path_Input, Max_Diff):#Todo STRING menos Max_Diff
    Plus = Solap_Location_Grouping(Path_Input, "plus", Max_Diff)
    Minus = Solap_Location_Grouping(Path_Input, "minus", Max_Diff)
    
    Plus_min = [min(x) for x in Plus[0]]
    Plus_max = [max(x) for x in Plus[1]]
    Minus_max = [max(x) for x in Minus[0]]
    Minus_min = [min(x) for x in Minus[1]]
    
    return(Plus_min, Plus_max, Minus_max, Minus_min)

#Solap_MinMax(Path_Input, Max_Diff)
    #Arg 0: STRING. Directorio del documento .csv a filtrar
    #Arg 1: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE

#5.4)Aqui, por ultimo, usando las funciones anteriores, eliminamos los solapamientos y los juntamos en un nuevo CSV con CSV_creator
def Solap_Main(Path_Input, Max_Diff, Writing_Path_Input):#Todo STRING menos Max_Diff
    Solap_Main_Matrix = []
    
    MinMax = Solap_MinMax(Path_Input, Max_Diff)#Aqui metemos los minimos y maximos anteriores en n array 3D. [0] --> Plus_min | [1] --> Plus_max | [2] --> Minus_max | [3] --> Minus_min
    Plus_Start_Matrix = []
    Plus_End_Matrix = []
    Minus_Start_Matrix = []
    Minus_End_Matrix = []
    
    #En esta parte se busca las secuencias que tengan tanto el minimo como el maximo, de tal forma, se obtiene el mayor alineamiento.
    with open (Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if "plus" in row[14]:
                if int(row[10]) in MinMax[0] and int(row[11]) in MinMax[1]:
                    Solap_Main_Matrix.append(row)
                    Plus_Start_Matrix.append(int(row[10]))
                    Plus_End_Matrix.append(int(row[11]))
            elif "minus" in row[14]:
                if int(row[10]) in MinMax[2] and int(row[11]) in MinMax[3]:
                    Solap_Main_Matrix.append(row)
                    Minus_Start_Matrix.append(int(row[10]))
                    Minus_End_Matrix.append(int(row[11]))
                    
    #Ahora, por si hay solapantes teniendo uno el minimo y otra secuencia el maximo, necesitaremos las dos. Por eso, revisamos sobre las coordenadas de lo anterior (para no repetir) y, buscamos esos solapantes.
    with open (Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if "plus" in row[14]:
                if int(row[10]) not in Plus_Start_Matrix and int(row[10]) in MinMax[0]:
                    Solap_Main_Matrix.append(row)

                if int(row[11]) not in Plus_End_Matrix and int(row[11]) in MinMax[1]:
                    Solap_Main_Matrix.append(row)
                
            elif "minus" in row[14]:
                if int(row[10]) not in Minus_Start_Matrix and int(row[10]) in MinMax[2]:
                    Solap_Main_Matrix.append(row)

                if int(row[11]) not in Minus_End_Matrix and int(row[11]) in MinMax[3]:
                    Solap_Main_Matrix.append(row)
                    
    CSV_Creator(Writing_Path_Input, Solap_Main_Matrix)

# Solap_Main(Path_Input, Max_Diff, Writing_Path_Input)

    #Arg 0: STRING. Directorio del documento .csv a filtrar
    #Arg 1: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    #Arg 2: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Recordar la extension .csv

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#______________________________MODULO GLOBAL: FILTERS DUPLICADOS_________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#1) Este filtro de duplicado ya funcionan a nivel del genoma entero. Se revisan las secuencias de las filas una,y por el cromosoma, sentido de la hebra y agurpaciones, se eliminan las duplicaciones.

#1.1) El primero es el filtro para Duplicados, funciona igual que el anterior pero es para todo un genoma y necesita tanto el CSV como el fasta.
def Genome_Pre_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, DNA_sense, Max_Diff): #Todo STRING menos Max_Diff
    
    Matrix_All_Genome = []
    Chromosome_Number = Chromosome_IDs(Genome_Fasta, Naming_Short)
    
    for Chromosome in Chromosome_Number: #Aqui nos metemos dentro del cromosoma a buscar.
        Location_Start = [] #Almacenaremos los START de ese cromosoma
        with open(Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter = ",")
            for row in reader:
                if Chromosome in row[1]: #CHROMOSOME FILTER
                    if DNA_sense in row[14]: #DNA_sense hay que poner PLUS y MINUS
                        Location_Start.append(int(row[10]))
            
            Matrix_Filter = [] 
            Position_Global = []##!!!!!MUY IMPORTANTE¡¡¡¡Con esto le decimos que no repita localizaciones, es vital
            with open(Path_Input, "r") as Main_File:
                reader = csv.reader(Main_File, delimiter = ",")
                for row in reader:
                    if Chromosome in row[1]: #CHROMOSOME FILTER
                        if DNA_sense in row[14] and int(row[10]) not in Position_Global:
                        #Arriba le hemos dicho que revise las localizaciones en el list GLOBAL, para que no se repita. Se hace aquí porque Position_RecPlus va cambiando constantemente.
                            Position_Rec = []
                            for Position in Location_Start:
                                if abs(int(row[10]) - Position) < Max_Diff:
                                    if Position not in Position_Rec:#Esta parte es para que en lo que dura una organizacion de Position_RecPlus, no se repitan ningun valor dentro de ella y por consiguiente, dentro del Global
                                        Position_Rec.append(Position)
                                        Position_Global.append(Position)

                            DNASeq_Filter = []
                            with open(Path_Input, "r") as Main_File:#Tengo que abrirlo de nuevo para empezar en la primera row siempre
                                reader = csv.reader(Main_File, delimiter = ",")
                                for row in reader:
                                    if Chromosome in row[1]: #CHROMOSOME FILTER
                                        if int(row[10]) in Position_Rec:#Comprobamos si la localizacion de START esta dentro de Position_RecPlus de este momento -recordar que va cambiando-. Si esta dentro de estas localizaciones, entonces nos centramos en mirar duplicaciones solo dentro de estas localizaciones
                                            if row[15] in DNASeq_Filter:
                                                continue#Si la secuencia esta dentro de nuestra base de datos, el CONTINUE salta el resto del codigo y vuelve al anterior loop FOR, evitando asi añadir duplicados
                                            else:#Al no estar la secuencia dentro de nuestra base de datos, se le añade a nuestra base de datos Y a la matrix final global.
                                                DNASeq_Filter.append(row[15])
                                                Matrix_Filter.append(row)
        Matrix_All_Genome += Matrix_Filter
    return(Matrix_All_Genome)

#Genome_Pre_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, DNA_sense, Max_Diff)

    #Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el FASTA del genoma entero
    #Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 3: STRING. Sentido de la hebra, puede ser "plus" o "minus"
    #Arg 4: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE



#1.2) Este es el principal:
def Genome_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input):#Todo STRING menos Max_Diff
    DNA_sense = ["plus", "minus"]
    
    Matrix_Main1 = Genome_Pre_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, DNA_sense[0], Max_Diff)
    
    Matrix_Main2 = Genome_Pre_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, DNA_sense[1], Max_Diff)
    
    Matrix_Main = Matrix_Main1 + Matrix_Main2
        
    CSV_Creator(Writing_Path_Input, Matrix_Main)
    
#Genome_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input)

    #Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    #Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 3: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    #Arg 4: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, Recordar la extension .csv
    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#____________________________MODULO GLOBAL: FILTERS SOLAPAMIENTOS________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#1)En este caso hacemos lo mismo que antes con los solapamientos pero a escala de todo el genoma. El principal cambio importante y MUY importante es el de la edicion del array Chromosome_Rows, el cual va hasta la primera funcion de Genome_Solap_Location_Filter

#1.1)Aqui metemos todos los START y END para tenerlos localizados.
def Genome_Solap_Location_Filter(Chromosome_Rows):#Todo STRING menos Chromosome_Rows
    Pos_Plus_Start = []
    Pos_Plus_End = []
    Pos_Minus_Start = []
    Pos_Minus_End = []

    for row in Chromosome_Rows:
            if "plus" in row[14]:#Concretamente esta parte
                Pos_Plus_Start.append(int(row[10]))
                Pos_Plus_End.append(int(row[11]))
            else:
                Pos_Minus_Start.append(int(row[10]))
                Pos_Minus_End.append(int(row[11]))
            
    return(Pos_Plus_Start, Pos_Plus_End, Pos_Minus_Start, Pos_Minus_End)

#Genome_Solap_Location_Filter(Chromosome_Rows)
    #Arg 0: Array con los datos a filtrar

#1.2)Aqui Metemos los START y END dentro de grupos por proximidad de coordenadas. Usamos la funcion anterior, a la que le impone el Chromosome_Rows
def Genome_Solap_Location_Grouping(Chromosome_Rows, DNA_sense, Max_Diff):#Todo STRING menos Number, Max_Diff y Chromosome_Rows
    #Number es: 0 para Plus:Start, 1 para Plus:End, 2 para Minus:Start, 3 para Minus:End
    Position_List = Genome_Solap_Location_Filter(Chromosome_Rows)
    if DNA_sense == "plus":
        Position_List_Main = [Position_List[0], Position_List[1]]
    elif DNA_sense == "minus":
        Position_List_Main = [Position_List[2], Position_List[3]]
        
    Matrix1 = []
    Matrix2 = []
    
    for List in Position_List_Main:
        if List == Position_List_Main[0]:
            Matrix = Matrix1
        elif List == Position_List_Main[1]:
            Matrix = Matrix2
        #De esta forma me aseguro que cada List vaya a una Matrix diferente y al final las junto
        
        for Position in List:
            Main_Statement = False
            for Group in Matrix:
                for Member in Group:
                    if abs(Member - Position) <= Max_Diff:
                        Group.append(Position)
                        Main_Statement = True
                        break#Romperia codigo y saldriamos del loop "for Member" para coontinuar a "if not Found_Matrix"
            if not Main_Statement:#Es lo mismo que si "If Found_Matrix == False". Si Group esta vacio entonces vendra aqui directamente a poner el numero en ([]), es decir, en un array 3D
                Matrix.append([Position])
                
    Matrix_Main = [Matrix1, Matrix2]
    return(Matrix_Main)

#Genome_Solap_Location_Grouping(Chromosome_Rows, DNA_sense, Max_Diff)
    #Arg 0: Array a filtrar.
    #Arg 1: STRING. Sentido de la hebra, puede ser "plus" o "minus"
    #Arg 2: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE

#1.3)Aqui obtenemos los minimos y maximos dependiendo de la hebra que tengamos. Se usa la funcion anterior y le impone el Genome_Rows
def Genome_Solap_MinMax(Chromosome_Rows, Max_Diff):#Todo STRING menos Max_Diff
    Plus = Genome_Solap_Location_Grouping(Chromosome_Rows, "plus", Max_Diff)
    Minus = Genome_Solap_Location_Grouping(Chromosome_Rows, "minus", Max_Diff)
    
    Plus_min = [min(x) for x in Plus[0]]
    Plus_max = [max(x) for x in Plus[1]]
    Minus_max = [max(x) for x in Minus[0]]
    Minus_min = [min(x) for x in Minus[1]]
    
    return(Plus_min, Plus_max, Minus_max, Minus_min)

#Genome_Solap_MinMax(Chromosome_Rows, Max_Diff)
    #Arg 0: Array a filtrar.
    #Arg 1: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE

#1.4)Modulo para solapantes, es MUY importante, ya que con el se filtran los solapantes parciales que no han sido filtrados por el modulo de solapantes totales de Genome_Solap_Main. Los analiza por parejas.
def Genome_Solap_by_Pairs(Rows_to_Filter):#Su argumento es un ARRAY 3D
    Rows_Final = []
    for First, Second in zip(*[iter(Rows_to_Filter)]*2):
        Two_Sequence_Rec = []
        Two_Sequence_Rec.append(First)
        Two_Sequence_Rec.append(Second)

        Sequence_START = []
        Sequence_END = []
        Homology1 = []
        E_Value1 = []
        Bit_Score1 = []
        for Sequence in Two_Sequence_Rec:

            Sequence_START.append(Sequence[10])
            Sequence_END.append(Sequence[11])
            Homology1.append(float(Sequence[2]))#Esta y las dos de abajo estan puestas por intentar mantener estos datos, pero en realidad no haria nada de falta
            E_Value1.append(float(Sequence[12]))
            Bit_Score1.append(float(Sequence[13]))

        Homology2 = str(round((Homology1[0] + Homology1[1])/2, 3))
        E_Value2 = (E_Value1[0] + E_Value1[1])/2
        E_Value2 = str("{:.2e}".format(E_Value2))
        Bit_Score2 = str(round((Bit_Score1[0] + Bit_Score1[1])/2, 1))

        if "plus" in First[14] and abs(int(First[10]) - int(Second[10])) <= 1000:#Este numero es VITAL
            Min_Start = min(Sequence_START)
            Max_End = max(Sequence_END)
            Seq_Length = str(int(Max_End) - int(Min_Start) + 1)

            Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + First[1] + " -range " + Min_Start + "-" + Max_End + " -strand plus -outfmt %s", shell = True, universal_newlines=True)#MUY IMPORTANTE EL SUBPROCESS
            Seq = Seq.strip()#Eliminar EoL caracteres

            New_Row = [First[0], First[1], Homology2, Seq_Length, First[4], First[5], "", "", "", "", str(Min_Start), str(Max_End), E_Value2, Bit_Score2, First[14], Seq]

            Rows_Final.append(New_Row)

        elif "minus" in First[14] and abs(int(First[10]) - int(Second[10])) <= 1000:
            Max_Start = max(Sequence_START)
            Min_End = min(Sequence_END)
            Seq_Length = str(int(Max_Start) - int(Min_End) + 1)

            Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + First[1] + " -range " + Min_End + "-" + Max_Start + " -strand minus -outfmt %s", shell = True, universal_newlines=True)#MUY IMPORTANTE EL SUBPROCESS
            Seq = Seq.strip()#Eliminar EoL caracteres

            New_Row = [First[0], First[1], Homology2, Seq_Length, First[4], First[5], "", "", "", "", str(Max_Start), str(Min_End), E_Value2, Bit_Score2, First[14], Seq]

            Rows_Final.append(New_Row)
    return(Rows_Final)

#Genome_Solap_by_Pairs(Rows_to_Filter)

    #Arg 0: array que se analizara por parejas obtenido de Genome_Solap_Main

#1.5)Este es el principal --> utiliza todos los de arriba y le lleva Chromosome_Rows
def Genome_Solap_Main(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input):#Todo STRING menos Max_Diff
    
    Genome_Solap_Main_Matrix = []
    
    Chromosome_Number = Chromosome_IDs(Genome_Fasta, Naming_Short)
    
    for Chromosome in Chromosome_Number:
        
        Solap_Main_Matrix = []
        
        Chromosome_Rows = []#ESENCIAL PARA LA PRIMERA DEFINICION A LA QUE SE LLAMA. Tiene que ir antes de MinMax
        with open (Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter = ",")
            for row in reader:
                if Chromosome in row[1]:#Chromosome FILTER
                    Chromosome_Rows.append(row)
        
        MinMax = Genome_Solap_MinMax(Chromosome_Rows, Max_Diff)#Aqui metemos los minimos y maximos anteriores en n array 3D. [0] --> Plus_min | [1] --> Plus_max | [2] --> Minus_max | [3] --> Minus_min
        Plus_Start_Matrix = []
        Plus_End_Matrix = []
        Minus_Start_Matrix = []
        Minus_End_Matrix = []
    
        #En esta parte se busca las secuencias que tengan tanto el minimo como el maximo, de tal forma, se obtiene el mayor alineamiento.
        with open (Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter = ",")
            for row in reader:
                if Chromosome in row[1]:#Chromosome FILTER
                    if "plus" in row[14]:
                        if int(row[10]) in MinMax[0] and int(row[11]) in MinMax[1]:
                            Solap_Main_Matrix.append(row)
                            Plus_Start_Matrix.append(int(row[10]))
                            Plus_End_Matrix.append(int(row[11]))
                    elif "minus" in row[14]:
                        if int(row[10]) in MinMax[2] and int(row[11]) in MinMax[3]:
                            Solap_Main_Matrix.append(row)
                            Minus_Start_Matrix.append(int(row[10]))
                            Minus_End_Matrix.append(int(row[11]))
                                              
        #Ahora, por si hay solapantes teniendo uno el minimo y otra secuencia el maximo, necesitaremos las dos. Por eso, revisamos sobre las coordenadas de lo anterior (para no repetir) y, buscamos esos solapantes.
        Solap_Segments = [] #Aqui meteria todos los segmentos solapados pequeñitos
        Solap_Segments_Plus_Start = []
        Solap_Segments_Plus_End = []
        Solap_Segments_Minus_Start = []
        Solap_Segments_Minus_End = []
        with open (Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter = ",")
            for row in reader:
                if Chromosome in row[1]:
                    if "plus" in row[14]:
                        if int(row[10]) not in Plus_Start_Matrix and int(row[10]) in MinMax[0]:
                            if int(row[10]) not in Solap_Segments_Plus_Start:
                                Solap_Segments.append(row)
                                Solap_Segments_Plus_Start.append(int(row[10]))
                                
                        if int(row[11]) not in Plus_End_Matrix and int(row[11]) in MinMax[1]:
                            if int(row[11]) not in Solap_Segments_Plus_End:
                                Solap_Segments.append(row)
                                Solap_Segments_Plus_End.append(int(row[11]))

                    elif "minus" in row[14]:
                        if int(row[10]) not in Minus_Start_Matrix and int(row[10]) in MinMax[2]:
                            if int(row[10]) not in Solap_Segments_Minus_Start:
                                Solap_Segments.append(row)
                                Solap_Segments_Minus_Start.append(int(row[10]))

                        if int(row[11]) not in Minus_End_Matrix and int(row[11]) in MinMax[3]:
                            if int(row[11]) not in Solap_Segments_Minus_End:
                                Solap_Segments.append(row)
                                Solap_Segments_Minus_End.append(int(row[11]))
                                
        Solap_by_Pairs_Definitive = Genome_Solap_by_Pairs(Solap_Segments)#Esto es CLAVE
        Solap_Main_Matrix += Solap_by_Pairs_Definitive
        
        Genome_Solap_Main_Matrix += Solap_Main_Matrix
                
    CSV_Creator(Writing_Path_Input, Genome_Solap_Main_Matrix)
    
#Genome_Solap_Main(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input)

    #Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    #Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos.
    #Arg 3: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    #Arg 4: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Recordar la extension .csv
    
    
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#_______________________________________GLOBAL FILTERS: ALl______________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#1)Este junta todos los anteriores de la primera prueba, cada uno escribe un CSV y lo va sobreescribiendo constantemente hasta el final.
def Global_Filters_Main(Path_Input, Writing_Path_Input, Genome_Fasta, Naming_Short, Max_Diff):
    
    Column = "length"
    Size_Filter = 100
    Filter_by_Column(Path_Input, Column, Size_Filter, Writing_Path_Input)
    
    Path_Input = Writing_Path_Input #Así le decimos que el archivo de entrada es el de salida del anterior, y que en el mismo, escriba los nuevos datos
    Dash_Filter(Path_Input, Writing_Path_Input)
    
    Genome_Duplicate_Filter(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input)
    
    Genome_Solap_Main(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input)
    #En este ultima ya imprime los resultados finales
    
#Global_Filters_Main(Path_Input, Writing_Path_Input, Genome_Fasta, Naming_Short, Max_Diff):

    #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 1: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado, Recordar la extension .csv
    #Arg 2: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    #Arg 3: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 4: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#________________________________________________________________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#---------------------------PARTE 2: AUTOMATIZACION DEL BLASTER---------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#________________________________________________________________________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#___________________________________MODULO: IDENTIFICADOR_______________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#1) Utilizando un archivo CSV a leer, nos obtiene las secuencias especificas para X cromosoma y las mueve a un nuevo archivo CSV, dentro de una carpeta especifica para el cromosoma
def Specific_Sequence_extractor(Path_Input, Chromosome_ID, Main_Folder_Path):
    
    Chr_X_Seqs = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
              if Chromosome_ID in row[1]:
                    Chr_X_Seqs.append(row)
    
    Folder_Path = Main_Folder_Path + "/" + Chromosome_ID
    Folder_Creator(Folder_Path)
    
    Writing_Path_Input = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + ".csv"
    
    CSV_Creator(Writing_Path_Input, Chr_X_Seqs)
    
    return(Folder_Path, Writing_Path_Input)##Es importante porque asi nos devuelven los nuevos directorios. No se las puedo añadir a variables globales, porque lamentablemente Python no funciona asi

#Specific_Sequence_extractor(Path_Input, Chromosome_ID, Main_Folder_Path)
    #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa. Lo utilizara para generar una subcarpeta con el nombre del cromosoma

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#2) Utiliza un archivo CSV, del cual lee las secuencias y las extiende hasta los 1000 nt, creando un archivo CSV resultante con esos datos
def Specific_Sequence_1000nt(Path_Input, Chromosome_ID, Main_Folder_Path):
    ChrX_1000nt = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if "plus" in row[14]:
                Seq_Length = int(row[11]) - int(row[10]) #Es mejor hacerlo con esto que con la row[3], porque esta puede en pocos casos haberse quedado igual aunque el resto haya cambiado mediante blastdbcmd
                Number_Add_Length = int((1000 - Seq_Length)/2)#Con esto tenemos lo que tenemos que añadir a las coordenadas para que las secuencias tengan 1000nt
                New_START = int(row[10]) - Number_Add_Length
                New_END = int(row[11]) + Number_Add_Length

                Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + row[1] + " -range " + str(New_START) + "-" + str(New_END) + " -strand plus -outfmt %s", shell = True, universal_newlines=True)#MUY IMPORTANTE EL SUBPROCESS
                Seq = Seq.strip()#Eliminar EoL caracteres

                New_Row = [row[0], row[1], "", str(len(Seq)), row[4], row[5], "", "", "", "", str(New_START), str(New_END), "", "", row[14], Seq]

                ChrX_1000nt.append(New_Row)

            elif "minus" in row[14]:
                Seq_Length = int(row[10]) - int(row[11]) #Al reves al ser minus. Podria hacaerlo sino en absoluto, pero bueno, he elegido esta opcion

                Number_Add_Length = int((1000 - Seq_Length)/2)
                New_START = int(row[10]) + Number_Add_Length #Al reves al ser minus
                New_END = int(row[11]) - Number_Add_Length #Al reves al ser minus

                Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + row[1] + " -range " + str(New_END) + "-" + str(New_START) + " -strand minus -outfmt %s", shell = True, universal_newlines=True)#MUY IMPORTANTE EL SUBPROCESS
                Seq = Seq.strip()#Eliminar EoL caracteres

                New_Row = [row[0], row[1], "", str(len(Seq)), row[4], row[5], "", "", "", "", str(New_START), str(New_END), "", "", row[14], Seq]

                ChrX_1000nt.append(New_Row)
                
    Writing_Path_Input = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_1000nt.csv"
    CSV_Creator(Writing_Path_Input, ChrX_1000nt)
    
    return(Writing_Path_Input)###Importante para saber el directorio de este archivo creado

#def Specific_Sequence_1000nt(Path_Input, Chromosome_ID, Main_Folder_Path)

    #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#3) Corrector se secuencias, obtendra las secuencias originales.
def Specific_Sequence_Corrected(Path_Input, Nucleotides1000_Directory, Main_Folder_Path, Chromosome_ID):
    Names = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if row[0] not in Names:
                Names.append(row[0])

    ChrX_Corrected = []
    for Query in Names:
        Start = []
        End = []
        Diference_End_minu_Start = 0 ##Si resulta ser mayor, se le añadira el valor mayor
        with open(Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter = ",")              
            for row in reader:
                if Query in row[1]:#Mejor hacerlo con este para realizar lo de las tablas del CSV que hice
                    if row[0] != Query:#Lo  row[0] != Query es para eliminar la secuencia que solapa consigo misma. Por otro lado, no hace falta diferenciar entre plus y minus. Porque hemos lanzado nuestras pequeñas secuencias a un blast propio, y en si se comportan como si todas fueran "plus"
                        Diference = int(row[11]) - int(row[10])#Por como esta el codigo ahora 11 siempre sera mayor que 10
                        if Diference > Diference_End_minu_Start:
                            
                            Diference_End_minu_Start = Diference
                            Start = []
                            End = []
                            Start.append(int(row[10]))
                            End.append(int(row[11]))
        
        #De nuevo, no hace falta diferenciar entre plus y minus por lo mismo de antes
        
        #Esta parte ya no se necesita despues de haber indicado lo de Diference, pero en otro momento la quito
        if len(Start) > 0 and len(End) > 0:#Asi creo que evito el hecho de que Seq5 no tenga homologia con nada.
            Min_Start = min(Start)
            Max_End = max(End)


            Correct_Seq = Query

            Number_for_Location = int(Correct_Seq[4]) - 1 #Asi cogemos con INT la cuarta posicion de Seq_X, siendo X un numero y le restamos 1, porque ya sabemos que en Python todo empieza en 0
            #Ahora filtramos el Correct_Seq para obtener un numero para poder filtrar la lista de csv de 4 x 1000 sin hacer blaster.

            Rows_by_Number = []#Esta parte la necesito para poder saber luego si estoy en la row correcta al hacer comparaciones, ya que puedo compararlo con Rows_by_Number[0] o [4] o [3], sin ir en orden. Esto se hace en el csv de 4 x 1000nt antes del blaster a ellos mismos
            with open(Nucleotides1000_Directory, "r") as Main_File:
                reader = csv.reader(Main_File, delimiter = ",")              
                for row in reader:
                    Rows_by_Number.append(row)

            with open(Nucleotides1000_Directory, "r") as Main_File:
                reader = csv.reader(Main_File, delimiter = ",")              
                for row in reader:
                    if row == Rows_by_Number[Number_for_Location]: #Asi me aseguro que estoy en la adecuada. Quizas es rizar el rizo pero no se me ocurre en el momento un paso mejor
                        if "plus" in row[14]:

                            x = 1000 - Max_End
                            New_START = int(row[10]) + Min_Start - 1
                            New_END = int(row[11]) - x

                            Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + row[1] + " -range " + str(New_START) + "-" + str(New_END) + " -strand plus -outfmt %s", shell = True, universal_newlines=True)#MUY IMPORTANTE EL SUBPROCESS
                            Seq = Seq.strip()#Eliminar EoL caracteres

                            New_Row = [Query, row[1], "", str(len(Seq)), row[4], row[5], "", "", "", "", str(New_START), str(New_END), "", "", row[14], Seq]

                            ChrX_Corrected.append(New_Row)

                        elif "minus" in row[14]:#Recordar que por como son las coordenadas de minus, que se definen segun las de la posicion plus en 3' --> 5'
                            x = 1000 - Max_End
                            New_START = int(row[10]) - Min_Start + 1
                            New_END = int(row[11]) + x

                            Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry " + row[1] + " -range " + str(New_END) + "-" + str(New_START) + " -strand minus -outfmt %s", shell = True, universal_newlines=True)#MUY IMPORTANTE EL SUBPROCESS
                            Seq = Seq.strip()#Eliminar EoL caracteres

                            New_Row = [Query, row[1], "", str(len(Seq)), row[4], row[5], "", "", "", "", str(New_START), str(New_END), "", "", row[14], Seq]

                            ChrX_Corrected.append(New_Row)
                            
        if len(Start) == 0 and len(End) == 0:#para casos en los que solo tenga homologia con el mismo, la secuencia se descarta, pero seria mejor cambiar este codigo para insertarla en los siguientes documentos pero no en la forma de 1000nt
            print("\nALERT: individual " + Query + " has no homology with no other seq, so it will not be added to the corrected Seqs")

    Writing_Path_Input = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_Corrected.csv"
    CSV_Creator(Writing_Path_Input, ChrX_Corrected)
    return(Writing_Path_Input)

#Specific_Sequence_Corrected(Path_Input, Nucleotides1000_Directory, Main_Folder_Path, Chromosome_ID)

    #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 1: Resultado del RETURN de la funcion Specific_Sequence_1000nt
    #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa
    #Arg 3: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#4)Este apartado simplemente nos generara un archivo con las subfamilias para ver como se clasifican
def Subfamily_Sorter(Path_Input, Corrected_Elements_Path, Writing_Path_Input):
    Names = []
    with open(Path_Input, "r") as Main_File:###Path_Input seria el documento del Blaster
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            if row[0] not in Names:
                Names.append(row[0])
                
    For_Subfamilies = []
    with open(Path_Input, "r") as Main_File:###Path_Input seria el documento del Blaster
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            New_Row = [row[0], row[1], row[2], "", ""]#importante las posiciones "" vacias, porque ahi añadiremos las longitudes para hacer el calculo del siguiente filtrado por longitud
            For_Subfamilies.append(New_Row)
            
    For_Subfamilies2 = []
    with open(Corrected_Elements_Path, "r") as Main_File:###Aqui tendriamos el documento de los elementos ya corregidos.
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            New_Row = [row[0], row[3]]
            For_Subfamilies2.append(New_Row)
                
    for row in For_Subfamilies2:
        for row2 in For_Subfamilies:
            if row[0] in row2[0]:
                row2[3] = row[1]

            if row[0] in row2[1]:
                row2[4] = row[1]
                                              
    #Aqui vamos a cambiar el % a valores de 1.
    For_Subfamilies3 = []
    for row in For_Subfamilies:
        if row[3] != "" and row[4] != "":#Para eliminar los que no hacen homologia con ninguno y tienen por tanto los valores vacios
            Subfamily_Filter = ((float(row[2]) * int(row[3]))/int(row[4]))/100

            if 0.85 <= Subfamily_Filter <= 1.15:
                New_Row = [row[0], row[1]]
                For_Subfamilies3.append(New_Row)
            
        if row[3] == "" and row[4] == "":##Para el que no hace homologia con ninguno salvo consigo mismo, al menos que en el archivo de subffamilias, aparezca, pero solo.
            Subfamily_Filter = 1.0
            New_Row = [row[0], row[1]]
            For_Subfamilies3.append(New_Row)
            
    
    For_Subfamilies4 = []
    Rec_For_Subfamilies = []
    for Sequence in Names:
        for Pair in For_Subfamilies3:
            if Sequence in Pair[0]:
                if Pair[1] not in Rec_For_Subfamilies:
                    Rec_For_Subfamilies.append(Pair[1])

        For_Subfamilies4.append(sorted(Rec_For_Subfamilies))
        Rec_For_Subfamilies = []
           
    #Ahora nos encargaremos de eliminar los duplicados:
    For_Subfamilies5 = []
    for Groups in For_Subfamilies4:
        if Groups not in For_Subfamilies5:
            For_Subfamilies5.append(Groups)
            
    CSV_Creator(Writing_Path_Input, For_Subfamilies5)

#Subfamily_Sorter(Path_Input, Corrected_Elements_Path, Writing_Path_Input)
    
        #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
        #Arg 1: RETURN del resultado de la funcion Specific_Sequence_Corrected
        #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#5) Funcion principal que utiliza todas las funciones anteriores descritas en esta parte, incluyendo filtraciones globales y guardado de datos en un archivo MIXER.csv

#Chromosome_ID: hay que indicarle en STRING la etiqueta del cromosoma, e.g., LinJ.01. Pero esto lo tiene que hacer un Script que invoque a Genome_Specific_Chromosome_Main

def Genome_Specific_Chromosome_Main(Path_Input, Chromosome_ID, Main_Folder_Path, Genome_Fasta, Naming_Short, Max_Diff):
    
    New_Directories = Specific_Sequence_extractor(Path_Input, Chromosome_ID, Main_Folder_Path)
    Folder_Path = New_Directories[0]#El directorio de la carpeta del cromosoma
    Last_Output = New_Directories[1]#El directorio del comando anterior. Estos dos pasos los hago por facilitar la lectura del codigo
    
    Nucleotides1000_Directory = Specific_Sequence_1000nt(Last_Output, Chromosome_ID, Main_Folder_Path)
    
    Fasta_Creator_Output = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_1000nt.fasta"
    Fasta_Creator(Nucleotides1000_Directory, Fasta_Creator_Output)
    
    BlastN_Dic(Fasta_Creator_Output)
    
    Blaster_Output = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_1000nt_Blaster.csv"
    BLASTN_Blaster(Fasta_Creator_Output,
                   Fasta_Creator_Output,
                   Blaster_Output,
                   "85")

    Filter_by_Column(Blaster_Output,
                     "length",
                     100,
                     Blaster_Output)
    
    Corrected_Sequences = Specific_Sequence_Corrected(Blaster_Output, Nucleotides1000_Directory, Main_Folder_Path, Chromosome_ID)
    
    Subfamilies_File_Path_Writing = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_Subfamily.csv"
    Subfamily_Sorter(Blaster_Output, Corrected_Sequences, Subfamilies_File_Path_Writing)
    
    Second_Fasta_Creator_Output = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_Corrected.fasta"
    Fasta_Creator(Corrected_Sequences, Second_Fasta_Creator_Output)
    
    Second_Blaster_Output = Main_Folder_Path + "/" + Chromosome_ID + "/" + Chromosome_ID + "_BLAST_MAIN.csv"
    BLASTN_Blaster(Second_Fasta_Creator_Output,
                   Genome_Fasta,
                   Second_Blaster_Output,
                   "60")
    
    Global_Filters_Main(Second_Blaster_Output,
                        Second_Blaster_Output,
                        Genome_Fasta,
                        Naming_Short,
                        Max_Diff)
    
    
    CSV_Mixer_Output = Main_Folder_Path + "/" + "MIXER.csv"
    if os.path.isfile(CSV_Mixer_Output) == False:##Cuando no existe, se crea
        CSV_Mixer(Path_Input, Second_Blaster_Output, CSV_Mixer_Output)#Para mezclar
    else:#Si existe ya el archivo porque ha sido creado, se cambia el Path_Input por CSV_Mixer_Output
        CSV_Mixer(CSV_Mixer_Output, Second_Blaster_Output, CSV_Mixer_Output)

#Genome_Specific_Chromosome_Main(Path_Input, Chromosome_ID, Main_Folder_Path, Genome_Fasta, Naming_Short, Max_Diff)

    #Arg 0: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 1: STRING. Identificacion del cromosoma, e.g., "LinJ.07"
    #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa
    #Arg 4: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    #Arg 5: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 6: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
        
        
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#_____________________________MODULO GLOBAL: BLASTER AUTOMATICO__________________________________
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------

#1) Aqui definiciones la funcion principal de absolutamente todo el programa, la que realizara todos los analisis por cada cromosoma de forma automatica y guardara los resultados de forma automatica tambien. Y dependiendo de los valores de Numbering y Maximun_Runs, se le podra ejecutar de forma iterativa un numero especifico de veces

def Repetitive_Blaster(Genome_Fasta, Path_Input, Folder_Path, Naming_Short, Max_Diff, Numbering, Maximun_Runs):

    print("\n\n" +
          "\n**********************************************************************************" +
          "\n************************************* RUN "+ str(Numbering) + " **************************************" +
          "\n**********************************************************************************" +
          "\n\n")
    
    Chr_IDs = []
    with open(Genome_Fasta, "r") as Main_File:
        for Line in Main_File:
            if ">" in Line:
                Chr_IDs.append(Line[1:8])#Así nos aseguramos coger lo que queremos

    CSV_IDs_All = []
    with open(Path_Input, "r") as Main_File:
        reader = csv.reader(Main_File, delimiter = ",")
        for row in reader:
            CSV_IDs_All.append(row[1])#Asi obtenemos todos los IDs de la segunda columna, habra muchos muchos repetidos, pero con la anterior podré saber cual esta y cual no.

    #Con esto tengo todos los cromosomas existentes en el documento que tienen mas de 1 representante, es decir, al menos 2 representantes.
    Chr_in_Objetive = []
    for Chromosome in Chr_IDs:
        Counter = 0
        for Objetive in CSV_IDs_All:
            if Chromosome in Objetive:
                Counter += 1
        if Counter > 1:
            Chr_in_Objetive.append(Chromosome)
    
            
    for Chromosome_ID in Chr_in_Objetive:
#         if Chromosome_ID != "LinJ.01": #Por si hay que eliminar de las busquedas un cromosoma en especial
        Genome_Specific_Chromosome_Main(Path_Input,
                                        Chromosome_ID,
                                        Folder_Path,
                                        Genome_Fasta,
                                        Naming_Short,
                                        Max_Diff)

    #Y cuando termine creando el archivo MIXER, lo que hago es purificarlo completamente
    Global_Filters_Main_Output = Folder_Path + "/MIXER.csv"
    Global_Filters_Main(Global_Filters_Main_Output,
                        Global_Filters_Main_Output,
                        Genome_Fasta,
                        Naming_Short,
                        Max_Diff)
    
    RUN_SAVER_Output = Folder_Path + "/RUNS/run_" + str(Numbering) + ".csv" 
    shutil.copyfile(Global_Filters_Main_Output, RUN_SAVER_Output)
    shutil.copyfile(Global_Filters_Main_Output, Path_Input)###Asi reseteo el Path input con el nuevo documento para lanzarlo todo de nuevo. El path input antiguo ya no existe porque ha sido sobre escrito (aunque se ha guardado en RUNS)
    os.remove(Global_Filters_Main_Output)###Eliminamos el Mixer, para que luego se cree de nuevo
    
    if Numbering == Maximun_Runs:
        print("\n\n\nEND of PROGRAM")
    if Numbering < Maximun_Runs:
        Numbering += 1
#         print("\n>>>>> Counter+1:", Numbering)
        Repetitive_Blaster(Genome_Fasta, Path_Input, Folder_Path, Naming_Short, Max_Diff, Numbering, Maximun_Runs)
    
#Repetitive_Blaster(Genome_Fasta, Path_Input, Folder_Path, Naming_Short, Max_Diff, Numbering, Maximun_Runs)

    #Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    #Arg 1: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos
    #Arg 2: STRING. Directorio de la carpeta en donde se disponen los resultados del programa. Lo utilizara para generar una subcarpeta con el nombre del cromosoma
    #Arg 3: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    #Arg 4: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    #Arg 5: INT. Numeracion de comienzo de ejecucion del programa, para si es un 0, las primeras runs empezaran en 0.
    #Arg 6: INT. NUmeracion del final de ejecucion del programa, si es un 10, cuando el Arg 5 se transforme automaticamente en 10, el programa parara.