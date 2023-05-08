import csv
import subprocess

from modules.filters import chromosome_IDs
from modules.files_manager import csv_creator
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# 1)En este caso hacemos lo mismo que antes con los solapamientos pero a escala de todo el genoma. El principal cambio importante y MUY importante es el de la edicion del array Chromosome_Rows, el cual va hasta la primera funcion de Genome_Solap_Location_Filter


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.1)Aqui metemos todos los START y END para tenerlos localizados.
def Genome_Solap_Location_Filter(Chromosome_Rows):  # Todo STRING menos Chromosome_Rows
    Pos_Plus_Start = []
    Pos_Plus_End = []
    Pos_Minus_Start = []
    Pos_Minus_End = []

    for row in Chromosome_Rows:
        if "plus" in row[14]:  # Concretamente esta parte
            Pos_Plus_Start.append(int(row[10]))
            Pos_Plus_End.append(int(row[11]))
        else:
            Pos_Minus_Start.append(int(row[10]))
            Pos_Minus_End.append(int(row[11]))

    return (Pos_Plus_Start, Pos_Plus_End, Pos_Minus_Start, Pos_Minus_End)

# Genome_Solap_Location_Filter(Chromosome_Rows)
    # Arg 0: Array con los datos a filtrar


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.2)Aqui Metemos los START y END dentro de grupos por proximidad de coordenadas. Usamos la funcion anterior, a la que le impone el Chromosome_Rows
def Genome_Solap_Location_Grouping(Chromosome_Rows, DNA_sense, Max_Diff):  # Todo STRING menos Number, Max_Diff y Chromosome_Rows
    # Number es: 0 para Plus:Start, 1 para Plus:End, 2 para Minus:Start, 3 para Minus:End
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
        # De esta forma me aseguro que cada List vaya a una Matrix diferente y al final las junto

        for Position in List:
            Main_Statement = False
            for Group in Matrix:
                for Member in Group:
                    if abs(Member - Position) <= Max_Diff:
                        Group.append(Position)
                        Main_Statement = True
                        break  # Romperia codigo y saldriamos del loop "for Member" para coontinuar a "if not Found_Matrix"
            if not Main_Statement:  # Es lo mismo que si "If Found_Matrix == False". Si Group esta vacio entonces vendra aqui directamente a poner el numero en ([]), es decir, en un array 3D
                Matrix.append([Position])

    Matrix_Main = [Matrix1, Matrix2]
    return(Matrix_Main)

# Genome_Solap_Location_Grouping(Chromosome_Rows, DNA_sense, Max_Diff)
    # Arg 0: Array a filtrar.
    # Arg 1: STRING. Sentido de la hebra, puede ser "plus" o "minus"
    # Arg 2: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.3)Aqui obtenemos los minimos y maximos dependiendo de la hebra que tengamos. Se usa la funcion anterior y le impone el Genome_Rows
def Genome_Solap_MinMax(Chromosome_Rows, Max_Diff):  # Todo STRING menos Max_Diff
    Plus = Genome_Solap_Location_Grouping(Chromosome_Rows, "plus", Max_Diff)
    Minus = Genome_Solap_Location_Grouping(Chromosome_Rows, "minus", Max_Diff)

    Plus_min = [min(x) for x in Plus[0]]
    Plus_max = [max(x) for x in Plus[1]]
    Minus_max = [max(x) for x in Minus[0]]
    Minus_min = [min(x) for x in Minus[1]]

    return(Plus_min, Plus_max, Minus_max, Minus_min)

# Genome_Solap_MinMax(Chromosome_Rows, Max_Diff)
    # Arg 0: Array a filtrar.
    # Arg 1: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.4)Modulo para solapantes, es MUY importante, ya que con el se filtran los solapantes parciales que no han sido filtrados por el modulo de solapantes totales de Genome_Solap_Main. Los analiza por parejas.
def Genome_Solap_by_Pairs(Rows_to_Filter):  # Su argumento es un ARRAY 3D
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
            Homology1.append(float(Sequence[2]))  # Esta y las dos de abajo estan puestas por intentar mantener estos datos, pero en realidad no haria nada de falta
            E_Value1.append(float(Sequence[12]))
            Bit_Score1.append(float(Sequence[13]))

        Homology2 = str(round((Homology1[0] + Homology1[1]) / 2, 3))
        E_Value2 = (E_Value1[0] + E_Value1[1]) / 2
        E_Value2 = str("{:.2e}".format(E_Value2))
        Bit_Score2 = str(round((Bit_Score1[0] + Bit_Score1[1]) / 2, 1))

        if "plus" in First[14] and abs(int(First[10]) - int(Second[10])) <= 1000:  # Este numero es VITAL
            Min_Start = min(Sequence_START)
            Max_End = max(Sequence_END)
            Seq_Length = str(int(Max_End) - int(Min_Start) + 1)

            Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry "
                                          + First[1] + " -range " + Min_Start + "-" + Max_End
                                          + " -strand plus -outfmt %s",
                                          shell=True,
                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
            Seq = Seq.strip()  # Eliminar EoL caracteres

            New_Row = [First[0], First[1], Homology2, Seq_Length, First[4], First[5], "", "", "", "", str(Min_Start), str(Max_End), E_Value2, Bit_Score2, First[14], Seq]

            Rows_Final.append(New_Row)

        elif "minus" in First[14] and abs(int(First[10]) - int(Second[10])) <= 1000:
            Max_Start = max(Sequence_START)
            Min_End = min(Sequence_END)
            Seq_Length = str(int(Max_Start) - int(Min_End) + 1)

            Seq = subprocess.check_output("blastdbcmd -db ./AA_Archivos/L_infantum_ALL_36Chr.fasta -entry "
                                          + First[1] + " -range " + Min_End + "-" + Max_Start
                                          + " -strand minus -outfmt %s",
                                          shell=True,
                                          universal_newlines=True)  # MUY IMPORTANTE EL SUBPROCESS
            Seq = Seq.strip()  # Eliminar EoL caracteres

            New_Row = [First[0], First[1], Homology2, Seq_Length, First[4], First[5], "", "", "", "", str(Max_Start), str(Min_End), E_Value2, Bit_Score2, First[14], Seq]

            Rows_Final.append(New_Row)
    return (Rows_Final)

# Genome_Solap_by_Pairs(Rows_to_Filter)

    # Arg 0: array que se analizara por parejas obtenido de Genome_Solap_Main


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# 1.5)Este es el principal --> utiliza todos los de arriba y le lleva Chromosome_Rows
def Genome_Solap_Main(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input):  # Todo STRING menos Max_Diff

    Genome_Solap_Main_Matrix = []

    Chromosome_Number = chromosome_IDs(Genome_Fasta, Naming_Short)

    for Chromosome in Chromosome_Number:

        Solap_Main_Matrix = []

        Chromosome_Rows = []  # ESENCIAL PARA LA PRIMERA DEFINICION A LA QUE SE LLAMA. Tiene que ir antes de MinMax
        with open(Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter=",")
            for row in reader:
                if Chromosome in row[1]:  # Chromosome FILTER
                    Chromosome_Rows.append(row)

        MinMax = Genome_Solap_MinMax(Chromosome_Rows, Max_Diff)  # Aqui metemos los minimos y maximos anteriores en n array 3D. [0] --> Plus_min | [1] --> Plus_max | [2] --> Minus_max | [3] --> Minus_min
        Plus_Start_Matrix = []
        Plus_End_Matrix = []
        Minus_Start_Matrix = []
        Minus_End_Matrix = []

        # En esta parte se busca las secuencias que tengan tanto el minimo como el maximo, de tal forma, se obtiene el mayor alineamiento.
        with open(Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter=",")
            for row in reader:
                if Chromosome in row[1]:  # Chromosome FILTER
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

        # Ahora, por si hay solapantes teniendo uno el minimo y otra secuencia el maximo, necesitaremos las dos. Por eso, revisamos sobre las coordenadas de lo anterior (para no repetir) y, buscamos esos solapantes.
        Solap_Segments = [] # Aqui meteria todos los segmentos solapados pequeñitos
        Solap_Segments_Plus_Start = []
        Solap_Segments_Plus_End = []
        Solap_Segments_Minus_Start = []
        Solap_Segments_Minus_End = []
        with open(Path_Input, "r") as Main_File:
            reader = csv.reader(Main_File, delimiter=",")
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

        Solap_by_Pairs_Definitive = Genome_Solap_by_Pairs(Solap_Segments)  # Esto es CLAVE
        Solap_Main_Matrix += Solap_by_Pairs_Definitive

        Genome_Solap_Main_Matrix += Solap_Main_Matrix

    csv_creator(Writing_Path_Input, Genome_Solap_Main_Matrix)

# Genome_Solap_Main(Genome_Fasta, Naming_Short, Path_Input, Max_Diff, Writing_Path_Input)

    # Arg 0: STRING. Directorio del archivo en formato fasta al que queremos leer la cantidad de cromosomas, es el fasta FASTA del genoma entero
    # Arg 1: STRING. Etiqueta para leer de identificacion y numeracion de cada cromosoma en el archivo CSV. Depende del propio archivo CSV. En el caso de L. infantum es "LinJ"
    # Arg 2: STRING. Directorio del archivo en formato CSV de donde leeremos y filtraremos los datos.
    # Arg 3: INT. Numeracion con la que le indicamos el maximo valor de proximidad para las diferentes secuancias cuando tienen que ser agrupadas. MUY IMPORTANTE
    # Arg 4: STRING. Directorio del archivo en formato CSV en donde guardaremos los resultados del filtrado. Recordar la extension .csv
