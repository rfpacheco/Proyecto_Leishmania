from extra.useful import fasta_extractor
from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.files_manager import folder_creator

repetitive_blaster("/home/viskuit/Desktop/Projects/Leishmania/SIDER_Data/Test_3/Leishmania_infantum_chr_1-5.fasta",
                   "/home/viskuit/Desktop/Projects/Leishmania/SIDER_Data/Test_3/First_Blaster.csv",
                   "/home/viskuit/Desktop/Projects/Leishmania/SIDER_Data/Test_3/Results/",
                   "LinJ",
                   600,
                   1,
                   2                  
)