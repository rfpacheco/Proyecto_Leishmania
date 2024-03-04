from extra.useful import fasta_extractor
from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.files_manager import folder_creator

repetitive_blaster("/home/viskuit/Desktop/Projects/Leishmania/SIDER_Data/Test_1/L_infantum_1-to-5Chr.fasta",
                   "/home/viskuit/Desktop/Projects/Leishmania/SIDER_Data/Test_1/First_Blaster.csv",
                   "/home/viskuit/Desktop/Projects/Leishmania/SIDER_Data/Test_1/Results100/",
                   "LinJ",
                   600,
                   1,
                   100                  
)