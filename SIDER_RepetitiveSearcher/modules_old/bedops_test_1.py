from extra.useful import fasta_extractor
from modules.blaster import blastn_dic, blastn_blaster, repetitive_blaster
from modules.files_manager import folder_creator

repetitive_blaster("/home/rfpacheco/Desktop/Projects/SIDER_Data/Z_BEDOPS_Test_1/1-5_chr.fasta",
                   "/home/rfpacheco/Desktop/Projects/SIDER_Data/Z_BEDOPS_Test_1/First_Blaster.csv",
                   "/home/rfpacheco/Desktop/Projects/SIDER_Data/Z_BEDOPS_Test_1/Results/",
                   "LinJ",
                   600,
                   1,
                   2                  
)