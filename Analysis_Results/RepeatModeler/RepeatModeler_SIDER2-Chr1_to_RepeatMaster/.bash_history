BuildDatabase SIDER2_Chr1_L_Infantum.fasta 
BuildDatabase -name SIDER2 SIDER2_Chr1_L_Infantum.fasta 
RepeatModeler -database SIDER2 -LTRStruct > run.out
RepeatModeler -h
RepeatModeler -database SIDER2 -LTRStruct -quick > run.out
exit
