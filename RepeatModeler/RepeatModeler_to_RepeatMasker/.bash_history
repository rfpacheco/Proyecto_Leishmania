BuildDatabase -name L_infantum L_infantum_ALL_36Chr.fasta 
RepeatModeler -database L_infantum -threads 12 -LTRStruct >& run.out &
ls
RepeatMasker -threads 2 -a -gff -no_is -lib L_infantum-families.fa L_infantum_ALL_36Chr.fasta &> RM.run.out &
tail -f RM.run.out 
RepeatMasker -threads 2 -a -gff -no_is -lib L_infantum-families.fa L_infantum_ALL_36Chr.fasta &> RM.run.out &
tail -f RM.run.out 
RepeatMasker -a -gff -no_is -lib L_infantum-families.fa L_infantum_ALL_36Chr.fasta &> RM.run.out &
