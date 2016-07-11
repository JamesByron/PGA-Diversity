set logscale x 10; set xtics 1,10,1e8

set term eps color enhanced lw 3
set termoption dashed
set key c r
set xlabel "generation"
set ylabel "phenotype diversity"
set output "phenodiv.eps"
plot [][] 'avgs-all.gp' index 9 using 1:9 w lines t "proportionate accuracy only", 'avgs-all.gp' index 0 using 1:9 w lines t "tournament accuracy only ", 'avgs-all.gp' index 1 using 1:9 w lines t "phenotype diversity", 'avgs-all.gp' index 5 using 1:9 w lines t "hamming diversity"

set ylabel "estimated hamming diversity"
set output "hammingdiv.eps"
plot [][] 'avgs-all.gp' index 9 using 1:13 w lines t "proportionate accuracy only", 'avgs-all.gp' index 0 using 1:13 w lines t "tournament accuracy only", 'avgs-all.gp' index 1 using 1:13 w lines t "phenotype diversity", 'avgs-all.gp' index 5 using 1:13 w lines t "hamming diversity"


set key b r
set ylabel "accuracy"
set output "divwgt1.0.eps"
plot [][0.05:.25] 'avgs-all.gp' index 0 using 1:2 w lines t "accuracy only", 'avgs-all.gp' index 1 using 1:2 w lines t "phenotype weight 1.0", 'avgs-all.gp' index 5 using 1:2 w lines t "hamming weight 1.0"

set output "divwgt0.75.eps"
plot [][0.05:.25] 'avgs-all.gp' index 9 using 1:2 w lines t "proportionate accuracy only", 'avgs-all.gp' index 0 using 1:2 w lines t "tournament accuracy only", 'avgs-all.gp' index 2 using 1:2 w lines t "phenotype weight 0.75", 'avgs-all.gp' index 6 using 1:2 w lines t "hamming weight 0.75"

set output "divwgt0.5.eps"
plot [][0.05:.25] 'avgs-all.gp' index 0 using 1:2 w lines t "accuracy only", 'avgs-all.gp' index 3 using 1:2 w lines t "phenotype weight 0.5", 'avgs-all.gp' index 7 using 1:2 w lines t "hamming weight 0.5"

set output "divwgt0.25.eps"
plot [][0.05:.25] 'avgs-all.gp' index 0 using 1:2 w lines t "accuracy only", 'avgs-all.gp' index 4 using 1:2 w lines t "phenotype weight 0.25", 'avgs-all.gp' index 8 using 1:2 w lines t "hamming weight 0.25"

