set output "random.svg"
set key top right
set pointsize 0.5
set term svg
set xlabel "steps"
set ylabel "E"
plot [4500:30000] "random.out" w l t "D_aux=16"
