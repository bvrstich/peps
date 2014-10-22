set output "energy.svg"
set key top right
set pointsize 0.5
set term svg
set xlabel "steps"
set ylabel "E"
plot [2500:2780] "energy_Dx=12.out" w l t "D_aux=12", "energy_Dx=16.out" w l t "D_aux=16", "energy_Dx=20.out" w l t "D_aux=20"
