set output "imag_tevol.svg"
set key right
set pointsize 0.5
set term svg
plot "t=0.01.out" u ($1):($2) w l t "t=0.01","t=0.001.out" u ($1/10):($2) w l t "t=0.001"
