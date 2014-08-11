set output "imag_tevol.svg"
set key off
set pointsize 0.5
set term svg
plot [3000:4000] "14x14.out" u ($1):($2) w l 
