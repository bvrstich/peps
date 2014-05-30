set output "imag_tevol.svg"
set key off
set pointsize 0.5
set term svg
plot [0:5000][-123.8:-122.5] "test" u ($1):($2) w l 
