set term postscript enhanced color "Helvetica" 16
set output "lsstest-absolute-error.eps"
set data style lines
set logscale y
plot  \
"lsstest-absolute-error.dat" title "absolute error"
