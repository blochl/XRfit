#!/usr/bin/gnuplot
#set terminal postscript eps color enhanced font "Verdana, 16"
set terminal pngcairo color enhanced font "Verdana, 10"

set style data lines
set style line 10 lt 1 lw 1 linecolor rgb "black"
set style line 20 lt 1 lw 1 linecolor rgb "red"
set style line 30 lt 1 lw 1 linecolor rgb "brown"
set style line 40 lt 1 lw 1 linecolor rgb "magenta"
set style line 50 lt 1 lw 1 linecolor rgb "cyan"
set style line 60 lt 1 lw 1 linecolor rgb "green"
#set noytic
set xlabel "2{/Symbol \\561}"
set ylabel "Intensity (log, normalized)"

###############################################
############## Edit - Basic Data ##############
set output 'ReplaceMe.png'
set fit logfile 'ReplaceMe.log'
set title "ReplaceTitle" font ",15"
set key top left
FILE = 'ReplaceMe.xye'
LAMBDA = 0.0000000 # Wavelength in A
XMIN = 10.0
XMAX = 11.0
DELTA = 0.01 # That's how much to show to left and right from the fit area on the plot.
################## End edit ###################
###############################################

fl(x,a,b,cl) = (a*cl**2)/((x - b)**2 + cl**2)          # Lorentzian
fg(x,a,b,cg) = a*exp(-((x - b)**2)/(2*cg**2))          # Gaussian
lin(x,m,l) = m*x + l                                   # Linear
g(x,min,max) = ( (min<=x && x<=max) ? 1.0 : (1/0) )    # Cutoff function
lr(n) = (tanh(n) + 1.0)/2.0                            # Lorentzianity limitation
Size(b,cl) = 0.1*LAMBDA/(cos(0.5*b*pi/180.0)*(2.0*abs(cl)*pi/180.0)) # result in nm
# Error estimaterd based on: sqrt(((dS/db)*Db)**2 + ((dS/dc)*Dc)**2)
SizeError(b,be,cl,cle) = sqrt(LAMBDA**2*be**2*sin(pi*b/360.0)**2/(1600.0*cl**2*cos(pi*b/360.0)**4) + 81*LAMBDA**2*cle**2/(pi**2*cl**4*cos(pi*b/360.0)**2))
Microstrain(b,cg) = (abs(cg)*pi/180.0)/(2*tan(0.5*b*pi/180.0))
# Error estimaterd based on: sqrt(((dM/db)*Db)**2 + ((dM/dc)*Dc)**2)
MicrostrainError(b,be,cg,cge) = sqrt(pi**4*cg**2*be**2*(tan(pi*b/360.0)**2 + 1)**2/(16796160000.0*tan(pi*b/360.0)**4) + pi**2*cge**2/(129600.0*tan(pi*b/360.0)**2))
f1(x) = lr(n1)*fl(x,a1,b1,cl1) + (1-lr(n1))*fg(x,a1,b1,cg1)
f2(x) = lr(n2)*fl(x,a2,b2,cl2) + (1-lr(n2))*fg(x,a2,b2,cg2)
f3(x) = lr(n3)*fl(x,a3,b3,cl3) + (1-lr(n3))*fg(x,a3,b3,cg3)
### Comment according to number of peaks: ####
f(x) = lin(x,m,l) + f1(x)# + f2(x)# + f3(x)
##############################################

stats FILE u (($1 >= XMIN && $1 <= XMAX) ? $2 : 1/0)
stats FILE u (($1 >= XMIN && $1 <= XMIN + DELTA) ? $2 : 1/0) name "smin"
stats FILE u (($1 >= XMAX - DELTA && $1 <= XMAX) ? $2 : 1/0) name "smax"
###############################################
############ Edit - Initial values ############
### First peak
cl1 = 0.032; cg1 = 0.032; # Broadening
a1 = 0.9;                 # Height
b1 = 10.500;              # Center 
### Second peak (COMMENT IF NOT NEEDED)
#cl2 = 0.025; cg2 = 0.025; # Broadening
#a2 = 0.18;                 # Height
#b2 = 10.482;             # Center
### Third peak (COMMENT IF NOT NEEDED)
#cl3 = 0.008; cg3 = 0.008; # Broadening
#a3 = 0.01;                 # Height
#b3 = 10.63;              # Center
### Lorentzianity: 4 => full Lorentzian, -4 => Full Gaussian
n1 = 0.7e-0;
#n2 = 0.1e-6;
#n3 = 0.1e-6;
### Linear background
LDOWN = 0.003 # Shift linear BG down by this value
################## End edit ###################
###############################################
# Linear background calculation
# y = m*x + l
m = (smax_mean/STATS_max - smin_mean/STATS_max)/(XMAX - XMIN);
l = (smin_mean/STATS_max - m*XMIN)-LDOWN;

set fit errorvariables
###############################################
### Touch this when number of peaks changes ###
## Comment before "fit" when just testing initial values ###
## Comment before the relevant row when fitting less than 3 peaks ###
fit [ XMIN : XMAX ] f(x) FILE u 1:($2/STATS_max):($3/STATS_max) via\
    cl1, cg1, a1, b1, n1, m, l\
#    , cl2, cg2, a2, b2, n2\
#    , cl3, cg3, a3, b3, n3
###############################################
###############################################

set xrange [ XMIN - DELTA : XMAX + DELTA ] noreverse nowriteback
set yrange [ * : 1.05 ] noreverse nowriteback
set logscale y # Comment if log scale in y is not needed.
set label sprintf("{/Symbol \\143}^2 = %.3f", FIT_STDFIT**2) at graph 0.28 , graph 0.94
set label sprintf("{/Symbol \\150} = %.3f", lr(n1)) at b1 , graph 0.22
set label sprintf("L = %.3f {/Symbol \\261} %.3f nm", Size(b1,cl1), SizeError(b1,b1_err,cl1,cl1_err)) at b1 , graph 0.17
set label sprintf("{/Symbol \\163} = %.3e {/Symbol \\261} %.3e", Microstrain(b1,cg1), MicrostrainError(b1,b1_err,cg1,cg1_err)) at b1 , graph 0.12
###############################################
### Touch this when number of peaks changes ###
## Second peak (COMMENT IF NOT NEEDED):
#set label sprintf("{/Symbol \\150} = %.3f", lr(n2)) at b2 , graph 0.21
#set label sprintf("L = %.3f nm", Size(b2,cl2)) at b2 , graph 0.16
#set label sprintf("{/Symbol \\163} = %.3e", Microstrain(b2,cg2)) at b2 , graph 0.11
## Third peak (COMMENT IF NOT NEEDED):
#set label sprintf("{/Symbol \\150} = %.3f", lr(n3)) at b3 , graph 0.21
#set label sprintf("L = %.3f nm", Size(b3,cl3)) at b3 , graph 0.16
#set label sprintf("{/Symbol \\163} = %.3e", Microstrain(b3,cg3)) at b3 , graph 0.11
###############################################
###############################################
set samples 1000
plot FILE  u 1:($2/STATS_max)\
     with points pointtype 2 pointsize 1 linecolor rgb "black"\
     title "Experiment",\
     f(x)*g(x,XMIN,XMAX) ls 20 title "Fitting"\
     , lin(x,m,l)*g(x,XMIN,XMAX) ls 30 title "Background" # Add "\" if more than one peak
#     , (lin(x,m,l) + f1(x))*g(x,XMIN,XMAX) ls 40 title "Peak 1"\
#     , (lin(x,m,l) + f2(x))*g(x,XMIN,XMAX) ls 50 title "Peak 2"\
#     , (lin(x,m,l) + f3(x))*g(x,XMIN,XMAX) ls 60 title "Peak 3"\
print "nu1 = ",lr(n1),"; +",lr(n1+n1_err)-lr(n1),lr(n1-n1_err)-lr(n1)
print "FWHM G1 = ",2*sqrt(2*log(2))*cg1
print "FWHM L1 = ",2*cl1
print "Size1 = ",Size(b1,cl1)," nm"
###############################################
### Touch this when number of peaks changes ###
## Second peak (COMMENT IF NOT NEEDED):
#print "nu2 = ",lr(n2),"; +",lr(n2+n2_err)-lr(n2),lr(n2-n2_err)-lr(n2)
#print "FWHM G2 =",2*sqrt(2*log(2))*cg2
#print "FWHM L2 =",2*cl2
#print "Size2 = ",Size(b2,cl2)," nm"
## Third peak (COMMENT IF NOT NEEDED):
#print "nu3 = ",lr(n3),"; +",lr(n3+n3_err)-lr(n3),lr(n3-n3_err)-lr(n3)
#print "FWHM G3 =",2*sqrt(2*log(2))*cg3
#print "FWHM L3 =",2*cl3
#print "Size3 = ",Size(b3,cl3)," nm"
###############################################
###############################################
