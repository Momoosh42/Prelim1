import Roots
import PyPlot
using PyPlot
using Roots
using Plots
ke = 39 #nuc/sec
L = 1000 #nucleotide length
ke_star = ke/L*3600 #1/h
R_polymerase = 5000 #polymerase/cell
mass_polymerase = 4.8e5 #Daltons
R_T = R_polymerase*mass_polymerase*1.66054e-24 * 1e8 #gDW/cell * 1e8 cells
Gene_copy_number = 2 #copies/cell
Average_Cell_Volume = 0.47 #um^3
G = Gene_copy_number*(1/(6.02e23))*(1e6)/(4.7e-16) #uMol
initiation_time = 22.0   #sec
k_I = 1/initiation_time*3600 #1/h
tau = ke_star/k_I
Transcription_Saturation = 0.0136 #uMol
half_life = 5 #minutes
deg = 1/(half_life/60)*log(2) #1/hr
double_time = 40 #minutes
dilution = log(2)/(double_time/60) #1/hr
W1 = 0; #unknown
W2 = 0; #unknown
K_inhibitor = 0 #unknown
n = 2.0 #hill coefficient
I = 1e3*[0.0, 0.0005, 0.005, 0.012, 0.053, 0.216, 1.0] #array of IPTG values mM
z = [19,21,41,67,86,93,93]          #array of lacZ values
print(I)
##############################################
#1a.
V = 1 #ml
Nc = 1 * 10^8 #cells/ml for OD = 0.1
mRNA_gDW = 2.3*10^-15 #gDW-mRNA/cell
Total_mRNA_copies = 1380 #molecules
mc = V*Nc*mRNA_gDW*z./(Total_mRNA_copies)  #gDW
print(mc)

#################################################
#1b,c,d
theta = (ke_star*R_T)/(dilution+deg)
K_x = theta*(G/((tau*Transcription_Saturation)+(tau+1)*G))

#When IPTG = 0 only value for u(I) should be leaky expression
#W1/(1+W1)
W1 = ((K_x-mc[1])/mc[1])^-1

#When I = 1 the system is at saturation therefore all of the inducer can be considered bound
#Then u(I) = (W1+W2)/(1+W1+W2)
println(K_x)
println(W1)
f(W2) = (K_x*(W1+W2)/(1+W1+W2))-mc[end]
W2 = fzero(f, 1)
q = (K_x*(W1+W2)/(1+W1+W2))
print(W2)

#Code continued in one_d.m
