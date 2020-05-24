%Continuation from 1d.jl for solving question 1c and d
%Values for mc,W1 and W2 were calculated in 1d.jl
%Values were copy and pasted into this script

clf
mc = [3.166666666666669e-9, 3.500000000000003e-9, 6.833333333333339e-9, 1.1166666666666676e-8, 1.4333333333333344e-8, 1.5500000000000013e-8, 1.5500000000000013e-8];
I = [0.0000001, 0.5, 5.0, 12.0, 53.0, 216.0, 1000.0];
semilogx(I,mc)
K = 8;
f_L = [];
n = 2;
K_x = 1.7041651782412267e-6;
W1 = 0.0018616516387985468;
W2 = 0.0073171957725891405;
%Loop calculates all the fractional binding values
for i=1:1:length(I)
        f_L(i)=((I(i)).^n)./(K^n+I(i).^n);
end
%Loop calculates all of the B values which will be compared to the actually
%measured data
for j = 1:1:length(I)
    mc_calc(j) = (K_x*(W1+W2*f_L(j))/(1+W1+W2*f_L(j)))
end
hold on
semilogx(I,mc_calc)
title('Calculated [lacZ] with actual [lacZ] vs. IPTG')
xlabel('[IPTG] (uMol)')
ylabel('lacZ mRNA (gDW)')
legend('actual data','calculated data','Location','northwest')