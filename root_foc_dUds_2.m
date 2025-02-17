function U_foc_s = root_foc_dUds_2(Std_s)

global Std_Ufoc_S;
global rep;
S=log(Std_Ufoc_S(rep)/(1-Std_Ufoc_S(rep)));
s=log(Std_s/(1-Std_s));

global r;
global nu_i;
global C;

global R1;
global R2;

a=(R1(2)-R2(2))/(exp(-R1(2)*(S-s))-exp(-R2(2)*(S-s)));
b=(R2(2)*exp(R2(2)*(S-s))-R1(2)*exp(R1(2)*(S-s)))/(exp(R2(2)*(S-s))-exp(R1(2)*(S-s)));

U_foc_s=a*(1+exp(-S))*(nu_i(2)+C(2)/r)+C(2)/r*(b*(1+exp(-s))-exp(-s));

end