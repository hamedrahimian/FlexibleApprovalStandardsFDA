function V_foc_S = root_foc_dVdS_2(Std_S)

global Std_Vfoc_s;
global rep;
s=log(Std_Vfoc_s(rep)/(1-Std_Vfoc_s(rep)));
S=log(Std_S/(1-Std_S));

global nu_G;
global nu_B;
global R1;
global R2;

f=(R1(2)*exp(-R1(2)*(S-s))-R2(2)*exp(-R2(2)*(S-s)))/(exp(-R1(2)*(S-s))-exp(-R2(2)*(S-s)));

V_foc_S=f*(nu_G(2)+exp(-S)*nu_B(2))-exp(-S)*nu_B(2);

end