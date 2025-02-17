function V_iso = V_iso_fun(Std_s)

global Std_S4V;
global pharma;
global nu_G;
global nu_B;
global c_G;
global c_B;
global sigma0_1;
global sigma0_2;
global R1;
global R2;
global rep1;
global rep2;
global v_iso_value;

s=log(Std_s/(1-Std_s));
S=log(Std_S4V/(1-Std_S4V));

i=pharma;
if i==1
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_1(rep1)-s))-exp(-R2(i)*(sigma0_1(rep1)-s)))/(exp(-R1(i)*(S-s))-exp(-R2(i)*(S-s)));
    psi_sigma0_G=(exp(R2(i)*(S-sigma0_1(rep1)))-exp(R1(i)*(S-sigma0_1(rep1))))/(exp(R2(i)*(S-s))-exp(R1(i)*(S-s)));
    V_iso = (exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*(Psi_sigma0_G*(nu_G(i)-exp(-S)*c_B(i))+psi_sigma0_G*(nu_B(i)-exp(-s)*c_G(i))))-v_iso_value(i);
else
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_2(rep2)-s))-exp(-R2(i)*(sigma0_2(rep2)-s)))/(exp(-R1(i)*(S-s))-exp(-R2(i)*(S-s)));
    psi_sigma0_G=(exp(R2(i)*(S-sigma0_2(rep1)))-exp(R1(i)*(S-sigma0_2(rep1))))/(exp(R2(i)*(S-s))-exp(R1(i)*(S-s)));
    V_iso = (exp(sigma0_2(rep1))/(1+exp(sigma0_2(rep1)))*(Psi_sigma0_G*(nu_G(i)-exp(-S)*c_B(i))+psi_sigma0_G*(nu_B(i)-exp(-s)*c_G(i))))-v_iso_value(i);
end

end