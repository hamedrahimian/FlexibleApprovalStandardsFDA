function U_iso = U_iso_fun(Std_S)

global Std_s4U;
global pharma;
global nu_i;
global sigma0_1;
global sigma0_2;
global R1;
global R2;
global rep1;
global rep2;
global u_iso_value;
global C;
global r;

S=log(Std_S/(1-Std_S));
s=log(Std_s4U/(1-Std_s4U));

i=pharma;
if i==1
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_1(rep1)-s))-exp(-R2(i)*(sigma0_1(rep1)-s)))/(exp(-R1(i)*(S-s))-exp(-R2(i)*(S-s)));
    Psi_sigma0_B=(exp(R2(i)*(sigma0_1(rep1)-s))-exp(R1(i)*(sigma0_1(rep1)-s)))/(exp(R2(i)*(S-s))-exp(R1(i)*(S-s)));
    
    psi_sigma0_B=(exp(-R1(i)*(S-sigma0_1(rep1)))-exp(-R2(i)*(S-sigma0_1(rep1))))/(exp(-R1(i)*(S-s))-exp(-R2(i)*(S-s)));
    psi_sigma0_G=(exp(R2(i)*(S-sigma0_1(rep1)))-exp(R1(i)*(S-sigma0_1(rep1))))/(exp(R2(i)*(S-s))-exp(R1(i)*(S-s)));

    Psi_sigma0=exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*Psi_sigma0_G+Psi_sigma0_B/(1+exp(sigma0_1(rep1)));
    psi_sigma0=exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*psi_sigma0_G+psi_sigma0_B/(1+exp(sigma0_1(rep1)));
else
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_2(rep2)-s))-exp(-R2(i)*(sigma0_2(rep2)-s)))/(exp(-R1(i)*(S-s))-exp(-R2(i)*(S-s)));
    Psi_sigma0_B=(exp(R2(i)*(sigma0_2(rep2)-s))-exp(R1(i)*(sigma0_2(rep2)-s)))/(exp(R2(i)*(S-s))-exp(R1(i)*(S-s)));

    psi_sigma0_B=(exp(-R1(i)*(S-sigma0_2(rep2)))-exp(-R2(i)*(S-sigma0_2(rep2))))/(exp(-R1(i)*(S-s))-exp(-R2(i)*(S-s)));
    psi_sigma0_G=(exp(R2(i)*(S-sigma0_2(rep2)))-exp(R1(i)*(S-sigma0_2(rep2))))/(exp(R2(i)*(S-s))-exp(R1(i)*(S-s)));

    Psi_sigma0=exp(sigma0_2(rep2))/(1+exp(sigma0_2(rep2)))*Psi_sigma0_G+Psi_sigma0_B/(1+exp(sigma0_2(rep2)));
    psi_sigma0=exp(sigma0_2(rep2))/(1+exp(sigma0_2(rep2)))*psi_sigma0_G+psi_sigma0_B/(1+exp(sigma0_2(rep2)));
end

U_iso = (Psi_sigma0*nu_i(i)-(1-psi_sigma0-Psi_sigma0)*C(i)/r)-u_iso_value(i);

%U_iso = (-C(i)/r+(S*Psi_sigma0_G+(1-S)*Psi_sigma0_B)*(nu_i(i)+C(i)/r)+(s*psi_sigma0_G+(1-s)*psi_sigma0_B)*C(i)/r)-u_iso_value(i);

end