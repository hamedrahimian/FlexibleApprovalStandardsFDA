function negative_U = neg_U_fun(Std_sS,i)

sS=log(Std_sS./(1-Std_sS));

global C;
global r;
global nu_i;
global sigma0_1;
global sigma0_2;
global R1;
global R2;
global rep1;
global rep2;

if i==1
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_1(rep1)-sS(i,1)))-exp(-R2(i)*(sigma0_1(rep1)-sS(i,1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
    Psi_sigma0_B=(exp(R2(i)*(sigma0_1(rep1)-sS(i,1)))-exp(R1(i)*(sigma0_1(rep1)-sS(i,1))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
    
    psi_sigma0_B=(exp(-R1(i)*(sS(i,2)-sigma0_1(rep1)))-exp(-R2(i)*(sS(i,2)-sigma0_1(rep1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
    psi_sigma0_G=(exp(R2(i)*(sS(i,2)-sigma0_1(rep1)))-exp(R1(i)*(sS(i,2)-sigma0_1(rep1))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
    
    Psi_sigma0=exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*Psi_sigma0_G+Psi_sigma0_B/(1+exp(sigma0_1(rep1)));
    psi_sigma0=exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*psi_sigma0_G+psi_sigma0_B/(1+exp(sigma0_1(rep1)));
else
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_2(rep2)-sS(i,1)))-exp(-R2(i)*(sigma0_2(rep2)-sS(i,1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
    Psi_sigma0_B=(exp(R2(i)*(sigma0_2(rep2)-sS(i,1)))-exp(R1(i)*(sigma0_2(rep2)-sS(i,1))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
    
    psi_sigma0_B=(exp(-R1(i)*(sS(i,2)-sigma0_2(rep2)))-exp(-R2(i)*(sS(i,2)-sigma0_2(rep2))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
    psi_sigma0_G=(exp(R2(i)*(sS(i,2)-sigma0_2(rep2)))-exp(R1(i)*(sS(i,2)-sigma0_2(rep2))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
   
    Psi_sigma0=exp(sigma0_2(rep2))/(1+exp(sigma0_2(rep2)))*Psi_sigma0_G+Psi_sigma0_B/(1+exp(sigma0_2(rep2)));
    psi_sigma0=exp(sigma0_2(rep2))/(1+exp(sigma0_2(rep2)))*psi_sigma0_G+psi_sigma0_B/(1+exp(sigma0_2(rep2)));
end

negative_U = -(Psi_sigma0*nu_i(i)-(1-psi_sigma0-Psi_sigma0)*C(i)/r);

end