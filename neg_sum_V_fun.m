function negative_sum_V = neg_sum_V_fun(Std_sS)

sS=log(Std_sS./(1-Std_sS));

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

V=zeros(1,2);
for i=1:2
    if i==1
        Psi_sigma0_G=(exp(-R1(i)*(sigma0_1(rep1)-sS(i,1)))-exp(-R2(i)*(sigma0_1(rep1)-sS(i,1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
        psi_sigma0_G=(exp(R2(i)*(sS(i,2)-sigma0_1(rep1)))-exp(R1(i)*(sS(i,2)-sigma0_1(rep1))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
        V(i) = exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*(Psi_sigma0_G*(nu_G(i)-exp(-sS(i,2))*c_B(i))+psi_sigma0_G*(exp(-sS(i,1))*nu_B(i)-c_G(i)));
    else
        Psi_sigma0_G=(exp(-R1(i)*(sigma0_2(rep2)-sS(i,1)))-exp(-R2(i)*(sigma0_2(rep2)-sS(i,1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
        psi_sigma0_G=(exp(R2(i)*(sS(i,2)-sigma0_2(rep2)))-exp(R1(i)*(sS(i,2)-sigma0_2(rep2))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
        V(i) = exp(sigma0_2(rep2))/(1+exp(sigma0_2(rep2)))*(Psi_sigma0_G*(nu_G(i)-exp(-sS(i,2))*c_B(i))+psi_sigma0_G*(exp(-sS(i,1))*nu_B(i)-c_G(i)));
    end
end
negative_sum_V=-sum(V);

end