function V_reward = V_fun_reward(Std_sS,i)

sS=log(Std_sS./(1-Std_sS));

global nu_G;
global nu_B;
global sigma0_1;
global sigma0_2;
global R1;
global R2;
global rep1;
global rep2;

if i==1
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_1(rep1)-sS(i,1)))-exp(-R2(i)*(sigma0_1(rep1)-sS(i,1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
    V_reward = Psi_sigma0_G*exp(sigma0_1(rep1))/(1+exp(sigma0_1(rep1)))*nu_G(i);
else
    Psi_sigma0_G=(exp(-R1(i)*(sigma0_2(rep2)-sS(i,1)))-exp(-R2(i)*(sigma0_2(rep2)-sS(i,1))))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
    V_reward = Psi_sigma0_G*exp(sigma0_2(rep2))/(1+exp(sigma0_2(rep2)))*nu_G(i);
end

end