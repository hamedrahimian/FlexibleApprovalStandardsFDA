clc;
clear all;
close all;

global C;
global r;
global R1;
global R2;
global nu_i;
global nu_G;
global nu_B;
global sigma0_1;
global sigma0_2;
global rep;
global rep1;
global rep2;
global Std_Ufoc_S;
global Std_s_i;
global fixed_S1;
global Std_s;


nu_G=[0.75,2];          
nu_B=[-0.5,-0.5];

%nu_i=[2,2];     %U^e_C<U^e_R
nu_i=[3,2];    %U^e_C>U^e_R
C=[15,20];

mu=[12,15];   

% nu_G=[0.75,1.25];          
% nu_B=[-0.5,-0.5];
% 
% nu_i=[3,2];
% C=[10,30];
% 
% mu=[12,15];       %S_C=S_R making the complete shift (common foc above rare)
% %mu=[12,24];      %S_C=S_R sort of making the shift (common foc conicide with rare)
% % mu=[8,20];       %S_C=S_R not making the shift (common foc below rare)


rho=[4,4];
mu_prime=2*mu.^2./rho.^2;
%mu_prime=12;
r=5;

R1=0.5*(1-sqrt(1+4*r./mu_prime));
R2=0.5*(1+sqrt(1+4*r./mu_prime));

rep=1;
qhat=nu_B./(nu_B-nu_G);

Std_Ufoc_S=qhat(1);
fun=@root_foc_dUds_1;
Std_s0 = Std_Ufoc_S/2;
Std_s_i(1) = fzero(fun,Std_s0);

Std_Ufoc_S=qhat(2);
fun=@root_foc_dUds_2;
Std_s0 = Std_Ufoc_S/2;
Std_s_i(2) = fzero(fun,Std_s0);

Std_sS_n=zeros(2,2);

fun=@root_nash_1;
Std_sS0 = [0.2,0.5];
Std_sS_n(1,:) = fsolve(fun,Std_sS0);

fun=@root_nash_2;
Std_sS0 = [0.2,0.5];
Std_sS_n(2,:) = fsolve(fun,Std_sS0);

q0_1=(ceil(Std_s_i(1)*100)/100:0.01:floor(Std_sS_n(1,2)*100)/100);
q0_2=(ceil(Std_s_i(2)*100)/100:0.01:floor(Std_sS_n(2,2)*100)/100);
n1=length(q0_1);
n2=length(q0_2);

q0_int=(max(ceil(Std_s_i*100)/100):0.01:min(floor(Std_sS_n(:,2)*100)/100));
[temp,n]=size(q0_int);

% sigma0_1=log(q0_1./(1-q0_1));
% sigma0_2=log(q0_2./(1-q0_2));

%% solve first order conditions dU/ds=0 and dV/dS=0

Std_Ufoc_S=(0.01:0.01:0.99);
Std_Ufoc_s=zeros(2,99);

for rep=1:99
    fun=@root_foc_dUds_1;
    Std_s0 = Std_Ufoc_S(rep)/2;
    Std_Ufoc_s(1,rep) = fzero(fun,Std_s0);
end

for rep=1:99
    fun=@root_foc_dUds_2;
    Std_s0 = Std_Ufoc_S(rep)/2;
    Std_Ufoc_s(2,rep) = fzero(fun,Std_s0);
end

%%

q0_C=0.5;
q0_R=0.5;
sigma0_1=log(q0_1./(1-q0_1));
sigma0_2=log(q0_2./(1-q0_2));
rep1=round((q0_C-min(q0_1)+0.01)*100);
rep2=round((q0_R-min(q0_2)+0.01)*100);

for i=1:length(q0_1)
    for j=1:length(q0_2)
        S_C=min(q0_1)+i/100-0.01;
        S_R=min(q0_2)+j/100-0.01;
        s_C=Std_Ufoc_s(1,floor(S_C*100));
        s_R=Std_Ufoc_s(2,floor(S_R*100));
        Std_sS=[s_C,S_C;s_R,S_R];
        U_RminusU_C(i,j)=neg_U_fun(Std_sS,1)-neg_U_fun(Std_sS,2);
    end
end


%%

figure;
meshc(q0_2,q0_1,U_RminusU_C,0,'LineWidth',4);
% zlim([-0.0001 0.0001]);
legend('U_R-U_C');
xlabel('S_R');
ylabel('S_C');
zlabel('U_R-U_C');
colorbar;

figure;
meshc(q0_2,q0_1,U_RminusU_C,'LineWidth',2);
zlim([-0.0001 0.0001]);
legend('U_R-U_C');
xlabel('S_R');
ylabel('S_C');
zlabel('U_R-U_C');
colorbar;

figure;
meshc(q0_2,q0_1,U_RminusU_C,'LineWidth',2);
% zlim([-0.0001 0.0001]);
legend('U_R-U_C');
xlabel('S_R');
ylabel('S_C');
zlabel('U_R-U_C');
colorbar;

