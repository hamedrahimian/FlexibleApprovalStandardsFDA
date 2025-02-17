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


nu_G=[0.75,1];          
nu_B=[-0.5,-0.5];

%nu_i=[2,2];     %U^e_C<U^e_R
nu_i=[3,2];    %U^e_C>U^e_R
C=[10,20];

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

%% Finding the feasible ranges

sigma0_1=log(q0_int./(1-q0_int));
sigma0_2=log(q0_int./(1-q0_int));
rep1=20;
rep2=50;
q0_1=q0_int(rep1);
q0_2=q0_int(rep2);
% rep2=rep1;

rep=100;
Std_Ufoc_S(rep)=Std_sS_n(1,2);
fun=@root_foc_dUds_1;
Std_s0 = Std_Ufoc_S(rep)/2;    
s1min=max(Std_s_i(1),Std_Ufoc_s(1,floor(q0_1*100)));
s1max=min(q0_1,fzero(fun,Std_s0));

Std_Ufoc_S(rep)=Std_sS_n(2,2);
fun=@root_foc_dUds_2;
Std_s0 = Std_Ufoc_S(rep)/2;  
s2min=max(Std_s_i(2),Std_Ufoc_s(2,floor(q0_2*100)));
s2max=min(q0_2,fzero(fun,Std_s0));

fun=@inv_of_root_foc_dUds_1;
Std_s=s1min;
Std_S0 = Std_s*5;
S1min = fzero(fun,Std_S0);
Std_s=s1max;
Std_S0 = Std_s*3;
S1max = fzero(fun,Std_S0);

fun=@inv_of_root_foc_dUds_2;
Std_s=s2min;
Std_S0 = Std_s*6;
S2min = fzero(fun,Std_S0);
Std_s=s2max;
Std_S0 = Std_s*3;
S2max = fzero(fun,Std_S0);

S1_set=ceil(S1min*100)/100:0.01:floor(S1max*100)/100;

%% Bilevel Optimization Problem for a fixed S1

Std_sS=ones(2,2,length(S1_set));

A = [];
b = [];
Aeq = [];
beq = [];

fun = @neg_sum_V_fun;
nonlcon = @fixed_S1_circlecon;
options = optimoptions(@fmincon,'Algorithm','interior-point');

for i=1:length(S1_set)
    fixed_S1=(S1_set(i));
    Std_sS0 = [0.2,0.3;0.2,0.3];
    lb = [Std_s_i(1),q0_int(rep1);Std_s_i(2),q0_int(rep2)];
    ub = [q0_int(rep1),Std_sS_n(1,2);q0_int(rep2),Std_sS_n(2,2)];
    [Std_sS(:,:,i),fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

for i=1:length(S1_set)
    temp_11(i)=Std_sS(1,1,i);
    temp_12(i)=Std_sS(1,2,i);
    temp_21(i)=Std_sS(2,1,i);
    temp_22(i)=Std_sS(2,2,i);
end

deltaSR_dev_deltaSC_first10=(Std_sS(2,2,11)-Std_sS(2,2,2))/(Std_sS(1,2,11)-Std_sS(1,2,2));
deltaSR_dev_deltaSC_first20=(Std_sS(2,2,20)-Std_sS(2,2,1))/(Std_sS(1,2,20)-Std_sS(1,2,1));
deltaSR_dev_deltaSC_first30=(Std_sS(2,2,30)-Std_sS(2,2,1))/(Std_sS(1,2,30)-Std_sS(1,2,1));
% deltaS2_dev_deltaS1_first40=(Std_sS(2,2,40)-Std_sS(2,2,1))/(Std_sS(1,2,40)-Std_sS(1,2,1));

%% Plotting U*1=U*2 in S1-S2 space

figure;
plot(temp_12,temp_22,'LineWidth',2);
xlabel('S_C');
ylabel('S_R');
% xlim([S1min S1max]);
% ylim([S2min S2max]);
legend('U_C=U_R');
