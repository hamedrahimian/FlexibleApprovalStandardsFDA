clc;
clear;
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
global Std_Vfoc_s;
global Std_s_i;

fixed_q0=0.4;
flag=0; %1 means R is changing

nu_B_C_list=-0.85:0.02:-0.15;
nu_G_R_list=0.25:0.02:1.25;


if flag==1
    itNum=length(nu_G_R_list);
else
    itNum=length(nu_B_C_list);
end

Std_sS=ones(2,2,itNum);
qhat=zeros(itNum,2);

for it=1:itNum

if flag==1
    nu_G=[0.75,nu_G_R_list(it)];          
    nu_B=[-0.5,-0.5];
else
    nu_G=[0.75,0.75];          
    nu_B=[nu_B_C_list(it),-0.5];
end

r=5;
nu_i=[2,2];
C=[10,10];

mu=[12,12];
rho=[4,4];
mu_prime=2*mu.^2./rho.^2;
%mu_prime=12;

R1=0.5*(1-sqrt(1+4*r./mu_prime));
R2=0.5*(1+sqrt(1+4*r./mu_prime));

rep=1;
qhat(it,:)=nu_B./(nu_B-nu_G);

Std_Ufoc_S=qhat(it,1);
fun=@root_foc_dUds_1;
Std_s0 = Std_Ufoc_S/2;
Std_s_i(1) = fzero(fun,Std_s0);

Std_Ufoc_S=qhat(it,2);
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

if min(q0_int)<fixed_q0 & max(q0_int)>fixed_q0
    q0_int=fixed_q0;
else
    dbstop();
end
[temp,n]=size(q0_int);


% %% Bilevel Optimization Problem

sigma0_1=log(fixed_q0./(1-fixed_q0));
sigma0_2=log(fixed_q0./(1-fixed_q0));

A = [];
b = [];
Aeq = [];
beq = [];

fun = @neg_sum_V_fun;
nonlcon = @circlecon;
options = optimoptions(@fmincon,'Algorithm','interior-point');
for rep1=1:n
    rep2=rep1;
    Std_sS0 = [0.2,0.3;0.2,0.3];
    lb = [Std_s_i(1),fixed_q0;Std_s_i(2),fixed_q0];
    ub = [fixed_q0,Std_sS_n(1,2);fixed_q0,Std_sS_n(2,2)];
    [Std_sS(:,:,it),fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

temp_11(it)=Std_sS(1,1,it);
temp_12(it)=Std_sS(1,2,it);
temp_21(it)=Std_sS(2,1,it);
temp_22(it)=Std_sS(2,2,it);

end

%% Plotting U* and V*

U_opt_val=zeros(2,itNum);
V_opt_val=zeros(2,itNum);
for pharma=1:2
    for it=1:itNum
        if flag==1
            nu_G=[0.75,nu_G_R_list(it)];
%             nu_B=[-0.5,-0.5];
        else
%             nu_G=[0.75,0.75];
            nu_B=[nu_B_C_list(it),-0.5];
        end
        U_opt_val(pharma,it)=-neg_U_fun(Std_sS(:,:,it),pharma);
        V_opt_val(pharma,it)=-neg_V_fun(Std_sS(:,:,it),pharma);
    end
end

if flag==1

figure;
plot(nu_G_R_list,temp_22,'LineWidth',2);
xlabel('nu^G');
ylabel('S^e');
legend('S^e(nu^G)');

figure;
plot(nu_G_R_list,temp_21,'LineWidth',2);
xlabel('nu^G');
ylabel('s^e');
legend('s^e(nu^G)');

figure;
plot(nu_G_R_list,U_opt_val(2,:),'LineWidth',2);
xlabel('nu^G');
ylabel('U*');
legend('U*(nu^G,S^e(nu^G))');

figure;
plot(qhat(:,2),V_opt_val(2,:),'LineWidth',2);
xlabel('qhat');
ylabel('V*');
legend('V*(qhat,S^e(qhat))');

figure;
plot(nu_G_R_list,V_opt_val(2,:),'LineWidth',2);
xlabel('nu^G');
ylabel('V*');
legend('V*(nu^G,S^e(nu^G))');

else
    
figure;
plot(nu_B_C_list,temp_12,'LineWidth',2);
xlabel('nu^B');
ylabel('S^e');
legend('S^e(nu^B)');

figure;
plot(nu_B_C_list,temp_11,'LineWidth',2);
xlabel('nu^B');
ylabel('s^e');
legend('s^e(nu^B)');

figure;
plot(nu_B_C_list,U_opt_val(1,:),'LineWidth',2);
xlabel('nu^B');
ylabel('U*');
legend('U*(nu^B,S^e(nu^B))');

figure;
plot(qhat(:,1),V_opt_val(1,:),'LineWidth',2);
xlabel('qhat');
ylabel('V*');
legend('V*(qhat,S^e(qhat))');

figure;
plot(nu_B_C_list,V_opt_val(1,:),'LineWidth',2);
xlabel('nu^B');
ylabel('V*');
legend('V*(nu^B,S^e(nu^B))');

end
