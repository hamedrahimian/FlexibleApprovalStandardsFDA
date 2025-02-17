clc;
clear all;
close all;

global C;
global c_e;
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

c_e=[10,10];

nu_G=[0.5,0.5];          
nu_B=[-0.5,-0.5];

nu_i=[1.7,1.7];
C=[15,20]; %30,100

mu=[12,17]; %30,20
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

sigma0_1=log(q0_1./(1-q0_1));
sigma0_2=log(q0_2./(1-q0_2));

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

Std_Vfoc_s=(0.01:0.01:0.99);
Std_Vfoc_S=zeros(2,99);

for rep=1:99
    fun=@root_foc_dVdS_1;
    if rep<25
        Std_S0 = Std_Vfoc_s(rep)+0.65;
    elseif rep<33
        Std_S0 = Std_Vfoc_s(rep)+0.6;
    elseif rep<40
        Std_S0 = Std_Vfoc_s(rep)+0.50;
    elseif rep<50
        Std_S0 = Std_Vfoc_s(rep)+0.35;
    else
        Std_S0 = Std_Vfoc_s(rep)+0.005;
    end
    Std_Vfoc_S(1,rep) = fzero(fun,Std_S0);
end
if mod(qhat(1)*100,1) == 0
    Std_Vfoc_S(1,qhat(1)*100)=qhat(1);
end

for rep=1:99
    fun=@root_foc_dVdS_2;
    if rep<25
        Std_S0 = Std_Vfoc_s(rep)+0.65;
    elseif rep<33
        Std_S0 = Std_Vfoc_s(rep)+0.6;
    elseif rep<40
        Std_S0 = Std_Vfoc_s(rep)+0.50;
    elseif rep<50
        Std_S0 = Std_Vfoc_s(rep)+0.35;
    else
        Std_S0 = Std_Vfoc_s(rep)+0.005;
    end
    Std_Vfoc_S(2,rep) = fzero(fun,Std_S0);
end
if mod(qhat(2)*100,1) == 0
    Std_Vfoc_S(2,qhat(2)*100)=qhat(2);
end

%% Plotting U for dU/ds=0

rep1=round((fixed_q0-min(q0_1)+0.01)*100);
rep2=round((fixed_q0-min(q0_2)+0.01)*100);
sigma0_1=log(q0_1./(1-q0_1));
sigma0_2=log(q0_2./(1-q0_2));

U_CofS=[];
for S_counter=round((fixed_q0+0.01)*100):99
    if Std_Ufoc_s(1,S_counter)<fixed_q0 
        fixed_Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        U_CofS=[U_CofS;[Std_Ufoc_S(S_counter),-neg_U_fun(fixed_Std_sS,1)]];
    end
end
U_RofS=[];
for S_counter=round((fixed_q0+0.01)*100):99
    if Std_Ufoc_s(2,S_counter)<fixed_q0
        fixed_Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        U_RofS=[U_RofS;[Std_Ufoc_S(S_counter),-neg_U_fun(fixed_Std_sS,2)]];
    end
end

figure;
plot(U_CofS(:,1),U_CofS(:,2),'LineWidth',2);
hold on;
plot(U_RofS(:,1),U_RofS(:,2),'LineWidth',2); 
xlabel('S');
ylabel('U');
% ylim([-10 20]);
legend('U_C(b_C(S),S)','U_R(b_R(S),S)');

%% Plotting V for dU/ds=0

V_CofS=[];
for S_counter=round((fixed_q0+0.01)*100):99
    if Std_Ufoc_s(1,S_counter)<fixed_q0
        fixed_Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        V_CofS=[V_CofS;[Std_Ufoc_S(S_counter),-neg_V_fun(fixed_Std_sS,1)]];
    end
end
V_RofS=[];
for S_counter=round((fixed_q0+0.01)*100):99
    if Std_Ufoc_s(2,S_counter)<fixed_q0
        fixed_Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        V_RofS=[V_RofS;[Std_Ufoc_S(S_counter),-neg_V_fun(fixed_Std_sS,2)]];
    end
end

figure;
plot(V_CofS(:,1),V_CofS(:,2),'LineWidth',2);
hold on;
% plot([0,V_CofS(:,1)],[Vmax,Vmax]);
% text(0,max(y)+0.05,'S^e_C');
plot(V_RofS(:,1),V_RofS(:,2),'LineWidth',2); 
xlabel('S');
ylabel('V');
legend('V_C(b_C(S),S)','V_R(b_R(S),S)');

%% Plotting Psi(sigma_0,G) for dU/ds=0

Psi_sigma0_G_CofS=[];
for S_counter=round((fixed_q0+0.01)*100):99
    if Std_Ufoc_s(1,S_counter)<fixed_q0
        fixed_Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        Psi_sigma0_G_CofS=[Psi_sigma0_G_CofS;[Std_Ufoc_S(S_counter),Psi_sigma0_G_fun(fixed_Std_sS,1)]];
    end
end
Psi_sigma0_G_RofS=[];
for S_counter=round((fixed_q0+0.01)*100):99
    if Std_Ufoc_s(2,S_counter)<fixed_q0
        fixed_Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        Psi_sigma0_G_RofS=[Psi_sigma0_G_RofS;[Std_Ufoc_S(S_counter),Psi_sigma0_G_fun(fixed_Std_sS,2)]];
    end
end

figure;
plot(Psi_sigma0_G_CofS(:,1),Psi_sigma0_G_CofS(:,2),'LineWidth',2);
hold on;
plot(Psi_sigma0_G_RofS(:,1),Psi_sigma0_G_RofS(:,2),'LineWidth',2); 
xlabel('S');
ylabel('Psi(sigma0,G)');
legend('Psi(sigma0,G)(b_C(S),S)','Psi(sigma0,G)(b_R(S),S)');

%% Bilevel Optimization Problem

sigma0_1=log(q0_int./(1-q0_int));
sigma0_2=log(q0_int./(1-q0_int));

A = [];
b = [];
Aeq = [];
beq = [];
fixed_Std_sS=ones(2,2,n);
flexible_Std_sS=ones(2,2,n);

fun = @neg_sum_V_fun;
fixed_nonlcon = @fixed_circlecon;
flexible_nonlcon = @flexible_circlecon;
options = optimoptions(@fmincon,'Algorithm','interior-point');
for rep1=1:n
    rep2=rep1;
    Std_sS0 = [0.2,0.5;0.2,0.5];
    lb = [Std_s_i(1),q0_int(rep1);Std_s_i(2),q0_int(rep2)];
    ub = [q0_int(rep1),Std_sS_n(1,2);q0_int(rep2),Std_sS_n(2,2)];
    [fixed_Std_sS(:,:,rep1),fixed_fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,fixed_nonlcon,options);
    [flexible_Std_sS(:,:,rep1),flexible_fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,flexible_nonlcon,options);
end

for i=1:n
    fixed_11(i)=fixed_Std_sS(1,1,i);
    fixed_12(i)=fixed_Std_sS(1,2,i);
    fixed_21(i)=fixed_Std_sS(2,1,i);
    fixed_22(i)=fixed_Std_sS(2,2,i);
    
    flexible_11(i)=flexible_Std_sS(1,1,i);
    flexible_12(i)=flexible_Std_sS(1,2,i);
    flexible_21(i)=flexible_Std_sS(2,1,i);
    flexible_22(i)=flexible_Std_sS(2,2,i);
end

%%

figure;
plot(Std_Ufoc_s(1,:),Std_Ufoc_S,'-','LineWidth',2);
hold on;
plot(Std_Ufoc_s(2,:),Std_Ufoc_S,'-','LineWidth',2);
plot(q0_int,fixed_11,':','LineWidth',2);
plot(q0_int,fixed_21,':','LineWidth',2);
plot(q0_int,fixed_12,':','LineWidth',2);
plot(q0_int,fixed_22,':','LineWidth',2);

plot(q0_int,flexible_11,'--','LineWidth',2);
plot(q0_int,flexible_21,'--','LineWidth',2);
plot(q0_int,flexible_12,'--','LineWidth',2);
plot(q0_int,flexible_22,'--','LineWidth',2);

legend('s^e_C=b_C(S^e_C)','s^e_R=b_R(S^e_R)','fixed s^e_C(q_0)','fixed s^e_R(q_0)','fixed S^e_C(q_0)','fixed S^e_R(q_0)','flexible s^e_C(q_0)','flexible s^e_R(q_0)','flexible S^e_C(q_0)','flexible S^e_R(q_0)');

xlim([0 1]);
ylim([0 1]);
xlabel('s , q_0');
ylabel('S');

%%% Plotting U* and V* as a function of q_0

fixed_U_opt_val=zeros(2,n);
fixed_V_opt_val=zeros(2,n);
flexible_U_opt_val=zeros(2,n);
flexible_V_opt_val=zeros(2,n);
for pharma=1:2
    for rep1=1:n
        rep2=rep1;
        fixed_U_opt_val(pharma,rep1)=-neg_U_fun(fixed_Std_sS(:,:,rep1),pharma);
        fixed_V_opt_val(pharma,rep1)=-neg_V_fun(fixed_Std_sS(:,:,rep1),pharma);
        
        flexible_U_opt_val(pharma,rep1)=-neg_U_fun(flexible_Std_sS(:,:,rep1),pharma);
        flexible_V_opt_val(pharma,rep1)=-neg_V_fun(flexible_Std_sS(:,:,rep1),pharma);
    end
end

figure;
plot(q0_int,fixed_U_opt_val(1,:),':r','LineWidth',2);
hold on;
plot(q0_int,fixed_U_opt_val(2,:),':b','LineWidth',2);
plot(q0_int,flexible_U_opt_val(1,:),'r','LineWidth',2);
plot(q0_int,flexible_U_opt_val(2,:),'b','LineWidth',2);
xlabel('q_0');
ylabel('U*');
legend('fixed U*_C(q_0,b_C(S^e_C),S^e_C)','fixed U*_R(q_0,b_R(S^e_R),S^e_R)','flexible U*_C(q_0,b_C(S^e_C),S^e_C)','flexible U*_R(q_0,b_R(S^e_R),S^e_R)');

% figure;
% plot(q0_int,fixed_V_opt_val(1,:),':r','LineWidth',2);
% hold on;
% plot(q0_int,fixed_V_opt_val(2,:),':b','LineWidth',2);
% plot(q0_int,flexible_V_opt_val(1,:),'r','LineWidth',2);
% plot(q0_int,flexible_V_opt_val(2,:),'b','LineWidth',2);
% xlabel('q_0');
% ylabel('V*');
% legend('fixed V*_C(q_0,b_C(S^e_C),S^e_C)','fixed V*_R(q_0,b_R(S^e_R),S^e_R)','flexible V*_C(q_0,b_C(S^e_C),S^e_C)','flexible V*_R(q_0,b_R(S^e_R),S^e_R)');

% figure;
% plot(q0_int,temp_12,'LineWidth',2);
% hold on;
% plot(q0_int,temp_22,'LineWidth',2);
% xlabel('q_0');
% ylabel('S^e');
% legend('S_C^e(q_0)','S_R^e(q_0)');