clc;
clear all;
close all;

global C;
global c_G;
global c_B;
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

c_G=[0,0];
c_B=[1.25,1];

nu_G=[1,1.25];          
nu_B=[0,0];

nu_i=[2.5,2];
C=[10,15];

mu=[9,9];
rho=[4,4];
mu_prime=2*mu.^2./rho.^2;
%mu_prime=12;
r=5;

R1=0.5*(1-sqrt(1+4*r./mu_prime));
R2=0.5*(1+sqrt(1+4*r./mu_prime));

rep=1;
qhat=(nu_B-c_B)./(nu_B-c_B-nu_G-c_G);

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

% q0=0.45;
% rep1=round((q0-min(q0_1)+0.01)*100);
% rep2=round((q0-min(q0_2)+0.01)*100);
% sigma0_1=log(q0_1./(1-q0_1));
% sigma0_2=log(q0_2./(1-q0_2));
% 
% U_CofS=[];
% for S_counter=round((q0+0.01)*100):99
%     if Std_Ufoc_s(1,S_counter)<q0 
%         Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
%         U_CofS=[U_CofS;[Std_Ufoc_S(S_counter),-neg_U_fun(Std_sS,1)]];
%     end
% end
% U_RofS=[];
% for S_counter=round((q0+0.01)*100):99
%     if Std_Ufoc_s(2,S_counter)<q0
%         Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
%         U_RofS=[U_RofS;[Std_Ufoc_S(S_counter),-neg_U_fun(Std_sS,2)]];
%     end
% end
% 
% figure;
% plot(U_CofS(:,1),U_CofS(:,2),'LineWidth',2);
% hold on;
% plot(U_RofS(:,1),U_RofS(:,2),'LineWidth',2); 
% xlabel('S');
% ylabel('U');
% % ylim([-10 20]);
% legend('U_C(b_C(S),S)','U_R(b_R(S),S)');
% 
% %% Plotting V for dU/ds=0
% 
% V_CofS=[];
% for S_counter=round((q0+0.01)*100):99
%     if Std_Ufoc_s(1,S_counter)<q0
%         Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
%         V_CofS=[V_CofS;[Std_Ufoc_S(S_counter),-neg_V_fun(Std_sS,1)]];
%     end
% end
% V_RofS=[];
% for S_counter=round((q0+0.01)*100):99
%     if Std_Ufoc_s(2,S_counter)<q0
%         Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
%         V_RofS=[V_RofS;[Std_Ufoc_S(S_counter),-neg_V_fun(Std_sS,2)]];
%     end
% end
% 
% figure;
% plot(V_CofS(:,1),V_CofS(:,2),'LineWidth',2);
% hold on;
% % plot([0,V_CofS(:,1)],[Vmax,Vmax]);
% % text(0,max(y)+0.05,'S^e_C');
% plot(V_RofS(:,1),V_RofS(:,2),'LineWidth',2); 
% xlabel('S');
% ylabel('V');
% legend('V_C(b_C(S),S)','V_R(b_R(S),S)');

%% Bilevel Optimization Problem

sigma0_1=log(q0_int./(1-q0_int));
sigma0_2=log(q0_int./(1-q0_int));

A = [];
b = [];
Aeq = [];
beq = [];
Std_sS=ones(2,2,n);

fun = @neg_sum_V_fun;
nonlcon = @flexible_circlecon;
options = optimoptions(@fmincon,'Algorithm','interior-point');
for rep1=1:n
    rep2=rep1;
    Std_sS0 = [0.2,0.3;0.2,0.3];
    lb = [Std_s_i(1),q0_int(rep1);Std_s_i(2),q0_int(rep2)];
    ub = [q0_int(rep1),Std_sS_n(1,2);q0_int(rep2),Std_sS_n(2,2)];
    [Std_sS(:,:,rep1),fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    temp_11(rep1)=Std_sS(1,1,rep1);
    temp_12(rep1)=Std_sS(1,2,rep1);
    temp_21(rep1)=Std_sS(2,1,rep1);
    temp_22(rep1)=Std_sS(2,2,rep1);
    Uhat_1(rep1)=-neg_U_fun(Std_sS(:,:,rep1),1);
    Uhat_2(rep1)=-neg_U_fun(Std_sS(:,:,rep1),2);
    Vhat_1(rep1)=-neg_V_fun(Std_sS(:,:,rep1),1);
    Vhat_2(rep1)=-neg_V_fun(Std_sS(:,:,rep1),2);
end

nonlcon = @fixed_circlecon;
for rep1=1:n
    rep2=rep1;
    Std_sS0 = [0.4,0.5;0.4,0.5];
    lb = [Std_s_i(1),q0_int(rep1);Std_s_i(2),q0_int(rep2)];
    ub = [q0_int(rep1),Std_sS_n(1,2);q0_int(rep2),Std_sS_n(2,2)];
    [Std_sS(:,:,rep1),fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    temp(rep1)=Std_sS(2,2,rep1);
    Ubar_1(rep1)=-neg_U_fun(Std_sS(:,:,rep1),1);
    Ubar_2(rep1)=-neg_U_fun(Std_sS(:,:,rep1),2);
    Vbar_1(rep1)=-neg_V_fun(Std_sS(:,:,rep1),1);
    Vbar_2(rep1)=-neg_V_fun(Std_sS(:,:,rep1),2);
end


% for i=1:n
%     temp_11(i)=Std_sS(1,1,i);
%     temp_12(i)=Std_sS(1,2,i);
%     temp_21(i)=Std_sS(2,1,i);
%     temp_22(i)=Std_sS(2,2,i);
% end

%%

figure;
plot(Std_Ufoc_s(1,:),Std_Ufoc_S,'-','LineWidth',2);
% hold on;
% plot(Std_Ufoc_s(2,:),Std_Ufoc_S,'-','LineWidth',2);
% plot(Std_Vfoc_s,Std_Vfoc_S(1,:),'-','LineWidth',2);
% plot(Std_Vfoc_s,Std_Vfoc_S(2,:),'-','LineWidth',2);
% plot(temp_11,temp_12,'*');
% plot(temp_21,temp_22,'*');
% plot(q0_int,temp_11,'b:','LineWidth',2);
hold on;
% plot(q0_int,temp_21,'r:','LineWidth',2);
plot(q0_int,temp_12,'b','LineWidth',3);
plot(q0_int,temp_22,'r:','LineWidth',3);
plot(q0_int,temp,'m--','LineWidth',3);
% '$b^1(\hat{A}^1)$','$b^2(\hat{A}^2)$',
h=legend('$\hat{A}^1$','$\hat{A}^2$','$\bar{A}$');
set(h,'interpreter','Latex','FontSize',14);
%'S_C^e=B(s_C)','S_R^e=B(s_R)',
xlim([0 1]);
ylim([0.45 0.85]);
xlabel('Initial belief of both markets (q_0^1=q_0^2)');
ylabel('Optimal approval standard');


figure;
plot(q0_int,Uhat_1,'b','LineWidth',3);
hold on;
plot(q0_int,Uhat_2,'r','LineWidth',3);
plot(q0_int,Ubar_1,'b:','LineWidth',3);
plot(q0_int,Ubar_2,'r:','LineWidth',3);
h=legend('$\hat{U}^1$','$\hat{U}^2$','$\bar{U}^1$','$\bar{U}^2$');
set(h,'interpreter','Latex','FontSize',14);
xlim([0 1]);
ylim([0 2.5]);
xlabel('Initial belief of both markets (q_0^1=q_0^2)');
ylabel('Pharmaceutical companies payoff');

figure;
plot(q0_int,Vhat_1,'b','LineWidth',3);
hold on;
plot(q0_int,Vhat_2,'r','LineWidth',3);
plot(q0_int,Vbar_1,'b:','LineWidth',3);
plot(q0_int,Vbar_2,'r:','LineWidth',3);
h=legend('$\hat{V}^1$','$\hat{V}^2$','$\bar{V}^1$','$\bar{V}^2$');
set(h,'interpreter','Latex','FontSize',14);
xlim([0 1]);
ylim([0 0.6]);
xlabel('Initial belief of both markets (q_0^1=q_0^2)');
ylabel('FDA payoff (social welfare)');

%% Plotting U* and V* as a function of q_0

% U_opt_val=zeros(2,n);
% for pharma=1:2
%     for rep1=1:n
%         rep2=rep1;
%         U_opt_val(pharma,rep1)=-neg_U_fun(Std_sS(:,:,rep1),pharma);
%         V_opt_val(pharma,rep1)=-neg_V_fun(Std_sS(:,:,rep1),pharma);
%     end
% end

% figure;
% plot(q0_int,U_opt_val(1,:),'LineWidth',2);
% hold on;
% plot(q0_int,U_opt_val(2,:),'LineWidth',2);
% xlabel('q_0');
% ylabel('U*');
% legend('U*_C(q_0,b_C(S^e_C),S^e_C)','U*_R(q_0,b_R(S^e_R),S^e_R)');

% figure;
% plot(q0_int,V_opt_val(1,:),'LineWidth',2);
% hold on;
% plot(q0_int,V_opt_val(2,:),'LineWidth',2);
% xlabel('q_0');
% ylabel('V*');
% legend('V*_C(q_0,b_C(S^e_C),S^e_C)','V*_R(q_0,b_R(S^e_R),S^e_R)');

% figure;
% plot(q0_int,temp_12,'LineWidth',2);
% hold on;
% plot(q0_int,temp_22,'LineWidth',2);
% xlabel('q_0');
% ylabel('S^e');
% legend('S_C^e(q_0)','S_R^e(q_0)');