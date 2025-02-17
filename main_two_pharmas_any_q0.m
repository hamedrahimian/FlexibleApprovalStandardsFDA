clc;
clear;
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

% global v_iso_value;
% global u_iso_value;
% global pharma;
% global Std_S4V;
% global Std_s4U;

% qhat=0.5
c_G=[0,0];
c_B=[1.25,1];

nu_G=[1,1.25];          
nu_B=[0,0];

% qhat=0.3
% c_G=[0.1,0.1];
% c_B=[0.5,0.5];
% 
% nu_G=[1,1];          
% nu_B=[0.2,0.2];

% qhat=0.3
% c_G=[23,23];
% c_B=[0.5,0.5];
% 
% nu_G=[1,1];          
% nu_B=[10,10];
% 
% % qhat=0.3
% c_G=[0.77,0.77];
% c_B=[0.5,0.5];
% 
% nu_G=[1,1];          
% nu_B=[0.4,0.4];

q0=[0.6,0.5];

r=5;
nu_i=[2.5,2];
C=[10,15];

mu=[9,9];
rho=[4,4];
mu_prime=2*mu.^2./rho.^2;
%mu_prime=12;

R1=0.5*(1-sqrt(1+4*r./mu_prime));
R2=0.5*(1+sqrt(1+4*r./mu_prime));

rep=1;
% qhat=c_B./(c_B-nu_G);
qhat=(nu_B-c_B)./(nu_B-c_B-nu_G+c_G);

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
% 

for rep=1:99
    fun=@root_foc_dVdS_2;
    if rep<25
        Std_S0 = Std_Vfoc_s(rep)+0.4;
    elseif rep<33
        Std_S0 = Std_Vfoc_s(rep)+0.5;
    else
        Std_S0 = Std_Vfoc_s(rep)+0.005;
    end
    Std_Vfoc_S(2,rep) = fzero(fun,Std_S0);
end
if mod(qhat(2)*100,1) == 0
    Std_Vfoc_S(2,qhat(2)*100)=qhat(2);
end

% Setting the q0 to plot U, V, and isopayoffs

rep1=round((q0(1)-min(q0_1)+0.01)*100);
rep2=round((q0(2)-min(q0_2)+0.01)*100);

%% solve U and V isopayoff

% v_iso_value=[0.0878,0.3164];
% v_iso_value=[0.05,0.17];
% Std_s4V=[];
% for pharma=1:2
%     for Std_S_pt=1:99
%         Std_S4V=Std_S_pt/100;
%         fun=@V_iso_fun;
%         Std_s0 = Std_S4V/2;
%         Std_s4V(pharma,Std_S_pt) = fzero(fun,Std_s0);
%     end
% end

% u_iso_value=[1,0.5];
% Std_S4U=[];
% 
% for pharma=1:2
%     for Std_s_pt=1:((q0*100)-0.01)
%         Std_s4U=Std_s_pt/100;
%         fun=@U_iso_fun;
%         if Std_s_pt<30
%             Std_S0 = Std_s4U+0.4;
%         else
%             Std_S0 = Std_s4U+0.2;
%         end
%         Std_S4U(pharma,Std_s_pt) = fzero(fun,Std_S0);
%     end
% end

% figure;
% plot(Std_Ufoc_s(1,:),Std_Ufoc_S,'-','LineWidth',2);
% hold on;
% plot(Std_Ufoc_s(2,:),Std_Ufoc_S,'-','LineWidth',2);
% % plot(0.01:0.01:q0-0.01,Std_S4U(1,:),'-','LineWidth',2);
% % plot(0.01:0.01:q0-0.01,Std_S4U(2,:),'-','LineWidth',2);
% % plot(Std_s4V(1,:),0.01:0.01:0.99,'-','LineWidth',2);
% % plot(Std_s4V(2,:),0.01:0.01:0.99,'-','LineWidth',2);
% legend('dU_C/ds=0','dU_R/ds=0','U_C isopayoff','U_R isopayoff','V_C isopayoff','V_R isopayoff');
% xlabel('s');
% ylabel('S');

%% Plotting U

% % s_pts=0.01:0.005:0.99;
% s_pts=0.07:0.005:0.4;
% S_pts=0.4:0.005:0.5;
% % S_pts=0.5*ones(1,length(s_pts));
% 
% s_ptsNum=length(s_pts);
% S_ptsNum=length(S_pts);
% U_C=zeros(s_ptsNum,S_ptsNum);
% d2U_Cds2=zeros(s_ptsNum,S_ptsNum);
% 
% sign_C=[];
% for s_counter=1:s_ptsNum
%     for S_counter=1:S_ptsNum
%         Std_sS=[s_pts(s_counter),S_pts(S_counter);s_pts(s_counter),S_pts(S_counter)];
%         U_C(s_counter,S_counter) = -neg_U_fun(Std_sS,1);
%         sign_C=[sign_C;signcheck(s_pts(s_counter),S_pts(S_counter))];
%     end
% end
% 
% figure;
% mesh(S_pts,s_pts,U_C,'LineWidth',2);
% grid on;
% hold on;
% plot(Std_Ufoc_S,Std_Ufoc_s(1,:),'-','LineWidth',2);
% % xlim([q0 1]);
% % ylim([0 q0]);
% xlabel('S');
% ylabel('s');
% zlabel('U_C');
% 
% s_pts=0:0.005:q0;
% % S_pts=q0+0.01:0.005:0.99;
% S_pts=0.5*ones(1,length(s_pts));
% 
% s_ptsNum=length(s_pts);
% S_ptsNum=length(S_pts);
% U_R=zeros(s_ptsNum,S_ptsNum);
% d2U_Rds2=zeros(s_ptsNum,S_ptsNum);
% 
% sign_R=[];
% for s_counter=1:s_ptsNum
%     for S_counter=1:S_ptsNum
%         Std_sS=[s_pts(s_counter),S_pts(S_counter);s_pts(s_counter),S_pts(S_counter)];
%         U_R(s_counter,S_counter) = -neg_U_fun(Std_sS,2);
%         sign_R=[sign_R;signcheck(s_pts(s_counter),S_pts(S_counter))];
%     end
% end
% 
% figure;
% mesh(S_pts,s_pts,U_R,'LineWidth',2);
% grid on;
% hold on;
% plot(Std_Ufoc_S,Std_Ufoc_s(2,:),'-','LineWidth',2);
% % xlim([q0 1]);
% % ylim([0 q0]);
% xlabel('S');
% ylabel('s');
% zlabel('U_R');
% 
% %% Plotting V
% 
% s_pts=0.01:0.005:q0-0.01;
% % s_pts=0.3*ones(1,length(s_pts));
% S_pts=q0+0.01:0.005:0.99;
% 
% s_ptsNum=length(s_pts);
% S_ptsNum=length(S_pts);
% V_C=zeros(s_ptsNum,S_ptsNum);
% 
% for s_counter=1:s_ptsNum
%     for S_counter=1:S_ptsNum
%         Std_sS=[s_pts(s_counter),S_pts(S_counter);s_pts(s_counter),S_pts(S_counter)];
%         V_C(s_counter,S_counter) = -neg_V_fun(Std_sS,1);
%     end
% end
% 
% figure;
% mesh(S_pts,s_pts,V_C,'LineWidth',2);
% hold on;
% plot(Std_Ufoc_S,Std_Ufoc_s(1,:),'-','LineWidth',2);
% plot(Std_Vfoc_S(1,:),Std_Vfoc_s,'-','LineWidth',2);
% xlim([q0 1]);
% ylim([0 q0]);
% % zlim([-5 5]);
% legend('V_C(sigma_0)','dU_C/ds=0','dV_R/dS=0');
% xlabel('S');
% ylabel('s');
% zlabel('V_C');
% 
% s_pts=0.01:0.005:q0-0.01;
% % s_pts=0.39*ones(1,length(s_pts));
% S_pts=q0+0.01:0.005:0.99;
% 
% s_ptsNum=length(s_pts);
% S_ptsNum=length(S_pts);
% V_R=zeros(s_ptsNum,S_ptsNum);
% 
% for s_counter=1:s_ptsNum
%     for S_counter=1:S_ptsNum
%         Std_sS=[s_pts(s_counter),S_pts(S_counter);s_pts(s_counter),S_pts(S_counter)];
%         V_R(s_counter,S_counter) = -neg_V_fun(Std_sS,2);
%     end
% end
% 
% figure;
% mesh(S_pts,s_pts,V_R,'LineWidth',2);
% hold on;
% plot(Std_Ufoc_S,Std_Ufoc_s(2,:),'-','LineWidth',2);
% plot(Std_Vfoc_S(2,:),Std_Vfoc_s,'-','LineWidth',2);
% xlim([q0 1]);
% ylim([0 q0]);
% % zlim([-5 5]);
% legend('V_R(sigma_0)','dU_R/ds=0','dV_R/dS=0');
% xlabel('S');
% ylabel('s');
% zlabel('V_R');

%%

figure;
plot(Std_Ufoc_S,Std_Ufoc_s(1,:),'LineWidth',2);
hold on;
plot(Std_Ufoc_S,Std_Ufoc_s(2,:),'--','LineWidth',2);
legend('m=1','m=2');
xlim([0 1]);
ylim([0 1]);
ylabel('r=b(A)');
xlabel('A');


%% Plotting V for dU/ds=0

V_CofS=[];
% C_signch=[];
for S_counter=ceil((max(q0(1),qhat(1)))*100):floor(Std_sS_n(1,2)*100)
    if Std_Ufoc_s(1,S_counter)<q0(1) & Std_Ufoc_s(1,S_counter)>Std_Ufoc_s(1,ceil(qhat(1)*100))
        Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        V_CofS=[V_CofS;[Std_Ufoc_S(S_counter),-neg_V_fun(Std_sS,1)]];
%         C_signch=[C_signch;fnu(Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter))];
    end
end
V_RofS=[];
% R_signch=[];
for S_counter=ceil((max(q0(2),qhat(2)))*100):floor(Std_sS_n(2,2)*100)
    if Std_Ufoc_s(2,S_counter)<q0(2) & Std_Ufoc_s(2,S_counter)>Std_Ufoc_s(2,ceil(qhat(2)*100))
        Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        V_RofS=[V_RofS;[Std_Ufoc_S(S_counter),-neg_V_fun(Std_sS,2)]];
%         R_signch=[R_signch;fnu_2(Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter))];
    end
end

figure;
plot(V_CofS(:,1),V_CofS(:,2),'LineWidth',2);
hold on;
plot(V_RofS(:,1),V_RofS(:,2),'--','LineWidth',2); 
xlabel('A');
ylabel('V(b(A),A)');
ylim([0.1 0.27]);
legend('m=1','m=2');

%% Plotting U for dU/ds=0

rep1=round((q0-min(q0_1)+0.01)*100);
rep2=round((q0-min(q0_2)+0.01)*100);

U_CofS=[];
for S_counter=ceil((max(q0(1),qhat(1)))*100):floor(Std_sS_n(1,2)*100)
    if Std_Ufoc_s(1,S_counter)<q0(1) & Std_Ufoc_s(1,S_counter)>Std_Ufoc_s(1,ceil(qhat(1)*100))
        Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        U_CofS=[U_CofS;[Std_Ufoc_S(S_counter),-neg_U_fun(Std_sS,1)]];
    end
end
U_RofS=[];
for S_counter=ceil((max(q0(2),qhat(2)))*100):floor(Std_sS_n(2,2)*100)
    if Std_Ufoc_s(2,S_counter)<q0(2) & Std_Ufoc_s(2,S_counter)>Std_Ufoc_s(2,ceil(qhat(2)*100))
        Std_sS=[Std_Ufoc_s(1,S_counter),Std_Ufoc_S(S_counter);Std_Ufoc_s(2,S_counter),Std_Ufoc_S(S_counter)];
        U_RofS=[U_RofS;[Std_Ufoc_S(S_counter),-neg_U_fun(Std_sS,2)]];
    end
end

figure;
plot(U_CofS(:,1),U_CofS(:,2),'LineWidth',2);
hold on;
plot(U_RofS(:,1),U_RofS(:,2),'--','LineWidth',2); 
xlabel('A');
ylabel('U(b(A),A)');
ylim([0.5 3.1]);
legend('m=1','m=2');

%% Bilevel Optimization Problem

% q0_1=q0_int;
% q0_2=q0_int;
% n1=n;
% n2=n;
% sigma0_1=log(q0_1./(1-q0_1));
% sigma0_2=log(q0_2./(1-q0_2));
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% Std_sS=ones(2,2,n1,n2);
% 
% fun = @neg_sum_V_fun;
% nonlcon = @circlecon;
% options = optimoptions(@fmincon,'Algorithm','interior-point');
% for rep1=1:n1
%     for rep2=1:n2
%         Std_sS0 = [0.2,0.3;0.2,0.3];
% %         if rep2>9 & rep2<16
% %             Std_sS0 = [0.3,0.4;0.3,0.4];
% %         end
%         lb = [Std_s_i(1),q0_1(rep1);Std_s_i(2),q0_2(rep2)];
%         ub = [q0_1(rep1),Std_sS_n(1,2);q0_2(rep2),Std_sS_n(2,2)];
%         [Std_sS(:,:,rep1,rep2),fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%     end
% end
% 
% temp11=zeros(n1,n2);
% temp12=zeros(n1,n2);
% temp21=zeros(n1,n2);
% temp22=zeros(n1,n2);
% for rep1=1:n1
%     for rep2=1:n2
%         temp_11(rep1,rep2)=Std_sS(1,1,rep1,rep2);
%         temp_12(rep1,rep2)=Std_sS(1,2,rep1,rep2);
%         temp_21(rep1,rep2)=Std_sS(2,1,rep1,rep2);
%         temp_22(rep1,rep2)=Std_sS(2,2,rep1,rep2);
%     end
% end
% 
% %%
% figure;
% mesh(q0_2,q0_1,temp_12);
% hold on;
% mesh(q0_2,q0_1,temp_22);
% xlabel('q_0^R');
% ylabel('q_0^C');
% zlabel('S_C & S_R');
% xlim([min(q0_2) max(q0_2)]);
% ylim([min(q0_1) max(q0_1)]);
% legend('S_C','S_R');
