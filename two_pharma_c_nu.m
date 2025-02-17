clc;
clear all;
%close all;

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
nu_B=[0,0];



mu=[12,12];
rho=[4,4];
mu_prime=2*mu.^2./rho.^2;
%mu_prime=12;
r=5;

R1=0.5*(1-sqrt(1+4*r./mu_prime));
R2=0.5*(1+sqrt(1+4*r./mu_prime));


c_list=12:1:18;
nu_list=1.4:0.1:2;

% c_list=0.25:0.1:0.75;
% nu_list=0.5:0.1:1;

c_length=length(c_list);
nu_length=length(nu_list);

Std_sS=ones(2,2,c_length,nu_length);

for c_counter=1:c_length
    for nu_counter=1:nu_length
        
        % lower
        c_B=[0.5,0.5];
        nu_G=[0.75,0.75];
        
        C=[15,c_list(c_counter)];
        nu_i=[1.7,nu_list(nu_counter)];
        
        % upper
%                 C=[15,15];
%                 nu_i=[1.7,1.7];
%         
%                 c_B=[0.5,c_list(c_counter)];
%                 nu_G=[0.75,nu_list(nu_counter)];
        
        
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
        
        % solve first order conditions dU/ds=0 and dV/dS=0
        
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
        
        q0=0.45;
        rep1=round((q0-min(q0_1)+0.01)*100);
        rep2=round((q0-min(q0_2)+0.01)*100);
        
        % Bilevel Optimization Problem
        
        % sigma0_1=log(q0_int./(1-q0_int));
        % sigma0_2=log(q0_int./(1-q0_int));
        %
        % rep1=round((q0-min(q0_1)+0.01)*100);
        % rep2=round((q0-min(q0_2)+0.01)*100);
        
        A = [];
        b = [];
        Aeq = [];
        beq = [];
        
        
        fun = @neg_sum_V_fun;
        nonlcon = @flexible_circlecon;
        options = optimoptions(@fmincon,'Algorithm','interior-point');
        
        
        Std_sS0 = [0.2,0.3;0.2,0.3];
        lb = [Std_s_i(1),q0_int(rep1);Std_s_i(2),q0_int(rep2)];
        ub = [q0_int(rep1),Std_sS_n(1,2);q0_int(rep2),Std_sS_n(2,2)];
        [Std_sS(:,:,c_counter,nu_counter),fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    end
end

%% Plotting U* and S* as a function of c and nu

temp_11=zeros(c_length,nu_length);
temp_12=zeros(c_length,nu_length);
temp_21=zeros(c_length,nu_length);
temp_22=zeros(c_length,nu_length);
for c_counter=1:c_length
    for nu_counter=1:nu_length
        temp_11(c_counter,nu_counter)=Std_sS(1,1,c_counter,nu_counter);
        temp_12(c_counter,nu_counter)=Std_sS(1,2,c_counter,nu_counter);
        temp_21(c_counter,nu_counter)=Std_sS(2,1,c_counter,nu_counter);
        temp_22(c_counter,nu_counter)=Std_sS(2,2,c_counter,nu_counter);
    end
end

U_1_opt_val=zeros(c_length,nu_length);
U_2_opt_val=zeros(c_length,nu_length);
for c_counter=1:c_length
    for nu_counter=1:nu_length
        U_1_opt_val(c_counter,nu_counter)=-neg_U_fun(Std_sS(:,:,c_counter,nu_counter),1);
        U_2_opt_val(c_counter,nu_counter)=-neg_U_fun(Std_sS(:,:,c_counter,nu_counter),2);
        %             V_opt_val(pharma,c_counter,nu_counter)=-neg_V_fun(Std_sS(:,:,c_counter,nu_counter),pharma);
    end
end

figure;
mesh(c_list,nu_list,U_1_opt_val,'LineWidth',2);
hold on;
mesh(c_list,nu_list,U_2_opt_val,'LineWidth',2);
xlabel('c');
ylabel('v');
legend('hat(U)^1(c,v)','hat(U)^2(c,v)');


figure;
mesh(c_list,nu_list,temp_12,'LineWidth',2);
hold on;
mesh(c_list,nu_list,temp_22,'LineWidth',2);
xlabel('c');
ylabel('v');
legend('S^1(c,v)','S^2(c,v)');

% figure;
% plot(q0_int,V_opt_val(1,:),'LineWidth',2);
% hold on;
% plot(q0_int,V_opt_val(2,:),'LineWidth',2);
% xlabel('q_0');
% ylabel('V*');
% legend('V*_C(q_0,b_C(S^e_C),S^e_C)','V*_R(q_0,b_R(S^e_R),S^e_R)');

