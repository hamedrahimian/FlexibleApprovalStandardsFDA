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

fixed_q0=0.45;
flag=0; % 1 means R is changing

nu_i_R_list=0.1:0.5:14.1;
c_C_list=1:5:500;

if flag==1
    itNum=length(nu_i_R_list);
else
    itNum=length(c_C_list);
end

Std_sS=ones(2,2,itNum);
qhat=zeros(itNum,2);

for it=1:itNum
    
    if flag==1
        nu_i=[2,nu_i_R_list(it)];
        C=[10,10];
    else
        nu_i=[1.7,1.7];
        C=[c_C_list(it),10];
    end
    
    nu_G=[0.5,0.5];
    nu_B=[-0.5,-0.5];
    
    r=5;
    
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
    if c_C_list(it)<5
        Std_s0 = Std_Ufoc_S/4;
    else
        Std_s0 = Std_Ufoc_S/2;
    end
    Std_s_i(1) = fzero(fun,Std_s0);
    
    Std_Ufoc_S=qhat(it,2);
    fun=@root_foc_dUds_2;
    Std_s0 = Std_Ufoc_S/2;
    Std_s_i(2) = fzero(fun,Std_s0);
    
    Std_sS_n=zeros(2,2);
    
    fun=@root_nash_1;
    Std_sS0 = [0.1,0.2];
    Std_sS_n(1,:) = fsolve(fun,Std_sS0);
    
    fun=@root_nash_2;
    Std_sS0 = [0.2,0.5];
    Std_sS_n(2,:) = fsolve(fun,Std_sS0);
    
    q0_1=(ceil(Std_s_i(1)*100)/100:0.01:floor(Std_sS_n(1,2)*100)/100);
    q0_2=(ceil(Std_s_i(2)*100)/100:0.01:floor(Std_sS_n(2,2)*100)/100);
    n1=length(q0_1);
    n2=length(q0_2);
    
    q0_int=(max(ceil(Std_s_i*100)/100):0.01:min(floor(Std_sS_n(:,2)*100)/100));
    
    if min(q0_int)<=fixed_q0 & max(q0_int)>=fixed_q0
%         q0_int=fixed_q0;
    else
        dbstop();
    end
    [temp,n]=size(fixed_q0);
    
%   %% Bilevel Optimization Problem
    
    sigma0_1=log(fixed_q0./(1-fixed_q0));
    sigma0_2=log(fixed_q0./(1-fixed_q0));
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    fun = @neg_sum_V_fun;
    nonlcon = @flexible_circlecon;
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
            nu_i=[2,nu_i_R_list(it)];
%             C=[10,10];
        else
%             nu_i=[3,3];
            C=[c_C_list(it),10];
        end
        U_opt_val(pharma,it)=-neg_U_fun(Std_sS(:,:,it),pharma);
        V_opt_val(pharma,it)=-neg_V_fun(Std_sS(:,:,it),pharma);
    end
end

if flag==1
    
    figure;
    plot(nu_i_R_list,temp_21,'LineWidth',2);
    xlabel('nu_i');
    ylabel('s^e');
%     legend('s^e(nu_i)');

    figure;
    plot(nu_i_R_list,temp_22,'LineWidth',2);
    xlabel('nu_i');
    ylabel('S^e');
%     legend('S^e(nu_i)');
    
    figure;
    plot(nu_i_R_list,U_opt_val(2,:),'LineWidth',2);
    xlabel('nu_i');
    ylabel('U*');
%     legend('U*(nu_i,S^e(nu_i))');
    
    figure;
    plot(nu_i_R_list,V_opt_val(2,:),'LineWidth',2);
    xlabel('nu_i');
    ylabel('V*');
%     legend('V*(nu_i,S^e(nu_i))');
    
else
    
    figure;
    plot(c_C_list,temp_11,'LineWidth',2);
    xlabel('C');
    ylabel('s^e');
%     legend('s^e(C)');
    
    figure;
    plot(c_C_list,temp_12,'LineWidth',2);
    xlabel('C');
    ylabel('S^e');
%     legend('S^e(C)');
    
    figure;
    plot(c_C_list,U_opt_val(1,:),'LineWidth',2);
    xlabel('C');
    ylabel('U*');
%     legend('U*(C,S^e(C))');
    
    figure;
    plot(c_C_list,V_opt_val(1,:),'LineWidth',2);
    xlabel('C');
    ylabel('V*');
%     legend('V*(C,S^e(C))');
    
end
