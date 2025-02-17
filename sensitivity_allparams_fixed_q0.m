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
global c_G;
global c_B;
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

fixed_q0_1=0.6;
fixed_q0_2=0.5;

flag=0; %0: c_i^2 is changing, 1: \nu_i^2 is changing, 3: c_e^2 is changing, and 4: \nu_e^2 is changing 

nu_i_list=0.2:0.1:10;
c_i_list=15;

nu_e_list=1.2:0.01:2;
c_e_list=1:0.1:5;

if flag==1
    itNum=length(nu_i_list);
elseif flag==0
    itNum=length(c_i_list);
elseif flag==3
    itNum=length(c_e_list);
elseif flag==4
    itNum=length(nu_e_list);
end

Std_sS=ones(2,2,itNum);
qhat=zeros(itNum,2);

for it=1:itNum
    
    if flag==1
        c_B=[1.25,1];
        nu_G=[1,1.25];
        nu_i=[2.5,nu_i_list(it)];
        C=[10,15];
    elseif flag==0
        c_B=[1.25,1];
        nu_G=[1,1.25];
        nu_i=[2.5,2];
        C=[10,c_i_list(it)];
    elseif flag==3
        c_B=[1.25,c_e_list(it)];
        nu_G=[1,1.25];
        nu_i=[2.5,2];
        C=[10,15];
    elseif flag==4
        c_B=[1.25,1];
        nu_G=[1,nu_e_list(it)];
        nu_i=[2.5,2];
        C=[10,15];
    end
    
    r=5;
    
    mu=[9,9];
    rho=[4,4];
    mu_prime=2*mu.^2./rho.^2;
    
    R1=0.5*(1-sqrt(1+4*r./mu_prime));
    R2=0.5*(1+sqrt(1+4*r./mu_prime));
    
    rep=1;
    qhat(it,:)=(nu_B-c_B)./(nu_B-c_B-nu_G-c_G);
    
    Std_Ufoc_S=qhat(it,1);
    fun=@root_foc_dUds_1;
%     if c_i_list(it)==0
%         Std_s0 = 0.2;
%     else
        Std_s0 = Std_Ufoc_S/2;
%     end
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
    
    if min(q0_int)<=fixed_q0_2 & max(q0_int)>=fixed_q0_2
%         q0_int=fixed_q0;
    else
%          dbstop();
    end
    n=length(fixed_q0_2);
    
% Bilevel Optimization Problem
    
    sigma0_1=log(fixed_q0_1./(1-fixed_q0_1));
    sigma0_2=log(fixed_q0_2./(1-fixed_q0_2));
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    fun = @neg_sum_V_fun;
    nonlcon = @flexible_circlecon;
    options = optimoptions(@fmincon,'Algorithm','interior-point');
    for rep1=1:n
        rep2=rep1;
        Std_sS0 = [0.3,0.5;0.3,0.5];
        lb = [Std_s_i(1),fixed_q0_1;Std_s_i(2),fixed_q0_2];
        ub = [fixed_q0_1,Std_sS_n(1,2);fixed_q0_2,Std_sS_n(2,2)];
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
            nu_i=[2.5,nu_i_list(it)];
        elseif flag==0
            C=[10,c_i_list(it)];
        elseif flag==3
            c_B=[1.25,c_e_list(it)];
        elseif flag==4
            nu_G=[1,nu_e_list(it)];
        end
        U_opt_val(pharma,it)=-neg_U_fun(Std_sS(:,:,it),pharma);
        V_opt_val(pharma,it)=-neg_V_fun(Std_sS(:,:,it),pharma);
    end
end

if flag==1
    
    figure;
    plot(nu_i_list,temp_21,'LineWidth',2);
    xlabel('\nu');
    ylabel('$$\hat{r}$$','Interpreter','Latex');

    figure;
    plot(nu_i_list,temp_22,'LineWidth',2);
    xlabel('\nu');
    ylabel('$$\hat{A}$$','Interpreter','Latex');
    
    figure;
    plot(nu_i_list,U_opt_val(2,:),'LineWidth',2);
    xlabel('\nu');
    ylabel('$$\hat{U}$$','Interpreter','Latex');
    
    figure;
    plot(nu_i_list,V_opt_val(2,:),'LineWidth',2);
    xlabel('\nu');
    ylabel('$$\hat{V}$$','Interpreter','Latex');
    
elseif flag==0
    
    figure;
    plot(c_i_list,temp_21,'LineWidth',2);
    xlabel('c');
    ylabel('$$\hat{r}$$','Interpreter','Latex');
    
    figure;
    plot(c_i_list,temp_22,'LineWidth',2);
    xlabel('c');
    ylabel('$$\hat{A}$$','Interpreter','Latex');
    
    figure;
    plot(c_i_list,U_opt_val(2,:),'LineWidth',2);
    xlabel('c');
    ylabel('$$\hat{U}$$','Interpreter','Latex');
    
    figure;
    plot(c_i_list,V_opt_val(2,:),'LineWidth',2);
    xlabel('c');
    ylabel('$$\hat{V}$$','Interpreter','Latex');
    
elseif flag==3
    
%     figure;
%     plot(c_e_list,temp_21,'LineWidth',2);
%     xlabel('c_e');
%     ylabel('$$\hat{r}$$','Interpreter','Latex');
%     
    figure;
    plot(c_e_list,temp_12,'LineWidth',2);
    hold on;
    plot(c_e_list,temp_22,'LineWidth',2);
    xlabel('c_e');
    ylabel('$$\hat{A}$$','Interpreter','Latex');
    legend("m=1","m=2");
    
    figure;
    plot(c_e_list,U_opt_val(1,:),'LineWidth',2);
    hold on;
    plot(c_e_list,U_opt_val(2,:),'LineWidth',2);
    xlabel('c_e');
    ylabel('$$\hat{U}$$','Interpreter','Latex');
    legend("m=1","m=2");
    
    figure;
    plot(c_e_list,V_opt_val(1,:),'LineWidth',2);
    hold on;
    plot(c_e_list,V_opt_val(2,:),'LineWidth',2);
    xlabel('c_e');
    ylabel('$$\hat{V}$$','Interpreter','Latex');
    legend("m=1","m=2");
    
elseif flag==4
    
%     figure;
%     plot(nu_e_list,temp_21,'LineWidth',2);
%     xlabel('\nu_e');
%     ylabel('$$\hat{r}$$','Interpreter','Latex');
% 
    figure;
    plot(nu_e_list,temp_12,'LineWidth',2);
    hold on;
    plot(nu_e_list,temp_22,'LineWidth',2);
    xlabel('\nu_e');
    ylabel('$$\hat{A}$$','Interpreter','Latex');
    legend("m=1","m=2");
    
    figure;
    plot(nu_e_list,U_opt_val(1,:),'LineWidth',2);
    hold on;
    plot(nu_e_list,U_opt_val(2,:),'LineWidth',2);
    xlabel('\nu_e');
    ylabel('$$\hat{U}$$','Interpreter','Latex');
    legend("m=1","m=2");
    
    figure;
    plot(nu_e_list,V_opt_val(1,:),'LineWidth',2);
    xlabel('\nu_e');
    hold on;
    plot(nu_e_list,V_opt_val(2,:),'LineWidth',2);
    ylabel('$$\hat{V}$$','Interpreter','Latex');
    legend("m=1","m=2");
    
end


%% solve first order conditions dU/ds=0 and dV/dS=0
C=[10,15]; %run with flag=0

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


figure;
plot(Std_Ufoc_S,Std_Ufoc_s(1,:),'-','LineWidth',2.5);
hold on;
plot(Std_Ufoc_S,Std_Ufoc_s(2,:),'--','LineWidth',2.5);
xlabel({'\textbf{\emph{r}}'},'Interpreter','latex');
ylabel({'\textbf{\emph{A}}'},'Interpreter','latex');
legend({'\textbf{\emph{m=1}}','\textbf{\emph{m=2}}'},'Interpreter','latex');

