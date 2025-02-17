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

q0_1=0.6;
q0_2=0.5;

C_R_list=12:0.1:18;
nu_i_R_list=1.8:0.01:2.2;

C_R_size=length(C_R_list);
nu_i_R_size=length(nu_i_R_list);

fixed_s=zeros(C_R_size,nu_i_R_size,2);
fixed_S=zeros(C_R_size,nu_i_R_size,2);

flexible_s=zeros(C_R_size,nu_i_R_size,2);
flexible_S=zeros(C_R_size,nu_i_R_size,2);

fixed_U_opt_val=zeros(C_R_size,nu_i_R_size,2);
fixed_V_opt_val=zeros(C_R_size,nu_i_R_size,2);

flexible_U_opt_val=zeros(C_R_size,nu_i_R_size,2);
flexible_V_opt_val=zeros(C_R_size,nu_i_R_size,2);

fixed_V_benefit=zeros(C_R_size,nu_i_R_size,2);
fixed_V_cost=zeros(C_R_size,nu_i_R_size,2);

flexible_V_benefit=zeros(C_R_size,nu_i_R_size,2);
flexible_V_cost=zeros(C_R_size,nu_i_R_size,2);

c_G=[0,0];
c_B=[4,0.04];	   % in million dollars - $100K for 1 DALY
%DALY (equal to one lost year of healthy life) in the US (age-standardized) is 1891/100,000 from drug use disorders. 
%19,850 patients have been prescibed Trikafta while 2,125,887 COPD patients have been prescribed Trelegy.
        
nu_G=[0.66,1.65]   % in million dollars -$100K for 1 QALY
%https://thorax.bmj.com/content/73/Suppl_4/A237.1        
%https://icer.org/wp-content/uploads/2020/08/ICER_CF_Evidence_Report_042720.pdf

nu_B=[0,0];

qhat=(nu_B-c_B)./(nu_B-c_B-nu_G-c_G);

mu=[9,9];
rho=[4,4];
mu_prime=2*mu.^2./rho.^2;

r=0.06;				%US 2022 inflation rate

R1=0.5*(1-sqrt(1+4*r./mu_prime));
R2=0.5*(1+sqrt(1+4*r./mu_prime));
        
c_base=8.7;                  % $8667 cost per patient - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9766529/
nu_i_base=40;                % 62.4 billion revenue projected till 2038  -  $https://www.pharmaceutical-technology.com/data-analysis/trelegy-ellipta-cagr-2-22-2038/#:~:text=The%20company%20reported%20revenues%20of,margin%20of%2022.8%25%20in%20FY2020 

c_rare=33.8;                 % $33,831 cost per patient  https://www.whatdotheyknow.com/request/791135/response/1896808/attach/html/3/5609%2040%20Trikafta%20trial%20WGH%20October%202021%20x2.doc.html
nu_i_rare=138.8             % $8.9 billion revenue in 2022 increasing %18 annually until 2038 - https://investors.vrtx.com/news-releases/news-release-details/vertex-reports-fourth-quarter-and-full-year-financial-2022

for C_R_counter=1:C_R_size
    for nu_i_R_counter=1:nu_i_R_size
        
        C=[c_base,C_R_list(C_R_counter)]; 
        nu_i=[nu_i_base,nu_i_R_list(nu_i_R_counter)];
               
        rep=1;
        
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
        
        q0_1_range=(ceil(Std_s_i(1)*100)/100:0.01:floor(Std_sS_n(1,2)*100)/100);
        q0_2_range=(ceil(Std_s_i(2)*100)/100:0.01:floor(Std_sS_n(2,2)*100)/100);
        n1=length(q0_1);
        n2=length(q0_2);
        
        q0_int=(max(ceil(Std_s_i*100)/100):0.01:min(floor(Std_sS_n(:,2)*100)/100));
        
        n=1;
        
        sigma0_1=log(q0_1./(1-q0_1));
        sigma0_2=log(q0_2./(1-q0_2));
        
% %% Bilevel Optimization Problem
        
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
            lb = [Std_s_i(1),q0_1;Std_s_i(2),q0_1];
            ub = [q0_2,Std_sS_n(1,2);q0_2,Std_sS_n(2,2)];
            rng default;
            [flexible_Std_sS(:,:,rep1),flexible_fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,flexible_nonlcon,options);
            [fixed_Std_sS(:,:,rep1),fixed_fval] = fmincon(fun,Std_sS0,A,b,Aeq,beq,lb,ub,fixed_nonlcon,options);
        end
        
        sig_num=4;
        for pharma=1:2
%            fixed_s(C_R_counter,nu_i_R_counter,pharma)=round(fixed_Std_sS(pharma,1,rep1),sig_num);
%            fixed_S(C_R_counter,nu_i_R_counter,pharma)=round(fixed_Std_sS(pharma,2,rep1),sig_num);
            
%            flexible_s(C_R_counter,nu_i_R_counter,pharma)=round(flexible_Std_sS(pharma,1,rep1),sig_num);
%            flexible_S(C_R_counter,nu_i_R_counter,pharma)=round(flexible_Std_sS(pharma,2,rep1),sig_num);
            
%            fixed_U_opt_val(C_R_counter,nu_i_R_counter,pharma)=round(-neg_U_fun(fixed_Std_sS(:,:,rep1),pharma),sig_num);
%            fixed_V_opt_val(C_R_counter,nu_i_R_counter,pharma)=round(-neg_V_fun(fixed_Std_sS(:,:,rep1),pharma),sig_num);
            
            flexible_U_opt_val(C_R_counter,nu_i_R_counter,pharma)=round(-neg_U_fun(flexible_Std_sS(:,:,rep1),pharma),sig_num);
            flexible_V_opt_val(C_R_counter,nu_i_R_counter,pharma)=round(-neg_V_fun(flexible_Std_sS(:,:,rep1),pharma),sig_num);
            
%            fixed_V_benefit(C_R_counter,nu_i_R_counter,pharma)=round(V_fun_reward(fixed_Std_sS(:,:,rep1),pharma),sig_num);
%            fixed_V_cost(C_R_counter,nu_i_R_counter,pharma)=round(V_fun_risk(fixed_Std_sS(:,:,rep1),pharma),sig_num);
            
%            flexible_V_benefit(C_R_counter,nu_i_R_counter,pharma)=round(V_fun_reward(flexible_Std_sS(:,:,rep1),pharma),sig_num);
%            flexible_V_cost(C_R_counter,nu_i_R_counter,pharma)=round(V_fun_risk(flexible_Std_sS(:,:,rep1),pharma),sig_num);
            
        end

    end
end


X=[];
Y=[];
figure;
hold on;
for C_R_counter=1:C_R_size
    for nu_i_R_counter=1:nu_i_R_size
         if fixed_U_opt_val(C_R_counter,nu_i_R_counter,2)==flexible_U_opt_val(C_R_counter,nu_i_R_counter,2)
             X=[X;C_R_list(C_R_counter)];
             Y=[Y;nu_i_R_list(nu_i_R_counter)];
         end    
         if C_R_list(C_R_counter)==(C(1)+min(C_R_list))/2 & nu_i_R_list(nu_i_R_counter)==(nu_i(1)+max(nu_i_R_list))/2
             plot(C_R_list(C_R_counter),nu_i_R_list(nu_i_R_counter),'b+','LineWidth',2);
         elseif C_R_list(C_R_counter)==(C(1)+max(C_R_list))/2 & nu_i_R_list(nu_i_R_counter)==(nu_i(1)+min(nu_i_R_list))/2
             plot(C_R_list(C_R_counter),nu_i_R_list(nu_i_R_counter),'rx','LineWidth',2);
         end
    end
end
plot(X,Y,'g','LineWidth',2);
plot(C(1),nu_i(1),'*','LineWidth',2);
xlabel('c');
ylabel('v');
%xlim([12.5 17.5]);