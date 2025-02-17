clc;
clear;

load ''equal_q0_joint_optimization_qhat_50_33.mat'';

for i=1:n
    temp_11(i)=Std_sS(1,1,i);
    temp_12(i)=Std_sS(1,2,i);
    temp_21(i)=Std_sS(2,1,i);
    temp_22(i)=Std_sS(2,2,i);
end

figure;
plot(Std_Ufoc_s(1,:),Std_Ufoc_S,'-','LineWidth',2);
hold on;
plot(Std_Ufoc_s(2,:),Std_Ufoc_S,'-','LineWidth',2);
% plot(Std_Vfoc_s,Std_Vfoc_S(1,:),'-','LineWidth',2);
% plot(Std_Vfoc_s,Std_Vfoc_S(2,:),'-','LineWidth',2);
plot(temp_11,temp_12,'*');
plot(temp_21,temp_22,'*');
plot(q0_1,temp_11,':','LineWidth',2);
plot(q0_2,temp_21,':','LineWidth',2);
plot(q0_1,temp_12,':','LineWidth',2);
plot(q0_2,temp_22,':','LineWidth',2);

% load ''equal_q0_joint_optimization_qhat_50_33_U1_equal_U2.mat'';
load ''equal_q0_separate_optimization_qhat_50_33.mat'';

for i=1:n
    temp_11(i)=Std_sS(1,1,i);
    temp_12(i)=Std_sS(1,2,i);
    temp_21(i)=Std_sS(2,1,i);
    temp_22(i)=Std_sS(2,2,i);
end

plot(temp_11,temp_12,'*');
plot(temp_21,temp_22,'*');
plot(q0_1,temp_11,':','LineWidth',2);
plot(q0_2,temp_21,':','LineWidth',2);
plot(q0_1,temp_12,':','LineWidth',2);
plot(q0_2,temp_22,':','LineWidth',2);
xlim([0 1]);
ylim([0 1]);
xlabel('s , q_0');
ylabel('S');

%%
% clear;
global rep1;
global rep2;
load ''equal_q0_joint_optimization_qhat_50_33.mat'';
%load ''equal_q0_joint_optimization_qhat_50_33_U1_equal_U2.mat'';

U_opt_val=zeros(2,n);
for pharma=1:2
    for rep1=1:n
        rep2=rep1;
        U_opt_val(pharma,rep1)=-neg_U_fun(Std_sS(:,:,rep1),pharma);
    end
end

figure;
plot(q0_1,U_opt_val(1,:),'LineWidth',2);
hold on;
plot(q0_1,U_opt_val(2,:),'LineWidth',2);
xlabel('q_0');
ylabel('U');
legend('U_C(q_0,b_C(S*_C),S*_C)','U_R(q_0,b_R(S*_R),S*_R)');

load ''equal_q0_separate_optimization_qhat_50_33.mat'';
% load ''equal_q0_joint_optimization_qhat_50_33_U1_equal_U2.mat'';

for pharma=1:2
    for rep1=1:n
        rep2=rep1;
        U_opt_val(pharma,rep1)=-neg_U_fun_3d(Std_sS(:,:,rep1),pharma);
    end
end

plot(q0_1,U_opt_val(1,:),'LineWidth',2);
plot(q0_1,U_opt_val(2,:),'LineWidth',2);
xlabel('q_0');
ylabel('U*');
legend('U_C(q_0,b_C(S*_C),S*_C)','U_R(q_0,b_R(S*_R),S*_R)');

%%
load ''separate.mat'';

