% close all;
clc;

%%

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_nu_e.mat'';
plot(nu_e_list,temp_12,'b','LineWidth',2.5);
hold on;
plot(nu_e_list,temp_22,'--b','LineWidth',2.5);
load ''fixed_nu_e.mat'';
plot(nu_e_list,temp_12,':r','LineWidth',2.5);
xlabel('\nu_e');
ylabel('Equilibrium Approval Standard');
legend({'$\hat{A}^1$','$\hat{A}^2$','$\bar{A}$'},'Interpreter','latex');

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_c_e.mat'';
plot(c_e_list,temp_12,'b','LineWidth',2.5);
hold on;
plot(c_e_list,temp_22,'--b','LineWidth',2.5);
load ''fixed_c_e.mat'';
plot(c_e_list,temp_12,':r','LineWidth',2.5);
xlabel('c_e');
ylabel("Equilibrium Approval Standard");
legend({'$\hat{A}^1$','$\hat{A}^2$','$\bar{A}$'},'Interpreter','latex');

%%

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_nu_e.mat'';
plot(nu_e_list,V_opt_val(1,:),'b','LineWidth',2.5);
hold on;
plot(nu_e_list,V_opt_val(2,:),'--b','LineWidth',2.5);
load ''fixed_nu_e.mat'';
plot(nu_e_list,V_opt_val(1,:),':r','LineWidth',2.5);
plot(nu_e_list,V_opt_val(2,:),'-.r','LineWidth',2.5);
% ylim([0 1.2]);
xlabel('\nu_e');
ylabel("FDA's Equilibrium Payoff");
legend({'$\hat{V}^1$','$\hat{V}^2$','$\bar{V}^1$','$\bar{V}^2$'},'Interpreter','latex');

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_c_e.mat'';
plot(c_e_list,V_opt_val(1,:),'b','LineWidth',2.5);
hold on;
plot(c_e_list,V_opt_val(2,:),'--b','LineWidth',2.5);
load ''fixed_c_e.mat'';
plot(c_e_list,V_opt_val(1,:),':r','LineWidth',2.5);
plot(c_e_list,V_opt_val(2,:),'-.r','LineWidth',2.5);
xlabel('c_e');
ylabel("FDA's Equilibrium Payoff");
legend({'$\hat{V}^1$','$\hat{V}^2$','$\bar{V}^1$','$\bar{V}^2$'},'Interpreter','latex');

%%

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_nu_e.mat'';
plot(nu_e_list,V_opt_val(1,:)+V_opt_val(2,:),'b','LineWidth',2.5);
hold on;
load ''fixed_nu_e.mat'';
plot(nu_e_list,V_opt_val(1,:)+V_opt_val(2,:),':r','LineWidth',2.5);
xlabel('\nu_e');
ylabel("FDA's Total Equilibrium Payoff");
legend({'$\hat{V}^1$+$\hat{V}^2$','$\bar{V}^1$+$\bar{V}^2$'},'Interpreter','latex');

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_c_e.mat'';
plot(c_e_list,V_opt_val(1,:)+V_opt_val(2,:),'b','LineWidth',2.5);
hold on;
load ''fixed_c_e.mat'';
plot(c_e_list,V_opt_val(1,:)+V_opt_val(2,:),':r','LineWidth',2.5);
xlabel('c_e');
ylabel("FDA's Total Equilibrium Payoff");
legend({'$\hat{V}^1$+$\hat{V}^2$','$\bar{V}^1$+$\bar{V}^2$'},'Interpreter','latex');

%%

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_nu_e.mat'';
plot(nu_e_list,U_opt_val(1,:),'b','LineWidth',2.5);
hold on;
plot(nu_e_list,U_opt_val(2,:),'--b','LineWidth',2.5);
load ''fixed_nu_e.mat'';
plot(nu_e_list,U_opt_val(1,:),':r','LineWidth',2.5);
plot(nu_e_list,U_opt_val(2,:),'-.r','LineWidth',2.5);
xlabel('\nu_e');
ylabel("Firms' Equilibrium Payoff");
legend({'$\hat{U}^1$','$\hat{U}^2$','$\bar{U}^1$','$\bar{U}^2$'},'Interpreter','latex');

figure('DefaultAxesFontSize',16);
set(gca,'FontWeight','bold');
load ''flexible_c_e.mat'';
plot(c_e_list,U_opt_val(1,:),'b','LineWidth',2.5);
hold on;
plot(c_e_list,U_opt_val(2,:),'--b','LineWidth',2.5);
load ''fixed_c_e.mat'';
plot(c_e_list,U_opt_val(1,:),':r','LineWidth',2.5);
plot(c_e_list,U_opt_val(2,:),'-.r','LineWidth',2.5);
xlabel('c_e');
ylabel("Firms' Equilibrium Payoff");
legend({'$\hat{U}^1$','$\hat{U}^2$','$\bar{U}^1$','$\bar{U}^2$'},'Interpreter','latex');

