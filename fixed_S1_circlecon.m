function [c,ceq] = fixed_S1_circlecon(Std_sS)

sS=log(Std_sS./(1-Std_sS));

global R1;
global R2;
global C;
global nu_i;
global r;
global fixed_S1;

S1=log(fixed_S1/(1-fixed_S1));

for i=1:2
a(i)=(R1(i)-R2(i))/(exp(-R1(i)*(sS(i,2)-sS(i,1)))-exp(-R2(i)*(sS(i,2)-sS(i,1))));
b(i)=(R2(i)*exp(R2(i)*(sS(i,2)-sS(i,1)))-R1(i)*exp(R1(i)*(sS(i,2)-sS(i,1))))/(exp(R2(i)*(sS(i,2)-sS(i,1)))-exp(R1(i)*(sS(i,2)-sS(i,1))));
ceq(i) = a(i)*(1+exp(-sS(i,2)))*(nu_i(i)+C(i)/r)+C(i)/r*(b(i)*(1+exp(-sS(i,1)))-exp(-sS(i,1)));
end

ceq(3) = sS(1,2)-S1;
ceq(4) = neg_U_fun(Std_sS,2)-neg_U_fun(Std_sS,1);

c=[];

end