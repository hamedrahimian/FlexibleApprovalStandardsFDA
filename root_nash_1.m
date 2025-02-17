function foc = root_nash_1(Std_sS)

s=log(Std_sS(1)/(1-Std_sS(1)));
S=log(Std_sS(2)/(1-Std_sS(2)));

global r;
global nu_i;
global nu_G;
global nu_B;
global c_G;
global c_B;
global C;

global R1;
global R2;

a=(R1(1)-R2(1))/(exp(-R1(1)*(S-s))-exp(-R2(1)*(S-s)));
b=(R2(1)*exp(R2(1)*(S-s))-R1(1)*exp(R1(1)*(S-s)))/(exp(R2(1)*(S-s))-exp(R1(1)*(S-s)));
f=(R1(1)*exp(-R1(1)*(S-s))-R2(1)*exp(-R2(1)*(S-s)))/(exp(-R1(1)*(S-s))-exp(-R2(1)*(S-s)));
g=-a/exp(S-s);

foc(1)=a*(1+exp(-S))*(nu_i(1)+C(1)/r)+C(1)/r*(b*(1+exp(-s))-exp(-s));
foc(2)=f*(nu_G(1)-exp(-S)*c_B(1))+exp(-S)*c_B(1)+g*(exp(-S)*nu_B(1)-c_G(1));

end