function foc = root_nash_2(Std_sS)

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

a=(R1(2)-R2(2))/(exp(-R1(2)*(S-s))-exp(-R2(2)*(S-s)));
b=(R2(2)*exp(R2(2)*(S-s))-R1(2)*exp(R1(2)*(S-s)))/(exp(R2(2)*(S-s))-exp(R1(2)*(S-s)));
f=(R1(2)*exp(-R1(2)*(S-s))-R2(2)*exp(-R2(2)*(S-s)))/(exp(-R1(2)*(S-s))-exp(-R2(2)*(S-s)));
g=-a/exp(S-s);

foc(1)=a*(1+exp(-S))*(nu_i(2)+C(2)/r)+C(2)/r*(b*(1+exp(-s))-exp(-s));
foc(2)=f*(nu_G(2)-exp(-S)*c_B(2))+exp(-S)*c_B(2)+g*(exp(-S)*nu_B(2)-c_G(2));

end