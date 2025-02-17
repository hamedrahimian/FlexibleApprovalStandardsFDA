function we = fnu(ss,SS)

s=log(ss/(1-ss));
S=log(SS/(1-SS));

global r;
global nu_i;
global nu_G;
global nu_B;
global c_G;
global c_B;
global C;

global R1;
global R2;

f=(R1(2)*exp(-R1(2)*(S-s))-R2(2)*exp(-R2(2)*(S-s)))/(exp(-R1(2)*(S-s))-exp(-R2(2)*(S-s)));

we=f*(nu_G(2)-c_B(2)*exp(-S))+c_B(2)*exp(-S);
% we=nu_G(2)-c_B(2)*exp(-S);

end