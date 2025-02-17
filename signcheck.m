function summ = signcheck(ss,SS)

s=log(ss/(1-ss));
S=log(SS/(1-SS));

global R1;
global R2;

a=(R1(1)-R2(1))/(exp(-R1(1)*(S-s))-exp(-R2(1)*(S-s)));
b=(R2(1)*exp(R2(1)*(S-s))-R1(1)*exp(R1(1)*(S-s)))/(exp(R2(1)*(S-s))-exp(R1(1)*(S-s)));
db_ds=a^2/exp(S-s);
summ=a*(1+exp(-S))+b^2+db_ds*(1+exp(-s))+(b-1)^2*exp(-s);
end