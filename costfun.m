function A_bar = costfun(A)

%A_common=8667*10^(-6);

r=0.06;

n=3;

m=1055;

A_bar=A*(((1+r/m)^(n*m)-1)/((r/m*(1+r/m)^(n*m))))/((exp(r*n)-1)/(r*exp(r*n)));

end