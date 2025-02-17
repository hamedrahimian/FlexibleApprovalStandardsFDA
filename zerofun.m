function s_i = s_i_fun(qhat,pharma)

global Std_Ufoc_S;
global rep;
rep=1;
Std_Ufoc_S=qhat;

Std_s0=0.0001;

if pharma==1
	fun=@root_foc_dUds_1;
elseif pharma==2
	fun=@root_foc_dUds_2;
end


while Std_s0<qhat
	s_i = fzero(fun,Std_s0);
	if (isreal(s_i) & s_i~=qhat) 
		break
	end
	Std_s0=Std_s0+0.005;
end

if Std_s0>qhat
	break
end	

end