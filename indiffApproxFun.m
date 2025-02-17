function [indiff_c,indiff_nu] = indiffApproxFun(Udiff,clist,nulist)

%To use this function run the following:
%[indiff_c,indiff_nu] = indiffApproxFun(fixed_U_opt_val-flexible_U_opt_val,C_R_list,nu_i_R_list)

c_length=length(clist);
nu_length=length(nulist);

indiff_c=clist;
indiff_nu=zeros(size(clist));

for i=1:c_length
	for j=2:(nu_length-1)
		if Udiff(i,j,2)<0 & Udiff(i,j+1,2)<0
			absSum=Udiff(i,j-1,2)+abs(Udiff(i,j,2));
			indiff_nu(i)=Udiff(i,j-1,2)/absSum*nulist(j-1)+abs(Udiff(i,j,2))/absSum*nulist(j);
			break
		end
	end
	
end

end

