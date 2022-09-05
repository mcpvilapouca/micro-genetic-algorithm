function [xr] = descod_ind(nvar,ind,dom,dim,p)

for i=1:nvar
k(i)=round(log2(dim(i,1)*10^p(i,1)));  %number of bits necessary for each var
end

fim=0;
ini=0;
for i=1:nvar
fim=fim+k(i);
ini=fim-k(i)+1;
ind1=ind(ini:fim);
int_bin=bin2dec(num2str(ind1));

if i==1
xr = bin_to_real(int_bin,dom(i,:),p(i,1));
else
xr=horzcat(xr,bin_to_real(int_bin,dom(i,:),p(i,1)));
end
end

end

