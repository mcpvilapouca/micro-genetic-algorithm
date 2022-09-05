function [bin,int_bin] =real_to_bin(x1,dom,p,nvar,dim)

%Converts real to binary
for i=1:nvar
%dim=dom(i,2)-dom(i,1)    %Domain dimension

k(i)=round(log2(dim(i,1)*10^p(i,1)));

int_bin(i,1)=round((x1(i,1)-dom(i,1))/(dim(i,1)/((2^k(i))-1)));

aux=int_bin(i,1);

aux1=dec2bin(aux)-'0';

if size(aux1,2)<k(i)
    var=k(i)-size(aux1,2);
    var1=zeros(1,var);
    aux1=horzcat(var1,aux1);
end

if i==1
    bin=aux1;
else
    bin=horzcat(bin,aux1);

end

end

