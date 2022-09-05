function [ind] = random_ind(nvar,dim,p)

%cno. of bits for each variable
tot_bits=0;
for i=1:nvar
k(i)=round(log2(dim(i,1)*10^p(i,1)));  %bits for each var
tot_bits=tot_bits+k(i);
end

%Creates a random individual
for i=1:nvar
aux=rand(1,k(i));
for j=1:size(aux,2)
    if aux(j)>0.5
        aux(j)=1;
    else
        aux(j)=0;
    end
end

if i==1
    ind=aux;
else
    ind=horzcat(ind,aux);
end

end

end

