function [conv_pop] = pop_conve(newpop,iFm)
%Check the convergence of each population
evalpop=newpop;

evalpop(iFm,:)=[];

iFmax=newpop(iFm,:);
count=0;

for j=1:size(evalpop,1)
for i=1:size(evalpop,2)
    if evalpop(j,i)~=iFmax(1,i)
        count=count+1;
    end
end
end
nbits=size(evalpop,1)*size(evalpop,2);
prob_conv=count/nbits;

if prob_conv<0.05
    conv_pop=1;
else
    conv_pop=0;
end

end

