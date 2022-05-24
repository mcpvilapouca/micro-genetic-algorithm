function [newpop] = new_generation(elitismo,Fm,pop,npop,nvar,dom,...
    dim,p,obj)

switch elitismo
    case 'y'
   [max_num,max_idx]=max(Fm(:,2));
   newpop(1,:)=pop(max_idx,:);
   nnewpop=1;
    case 'n'
   nnewpop=0;
end

nt=2;   %dimens�o do torneio

while nnewpop<npop
Ft=Fm;
for i=1:nt

%Tournment
 Fm_rand=Ft(randperm(end),:);       %randomize the order
 [max_n,max_i]=max(Fm_rand(1:2,2));  %Find the fittest in the first individuals
 m=Fm_rand(max_i,1);
 parent(i,:)=pop(m,:);             %Define the father as the fittest
 I = find(Ft(:,1)==m);
 Ft(I,:)=[];                        %Eliminates from the tournment the one already selected
end

cross_rand=rand(1,size(parent,2));
%Loop a cada gene para gerar filho
for j=1:size(parent,2)
    if cross_rand(1,j)>0.5
        child(1,j)=parent(2,j);
    else
        child(1,j)=parent(1,j);
    end
end

nnewpop=nnewpop+1;
newpop(nnewpop,:)=child;    %Alocates the child to the new population
end

[Fm_new] = merit(npop,newpop,nvar,dom,dim,p,obj);

disp_res(:,1:size(newpop,2))=newpop;
for i=1:npop
 ind=newpop(i,:);
 [xr] = descod_ind(nvar,ind,dom,dim,p)
disp_xr(i,:)=xr;
space(i,:)='    ';
end

T=[num2str(newpop),space,num2str(disp_xr),space,num2str(Fm_new(:,2))];
disp(T)
end

