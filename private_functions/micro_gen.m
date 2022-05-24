function [xfinal,Fo,maxFm,restart,nger,meanFm] = micro_gen(obj,...
    elitismo,nvar,dim,dom,p,npop,ger_max,conv_iter,tol,nrestart,...
    h0,initialpop)

%calculates the number of required bits for each parameter
tot_bits=0;
for i=1:nvar
k(i)=round(log2(dim(i,1)*10^p));  %bits required for each var
tot_bits=tot_bits+k(i);           %total bits of an individual
end

switch h0
case 'yes'
%Randomly start the initial population
[pop,space,disp_xri] = random_pop(nvar,dim,p,npop,dom);
case 'no'
[pop,space,disp_xri] = random_pop(nvar,dim,p,npop,dom);
[bin,int_bin] =real_to_bin(initialpop,dom,p,nvar,dim);
pop(1,:)=bin;
disp_xri(1,:)=initialpop';
end

ger=0;
%Obtain the merit of the initial population
[Fm0,Fev0] = merit(npop,pop,nvar,dom,dim,p,obj);
%Display to the Command Window
disp(strcat('Generation ',num2str(ger)))
T0=[num2str(pop),space,num2str(disp_xri),space,num2str(Fm0(:,2))];
disp(T0)

restart=0;

while ger<ger_max

if ger==0
Fm=Fm0;
else
%Evaluates the population merit
[Fm,Fev] = merit(npop,pop,nvar,dom,dim,p,obj);
end

ger=ger+1;
disp(strcat('Generation ',num2str(ger)))

%Defines a new generation
[newpop] = new_generation(elitismo,Fm,pop,npop,nvar,dom,dim,p,obj);

%Evaluate the merit of the new population
[Fm_new,Fev_new] = merit(npop,pop,nvar,dom,dim,p,obj);

%saves the generation, the median merit and the maximum merit of the
%population
nger(ger,1)=ger;
meanFm(ger,1)=mean(Fm_new(:,2));
[maxFm(ger,1),iFm]=max(Fm_new(:,2));

%saves the value of the objective function for the parameters that had
%maximum merit of the population

iobj=newpop(iFm,:);
[xobj] = descod_ind(nvar,iobj,dom,dim,p);
[Fobj] = Fev_new(iFm,2);
Fo(ger,1)=Fobj;

%Check global convergence
if (ger>1) && (abs(maxFm(ger,1)-maxFm(ger-1,1))<tol)
   n_it=n_it+1;
else
    n_it=0;
end
if (n_it>=conv_iter) && (restart>nrestart)
    disp(' ')
    disp(strcat('Converged at the generation= ',num2str(ger)))
    break
end
%Check convergence of the population
[conv_pop] = pop_conve(newpop,iFm);

%If the population converges, generate a new random population
%keeping only the individual with the maximum merit

if conv_pop==1&&ger<ger_max
    disp('*******************Restart micro-population*******************')
    restart=restart+1;
    pop=zeros(npop,size(newpop,2));

    %Evaluate the merit of the population to keep the most fit
    [Fm_new,Fev_new] = merit(npop,pop,nvar,dom,dim,p,obj);
    [max_num,max_new]=max(Fm_new(:,2));
     pop1(1,:)=newpop(max_new,:);
    [xrfit] = descod_ind(nvar,pop1,dom,dim,p);

     %Randomly generate the remainining population elements
     npop_cont=npop-1;
     [pop_cont,spacer,xrr] = random_pop(nvar,dim,p,npop_cont,dom);
     pop=vertcat(pop1,pop_cont);
     disp_xrr=vertcat(xrfit,xrr);

     %Evaluate the merit of the restarted population
     [Fm_restart,Fev_restart] = merit(npop,pop,nvar,dom,dim,p,obj);

     %Display of the new generation at the command window
     ger=ger+1;
     disp(strcat('Generation ',num2str(ger)))
     Tr=[num2str(pop),space,num2str(disp_xrr),space,num2str(Fm_restart(:,2))];
     disp(Tr)

     %Saves the generation the average merit and the maximum merit of  the pop
     nger(ger,1)=ger;
     meanFm(ger,1)=mean(Fm_restart(:,2));
     [maxFm(ger,1),iFm]=max(Fm_restart(:,2));

     %Saves the value of the objective function for the parametes with the
     %highest merit of the generation
     iobj=newpop(iFm,:);
     [xobj] = descod_ind(nvar,iobj,dom,dim,p);
     [Fobj] = Fev_restart(iFm,2);
     Fo(ger,1)=Fobj;
else
    pop=newpop;
end
%Writes the following to a txt file, for intermediate and easy verification
%of the optimization process:
%number of the generation
%parameter 1
%parameter 2
%...
%parameter n
%Objective function

fid=fopen('ger_var_Fob.txt','wt');
fprintf(fid, '%f \n', ger, xobj, Fobj);
fclose(fid)

end

%Evaluates the merit of the final population
[Fm_new] = merit(npop,pop,nvar,dom,dim,p,obj);

%Gets the variable with the highest fit
[max_num,max_new]=max(Fm_new(:,2));

%Decode the variables of highest fit from binary to real
final=newpop(max_new,:);
[xfinal] = descod_ind(nvar,final,dom,dim,p);



end
