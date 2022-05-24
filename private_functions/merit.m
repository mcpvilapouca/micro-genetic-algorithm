function [Fm,Fev] = merit(npop,pop,nvar,dom,dim,p,obj,filename,umatname)
%Evaluate the merit of the population

for i=1:npop
 ind=pop(i,:);
 xr(i,:) = descod_ind(nvar,ind,dom,dim,p);
end

%run abaqus with xr for all elements of the population at the same time
% in each folder

for i=1:npop
  %go to the folder of the pop
  cd([pwd,'/xr',num2str(i)]);

 %write the variables to a file
fid1=fopen('par.txt','wt');
fprintf(fid1, '%f \n', [xr(i,:)]');
fclose(fid1);

%call python file to replace parameters
system('python3 subs_param.py');

%call abaqus to run the analysis
if i==npop
system('abaqus job=filename user=umatname ask_delete=OFF cpus=1 -interactive');
else
system('abaqus job=filename user=umatname ask_delete=OFF cpus=1');
end

cd ../
end

for i=1:npop
    cd([pwd,'/xr',num2str(i)]);
    while exist(filename+".lck", 'file')
    end
 cd ../
end

for i=1:npop
 cd([pwd,'/xr',num2str(i)]);
 [Fe] = eval_fun(filename);
 Fm(i,1)=i;
 Fev(i,1)=i;
 switch obj
   case 'max'
      Fm(i,2)=Fe;
      Fev(i,2)=Fe;
   case 'min'
      Fm(i,2)=1/(1+(Fe));
      Fev(i,2)=Fe;
 end
   cd ../
end

end

