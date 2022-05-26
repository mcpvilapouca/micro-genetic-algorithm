%*************************************************************************%
%   Micro-genetic algorithm for optimization
%*************************************************************************%
%         Date               Programmer
%       ========        ====================
%      23/01/2018          Maria Vila Pouca
%*************************************************************************%
%*************************************************************************%
% initialization
clearvars; close all; clc; clear all
addpath([pwd,'/private_functions']);
%*************************************************************************%
%--To be defined by user-------------------------------------------------
%
%Note: change filename (name of ABAQUS input file or job) in the merit.m
% change name of umat if using in the merit.m
%default: cube_umat.input
%umat: umat_hgo_visco.std-o
%
p=2;              %no. decimal points
nvar=7;           %no. variables to optimize
npop=5;           %population dimension
dom1=[0.01,1.0];   %C10    %variable domain
dom2=[0.01,3.0];   %K11
dom3=[0.01,50];   %K12
dom4=[0.0,2.0];   %BETA1
dom5=[0.01,100];   %TAU1
dom6=[0.0,2];      %BETA2
dom7=[0.01,1500];   %TAU2

%Random population? h0='no' or h0='yes'
h0='no';
%if h0='no', define initial population
initialpop=[0.21;2.16;29.69;0.44;78.33;0.1;1234.86];

%how many times the micro-genetic algorithm will be performed
nloops=1;

%Absolute Convergency
ger_max=500;   %maximum number of generations
conv_iter=3000; %maximum number of converged iterations
tol=1E-10;      %tolerance
nrestart=40;    %minimum restarts number

%Aim of the analysis
%   maximize: 'max'
%   minimize: 'min'
obj='min';
%--------------------------------------------------------------------------
%Elitism option: 'y'/'n'
elitismo='y';

%add your domains
dom=vertcat(dom1,dom2,dom3,dom4,dom5,dom6,dom7);

%**************************************************************************
%Private Functions*********************************************************
%**************************************************************************
%
%
%
for i=1:nvar
dim(i,1)=dom(i,2)-dom(i,1);
end

%Loop to the micro-geneetic algorithm
for i=1:nloops
   [xfinal,Fo,maxFm,restart,nger,meanFm] = micro_gen(obj,...
       elitismo,nvar,dim,dom,p,npop,ger_max,conv_iter,tol,nrestart,...
       h0,initialpop);

   %Save the final value of each parameter and the merit function
   % in a matrix, where the vals of each liine match one evaluation of
   %the micro-genetic

   FFO(:,i)=Fo(:,1);
   FFmax(:,i)=maxFm(:,1);
   FFmean(:,i)=meanFm;
   rr(1,i)=restart;
end


%Average of the final values of the nloops

FFOfinal=sum(FFO,2)/size(FFO,2);
FFmaxfinal=sum(FFmax,2)/size(FFmax,2);
FFmeanfinal=sum(FFmean,2)/size(FFmean,2);

rrmean=mean(rr);

%Save in a text file the variable that have the maximum score
%in every generation and the value of the objective function

fid=fopen('Fm_max.txt','wt');
fprintf(fid, '%f \n', FFmaxfinal(:,1));
fclose(fid);
fid=fopen('Fobj.txt','wt');
fprintf(fid, '%f \n', FFOfinal(:,1));
fclose(fid);
fid=fopen('var.txt','wt');
fprintf(fid, '%f \n', xfinal);
fclose(fid);

% %**************************************************************************
% %**************************************************************************
exit
