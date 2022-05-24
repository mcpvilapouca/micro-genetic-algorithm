function [pop,space,disp_xri] = random_pop(nvar,dim,p,npop,dom)
%Generates the population randomly
for i=1:npop
    [ind] = random_ind(nvar,dim,p);
    [xr] = descod_ind(nvar,ind,dom,dim,p);
    if i==1
    pop=ind;
    else
    pop=vertcat(pop,ind);
    end
    space(i,:)='    ';

    disp_xri(i,:)=xr;
end
end

