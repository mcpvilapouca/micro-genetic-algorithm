function [Fe] = eval_fun(filename)

errf1=0;

%Data from Experimental Curve
        d2=dlmread('experimental.txt');
        xx_ref=d2(:,1);
        yy_ref=d2(:,2);

if exist(filename+".lck", 'file')
    delete filename+".lck"
    errf1=errf1+1;
else
    %check if the job was completed successfuly
    fid=fopen(filename+'.sta');
    while ~feof(fid)
    tline = fgetl(fid);
    IndexC = strfind(tline,'NOT BEEN COMPLETED');
    end
    fclose(fid);
    cond=isempty(IndexC);
     if cond==1
        %Run python script that generate the report
        system('python3 get_num.py');

        if exist("report.rpt", 'file')
        d1=dlmread('report.rpt');
        xx=d1(:,1);
        yy=d1(:,2);
        else
        errf1=errf1+1;
        end

  %check if there is NaN in the results
        if sum(isnan(d1(:,2)))>=1
            errf1=errf1+1;
        else
            errf1=errf1+0;
        end

	else
	errf1=errf1+1;
     end
end

%check if there is the same number of results. If not it's probably
%the error "primary variable is not available", from ABAQUS
 if errf1==0
        if size(yy,1)~=size(xx_ref,1)
            errf1=errf1+1;
        else
            errf1=errf1+0;
        end
 end


  if errf1==0
        %Generate function
        sum1=0;

        for i=1:size(xx_ref,1)
            sum1=sum1+((yy_ref(i,1)-yy(i,1)))^2;
        end
        se1=sqrt(sum1);
        [Fe]=se1;
     else
        [Fe]=1000;
  end


end
