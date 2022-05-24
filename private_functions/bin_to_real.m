function [xr] = bin_to_real(int_bin,dom,p)

dim=dom(1,2)-dom(1,1);   %Domain dimension

n=round(log2(dim*10^p));  %No. of bits

%Real decoded number
xr=round(int_bin*((dim/((2^n)-1)))+dom(1,1),p);

end

