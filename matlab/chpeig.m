function [P,D]=chpeig(A)
if isa(A,'double')
    A=cmulti(A);
elseif isa(A,'rmulti')
    A=cmulti(A);
end
if isa(A,'cmulti')
    [Pdata,Ddata]=cmulti_data(get_prec,'one','hpeig',A.data);
    P=cmulti(-1,Pdata);
    D=cmulti(-1,Ddata);
else
    error('chpeig');
end