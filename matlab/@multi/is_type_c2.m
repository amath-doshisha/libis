function tf=is_type_c2(A)
if isa(A,'multi')
    tf=prod(isfield(A.data,{'C0r_digits','C0r_exp','C0r_prec','C0r_sign','C0i_digits','C0i_exp','C0i_prec','C0i_sign','C1r_digits','C1r_exp','C1r_prec','C1r_sign','C1i_digits','C1i_exp','C1i_prec','C1i_sign'}));
else
    tf=0;
end