function tf=isreal(A)
tf=logical(prod(isfield(A.data,{'r_digits','r_exp','r_prec','r_sign'})));