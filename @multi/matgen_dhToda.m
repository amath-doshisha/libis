function A=matgen_dhToda(N,M,lambda,c)
A=multi('matgen_dhToda',N,M,multi(lambda).data,multi(c).data);

%A=multi();
%A.data=multi_mex('matgen_dhToda',get_default_prec(),double(get_auto_prec_mode()),N,M,multi(lambda).data,multi(c).data);

