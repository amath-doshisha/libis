function tf=is_type_c1(A)
if isa(A,'multi')
tf=prod(isfield(A.data,{'cr_digits','cr_exp','cr_prec','cr_sign','ci_digits','ci_exp','ci_prec','ci_sign'}));
else
    tf=0; % ���͂�multi�^�ȊO�Ȃ炻�������Ă΂�Ȃ�
end