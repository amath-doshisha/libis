function tf=is_type_r(A)
if isa(A,'multi')
tf=prod(isfield(A.data,{'r_digits','r_exp','r_prec','r_sign'}));
else
    tf=0; % “ü—Í‚ªmultiŒ^ˆÈŠO‚È‚ç‚»‚à‚»‚àŒÄ‚Î‚ê‚È‚¢
end