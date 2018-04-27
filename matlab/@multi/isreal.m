function tf=isreal(A)
if isa(A,'multi')
    tf=prod(isfield(A.data,{'r_digits','r_exp','r_prec','r_sign'}),'native') || prod(isfield(A.data,{'R0_digits','R0_exp','R0_prec','R0_sign','R1_digits','R1_exp','R1_prec','R1_sign'}),'native');
else
    tf=0; % “ü—Í‚ªmultiŒ^ˆÈŠO‚È‚ç‚»‚à‚»‚àŒÄ‚Î‚ê‚È‚¢
end
