function Rmax=TEST_solve_multi

set_default_prec(128);

n=5;
disp('Solve two equations A*X(:,1)=B(:,1), A*X(:,2)=B(:,2) simultanuously.');
A=rand(n,n,'multi')
B=[ones(n,1,'multi') rand(n,1,'multi')]
disp('Solutions are X(:,1) and X(:,1).');
X=A\B
disp('Residuals of solutions are');
R=A*X-B
disp('Maximus norms of residuals of solutions are');
Rmax=max(abs(R))