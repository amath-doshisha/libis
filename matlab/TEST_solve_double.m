function Rmax=TEST_solve_double

n=5;
disp('Solve two equations A*X(:,1)=B(:,1), A*X(:,2)=B(:,2) simultanuously.');
B=[ones(n,1) rand(n,1)]
A=rand(n,n)
disp('Solutions are X(:,1) and X(:,1).');
X=A\B
disp('Residuals of solutions are');
R=A*X-B
disp('Maximus norms of residuals of solutions are');
Rmax=max(abs(R))