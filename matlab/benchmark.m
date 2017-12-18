clear all
close all

N=500;prec=53;
% ŽŽ‰^“] 1‰ñ
tic
set_default_prec(77);    
a=2*rand(77,77)-1;
b=2*rand(77,1)-1;   
a=imulti(a);
b=multi(b);   
toc

% % imulti
% set_default_prec(prec);
% A=2*rand(N,N,N)-1;
% A=imulti(A);
% tic
% %A.*A.*A.*A.*A.*A.*A.*A.*A.*A;
% %A-A-A-A-A-A-A-A-A-A;
% %A./A./A./A./A./A./A./A./A./A;
% t1=toc

% 
% multi
% set_default_prec(prec);
% B=2*rand(N,N)-1;
% B=multi(B);
% tic
% %B.*B.*B.*B.*B.*B.*B.*B.*B.*B;
% %B-B-B-B-B-B-B-B-B-B;
% %B./B./B./B./B./B./B./B./B./B;
% down(B);
% B+B+B+B+B+B+B+B+B+B;
% t2=toc



% %multi
% set_default_prec(prec);
% B=2*rand(N,N)-1;
% B=multi(B);
% tic
% B=abs(B);
% t2=toc

set_default_prec(prec);
A=2*rand(N,N)-1;
A=multi(A);
B=2*rand(N,1)-1;
B=multi(B);
tic
A\B;
t1=toc
%disp('ans=');disp(num2str(ans,'%+.7e'))