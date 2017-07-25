%function sample5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%n=10;A=rand(n)
A=toeplitz([1 2 0 1 0],[1 2 0 0 0])
lambda=sort(eig(A))

set_default_prec(512)
A=multi(A);
lambda2=sort(double(eig(A)))

plot(real(lambda),imag(lambda),'o',real(lambda2),imag(lambda2),'x')
axis square
grid on




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all
% close all
% 
% set_default_prec(512)
% A=zeros(10,10,'multi');
% A(:)=1:100;
% A
% A(44)
% A(7,:)
% A(:,5)
% A(1:2:10,2:2:10)
% i=1;
% while i<=100
%     A(i)=-1*A(i)
%     i=i+1;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
%close all
%set_default_prec(256);
%
%A=[1 -1; -1 1]
%A=multi(A);
%
%disp(num2str(A))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all
%close all
%
%N=100;
%prec=1024;
%set_default_prec(prec);
%%A=eye(N,N);
%A=rand(N,N)>0.5
%
%% i=1;
%% while i<=N
%%     j=1;
%%     while j<=N
%%         if i==j-1
%%             A(i,j)=1;
%%         elseif i-1==j
%%             A(i,j)=-3;
%%         end
%%         j=j+1;
%%     end
%%     i=i+1;
%% end
%
%
%
%disp('Generated matrix: ');
%disp(num2str(A,'%.2f     '));
%disp('Compute eigenvalues of generated matrix A on double precision: ');
%l=sort(eig(double(A)));
%disp(num2str(l,'%.40f'))
%figure(1);plot(real(l),imag(l),'o')
%
%disp(sprintf('Compute eigenvalues of A on %d-bit:',prec));
%D=eig(multi(A));
%disp(num2str(D,'%.40f'))
%D=double(D);
%figure(2);plot(real(D),imag(D),'o',real(l),imag(l),'x')
