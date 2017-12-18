clear all
close all


j=1;k=1;kmax=1;s1=0;s2=0;
pmin=1320;pmax=1320;pplus=1;N=1000;Nplus=10;
n=(pmin:pplus:pmax);
[s1,s2]=size(n);

prec=pmin;
times=zeros(kmax,s2);
rad1=zeros(kmax,s2);
while(k<=kmax)
    while(j<=s2)
        tic;
        set_default_prec(prec);    
        A=2*rand(N,N)-1;
        b=2*rand(N,1)-1;   
        A=imulti(A);
        b=imulti(b);   
        x=A\b;
%disp('A=');disp(num2str(A,'%+.7e'))
%disp('b=');disp(num2str(b,'%+.7e'))
%disp('x=');disp(num2str(x,'%+.7e'))
        prec
%disp('D=');disp(num2str(D,'%+.7e'))
            if max(max(isnan(double(x))))==1
                rad1(k,j)=1;
            else
                rad1(k,j)=max(rad(x));
            end
        times(k,j)=toc;
        prec=prec+pplus;
        j=j+1;
    end
j=1;
k=k+1;N=N+Nplus;
prec=pmin;
k
end
k=1;
while(k<=kmax)   
    figure(1);
    semilogy(n,rad1(k,:));
    hold on
    figure(2);
    plot(n,times(k,:));
    hold on
%     axis([nmin nmax 0.1 500])
k=k+1;
end
% k=1;
% while(k<=kmax)   
%     figure(3);
%     semilogy(n,times(k,:));
%     hold on
% %     axis([nmin nmax 0.1 500])
% k=k+1;
% end
% k1=(1:13);
% figure(4);semilogy(k1,times(:,61),'-s');
% figure(5);plot(k1,times(:,61),'G-s');