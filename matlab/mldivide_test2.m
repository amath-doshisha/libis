clear all
close all


j=1;k=1;kmax=4;
%N=100;Nplus=100;
p=[825,950,1075,1200,1300];
%p=[175,200,250,300,700,750,775,1500,2000];
%p=[53,64,128,256,512,1024,2048,4096,8192,16384];
n=[600,700,800,900];
[s1,s2]=size(p);

prec=p(1,1);
times=zeros(kmax,s2);
rad1=zeros(kmax,s2);
while(k<=kmax)
    while(j<=s2)
        tic;
        prec=p(1,j);
        N=n(1,k);
        set_default_prec(prec);    
        A=2*rand(N,N)-1;
        b=2*rand(N,1)-1;   
        A=imulti(A);
        b=imulti(b);   
        x=A\b;
        prec
            if max(max(isnan(double(x))))==1
                rad1(k,j)=NaN;
            else
                if max(rad(x))>1
                    rad1(k,j)=NaN;
                else
                    rad1(k,j)=max(rad(x));
                end
            end
        times(k,j)=toc;
        j=j+1;
    end
j=1;
k=k+1;
prec=p(1,1);
k
end
k=1;
while(k<=kmax)   
    figure(1);
    semilogy(p,rad1(k,:),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');  
    grid on
    hold on
    axis([0 p(1,s2) 10^(-50) 10^(0)]);
    figure(2);
    semilogy(p,times(k,:),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
    grid on
    hold on
k=k+1;
end
l=1;lmax=s2;
while(l<=lmax)
    figure(3);
    loglog(n,times(:,l),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
    grid on
    hold on
    l=l+1;
end

