clear all
close all
%load mldivide_teat3.mat

data=NaN(1,4);

j=1;k=1;kmax=1;p=0;
pmin=50;pmax=100;pplus=4;N=100;Nplus=1;

prec=pmin;
data(end+1,1)=N;
data(end,2)=pmin;
while(k<=kmax)
    while(p<=pmax)
        tic;
        set_default_prec(prec);    
        A=2*rand(N,N)-1;
        b=2*rand(N,1)-1;   
        A=imulti(A);
        b=imulti(b);   
        x=A\b;
        prec
            if max(max(isnan(double(x))))==1
                data(end,4)=NaN;
            else
                data(end,4)=max(rad(x));
            end
        data(end,3)=toc;
        prec=prec+pplus;
        data(end+1,2)=prec;
        j=j+1;
    end
j=1;
k=k+1;N=N+Nplus;
prec=pmin;
k
data(end+1,1)=N;
data(end,2)=prec;
end













save 'mldivide_test3.mat' data

% data =[10   100   128;
%     20   120   300;
%     15   120   300
%     ]
% 
% m=size(data,1);
% i=1;
% while i<=m
%     j=i+1;
%     while j<=m
%         if data(j,1)<data(i,1)
%             a=data(j,:);
%             data(j,:)=data(i,:);
%             data(i,:)=a;
%         end
%         
%         j=j+1;
%     end
%     i=i+1;
% end
% 
% data
% 