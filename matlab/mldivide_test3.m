clear all
close all
% ���{����ԉ��Z
%�؂�ւ���
load mldivide_test3.mat
%data=zeros(1,4);

j=1;
% N�ω�
k=1;kmax=91;N=100;Nplus=10;
% prec�ω�
pmin=53;pmax=1696;pplus=1;ptime=2;
prec=pmin;
% ���^�] 1��
tic
set_default_prec(77);    
A=2*rand(77,77)-1;
b=2*rand(77,1)-1;   
A=imulti(A);
b=imulti(b);   
x=A\b;
toc
%�����܂�
while(1)
    while(1)
        [s1,s2]=size(data);
        Ndata=NaN(s1,1);Precdata=NaN(s1,1);
        Ndata=data(:,1)==N;Precdata=data(:,2)==prec;
        if max(Ndata.*Precdata)==1 
            % �؂�ւ���
            %prec=prec+pplus;
            prec=prec*ptime;
            if prec > pmax
                    break
            end
        else
            prec
            data(end+1,3)=NaN;                   
            set_default_prec(prec);    
            A=2*rand(N,N)-1;
            b=2*rand(N,1)-1;   
            A=imulti(A);
            b=imulti(b);  
            tic;  
            x=A\b;
            data(end,3)=toc;
                if max(max(isnan(double(x))))==1
                    data(end,4)=NaN;
                    data(end,5)=NaN;
                else
                    data(end,4)=max(rad(x));
                    data(end,5)=log10(max(rad(x)));                
                end         
            data(end,1)=N; 
            data(end,2)=prec;           
            data=sortrows(data);
            
            save 'mldivide_test3.mat' data
           
            % �؂�ւ���
            %prec=prec+pplus;
            prec=prec*ptime;
                if prec > pmax
                    break
                end          
        end
    end
j=1;
k=k+1
N=N+Nplus
prec=pmin;
    if k >kmax
        break
    end
end