clear all
close all
% ���{�����Z
%�؂�ւ���
load mldivide_test3_multi.mat
%mdata=zeros(1,4);

j=1;
% N�ω�
k=1;kmax=10;N=100;Nplus=100;
% prec�ω�
pmin=53;pmax=1696;pplus=1;ptime=2;
prec=pmin;
% ���^�] 1��
tic
set_default_prec(77);    
A=2*rand(77,77)-1;
b=2*rand(77,1)-1;   
A=multi(A);
b=multi(b);   
x=A\b;
toc
%�����܂�
while(1)
    while(1)
        [s1,s2]=size(mdata);
        Ndata=NaN(s1,1);Precdata=NaN(s1,1);
        Ndata=mdata(:,1)==N;Precdata=mdata(:,2)==prec;
        if max(Ndata.*Precdata)==1 
            % �؂�ւ���
            %prec=prec+pplus;
            prec=prec*ptime;
            if prec > pmax
                    break
            end
        else
            prec
            mdata(end+1,3)=NaN;              
            set_default_prec(prec);    
            A=2*rand(N,N)-1;
            b=2*rand(N,1)-1;   
            A=multi(A);
            b=multi(b); 
            tic;  
            x=A\b;
            mdata(end,3)=toc;
%                 if max(max(isnan(double(x))))==1
%                     mdata(end,4)=NaN;
%                     mdata(end,5)=NaN;
%                 else
%                     mdata(end,4)=max(rad(x));
%                     mdata(end,5)=log10(max(rad(x)));                
%                 end         
            mdata(end,1)=N; 
            mdata(end,2)=prec;           
            mdata=sortrows(mdata);
            
            save 'mldivide_test3_multi.mat' mdata
           
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