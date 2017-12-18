clear all
close all
% グラフ作成
load mldivide_test3.mat
load mldivide_test3_multi.mat
% data=(N,prec,time,rad,log10(rad),minprec,16prec)

% 値採用ルール
% 必要最小prec　半径にNaN以外が2つ続くところ
% 必要最小prec(○倍精度)　log10半径が２つ続けて×以下になるところ （倍精度なら-16，4倍なら-35)

% タグ付け 半径にNaNが発生しない必要最小prec
j=2;N=0;count=0;tag=0;
[s1,s2]=size(data);
while(1)
    k=j;
    N=data(j,1);
    while(data(k,1)==N)
        if(isnan(data(k,5)))
            count=0;
            data(k,6)=NaN;
        else
            count=count+1;
            if(tag==0)
                if(count==2)
                    data(k-1,6)=1;
                    data(k,6)=NaN;
                    tag=1;
                else
                    data(k,6)=NaN;
                end
            else
                data(k,6)=NaN;
            end
        end
        k=k+1;
        if(k>s1)
            break
        end
   end
   j=k;count=0;tag=0;
   if(j>s1)
       break
   end
end

% タグ付け 半径<=10^(-16)必要最小prec
j=2;N=0;count=0;tag=0;
[s1,s2]=size(data);
while(1)
    k=j;
    N=data(j,1);
    while(data(k,1)==N)
       
        if(data(k,5)<=-16)
            count=count+1;
            if(tag==0)
                if(count==2)
                    data(k-1,7)=1;
                    data(k,7)=NaN;
                    tag=1;
                else
                    data(k,7)=NaN;
                end
            else
                data(k,7)=NaN;
            end
        else
            count=0;
            data(k,7)=NaN;
        end
        k=k+1;
        if(k>s1)
            break
        end
   end
   j=k;count=0;tag=0;
   if(j>s1)
       break
   end
end
    
% plot(prec,log10rad)
N1=10;
N2=100;
N3=1000;
x1=[10,1696];y1=[6.1786,-505.3538];
x2=[95,1696];y2=[13.2036,-472.3398];
x3=[806,1696];y3=[125.6262,-144.3998];
data1=data(data(:,1)==N1 & (data(:,2)==10 | data(:,2)==53 | data(:,2)==106 | data(:,2)==212 | data(:,2)==424 | data(:,2)==848 | data(:,2)==1696),:);
data2=data(data(:,1)==N2 & (data(:,2)==95 | data(:,2)==53 | data(:,2)==106 | data(:,2)==212 | data(:,2)==424 | data(:,2)==848 | data(:,2)==1696),:);
data3=data(data(:,1)==N3 & (data(:,2)==53 | data(:,2)==106 | data(:,2)==212 | data(:,2)==424 | data(:,2)==848 | data(:,2)==1696 | data(:,2)==806 | data(:,2)==1272),:);
figure(1);
%plot(data1(:,2),data1(:,5),'o','LineWidth',1,'MarkerSize',15,'MarkerFaceColor','w'); 

%plot(data2(:,2),data2(:,5),'o','LineWidth',1,'MarkerSize',15,'MarkerFaceColor','w'); 
%plot(data3(:,2),data3(:,5),'o','LineWidth',1,'MarkerSize',15,'MarkerFaceColor','w'); 
plot(x1,y1,'-ok','LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w','MarkerEdgeColor','b');grid on;
hold on;plot(x2,y2,'-ok','LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w','MarkerEdgeColor','r');
plot(x3,y3,'-ok','LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w','MarkerEdgeColor','g');
axis([0 1696 -600 400]);

% semilogy(prec,time)
figure(2);
semilogy(data1(:,2),data1(:,3),'-o','LineWidth',3,'MarkerSize',15,'MarkerFaceColor','w');
grid on
hold on
semilogy(data2(:,2),data2(:,3),'-o','LineWidth',3,'MarkerSize',15,'MarkerFaceColor','w');
semilogy(data3(:,2),data3(:,3),'-o','LineWidth',3,'MarkerSize',15,'MarkerFaceColor','w');
axis([0 1696 10^(-3) 10^3])


% dataをprecごとに取り出す
% 
prec1=53;
prec2=424;
prec3=1696;
data1=data(data(:,2)==prec1,:);
data2=data(data(:,2)==prec2,:);
data3=data(data(:,2)==prec3,:);
figure(3);
loglog(data1(:,1),data1(:,3),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
grid on
hold on
loglog(data2(:,1),data2(:,3),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
loglog(data3(:,1),data3(:,3),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');


%
prec1=424;
prec2=848;
prec3=1696;
data1=data(data(:,2)==prec1,:);
data2=data(data(:,2)==prec2,:);
data3=data(data(:,2)==prec3,:);
figure(6);
plot(data1(:,1),data1(:,5),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
grid on
hold on
plot(data2(:,1),data2(:,5),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
plot(data3(:,1),data3(:,5),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
axis([0 1000 -550 200])
%

% 計算必要最小prec (N,prec)
% 半径(10^-16)のときの(N,prec)
data1=data(data(:,6)==1 & (data(:,1)<=150 | data(:,1)==200 | data(:,1)==300| data(:,1)==400| data(:,1)==500 | data(:,1)==1000),:);
data2=data(data(:,7)==1 & (data(:,1)<=150 | data(:,1)==200 | data(:,1)==300| data(:,1)==400| data(:,1)==500 | data(:,1)==1000),:);
figure(4)
plot(data1(:,1),data1(:,2),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
grid on
hold on
plot(data2(:,1),data2(:,2),'-o','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
axis([0 1000 0 1400])


% loglog(N,time) 計算必要最小precで区間と多倍長を比較

data3=data(data(:,2)==10 &   (data(:,1)<=50| data(:,1)==100 | data(:,1)==200| data(:,1)==300| data(:,1)==400| data(:,1)==500| data(:,1)==600 |data(:,1)==700|data(:,1)==800|data(:,1)==900| data(:,1)==1000),:); % prec=10 固定
data4=data(data(:,6)==1 &    (data(:,1)<=50| data(:,1)==100 | data(:,1)==200| data(:,1)==300| data(:,1)==400| data(:,1)==500| data(:,1)==600 |data(:,1)==700|data(:,1)==800|data(:,1)==900| data(:,1)==1000),:); % 必要最小prec 変動
data1=mdata(mdata(:,2)==10 & (mdata(:,1)<=50|mdata(:,1)==100| mdata(:,1)==200| mdata(:,1)==300|mdata(:,1)==400|mdata(:,1)==500|mdata(:,1)==600|mdata(:,1)==700|mdata(:,1)==800|mdata(:,1)==900| mdata(:,1)==1000),:); % prec=10 固定
data2=zeros(1,4);
j=1;k=1;[s1,s2]=size(data4);
while(1)
    N=data4(j,1);prec=data4(j,2);
    [ss1,ss2]=size(mdata(mdata(:,1)==N & mdata(:,2)==prec,:));
    if ss1>0
        data2(k,:)=mdata(mdata(:,1)==N & mdata(:,2)==prec,:);
        k=k+1;
    end
    j=j+1;
    if j>s1
        break;
    end
end

figure(5);
loglog(data1(:,1),data1(:,3),'-bo','LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w');
grid on
hold on
loglog(data2(:,1),data2(:,3),'-ro','LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w');
loglog(data3(:,1),data3(:,3),'-o','color',[0,0.7,1],'LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w');
loglog(data4(:,1),data4(:,3),'-o','color',[1,0,0.7],'LineWidth',3,'MarkerSize',12,'MarkerFaceColor','w');
axis([9 1000 0.001 1000])


