function [A,wr,wi]=eig_hqr(A)

% define in C
% [int n, float A, float wr, float wi ]
% int nn,m,l,k,j,its,i,mmin;
% float z,y,x,w,v,u,t,s,r,q,p,anorm;

%void nrerror(char error_text[])
%/* Numerical Recipes standard error handler */
%{
%	fprintf(stderr,"Numerical Recipes run-time error...\n");
%	fprintf(stderr,"%s\n",error_text);
%	fprintf(stderr,"...now exiting to system...\n");
%	exit(1);
%}

% check size
[m,n]=size(A);
if m~=n
    disp('Error');
end

% ƒwƒbƒZƒ“ƒxƒ‹ƒOŒ^‚É•ÏŠ·
[B,H0,alpha0]=eig_mat_hessenberg(A);
A=B;
% disp('Hessenberg matrix:');
% disp(A);

% start
anorm=abs(A(1,1));

i=2;
while i<=n
    j=i-1;
    while j<=n
        anorm=anorm+abs(A(i,j));
    j=j+1;
    end
i=i+1;    
end
%debag---------
%disp(anorm);
%--------------

nn=n;
t=0.0;
while nn>=1
    its=0;
    l=-1;
    while l<(nn-1)         %do{ 
    l=nn;
    while l>=2
        s=abs(A(l-1,l-1))+abs(A(l,l));
        if s==0.0
           s=anorm;
        end
        if (abs(A(l,l-1))+s)==s
            %---------
%            disp('s=');disp(s);
            %---------
%            disp('A(l-1,l-1)=');  disp(A(l-1,l-1));
%           disp('A(l,l)=');         disp(A(l,l));
%            disp('A(l,l-1)=');      disp(A(l,l-1));
            %---------
            break;
        end
        l=l-1;
    end
    x=A(nn,nn);
    if l==nn
        %debag---------
        %disp('l=');disp(l);disp('nn=');disp(nn);
        %--------------
        wr(nn)=x+t;
        wi(nn)=0.0;
        nn=nn-1;
    else
        y=A(nn-1,nn-1);
        w=A(nn,nn-1)*A(nn-1,nn);
        if l==(nn-1)
            %debag---------
            %disp('here');
            %--------------
            p=0.5*(y-x);
            q=p*p+w;
            z=sqrt(abs(q));
            x=x+t;
            if q>=0.0
                z=p+eig_sign(z,p);
                wr(nn-1)=x+z;
                wr(nn)=x+z;
                if z~=0
                    wr(nn)=x-w/z;
                end
                wi(nn-1)=0.0;
                wi(nn)=0.0;
            else
                wr(nn-1)=x+p;
                wr(nn)=x+p;
                wi(nn-1)=-z;
                wi(nn)=z;
            end
            nn=nn-2;
        else
            if its==30
                fprintf('Too many iteration in HQR');
                exit(1);
            end
            if its==10 || its==20
                t=t+x;
                i=1;
                while i<=nn
                    A(i,i)=A(i,i)-x;
                    i=i+1;
                end
                s=abs(A(nn,nn-1))+abs(A(nn-1,nn-2));
                y=0.75*s;
                x=0.75*s;
                w=-0.4375*s*s;
            end
            its=its+1;
            %debag--------
            %disp('its=');disp(its);
            %-------------
            m=nn-2;
            while m>=l
                z=A(m,m);
                r=x-z;
                s=y-z;
                p=(r*s-w)/A(m+1,m)+A(m,m+1);
                q=A(m+1,m+1)-z-r-s;
                r=A(m+2,m+1);
                s=abs(p)+abs(q)+abs(r);
                p=p/s;
                q=q/s;
                r=r/s;
                if m==l
                   break; 
                end
                u=abs(A(m,m-1))*(abs(q)+abs(r));
                v=abs(p)*(abs(A(m-1,m-1))+abs(z)+abs(A(m+1,m+1)));
                if (u+v)==v
                    break;
                end
            m=m-1;    
            end
            i=m+2;
            while i<=nn
               A(i,i-2)=0.0;
               if i~=(m+2)
                  A(i,i-3)=0.0; 
               end
            i=i+1;   
            end
            k=m;
            while k<=nn-1
                if k~=m
                  p=A(k,k-1);
                  q=A(k+1,k-1);
                  r=0.0;
                  if k~=(nn-1)
                      r=A(k+2,k-1);
                  end
                  x=abs(p)+abs(q)+abs(r);
                  if x~=0.0
                      p=p/x;
                      q=q/x;
                      r=r/x;
                  end
               end
               s=eig_sign(sqrt(p*p+q*q+r*r),p);
               if s~=0.0
                   if k==m
                       if l~=m
                          A(k,k-1)=-A(k,k-1);
                       end
                   else
                       A(k,k-1)=-s*x;
                   end
                   p=p+s;
                   x=p/s;
                   y=q/s;
                   z=r/s;
                   q=q/p;
                   r=r/p;
                   j=k;
                   while j<=nn
                       p=A(k,j)+q*A(k+1,j);
                       if k~=(nn-1)
                          p=p+r*A(k+2,j);
                          A(k+2,j)=A(k+2,j)-p*z;
                       end
                       A(k+1,j)=A(k+1,j)-p*y;
                       A(k,j)=A(k,j)-p*x;
                       j=j+1;
                       end
                       
                       %mmin = nn < k+3 ? nn : k+3; —ñ‚ÌC³
                       if nn<k+3
                           mmin=nn;
                       else
                           mmin=k+3;
                       end
                       
                       i=l;
                       while i<=mmin
                           p=x*A(i,k)+y*A(i,k+1);
                           if k~=(nn-1)
                              p=p+z*A(i,k+2);
                              A(i,k+2)=A(i,k+2)-p*r;
                           end
                           A(i,k+1)=A(i,k+1)-p*q;
                           A(i,k)=A(i,k)-p;
                           i=i+1;
                       end
               end
               k=k+1;
            end
        end
        end
    end          %}while l<(nn-1)
    %disp('ok-------------------');
end

