clear all
 clc
format long 
%N=100
 
   M=8; N=15

A=rand(N,N)
%A=A'*A;
% for k=1:N
%     A(k,k)=0
% end
A=(sqrt((N)/trace(A'*A)))*A; 

I=eye(N);
 

B=rand(N,M); 
B=I(:,1:M)
% 
% B=[1 0 0 0 0 0; 0 0 0 0 0 1]';
%B=[0 1]';

 %B=[1 0; 0 0; 0 1; 0 0];
I1=eye(M);
tf=1


% WB0=zeros(N,N);
% ot=0.025;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
% end
%  C0=pinv(WB0); 
 
pf=expm(A*tf);

%xf=pf*x0;


% syms s 
% AA=s*I-A;
% invAA=inv(A);
% p=ilaplace(invAA,t);
% p=simplify(p);




v=0.01;


iterations=100
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  



 
for ii=1:iterations
   
  
    pf=expm(A*tf);
    
 WB0=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=pinv(WB0); 

   
Cost1(ii)=trace(C*pf*pf');

AA(:,:,ii)=A;

   
V=eig(A);

BB=fliplr(vander(V));

% if rank(diag(V)) < N;
% 
%  break
% else  


    
CC2=diag(expm(diag(V)*tf));
%alpha_tf=pinv(BB)*CC2;
alpha_tf=BB\CC2;


Cpf=C*pf;

F111=zeros(N,N);

ot=0.02;
for kk=1:tf/ot
    
CC1=diag(expm(diag(V)*ot*kk));
%alpha_t=pinv(BB)*CC1;
alpha_t=BB\CC1;

CB=Cpf*pf'*C*expm(A*ot*kk)*B*B';


F11=zeros(N,N);
for i=1:N-1
    
F1=zeros(N,N);
    for k=1:i
        
        F1=A'^(k-1)*CB*A'^(i-k)+F1;

    end

 F11=-2*alpha_t (i+1)*F1+F11;
end

F111=real(F11)*ot+F111;
  
end 

    
    
  F22=zeros(N,N);  
for i=1:N-1
  F2=zeros(N,N);
    for k=1:i
        
        F2=A'^(k-1)*Cpf*A'^(i-k)+F2;
    end

F22=2*alpha_tf(i+1)*F2+F22;
end
   
 F=F111+real(F22);   
      
% F=F/norm(F);
% 
% F3=(trace(A'*A)-N)*4*A; % gradient of norm function
% 
% F3=F3/norm(F3);
% 
% F0=reshape((eye(N*N)-F3(:)*pinv(F3(:)))*F(:),N,N); 


 % F=F'+F-diag(diag(F));
%  for kkk=1:N
%     F(kkk,kkk)=0;
% end

A=A-v*F/norm(F);
A=(sqrt((N)/trace(A'*A)))*A; 


ii  
end

% for i=1:length(A)
%     for j=1:length(A)
%        if abs(A(i,j))<1.5*sum(sum(A))/(N*N)
%        A(i,j)=0;
%        end
%     
%     end
% end


  plot((Cost1),'r-*')
  
  figure(2)
  plot(log10(Cost1),'r-*')
  
  
  
  
 %subs(x,'t',tf)

AL=zeros(N,N);
for kkkk=1:N
    
   AL= alpha_tf(kkkk)*A^(kkkk-1)+AL;
end

AAAA=expm(A)
AAAA-AL