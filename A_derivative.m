clear all
 clc
format long 
%N=100
 
   M=2; N=4



A=rand(N,N)
I=eye(N);


%B=rand(N,M); 
B=I(:,1:M)

I1=eye(M);
tf=1


% WB0=zeros(N,N);
% ot=0.025;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
% end
%  C0=pinv(WB0); 
 
 

syms t 

p=expm(A*t);

pf=expm(A*tf);

%xf=pf*x0;


% syms s 
% AA=s*I-A;
% invAA=inv(A);
% p=ilaplace(invAA,t);
% p=simplify(p);




v=0.0005;


iterations=10
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  



 
for ii=1:iterations
   
    
    p=expm(A*t);
    pf=expm(A*tf);
    
  WB=double(int(p*B*B'*p',t,0,tf));
 C=pinv(WB);  
 
%  WB0=zeros(N,N);
% ot=0.025;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
% end
%  C=pinv(WB0); 
 
Cost1(ii)=trace(C*pf*pf');


V=eig(A);

BB=fliplr(vander(V));

% if rank(diag(V)) < N;
% 
%  break
% else   
CC1=diag(expm(diag(V)*t));
alpha_t=inv(BB)*CC1;

    
CC2=diag(expm(diag(V)*tf));
alpha_tf=inv(BB)*CC2;



Cpf=C*pf;
CB=Cpf*pf'*C*p*B*B';

F11=zeros(N,N);
F22=zeros(N,N);


for i=1:N-1
    
    
F1=zeros(N,N);
F2=zeros(N,N);


    for k=1:i
        
        F1=A'^(k-1)*CB*A'^(i-k)+F1;
        
        F2=A'^(k-1)*Cpf*A'^(i-k)+F2;
    end
    
   F11=-2*alpha_t(i+1)*F1+F11;
   
   F22=2*alpha_tf(i+1)*F2+F22;
   
end 


F=int(F11,t,0,tf)+F22;
 F=double(F); 
   
       
% F111=zeros(N,N);
%    ot=0.025; 
%    for kk=1:tf/ot
%      F111=subs(F11,t,ot*kk)*ot+F111;
%    end
%    
%    F=F111+ double(F22);
     
%F=F/norm(F);

%F3=(trace(A'*A)-N)*4*A; % gradient of norm function

%F3=F3/norm(F3);

%F0=reshape((eye(N*N)-F3(:)*pinv(F3(:)))*F(:),N,N); 




A=A-v*F/norm(F);


A=(sqrt((N+0.01)/trace(A'*A)))*A; 





ii  
end

  plot(Cost1,'r-*')
  
%  subs(x,'t',tf)

AL=zeros(N,N);
for kk=1:N
    
   AL= alpha_tf(kk)*A^(kk-1)+AL
end