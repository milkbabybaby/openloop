clear all
 clc
format long 
%N=100
 
   M=2; N=3
 A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];  % stem

N=length(A);

B=rand(N,M);
x0=rand(N,1);
 
 

tf=1


% WB0=zeros(N,N);
% ot=0.025;
% for k=1:tf/ot
%     WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
% end
%  C0=pinv(WB0); 
 
 

syms t s real
p=simple(expm(A*t));

pf=expm(A*tf);
q=simple(expm(A*(t-s)));
qf=simple(expm(A*(tf-s)));
qft=simple(expm(A*(tf-t)));

xf=pf*x0;







v=0.03;



iterations=10
Cost1=zeros(iterations,1);

   
for ii=1:iterations
    
    
WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 

u=-B'*qft'*C*xf;
qu=int(qft*B*u,t,0,tf);

Cost0=int(u'*u,0,tf);
%F=int(int(q'*x*x0'*pf*C*qf*B,s,0,t),t,0,tf);
F=-2*int(qft'*C*xf*u',t,0,tf)+2*int(p'*C*qu*xf'*C*p*B,t,0,tf)+2*int(p'*C*xf*qu'*C*p*B,t,0,tf);

 
F=double(F);
B=B-v*F/norm(F);

Cost1(ii)=int(u'*u,0,tf);
ii  
end

  plot(Cost1,'r-*')
  
%  subs(x,'t',tf)