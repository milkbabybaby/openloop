clear all
 clc
format long 
%N=100
 
   M=2; N=3
 A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];  % stem

N=length(A);


B=rand(N,M);
I=eye(N);
I1=eye(M);
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







v=0.0005;
alpha=0.1


iterations=20
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  
for ii=1:iterations
    
    
WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 

us=-B'*qf'*C*xf;
u=subs(us,'s',t);
x=simple(p*x0+int(q*B*us, s, 0, t));


qu=int(qft*B*u,t,0,tf);

pq=int(int(x'*q*B*B'*qf',s,0,t),t,0,tf);

%F=int(int(q'*x*x0'*pf*C*qf*B,s,0,t),t,0,tf);
F1=-2*int(int(q'*x*xf'*C*qf*B,s,0,t),t,0,tf)-2*int(int(qf'*C*xf*x'*q*B,s,0,t),t,0,tf)+2*int(p'*C*pq'*xf'*C*p*B,t,0,tf)+2*int(p'*C*xf*pq*C*p*B,t,0,tf);

F1=double(F1);

F2=-2*int(qft'*C*xf*u',t,0,tf)+2*int(p'*C*qu*xf'*C*p*B,t,0,tf)+2*int(p'*C*xf*qu'*C*p*B,t,0,tf);
F2=double(F2);

F=F1+alpha*F2;

B1=B-v*(I-B*B')*F/norm(F);


B=sqrt(trace(B1'*B1)/trace(B1'*B1*B1'*B1))*B1;

Cost1(ii)=int(x'*x,0,tf)+alpha*int(u'*u,0,tf);
Cost2(ii)=int(x'*x,0,tf);
Cost3(ii)=int(u'*u,0,tf);
Cost4(ii)=trace((B'*B-I1)'*(B'*B-I1));  
ii  
end

  plot(Cost1,'r-*')
  
%  subs(x,'t',tf)