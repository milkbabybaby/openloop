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
%xf=pf*x0;







v=0.02;
alpha=0.1


iterations=80
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1);  

WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 

us=-B'*qf'*C*pf;
u=subs(us,'s',t);
x=simple(p+int(q*B*us, s, 0, t));


 
for ii=1:iterations
    
    



%qu=int(qft*B*u,t,0,tf);

pq=int(int(x'*q*B*B'*qf',s,0,t),t,0,tf);

%F=int(int(q'*x*x0'*pf*C*qf*B,s,0,t),t,0,tf);
F1=-2*int(int(q'*x*pf'*C*qf*B,s,0,t),t,0,tf)-2*int(int(qf'*C*pf*x'*q*B,s,0,t),t,0,tf)+2*int(p'*C*pq'*pf'*C*p*B,t,0,tf)+2*int(p'*C*pf*pq*C*p*B,t,0,tf);

F1=double(F1);

% F2=2*int(p'*C*qu*pf'*C*p*B,t,0,tf);
% F2=double(F2);

F2=-2*int(p'*C*pf*pf'*C*p*B,t,0,tf);
F2=double(F2);

F=F1+alpha*F2;


F3=(2*trace(B'*B)-2*M)*2*B; % gradient of norm function

F3=F3/norm(F3);

F0=reshape((eye(N*M)-F3(:)*pinv(F3(:)))*F(:),N,M); 


B=B-v*F0/norm(F0);


B=(sqrt(1*(M+0.01)/trace(B'*B)))*B; 

WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 

us=-B'*qf'*C*pf;
u=subs(us,'s',t);
x=simple(p+int(q*B*us, s, 0, t));


Cost1(ii)=int(trace(x*x'),0,tf)+alpha*int(trace(u*u'),0,tf);
Cost2(ii)=int(trace(x*x'),0,tf);
Cost3(ii)=int(trace(u*u'),0,tf);
Cost4(ii)=trace((B'*B-I1)'*(B'*B-I1));  
ii  
end

  plot(Cost1,'r-*')
  
%  subs(x,'t',tf)