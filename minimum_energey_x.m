clear all
 clc
format long 
%N=100
 
   M=2; N=4
 A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];  % stem
%A(1,N)=1;
 
% 
%   A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)]; 
%    A(1,6)=1;   %circle
%  

%   A=zeros(6,6);
%   A(2,1)=1;
%     A(3,2)=1;
%         A(4,3)=1;
%  A(5,1)=1;
%   A(6,5)=1;
%  
%  A2=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
%  A2(1,3)=1;
% 
%  A=[A zeros(3,3); zeros(3,3) A2];
%  
%  A(6,2)=1;
%  M=1; N=6

% A good example to show the minimization of  the maximum matching paths
% corresponds to the optimal energy control.
 
%A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];
% A=csvread('E:\network\Electronic networks\s208_st.csvadjacencymatrix.csv');
% A=[0 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 1 0 0 0 1 1
% 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0
% 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 1 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 1 0
% 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 1 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0
% 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0
% 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1
% 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0];
% %   A=A';
%    N=size(A);
%    M=8
%  M=29
%A(N,1)=1;

% A=[0 1 0 0 0 0 0 0  0 0 0 0;
%    0 0 1 0 0 0 1 0  0 0 0  0;
%    0 0 0 1 0 0 0 0  0 0 0  0;
%    0 0 0 0 1 0 0 0  0 0 0  0;
%    0 0 0 0 0 1 0 0  0 0 0  0;
%    0 0 0 0 0 0 0 0  0 0 0  0;
%    0 0 0 0 0 0 0 1  0 0 0  0;
%    0 0 0 0 0 0 0 0  1 0 0  0;
%    0 0 0 0 0 0 0 0  0 1 0  0;
%    0 0 0 0 0 0 0 0  0 0 1  0;
%    0 0 0 0 0 0 0 0  0 0 0  1
%    0 0 0 0 0 0 0 0  0 0 0 0];
% 
% 
% 
%  A=A';
 
% A = zeros(N,N);
% A(2,1) = 1;
% A(3,2) = 1;
% A(4,3) = 1;
% A(5,4) = 1;
% A(6,5) = 1;
% A(7,6) = 1;
% A(8,7) = 1;
% A(9,3) = 1;
% A(10,9) = 1;
% A(11,10) = 1;
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

xf=pf*x0;





v=0.01;



iterations=20
Cost1=zeros(iterations,1);

   
for ii=1:iterations
    
    
WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 

    
pq0=int(q*B*B'*qf',s,0,t);  
x=simple(p*x0-pq0*C*xf);

% 
% u=-B'*qf'*C*xf;
% xx=simple(p*x0+int(q*B*u, s, 0, t));

pq=double(int(x'*pq0,t,0,tf));

%F=int(int(q'*x*x0'*pf*C*qf*B,s,0,t),t,0,tf);
F=-2*int(int(q'*x*xf'*C*qf*B,s,0,t),t,0,tf)-2*int(int(qf'*C*xf*x'*q*B,s,0,t),t,0,tf)+2*int(p'*C*pq'*xf'*C*p*B,t,0,tf)+2*int(p'*C*xf*pq*C*p*B,t,0,tf);

 
F=double(F);
B=B-v*F/norm(F);

Cost1(ii)=int(x'*x,t,0,tf);
ii  
end

  plot(Cost1,'r-*')

