clear all
close all
 clc
format long 
%N=100
 
   M=3 ; N=8
 A=[zeros(1,N-1)  0; eye(N-1) zeros(N-1,1)];  % stem

N=length(A);
I=eye(N);

I1=eye(M);


B=rand(N,M);
B=orth(B);
%B=I(:,1:M)
%B=[-0.666549601221895,0.258757856602330;-0.529674783868792,-0.240268472171118;-0.300205473874438,0.714625909757013;-0.430154537230872,-0.603842067126071]



 x0=ones(N,1);
% x0=[1,2,3,4]'
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







v=0.02
alpha=0.6



iterations=60
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1); 
Cost20=zeros(iterations,1);
Cost30=zeros(iterations,1);










for ii=1:iterations
    
    
        
 WB0=zeros(N,N);
ot=0.02;
for k=1:tf/ot
    WB0=WB0+expm(A*(ot*k))*B*B'*expm(A'*(ot*k))*ot;
end
 C=pinv(WB0); 

 ef=expm(A*tf);
 
     F1=zeros(N,N);
  for kg=1:tf/ot
          F11=zeros(N,N);
          F111=zeros(N,N);
              for tao=ot:ot:ot*kg
                  
                  F11=F11+expm(A*(ot*kg-tao))*B*B'*expm(A'*(tf-tao))*ot;
                  F111=F111+expm(A'*(ot*kg-tao))*(expm(A*ot*kg)-F11*C*ef)*ef'*C*expm(A*(tf-tao))*B*ot;
              end 
         
              F1=F1+F111*ot;
  end 
                  
                  

  
   F2=zeros(N,N);
  for kg=1:tf/ot
          F22=zeros(N,N);
          F222=zeros(N,N);
              for tao=ot:ot:ot*kg
                  
                  F22=F22+expm(A*(tf-tao))*B*B'*expm(A'*(ot*kg-tao))*ot;
                  F222=F222+expm(A'*(tf-tao))*(expm(A'*ot*kg)-ef'*C*F22)*expm(A*(ot*kg-tao))*B*ot;
              end 
         
              F2=F2+F22*ot;
              
 
  end 
                  
                  
 
    
   F3=zeros(N,N);
  for kg=1:tf/ot
         
          F3333=zeros(N,N);
         for  kkg=1:tf/ot
              F33=zeros(N,N);
               F333=zeros(N,N);
              for tao=ot:ot:ot*kkg
                  
                  F33=F33+expm(A*(ot*kkg-tao))*B*B'*expm(A'*(tf-tao))*ot;
                  F333=F333+expm(A*(tf-tao))*B*B'expm(A'*(ot*kkg-tao))*(expm(A*ot*kkg)-F33*C*ef)*ef'*ot;
              end 
         
              F3333=F3333+F3333*ot;
         end 
      F3=F3+ expm(A'*kg*ot)*C*F3333*C*expm(A*kg*ot)*B*ot;     
 
  end
  
  
              
 
     F4=zeros(N,N);
  for kg=1:tf/ot
         
          F4444=zeros(N,N);
         for  kkg=1:tf/ot
              F44=zeros(N,N);
               F444=zeros(N,N);
              for tao=ot:ot:ot*kkg
                  
                  F44=F44+expm(A*(tf-tao))*B*B'*expm(A'*(ot*kg-tao))*ot;
                  F444=F444+ef*(expm(A*ot*kkg)-expm(A'*tf)*C*F44)*expm(A*(kkg*ot-tao))*B*B'*expm(A'*(tf-tao))*ot;
              end 
         
              F4444=F4444+F4444*ot;
         end 
      F4=F4+ expm(A'*kg*ot)*C*F4444*C*expm(A*kg*ot)*B*ot;     
 
  end
  
  
  
              
F01=-2*F1-2*F2+2*F3+2*F4;             
              
              
              
              
 
 
 F02=zeros(N,M); 
 for k1=1:tf/ot
 F02=F02-(expm(A'*(ot*k1))*C*Xf*C*expm(A*(ot*k1)))*B*ot;
end
 
     
              


F=(1-alpha)*F01+alpha*F02;













B1=B-v*(I-B*B')*F/norm(F);


B=sqrt(trace(B1'*B1)/trace(B1'*B1*B1'*B1))*B1;


WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 
us=-B'*qf'*C*xf;
u=subs(us,'s',t);
x=simple(p*x0+int(q*B*us, s, 0, t));

Cost1(ii)=(1-alpha)*int(x'*x,0,tf)+alpha*trace(C*ef*ef');
Cost2(ii)=(1-alpha)*int(x'*x,0,tf);
Cost3(ii)=alpha*int(u'*u,0,tf);
Cost20(ii)=int(x'*x,0,tf);
Cost30(ii)=int(u'*u,0,tf);
Cost4(ii)=trace((B'*B-I1)'*(B'*B-I1));  
ii  
end



