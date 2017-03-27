clear all
close all
 clc
format long 
%N=100
 
   M=2 ; N=6
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







v=0.008
alpha=0.6



iterations=150
Cost1=zeros(iterations,1);
Cost2=zeros(iterations,1);
Cost3=zeros(iterations,1);
Cost4=zeros(iterations,1); 
Cost20=zeros(iterations,1);
Cost30=zeros(iterations,1);
WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 

us=-B'*qf'*C*xf;
u=subs(us,'s',t);
x=simple(p*x0+int(q*B*us, s, 0, t));
x00=x;
u00=u;



for ii=1:150
 % us=-B'*qf'*C*xf;
% u=subs(us,'s',t);
% x=simple(p*x0+int(q*B*us, s, 0, t));


qu=int(qft*B*u,t,0,tf);

pq=int(int(x'*q*B*B'*qf',s,0,t),t,0,tf);

%F=int(int(q'*x*x0'*pf*C*qf*B,s,0,t),t,0,tf);
F1=-2*int(int(q'*x*xf'*C*qf*B,s,0,t),t,0,tf)-2*int(int(qf'*C*xf*x'*q*B,s,0,t),t,0,tf)+2*int(p'*C*pq'*xf'*C*p*B,t,0,tf)+2*int(p'*C*xf*pq*C*p*B,t,0,tf);

F1=double(F1);

F2=-2*int(qft'*C*xf*u',t,0,tf)+2*int(p'*C*qu*xf'*C*p*B,t,0,tf)+2*int(p'*C*xf*qu'*C*p*B,t,0,tf);
F2=double(F2);

F=(1-alpha)*F1+alpha*F2;

B1=B-v*(I-B*pinv(B))*F/norm(F);


B=sqrt(trace(B1'*B1)/trace(B1'*B1*B1'*B1))*B1;


WB=double(int(p*B*B'*p',t,0,tf));
C=pinv(WB); 
us=-B'*qf'*C*xf;
u=subs(us,'s',t);
x=simple(p*x0+int(q*B*us, s, 0, t));

Cost1(ii)=(1-alpha)*int(x'*x,0,tf)+alpha*int(u'*u,0,tf);
Cost2(ii)=(1-alpha)*int(x'*x,0,tf);
Cost3(ii)=alpha*int(u'*u,0,tf);
Cost20(ii)=int(x'*x,0,tf);
Cost30(ii)=int(u'*u,0,tf);
Cost4(ii)=trace((B'*B-I1)'*(B'*B-I1));  
ii  
end

  figure
  plot(Cost1,'r-*','LineWidth', 1.5,'MarkerSize',8)
  hold on 
    plot(Cost30,'gx-.','LineWidth', 1.5,'MarkerSize',8)
   
     hold on 
    legend('E(t_f,B)', 'E_u(t_f,B)')
      set(gca, 'LineWidth', 1.5);
     xlabel('Iteration number ','FontName','Times New Roman','FontWeight','bold');
 ylabel( 'Energy cost','FontName','Times New Roman','FontWeight','bold');
  set(gca,'FontName','Times New Roman','FontWeight','bold')
  
    %export_fig cost1_alpha_10.eps -painters -transparent
    
    
   figure
      plot(Cost20,'b-o','LineWidth', 1.5,'MarkerSize',8)  
   
    %myzoom([0.3,0.4,0.5,0.4],[0,20,0,50])
     
    legend( 'E_x(t_f, B)')
    set(gca, 'LineWidth', 1.5);
     xlabel('Iteration number ','FontName','Times New Roman','FontWeight','bold');
 ylabel('Energy cost','FontName','Times New Roman','FontWeight','bold');
  set(gca,'FontName','Times New Roman','FontWeight','bold')
    %export_fig cost2_alpha_10.eps -painters -transparent
% % eval(subs(x,'t',tf))
% xxx=[];
% for iii=0.001:0.001:tf
%     xx=eval(subs(x,'t',iii));
%     xxx=[xxx,xx];
% end
% xxx0=[];
% for iii=0.001:0.001:tf
%     xx0=eval(subs(x00,'t',iii));
%     xxx0=[xxx0,xx0];
% end
%    
%  for ii=1:N
% x0=plot(xxx0(ii,:),'-');
% hold on
%  end
%  hold on
%  for ii=1:N
% x2=plot(xxx(ii,:),'r-');
% hold on
%  end
% 
%  legend( [x0,x2],'x0',' x*');
 
 %figure 2
fg=figure
 for ii=1:N
ez1=ezplot(x00(ii),[0, tf]);
hold on 
axis auto
 end
 
 
 hold on
 for ii=1:N
ez2(ii)=ezplot(x(ii),[0, tf]);
hold on 
axis auto
 end
% set(findobj(gca,'color','b'),'color','r')
% legend('y1=t^2','y2=t^2-16')
delete(findall(findall(gcf,'Type','axe'),'Type','text'))
legend([ez1,ez2(1)],' x_0', 'x^*',0)
for ii=1:N
set(ez2(ii),'color',[1 0 0])

end
 
 xlabel('t','FontName','Times New Roman','FontWeight','bold');
 ylabel('The state of x ','FontName','Times New Roman','FontWeight','bold');
  set(gca, 'LineWidth', 1.5);
  set(gca,'FontName','Times New Roman','FontWeight','bold')
 %export_fig x_alpha_10.eps -painters -transparent
 
 
 
 
fg=figure
 for ii=1:M
ez1=ezplot(u00(ii),[0, tf]);
hold on 
axis auto
 end
 
 
 hold on
 for ii=1:M
ez2(ii)=ezplot(u(ii),[0, tf]);
hold on 
axis auto
 end
% set(findobj(gca,'color','b'),'color','r')
% legend('y1=t^2','y2=t^2-16')
delete(findall(findall(gcf,'Type','axe'),'Type','text'))
legend([ez1,ez2(1)],'u_0', 'u^*',0)
for ii=1:M
set(ez2(ii),'color',[1 0 0])

end
 set(gca, 'LineWidth', 1.5);
  set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('t','FontName','Times New Roman','FontWeight','bold');
 ylabel('The input u ','FontName','Times New Roman','FontWeight','bold');
 %export_fig u_alpha_10.eps -painters -transparent
 
 
 

%  semilogy(Cost20) 
%  hold on  
%  semilogy(Cost30)  
%  hold on  
%  semilogy(Cost1)