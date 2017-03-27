close all
load('00.mat');
tf=1 
N=4;
M=2;
figure
  plot(Cost1,'r-*','LineWidth', 1.5,'MarkerSize',8)
  hold on 
    
     plot(Cost20,'b-o','LineWidth', 1.5,'MarkerSize',8)  
     hold on 
    legend('E (t_f,B)', 'E_x(t_f,B)')
      set(gca, 'LineWidth', 1.5);
     xlabel('Iteration number ','FontName','Times New Roman','FontWeight','bold');
 ylabel( 'Control cost','FontName','Times New Roman','FontWeight','bold');
  set(gca,'FontName','Times New Roman','FontWeight','bold')
  
    export_fig cost1_alpha_0.eps -painters -transparent
    
    
   figure
    plot(Cost30,'gx-.','LineWidth', 1.5,'MarkerSize',8)
   
    %myzoom([0.3,0.4,0.5,0.4],[0,20,0,50])
     
    legend( 'E_u(t_f,B)')
    set(gca, 'LineWidth', 1.5);
     xlabel('Iteration number ','FontName','Times New Roman','FontWeight','bold');
 ylabel('Control cost','FontName','Times New Roman','FontWeight','bold');
  set(gca,'FontName','Times New Roman','FontWeight','bold')
    export_fig cost2_alpha_0.eps -painters -transparent
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
ez1(ii)=ezplot(x00(ii),[0, 1]);
hold on 
axis auto
 end
 
 
 hold on
 for ii=1:N
ez2(ii)=ezplot(x(ii),[0, 1]);
hold on 
axis auto
 end
% set(findobj(gca,'color','b'),'color','r')
% legend('y1=t^2','y2=t^2-16')
delete(findall(findall(gcf,'Type','axe'),'Type','text'))
legend([ez1(1),ez2(1)],'x(t,B_0)', 'x(t,B^*)',4)
for ii=1:N
set(ez2(ii),'color',[1 0 0])
set(ez1(ii), 'LineWidth', 1.5,'LineStyle','--');
%set(ez2(ii), 'LineWidth', 1.5);

end
 
 xlabel('t','FontName','Times New Roman','FontWeight','bold');
 ylabel('The state of x ','FontName','Times New Roman','FontWeight','bold');
  set(gca, 'LineWidth', 1.5);
  set(gca,'FontName','Times New Roman','FontWeight','bold')
 export_fig x_alpha_0.eps -painters -transparent
 
 
 
 
fg=figure
 for ii=1:M
ez1(ii)=ezplot(u00(ii),[0, 1]);
hold on 
axis auto
 end
 
 
 hold on
 for ii=1:M
ez2(ii)=ezplot(u(ii),[0, 1]);
hold on 
axis auto
 end
% set(findobj(gca,'color','b'),'color','r')
% legend('y1=t^2','y2=t^2-16')
delete(findall(findall(gcf,'Type','axe'),'Type','text'))
legend([ez1(1),ez2(1)],'u(t,B_0)', 'u(t,B^*)',4)
for ii=1:M
set(ez2(ii),'color',[1 0 0])
 %set(ez2(ii), 'LineWidth', 1.5);
 set(ez1(ii), 'LineWidth', 1.5,'LineStyle','--');
end
  set(gca, 'LineWidth', 1.5);
  set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('t','FontName','Times New Roman','FontWeight','bold');
 ylabel('The input u ','FontName','Times New Roman','FontWeight','bold');
 export_fig u_alpha_0.eps -painters -transparent
 