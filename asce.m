load('matlab2.mat');

  figure
  plot(Cost1,'r-*','LineWidth', 1.5,'MarkerSize',6)
  hold on 
    plot(Cost30,'gx-.','LineWidth', 1.5,'MarkerSize',6)
   
     hold on 
    legend('E(t_f,B)', 'E_u(t_f,B)')
      set(gca, 'LineWidth', 1.5);
     xlabel('Iteration number ','FontName','Times New Roman','FontWeight','bold');
 ylabel( 'Control cost','FontName','Times New Roman','FontWeight','bold');
  set(gca,'FontName','Times New Roman','FontWeight','bold')
  
    export_fig cost_u_ascending.eps -painters -transparent
    
    
   figure
      plot(Cost20,'b-o','LineWidth', 1.5,'MarkerSize',6)  
   
    %myzoom([0.3,0.4,0.5,0.4],[0,20,0,50])
     
    legend( 'E_x(t_f, B)')
    set(gca, 'LineWidth', 1.5);
     xlabel('Iteration number ','FontName','Times New Roman','FontWeight','bold');
 ylabel('Control cost','FontName','Times New Roman','FontWeight','bold');
  set(gca,'FontName','Times New Roman','FontWeight','bold')
    export_fig cost_x_ascending.eps -painters -transparent