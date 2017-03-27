x1=[0,      0.05 , 0.1 ,    0.2 ,  0.3 ,   0.4 ,    0.5,    0.6 ,  0.7 ,   0.8  ,  0.9, 0.95 ,   1  ];
  y1= [ 0.4216 ,   0.6061  ,  0.7253 ,  0.8855,    0.9928,   1.0668,   1.1138,   1.1378 ,  1.1258,    1.1448,  1.2144,  1.1584,   1.2412];
 	y2=[	16.4983,  10.7648,    9.4574,   8.4382 ,   8.1078,     7.9298 ,    7.8888 ,  7.8813 ,  7.8683 ,   7.8392 ,   7.8049 ,  7.7589 , 7.8178];
        
    %Coefficients (with 95% confidence bounds):
       p1 =       7.472  %(7.383, 7.561)
       p2 =      0.4555  %(0.4119, 0.4991)
       q1 =     0.02759  %(0.02486, 0.03033)
        
        
        A=polyfit(x1,y1,1);
z=polyval(A,x1);
l1=plot(x1,y1,'r*')
hold on
plot(x1,z,'k')
hold on


%  A2=polyfit(x,y2,4);
 %z2=polyval(A2,x);

 l2=plot(x1,y2,'b*')
 hold on
 
     

         x=0:0.001:1;
  f= (p1*x + p2) ./ (x + q1);
 plot(x,f,'g')
 
 
 legend([l1,l2],'E_x(t_f,B^*)', 'E_u(t_f,B^*)')
  set(gca, 'LineWidth', 1.5);
    set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('\alpha','FontName','Times New Roman','FontWeight','bold');
 ylabel('Control cost','FontName','Times New Roman','FontWeight','bold');
 export_fig cost_fitting_2.eps -painters -transparent 