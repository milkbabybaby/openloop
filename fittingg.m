clear;
close all
x1=[	 0   , 0.05  , 0.1   ,  0.2 , 0.3  , 0.4   , 0.5 ,  0.6 , 0.7  , 0.8 ,  0.9  , 0.95 , 1  ];
		
y1=[		0.4427   , 0.5132 ,0.6140   , 0.6422  , 0.7641 , 0.8314, 0.9315 ,  1.0318 , 1.1071,1.1375 , 1.2584, 1.2709,1.3847];

 
y2=[9.1607 , 7.2288 , 6.2677 ,6.1101  , 5.5602  ,5.3554, 5.0669 , 4.9806, 4.7776,4.7544 ,4.7083 , 4.6735, 4.6492 ];

 %syms x
     
%        a =       3.105 ; 
%        b =      -18.09 ; 
%        c =       6.049;
%        d =     -0.2912; 
%         f= a*exp(b*x) + c*exp(d*x);


       p1 =       4.137  %(3.957, 4.317)
       p2 =        1.04  %(0.846, 1.235)
       q1 =      0.1136  %0.09076, 0.1364)
         

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
 export_fig cost_fitting_1.eps -painters -transparent 