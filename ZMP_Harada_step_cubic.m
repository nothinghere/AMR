
% ZMP ---> CoM
% Harada solution
% 2 steps
% ZMP = 5 segments: Linear/Cubic/Linear/Cubic/linear
%
% Last update: 22/06/2014
%


clear all
close all

rossomattone = [ 1.0 0.3 0.3 ];
bluoceano = [ 0.0 0.5 0.9 ];
gialloocra = [ 0.75 0.75 0.0 ];

% Cost-Lin-Cost parameters
% aij = a_i^{(j)} Harada's notation

tf = 4;

t = 0:0.001:4;

g = 9.8;
zG = 0.8; % Center of Mass constant height
zzmp=0;

Tc = sqrt(g/(zG-zzmp));
% Steps time duration

stepss1 = 0.7; % First Single support step duration (1st segment)
stepds2 = 0.1; % First Double support step duration (2nd segment)
stepss3 = 0.7; % Second Single support step duration (3rd segment)
stepds4 = 0.1; % Second Double support step duration (4th segment)
stepss5 = 0.7; % Third Single support step duration (5th segment)


steplength1 = 0.3; % First step length
steplength2 = 0.3; % Second step length


t0 = 0.0; % Initial time
t1 = t0 + stepss1;
t2 = t1 + stepds2;
t3 = t2 + stepss3;
t4 = t3 + stepds4;
t5 = t4 + stepss5 + 20; % These 20 seconds are added since calculus is anticipative

% Terms and coefficients used in the desired ZMP trajectory
T10 = t1 - t0;
T21 = t2 - t1;
T32 = t3 - t2;
T43 = t4 - t3;
T54 = t5 - t4;

a01 = 0;
a11=0;
%a11=0.05;



a02 = a01+a11*T10;
a12 = a01;
a13=a11;

a03 = steplength1;
a04 = a03+a13*T32;
a05 = a03 + steplength2;

a22 = 3/T21^2*(a03 - a02);

a32 = 2/T21^3*(a02 - a03);

a14 = a12;
a24 = 3/T43^2*(a05 - a04);
a34 = 2/T43^3*(a04 - a05);

a15=a11;



Zmp = (a01 + a11*(t - t0)).*(stepfun(t,t0)-stepfun(t,t1)) +...
        (a02 + a12*(t - t1) + a22*(t-t1).^2 + a32*(t-t1).^3).*...
        (stepfun(t,t1)-stepfun(t,t2)) + ...
        (a03 + a13*(t - t2)).*(stepfun(t,t2)-stepfun(t,t3)) + ...
        (a04 + a14*(t - t3) + a24*(t-t3).^2 + a34*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4)) + ...
        (a05 + a15*(t - t4)).*(stepfun(t,t4)-stepfun(t,t5));

       
       
%i=0..1; exposants n=1
%j=0..3; temps m=3

%Aij

A11=a11;
A01=a01;

A32=a32;
A22= a22;
A12= a12+1/Tc^2*(1+1)*(1+2)*A32;
A02= a02+1/Tc^2*(0+1)*(0+2)*A22;



A13=a13;
A03=a03;


A34=a34;
A24= a24;
A14= a14+1/Tc^2*(1+1)*(1+2)*A34;
A04= a04+1/Tc^2*(0+1)*(0+2)*A24;

A15=a15;
A05=a05;

xG1=a01;
xG2=a02;
xG3=a03;
xG4=a04;
xG5=a05;



o=[1 0 0 0 ];
l1 =li(t0,t1,Tc);
l2=li(t1,t2,Tc);
l3=li(t2,t3,Tc);
l4=li(t3,t4,Tc);


Mm1= [0 0 cosh(Tc*(t5-t4)) sinh(Tc*(t5-t4))];


l=[o,zeros(1,6);
    l1,zeros(2,6);
    zeros(2,2),l2,zeros(2,4);
    zeros(2,4),l3,zeros(2,2);
    zeros(2,6),l4;
    zeros(1,6),Mm1];
             
 
D=[xG1-A01;A02-A01*(t1-t0)^0-A11*(t1-t0)^1;
          A12-A11*1*(t1-t0)^0;
          
          A03-A02*(t2-t1)^0-A12*(t2-t1)^1-A22*(t2-t1)^2-A32*(t2-t1)^3;
          A13-A12*1*(t2-t1)^0-A22*2*(t2-t1)^1-A32*3*(t2-t1)^2;
          
          A04-A03*(t3-t2)^0-A13*(t3-t2)^1;
          A14-A13*1*(t3-t2)^0;
          
          A05-A04*(t4-t3)^0-A14*(t4-t3)^1-A24*(t4-t3)^2-A34*(t4-t3)^3;
          A15-A14*1*(t4-t3)^0-A24*2*(t4-t3)^1-A34*3*(t4-t3)^2;
         xG5-A05*(t5-t4)^0-A15*(t5-t4)^1];
 
     
     
VW=l\D;




Xc = (VW(1,1)*cosh(Tc*(t-t0))+VW(2,1)*sinh(Tc*(t-t0))+A01*(t-t0).^0+A11*(t-t0).^1).*...
     (stepfun(t,t0)-stepfun(t,t1))+...
    (VW(3,1)*cosh(Tc*(t-t1))+VW(4,1)*sinh(Tc*(t-t1))+A02*(t-t1).^0+A12*(t-t1).^1+A22*(t-t1).^2+A32*(t-t1).^3).*...
     (stepfun(t,t1)-stepfun(t,t2))+...
    (VW(5,1)*cosh(Tc*(t-t2))+VW(6,1)*sinh(Tc*(t-t2))+A03*(t-t2).^0+A13*(t-t2).^1).*...
     (stepfun(t,t2)-stepfun(t,t3))+...
     (VW(7,1)*cosh(Tc*(t-t3))+VW(8,1)*sinh(Tc*(t-t3))+A04*(t-t3).^0+A14*(t-t3).^1+A24*(t-t3).^2+A34*(t-t3).^3).*...
     (stepfun(t,t3)-stepfun(t,t4))+...
     (VW(9,1)*cosh(Tc*(t-t4))+VW(10,1)*sinh(Tc*(t-t4))+A05*(t-t4).^0+A15*(t-t4).^1).*...
     (stepfun(t,t4)-stepfun(t,t5));
%  
%    plot(t,Zmp);
%     hold on    
%       plot(t,Xc,'r');grid
      


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This second part is only needed if a step change occurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Terms and coefficients used in the desired ZMP trajectory
% New terms since t4 has changed in t4n
% We first compute the new ZMP deriving from the new t4n
% and then add the variation (it's just a choice)

t4n = t4 + 0.2; % new t4n
T43n = t4n - t3; % and therefore new time interval
a05n=a05+a15*(t4n-t4)+steplength2/2;


A15n=a15;
A05n=a05n;

xG5n=A05n;
l0= [1, 0,-2/Tc^2,0,0;                                      %xzmp1
    1,T43n,T43n^2-2/Tc^2,T43n^3-6*T43n/Tc^2,T43n^4-12*T43n^2/Tc^2;      %xzmp2
    0,1,0,-6/Tc^2,0;                                      %xzmpp1
    0,1,2*T43n,3*T43n^2-6/Tc^2,4*T43n^3-2*12*T43n/Tc^2;   %xzmpp2
    
    ];%connection to the final value
l11=[1,0,0,0,0;
    0,1,0,0,0;    
    1,T43n,T43n^2,T43n^3,T43n^4;
    0,1,2*T43n,3*T43n^2,4*T43n^3];
l12=[1,0,0,0;
    0,Tc,0,0;
    cosh(Tc*T43n),sinh(Tc*T43n),-1,0;
    Tc*sinh(Tc*T43n),Tc*cosh(Tc*T43n),0,-Tc];


Mm1n= [0 0 cosh(Tc*(t5-t4n)) sinh(Tc*(t5-t4n))];

ln=[l0,zeros(4,4);
    l11,l12;
    zeros(1,5),Mm1n];

             
%connection at t3;
ZmpT3=Zmp(round(t3/0.001)+1);
ZmpT4n=A05n;

xG3p = Tc*VW(5,1)*sinh(Tc*(t3-t2))+Tc*VW(6,1)*cosh(Tc*(t3-t2))+A13*(t3-t2).^0;

xG3=Xc(round(t3/0.001)+1);
Dn=[ZmpT3;
    ZmpT4n;
    a13; %continuity velocity new Zmp initial
    a15; %continuity velocity new Zmp  final
    xG3;xG3p;A05n;A15n; 
          xG5n-A05n*(t5-t4n)^0-A15n*(t5-t4n)^1];

    
  
VWn=ln\Dn;
 

A04n=VWn(1,1);
A14n=VWn(2,1);
A24n=VWn(3,1);
A34n=VWn(4,1);
A44n=VWn(5,1);
Xcn = (VW(1,1)*cosh(Tc*(t-t0))+VW(2,1)*sinh(Tc*(t-t0))+A01*(t-t0).^0+A11*(t-t0).^1).*...
     (stepfun(t,t0)-stepfun(t,t1))+...
    (VW(3,1)*cosh(Tc*(t-t1))+VW(4,1)*sinh(Tc*(t-t1))+A02*(t-t1).^0+A12*(t-t1).^1+A22*(t-t1).^2+A32*(t-t1).^3).*...
     (stepfun(t,t1)-stepfun(t,t2))+...
    (VW(5,1)*cosh(Tc*(t-t2))+VW(6,1)*sinh(Tc*(t-t2))+A03*(t-t2).^0+A13*(t-t2).^1).*...
     (stepfun(t,t2)-stepfun(t,t3))+...
     (VWn(6,1)*cosh(Tc*(t-t3))+VWn(7,1)*sinh(Tc*(t-t3))+A04n*(t-t3).^0+A14n*(t-t3).^1+A24n*(t-t3).^2+A34n*(t-t3).^3+A44n*(t-t3).^4).*...
      (stepfun(t,t3)-stepfun(t,t4n))+...
       (VWn(8,1)*cosh(Tc*(t-t4n))+VWn(9,1)*sinh(Tc*(t-t4n))+A05n*(t-t4n).^0+A15n*(t-t4n).^1).*...
       (stepfun(t,t4n)-stepfun(t,t5));  
    
   
a44n=A44n;
a34n=A34n;
a24n=A24n-1/Tc^2*(2+1)*(2+2)*A44n;
a14n=A14n-1/Tc^2*(1+1)*(1+2)*A34n;
a04n=A04n-1/Tc^2*(0+1)*(0+2)*A24n;
% New desired ZMP trajectory
 Zmpn = (a01 + a11*(t - t0)).*(stepfun(t,t0)-stepfun(t,t1)) +...
        (a02 + a12*(t - t1) + a22*(t-t1).^2+a32*(t-t1).^3).*...
        (stepfun(t,t1)-stepfun(t,t2)) + ...
        (a03 + a13*(t - t2)).*(stepfun(t,t2)-stepfun(t,t3)) + ...
        (a04n + a14n*(t - t3) + a24n*(t-t3).^2+a34n*(t-t3).^3+a44n*(t-t3).^4 ).*...
        (stepfun(t,t3)-stepfun(t,t4n)) + ...
        (a05n + a15*(t - t4n)).*(stepfun(t,t4n)-stepfun(t,t5));

    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
A0 = plot(t,0, t, 0, t, Zmp, t,Xc, t,Zmpn, t, Xcn);

A1 = line([t0,t0],[0,0.8]);
A2 = line([t1,t1],[0,0.8]);
A3 = line([t2,t2],[0,0.8]);
A4 = line([t3,t3],[0,0.8]);
A5 = line([t4,t4],[0,0.8]);
%A6 = line([tvar,tvar],[0,0.8]);
A7 = line([t4n,t4n],[0,1]);


grid
set(A1(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A2(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A3(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A4(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A5(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
%set(A6(1),'Color',rossomattone,'LineWidth',0.5,'LineStyle','-.');
set(A7(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');


set(A0(1),'Color',bluoceano,'LineWidth',1.2);
set(A0(2),'Color',bluoceano,'LineWidth',1.2,'LineStyle','--'); % Delta ZMP
%set(A0(5),'Color',rossomattone,'LineWidth',1.2,'LineStyle','--'); % Delta Xc
%set(A0(4),'Color',bluoceano,'LineWidth',1.2);
%set(A0(3),'Color',rossomattone,'LineWidth',1.2);
%set(A0(6),'Color',rossomattone,'LineWidth',1.2);

xlabel('t (sec)', 'FontName','cmr','FontSize',12)


     
     
