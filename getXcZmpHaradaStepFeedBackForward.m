function [ Xchn,Zmpn,Xch,Zmp,y,ydot,xcdotdes, t] = getXcZmpHaradaStepFeedBackForward(steplengthChange,steptimeChange, xG1, DesEig, x0)
%UNTITLED2 Summary of this function goes here 
%   Detailed explanation goes here

tF = 4;

t = 0:0.001:tF;

g = 9.8;
zG = 0.8; % Center of Mass constant height
zzmp=0;

Tc = sqrt(g/(zG-zzmp));
% Steps time duration

T10 = 0.7; % First Single support step duration (1st segment)
T21 = 0.1; % First Double support step duration (2nd segment)
T32 = 0.7; % Second Single support step duration (3rd segment)
T43 = 0.1; % Second Double support step duration (4th segment)
T54 = 0.7 + 5; % Third Single support step duration (5th segment) 
                % These 20 seconds are added since calculus is anticipative


steplength1 = 0.3; % First step length
steplength2 = 0.3; % Second step length


t0 = 0.0; % Initial time
t1 = t0 + T10;
t2 = t1 + T21;
t3 = t2 + T32;
t4 = t3 + T43;
t5 = t4 + T54 ;



% Desired CoM 

%xG1=0;     %initial desired value for CoM 
xG5=steplength2+steplength1;      %final desired value for CoM




% Terms and coefficients used in the desired ZMP trajectory

a01 = 0; %constant term
%a11=0;  %coef to get constant trajectory
a11=0.05;   %coef to get Linear trajectory
a15=0;             %linear term segment 5



a02 = a01+a11*T10; %constant term segment 2
a12 = a01;        %linear term segment 2
a13=a11;          %linear term segment 3

a03 = steplength1; %constant term segment 3
a04 = a03+a13*T32;  %constant term segment 4
a05 = a03 + steplength2;  %constant term segment 5

a22 = 3/T21^2*(a03 - a02); %quadric term segment 2

a32 = 2/T21^3*(a02 - a03); %cubic term segment 2

a14 = a12;                   %linear term segment 4
a24 = 3/T43^2*(a05 - a04);   %quadric term segment 4
a34 = 2/T43^3*(a04 - a05);   %cubic term segment 4



Zmp = (a01 + a11*(t - t0)).*(stepfun(t,t0)-stepfun(t,t1)) +...
        (a02 + a12*(t - t1) + a22*(t-t1).^2 + a32*(t-t1).^3).*...
        (stepfun(t,t1)-stepfun(t,t2)) + ...
        (a03 + a13*(t - t2)).*(stepfun(t,t2)-stepfun(t,t3)) + ...
        (a04 + a14*(t - t3) + a24*(t-t3).^2 + a34*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4)) + ...
        (a05 + a15*(t - t4)).*(stepfun(t,t4)-stepfun(t,t5));

       
       
% Terms and coefficients used in the CoM trajectory

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


%Matrices as defined in Harada 2004
% L coefficient matrix construction
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
             
 
D=[xG1-A01;                                      % Init CoM Pos
   A02-A01*(t1-t0)^0-A11*(t1-t0)^1;              % Pos Cont in t1
   A12-A11*1*(t1-t0)^0;                          % Vel Cont in t1
          
   A03-A02*(t2-t1)^0-A12*(t2-t1)^1-A22*(t2-t1)^2-A32*(t2-t1)^3; % Pos Cont in t2
   A13-A12*1*(t2-t1)^0-A22*2*(t2-t1)^1-A32*3*(t2-t1)^2;         % Vel Cont in t2
          
   A04-A03*(t3-t2)^0-A13*(t3-t2)^1;         % Pos Cont in t3
   A14-A13*1*(t3-t2)^0;                     % Vel Cont in t3
          
   A05-A04*(t4-t3)^0-A14*(t4-t3)^1-A24*(t4-t3)^2-A34*(t4-t3)^3;     % Pos Cont in t4
   A15-A14*1*(t4-t3)^0-A24*2*(t4-t3)^1-A34*3*(t4-t3)^2;             % Vel Cont in t4
   xG5-A05*(t5-t4)^0-A15*(t5-t4)^1];                                % Pos Cont in t5
 
     
     
VW=l\D; %=(V1,W1,V2,W2,V3,W3,V4,W4,V5,W5)'





Xch = (VW(1,1)*cosh(Tc*(t-t0))+VW(2,1)*sinh(Tc*(t-t0))+A01*(t-t0).^0+A11*(t-t0).^1).*...
     (stepfun(t,t0)-stepfun(t,t1))+...
    (VW(3,1)*cosh(Tc*(t-t1))+VW(4,1)*sinh(Tc*(t-t1))+A02*(t-t1).^0+A12*(t-t1).^1+A22*(t-t1).^2+A32*(t-t1).^3).*...
     (stepfun(t,t1)-stepfun(t,t2))+...
    (VW(5,1)*cosh(Tc*(t-t2))+VW(6,1)*sinh(Tc*(t-t2))+A03*(t-t2).^0+A13*(t-t2).^1).*...
     (stepfun(t,t2)-stepfun(t,t3))+...
     (VW(7,1)*cosh(Tc*(t-t3))+VW(8,1)*sinh(Tc*(t-t3))+A04*(t-t3).^0+A14*(t-t3).^1+A24*(t-t3).^2+A34*(t-t3).^3).*...
     (stepfun(t,t3)-stepfun(t,t4))+...
     (VW(9,1)*cosh(Tc*(t-t4))+VW(10,1)*sinh(Tc*(t-t4))+A05*(t-t4).^0+A15*(t-t4).^1).*...
     (stepfun(t,t4)-stepfun(t,t5));

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [0, 1; Tc^2, 0];
B = [0; -Tc^2];
C = [1, 0];
D = 0;
 
 
num=[0.01,1];
den=[1,0];
SysD = tf(den,num);
 

xcdes = Xch;
 
xcdotdes = lsim(SysD,xcdes,t);
 
SysA = ss(A,B,C,D);
% DesEig = [-5, -10];
F = -place(A,B,DesEig);
SysAs = ss(A+B*F,B,C,D);  % Closed-loop system
 
 
% x0 = [0, 0];
 
 
xd = [xcdes; xcdotdes'];
  
u = Zmp - F*xd;
y = lsim(SysAs, u,t,x0);
ydot=lsim(SysD,y,t);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This second part is only needed if a step change occurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Terms and coefficients used in the desired ZMP trajectory
% New terms since t4 has changed in t4n
% We first compute the new ZMP deriving from the new t4n
% and then add the variation (it's just a choice)

t4n = t4 + steptimeChange; % new t4n
T43n = t4n - t3; % and therefore new time interval
a05n=a05+a15*(t4n-t4)+steplengthChange;

xG5n=a05n;% new final value of CoM


A15n=a15;
A05n=a05n;


l0= [1, 0,-2/Tc^2,0,0;                                                  %xzmpT3
    1,T43n,T43n^2-2/Tc^2,T43n^3-6*T43n/Tc^2,T43n^4-12*T43n^2/Tc^2;      %xzmpT4n
    0,1,0,-6/Tc^2,0;                                                    %xzmpdotT3
    0,1,2*T43n,3*T43n^2-6/Tc^2,4*T43n^3-2*12*T43n/Tc^2;                 %xzmpdotT4n
    
    ];
l11=[1,0,0,0,0;
    0,1,0,0,0;    
    1,T43n,T43n^2,T43n^3,T43n^4; 
    0,1,2*T43n,3*T43n^2,4*T43n^3];
l12=[1,0,0,0;
    0,Tc,0,0;
    cosh(Tc*T43n),sinh(Tc*T43n),-1,0;
    Tc*sinh(Tc*T43n),Tc*cosh(Tc*T43n),0,-Tc];


Mm1n= [0 0 cosh(Tc*(t5-t4n)) sinh(Tc*(t5-t4n))]; 

ln=[l0,zeros(4,4); % continuity position and velocity at t3 and t4n
    l11,l12;       % coef V W for new CoM
    zeros(1,5),Mm1n]; % coef V W for new CoM

             
%connection at t3 and t4n;
ZmpT3=a03 + a13*(t3 - t2);   % position new Zmp at t3
ZmpT4n=a05n;                 % position new Zmp at t4n
ZmpdotT3=a13;                % velocity new Zmp at t3
ZmpdotT4n=a15;               % velocity new Zmp at t4n


xG3dot = Tc*VW(5,1)*sinh(Tc*(t3-t2))+Tc*VW(6,1)*cosh(Tc*(t3-t2))+A13*(t3-t2).^0;
xG3=VW(5,1)*cosh(Tc*(t3-t2))+VW(6,1)*sinh(Tc*(t3-t2))+A03*(t3-t2).^0+A13*(t3-t2).^1;




Dn=[ZmpT3;     %continuity position new Zmp at t3
    ZmpT4n;    %continuity position new Zmp at t4n
    ZmpdotT3;  %continuity velocity new Zmp at t3
    ZmpdotT4n; %continuity velocity new Zmp at t4n
    xG3;       %continuity position new CoM at t3
    xG3dot;    %continuity velocity new CoM at t3
    A05n;      %continuity position new CoM at t4n
    A15;       %continuity velocity new CoM at t4n
    xG5n-A05n*(t5-t4n)^0-A15n*(t5-t4n)^1]; %final value of new CoM at t5

VWn=ln\Dn; % (A04n,A14n,A24n,A34n,A44n,V4n,W4n,V5n,W5n)

A04n=VWn(1,1);
A14n=VWn(2,1);
A24n=VWn(3,1);
A34n=VWn(4,1);
A44n=VWn(5,1);

%new CoM
Xchn = (VW(1,1)*cosh(Tc*(t-t0))+VW(2,1)*sinh(Tc*(t-t0))+A01*(t-t0).^0+A11*(t-t0).^1).*...
     (stepfun(t,t0)-stepfun(t,t1))+...
    (VW(3,1)*cosh(Tc*(t-t1))+VW(4,1)*sinh(Tc*(t-t1))+A02*(t-t1).^0+A12*(t-t1).^1+A22*(t-t1).^2+A32*(t-t1).^3).*...
     (stepfun(t,t1)-stepfun(t,t2))+...
    (VW(5,1)*cosh(Tc*(t-t2))+VW(6,1)*sinh(Tc*(t-t2))+A03*(t-t2).^0+A13*(t-t2).^1).*...
     (stepfun(t,t2)-stepfun(t,t3))+...
     (VWn(6,1)*cosh(Tc*(t-t3))+VWn(7,1)*sinh(Tc*(t-t3))+A04n*(t-t3).^0+A14n*(t-t3).^1+A24n*(t-t3).^2+A34n*(t-t3).^3+A44n*(t-t3).^4).*...
      (stepfun(t,t3)-stepfun(t,t4n))+...
       (VWn(8,1)*cosh(Tc*(t-t4n))+VWn(9,1)*sinh(Tc*(t-t4n))+A05n*(t-t4n).^0+A15n*(t-t4n).^1).*...
       (stepfun(t,t4n)-stepfun(t,t5));  
    
%Coef of ZMP for new branch between t3 and t4n 
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


    
end

