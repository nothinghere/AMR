
% ZMP ---> CoM
% Bounded solution
% 2 steps
% ZMP = 5 segments: Constant/Cubic/Constant/Cubic/Constant
% 
% Last update: 22/06/2014
% 


clear all


rossomattone = [ 1.0 0.3 0.3 ];
bluoceano    = [ 0.0 0.5 0.9 ];
gialloocra   = [ 0.75 0.75 0.0 ];

% Cost-Lin-Cost parameters
% aij = a_i^{(j)} Harada's notation

tf = 4;

t = 0:0.001:4;

g = 9.8;
zh = 0.8;   % Center of Mass constant height

om = sqrt(g/zh);

% Steps time duration

stepss1 = 0.7; % First Single support step duration (1st segment)
stepds2 = 0.1; % First Double support step duration (2nd segment)
stepss3 = 0.7; % Second Single support step duration (3rd segment)
stepds4 = 0.1; % Second Double support step duration (4th segment)
stepss5 = 0.7; % Third Single support step duration (5th segment)


steplength1 = 0.3; % First step length
steplength2 = 0.3; % Second step length


t0 = 0.0;          % Initial time
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

a01 = 0;
a11=0.05;
a03 = steplength1;


a05 = a03 + steplength2;

a02 = a01+a11*T10;
a12 = a01;
a22 = 3/T21^2*(a03 - a02);

a13=a11;
a32 = 2/T21^3*(a02 - a03);

a04 = a03+a13*T32;
a14 = a12;
a24 = 3/T43^2*(a05 - a04);
a34 = 2/T43^3*(a04 - a05);
a15=a11;
% Desired ZMP trajectory
Zmp = (a01 + a11*(t - t0)).*(stepfun(t,t0)-stepfun(t,t1)) +...
        (a02 + a12*(t - t1) + a22*(t-t1).^2 + a32*(t-t1).^3).*...
        (stepfun(t,t1)-stepfun(t,t2)) + ...
        (a03 + a13*(t - t2)).*(stepfun(t,t2)-stepfun(t,t3)) + ...
        (a04 + a14*(t - t3) + a24*(t-t3).^2 + a34*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4)) + ...
        (a05 + a15*(t - t4)).*(stepfun(t,t4)-stepfun(t,t5));
       % a05*(stepfun(t,t4)-stepfun(t,t5));



% Constant segment from t0 to t1

Xc1 = x_constant( t0, t1, a01, om, t);

% Second segment from t1 to t2 --- Spline
% various terms are computed separately and added at the end

% Cubic terms

Cubic2 =x_cubic( t1,t2,T21,om,t );
% Quadratic terms

 Quad2 =   x_quadric( t1,t2,T21,om,t );
% Linear terms
 Lin2 =x_linear( t1,t2,T21,om,t );
% Constant term

Cost2 = x_constant(t1,t2,a02,om,t);

% Final Spline from t1 to t2 (second segment)
Xc2 = a32*Cubic2 + a22*Quad2 + a12*Lin2 + Cost2;

% Third segment from t2 to t3: Constant
Xc3  =x_constant(t2,t3,a03,om,t);

% Fourth segment from t3 to t4 --- Spline

% Cubic terms

Cubic4 = x_cubic(t3,t4,T43,om,t);

% Quadratic terms   

Quad4 = x_quadric(t3,t4,T43,om,t);


% Linear terms

Lin4 = x_linear(t3,t4,T43,om,t);

% Linear terms
Cost4 =x_constant(t3,t4,1,om,t);
         

% Final spline from t3 to t4 --- Fourth segment 
Xc4 = a34*Cubic4 + a24*Quad4 + a14*Lin4 + a04*Cost4;

         

% Fifth segment: Constant from t4 to t5
Xc5  = x_constant(t4,t5,a05,om,t);

         
% Xcstar final

Xc = Xc1 + Xc2 + Xc3 + Xc4 + Xc5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% This second part is only needed if a step change occurs         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Terms and coefficients used in the desired ZMP trajectory
% New terms since t4 has changed in t4n
% We first compute the new ZMP deriving from the new t4n
% and then add the variation (it's just a choice)

t4n = t4 + 0.2;          % new t4n
T43n = t4n - t3;         % and therefore new time interval

a24 = 3/T43n^2*(a05 - a04);
a34 = 2/T43n^3*(a04 - a05);


% New desired ZMP trajectory
Zmpn = (a01 + a11*(t - t0)).*(stepfun(t,t0)-stepfun(t,t1)) +...
        (a02 + a12*(t - t1) + a22*(t-t1).^2 + a32*(t-t1).^3).*...
        (stepfun(t,t1)-stepfun(t,t2)) + ...
        (a03 + a13*(t - t2)).*(stepfun(t,t2)-stepfun(t,t3)) + ...
        (a04 + a14*(t - t3) + a24*(t-t3).^2 + a34*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4n)) + ...        
        (a05 + a15*(t - t4n)).*(stepfun(t,t4n)-stepfun(t,t5));
        %a05*(stepfun(t,t4n)-stepfun(t,t5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% ZMP variation in tvar
% New trajectory with 3 segments
% from tvar to t3: constant = a03v
% from t3 to t4n: spline (t4n may be different from t4)
% from t4n on, a03v + varstep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 

% Here the change occurs at half time of the third segment
% it can obviously be different

tvar = t2 + 0.5*stepss3; % change at half of second step
varstep = steplength2/2; % variation amplitude, just an example


a03v = 0;
a04v = a03v;
a14v = 0;
a24v = 3/T43n^2*varstep;
a34v = -2/T43n^3*varstep;
a05v = a03v + varstep;

ZmpVar = a03v*(stepfun(t,tvar)-stepfun(t,t3)) + ...
        (a04v + a14v*(t - t3) + a24v*(t-t3).^2 + a34v*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4n)) + ...
        a05v*stepfun(t,t4n);

% CoM corresponding variation

% Constant term from tvar to t3 (is equal to 0)
Xcv1  = x_constant(tvar,t3,a03v,om,t);
% Spline segment from t3 to t4 
% Cubic terms

Cubic4v = x_cubic(t3,t4n,T43n,om,t);

% Quadratic terms   

Quad4v = x_quadric(t3,t4n,T43n,om,t);

% Linear terms

Lin4v = x_linear(t3,t4n,T43n,om,t);

% Constant term
Cost4v = x_constant(t2,t4n,1,om,t);


Xcv2 = a34v*Cubic4v + a24v*Quad4v + a14v*Lin4v + a04v*Cost4v;

% Constant segment from t4n to t5
Xcv3 = x_constant(t4n,t5,a05v,om,t);

XcVar = Xcv1 + Xcv2 + Xcv3;


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(2)
A0 = plot(t,Zmp); 
%     t, ZmpVar, t, Zmpn+ZmpVar, t,Xc, t,XcVar, t, Xc+XcVar

A1 = line([t0,t0],[0,0.8]);
A2 = line([t1,t1],[0,0.8]);
A3 = line([t2,t2],[0,0.8]);
A4 = line([t3,t3],[0,0.8]);
A5 = line([t4,t4],[0,0.8]);

A6 = line([tvar,tvar],[0,0.8]);
A7 = line([t4n,t4n],[0,0.8]);


grid
set(A1(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A2(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A3(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A4(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A5(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
set(A6(1),'Color',rossomattone,'LineWidth',0.5,'LineStyle','-.');
set(A7(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');


set(A0(1),'Color',gialloocra,'LineWidth',1.2);
set(A0(2),'Color',bluoceano,'LineWidth',1.2,'LineStyle','--'); % Delta ZMP
set(A0(5),'Color',rossomattone,'LineWidth',1.2,'LineStyle','--'); % Delta Xc
set(A0(4),'Color',bluoceano,'LineWidth',1.2);
set(A0(3),'Color',rossomattone,'LineWidth',1.2);
set(A0(6),'Color',rossomattone,'LineWidth',1.2);

xlabel('t (sec)', 'FontName','cmr','FontSize',12)


