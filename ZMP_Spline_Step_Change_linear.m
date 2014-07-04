
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
green = [0 1 0];

% Cost-Lin-Cost parameters
% aij = a_i^{(j)} Harada's notation

% tf = 4;

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
T54 = t5 - t4;
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

%Xc1 = x_constant( t0, t1, a01, om, t);


Xc1 = a11*x_linear( t0, t1, T10,om, t)+x_constant( t0, t1, a01, om, t);

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
%Xc3  =x_constant(t2,t3,a03,om,t);
Xc3  =a13*x_linear(t2,t3,T32,om,t)+x_constant(t2,t3,a03,om,t);
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
%Xc5  = x_constant(t4,t5,a05,om,t);
Xc5  = a15*x_linear(t4,t5,T54,om,t)+ x_constant(t4,t5,a05,om,t);

         
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
a05n=a05+a15*(t4n-t4);
a24n = 3/T43n^2*(a05n - a04);
a34n = 2/T43n^3*(a04 - a05n);
% 
% a24 = 3/T43n^2*(a05 - a04);
% a34 = 2/T43n^3*(a04 - a05);

% New desired ZMP trajectory
Zmpn = (a01 + a11*(t - t0)).*(stepfun(t,t0)-stepfun(t,t1)) +...
        (a02 + a12*(t - t1) + a22*(t-t1).^2 + a32*(t-t1).^3).*...
        (stepfun(t,t1)-stepfun(t,t2)) + ...
        (a03 + a13*(t - t2)).*(stepfun(t,t2)-stepfun(t,t3)) + ...
        (a04 + a14*(t - t3) + a24n*(t-t3).^2 + a34n*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4n)) + ...        
        (a05n + a15*(t - t4n)).*(stepfun(t,t4n)-stepfun(t,t5));
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
a14v = 0;

a04v = a03v+a14v*(tvar-t2);
a24v = 3/T43n^2*varstep;
a34v = -2/T43n^3*varstep;
a05v = a04v + varstep;

ZmpVar = a03v*(stepfun(t,tvar)-stepfun(t,t3)) + ...
        (a04v + a14v*(t - t3) + a24v*(t-t3).^2 + a34v*(t-t3).^3).*...
        (stepfun(t,t3)-stepfun(t,t4n)) + ...
        (a05v)*stepfun(t,t4n);

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
Cost4v = x_constant(t3,t4n,1,om,t);


Xcv2 = a34v*Cubic4v + a24v*Quad4v + a14v*Lin4v + a04v*Cost4v;

% Constant segment from t4n to t5
Xcv3 = x_constant(t4n,t5,a05v,om,t);

XcVar = Xcv1 + Xcv2 + Xcv3;


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(2)
% A0 = plot(t,Zmp, t, ZmpVar, t, Zmpn+ZmpVar, t,Xc, t,XcVar, t, Xc+XcVar);
% A1 = line([t0,t0],[0,0.8]);
% A2 = line([t1,t1],[0,0.8]);
% A3 = line([t2,t2],[0,0.8]);
% A4 = line([t3,t3],[0,0.8]);
% A5 = line([t4,t4],[0,0.8]);
% 
% A6 = line([tvar,tvar],[0,0.8]);
% A7 = line([t4n,t4n],[0,0.8]);
% 
% 
% grid
% set(A1(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
% set(A2(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
% set(A3(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
% set(A4(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
% set(A5(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
% set(A6(1),'Color',rossomattone,'LineWidth',0.5,'LineStyle','-.');
% set(A7(1),'Color',bluoceano,'LineWidth',1, 'LineStyle',':');
% 
% 
% set(A0(1),'Color',bluoceano,'LineWidth',1.2);
% set(A0(2),'Color',bluoceano,'LineWidth',1.2,'LineStyle','--'); % Delta ZMP
% set(A0(5),'Color',rossomattone,'LineWidth',1.2,'LineStyle','--'); % Delta Xc
% set(A0(4),'Color',bluoceano,'LineWidth',1.2);
% set(A0(3),'Color',rossomattone,'LineWidth',1.2);
% set(A0(6),'Color',rossomattone,'LineWidth',1.2);
% 
% xlabel('t (sec)', 'FontName','cmr','FontSize',12)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Feedforward/feedback implementation%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [0, 1; om^2, 0];
B = [0; -om^2];
C = [1, 0];
D = 0;
 
 
 
SysD = tf([1, 0],[0.01, 1]);
 
Deltaxc = XcVar.*stepfun(t,tvar); % Azzerato il contributo fino a t*
xcdes = Xc+Deltaxc;
 
xcdotdes = lsim(SysD,xcdes,t);
 
SysA = ss(A,B,C,D);
DesEig = [-5, -10];
F = -place(A,B,DesEig);
SysAs = ss(A+B*F,B,C,D);  % Closed-loop system
 
 
x0 = [0, 0];

% try with diff x0 values

x0_1 = [0.1, 0.1];
x0_2 = [0.2, 0.2];
x0_3 = [0.3, 0.3];
x0_4 = [0.4, 0.4];
x0_5 = [0.5, 0.5];
x0_6 = [0.6, 0.6];
x0_7 = [0.7, 0.7];
x0_8 = [0.8, 0.8];
x0_9 = [0.9, 0.9];
 
xd = [xcdes; xcdotdes'];
 
% figure('Name', 'xc vs vc_des') 
 
u = Zmp + ZmpVar .* stepfun(t,tvar) - F*xd;
y = lsim(SysAs, u,t,x0);
xcdot = lsim(SysD,y,t);

% err_01 = abs(lsim(SysAs, u, t, x0_1)' - xcdes);
% err_02 = abs(lsim(SysAs, u, t, x0_2)' - xcdes);
% err_03 = abs(lsim(SysAs, u, t, x0_3)' - xcdes);
% err_04 = abs(lsim(SysAs, u, t, x0_4)' - xcdes);
% err_05 = abs(lsim(SysAs, u, t, x0_5)' - xcdes);
% err_06 = abs(lsim(SysAs, u, t, x0_6)' - xcdes);
% err_07 = abs(lsim(SysAs, u, t, x0_7)' - xcdes);
% err_08 = abs(lsim(SysAs, u, t, x0_8)' - xcdes);
% err_09 = abs(lsim(SysAs, u, t, x0_9)' - xcdes);

% try with diff eigenvalues
% 
DesEig_1 = [-1, -2];
DesEig_2 = [-2, -4];
DesEig_3 = [-3, -6];
DesEig_4 = [-4, -8];
DesEig_5 = [-5, -10];
DesEig_6 = [-6, -12];
% DesEig_7 = [-2, -1];
% DesEig_8 = [-4, -2];
% DesEig_9 = [-6, -3];
% DesEig_10 = [-8, -4];
% DesEig_11 = [-10, -5];
% DesEig_12 = [-12, -6];

err_01 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_1);
err_02 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_2);
err_03 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_3);
err_04 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_4);
err_05 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_5);
err_06 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_6);
err_07 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_7);
err_08 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_8);
err_09 = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0_9);

err_11 = compute_error(A,B,C,D, DesEig_1, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
err_12 = compute_error(A,B,C,D, DesEig_2, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
err_13 = compute_error(A,B,C,D, DesEig_3, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
err_14 = compute_error(A,B,C,D, DesEig_4, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
err_15 = compute_error(A,B,C,D, DesEig_5, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
err_16 = compute_error(A,B,C,D, DesEig_6, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
% err_17 = compute_error(A,B,C,D, DesEig_7, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
% err_18 = compute_error(A,B,C,D, DesEig_8, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
% err_19 = compute_error(A,B,C,D, DesEig_9, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
% err_110 = compute_error(A,B,C,D, DesEig_10, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
% err_111 = compute_error(A,B,C,D, DesEig_11, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);
% err_112 = compute_error(A,B,C,D, DesEig_12, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0);

% errors_vectors = [y' - xcdes; xcdot' - xcdotdes']';

% error = sqrt(sum(abs(errors_vectors).^2,2));

% C1 = plot(t, xcdes, t, y, t, error);

% set(C1(1),'Color',rossomattone,'LineWidth',1.2);
% set(C1(2),'Color',bluoceano,'LineWidth',1.2);
% set(C1(3), 'Color', rossomattone, 'LineWidth', 1.2);
% set(C1(4), 'Color', bluoceano, 'LineWidth', 1.2);
% set(C1(5), 'Color', rossoma, 'LineWidth', 1.2);
% legend('x_c^{des}', 'x_c', 'error', 'Location','Best')

figure('Name', 'x0 vs corresponding errors')

C2 = plot(t, err_01,t, err_02,t, err_03, t, err_04, t, err_05, t, err_06, t, err_07, t, err_08, t, err_09);
legend('[0.2, 0.2]', '[0.3, 0.3]', '[0.4, 0.4]', '[0.5, 0.5]', '[0.6, 0.6]', '[0.7, 0.7]', '[0.8, 0.8]', '[0.9, 0.9]');

figure('Name', 'eigenvalues vs corresponding errors')

C3 = plot(t, err_11, t, err_12, t, err_13, t, err_14, t, err_15, t, err_16);

legend('[-1, -2]','[-2, -4]','[-3, -6]','[-4, -8]','[-5, -10]','[-6, -12]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%