
g = 9.8;
zh = 0.8;   % Center of Mass constant height

om = sqrt(g/zh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [0, 1; om^2, 0];
B = [0; -om^2];
C = [1, 0];
D = 0;
 
SysD = tf([1, 0],[0.01, 1]);
 
Deltaxc = xcVar.*stepfun(t,tvar); % Azzerato il contributo fino a t*
xcdes = xc+Deltaxc;
 
xcdotdes = lsim(SysD,xcdes,t);
 
SysA = ss(A,B,C,D);
DesEig = [-5, -10];
F = -place(A,B,DesEig);
SysAs = ss(A+B*F,B,C,D);  % Closed-loop system
 
 
x0 = [0, 0];
 
 
xd = [xcdes; xcdotdes'];
 
figure(2) 
 
u = xzmp + zmpVar.*stepfun(t,tvar) - F*xd;
y = lsim(SysAs, u,t,x0);
C1 = plot(t,xcdes,t,y), grid
set(C1(1),'Color',rossomattone,'LineWidth',1.2);
set(C1(2),'Color',bluoceano,'LineWidth',1.2);
legend('x_c^{des}', 'x_c','Location','Best')
 
 
figure(3)
plot(t, y' - xcdes), grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%