function [ error ] = compute_error(A,B,C,D, DesEig, Zmp, ZmpVar, t, tvar, xcdes, xcdotdes, x0)

xd = [xcdes; xcdotdes'];

SysD = tf([1, 0],[0.01, 1]);
SysA = ss(A,B,C,D);
F = -place(A,B,DesEig);
SysAs = ss(A+B*F,B,C,D);

u = Zmp + ZmpVar .* stepfun(t,tvar) - F*xd;

y = lsim(SysAs, u,t,x0);
xcdot = lsim(SysD,y',t);

errors_vectors = [y' - xcdes; xcdot' - xcdotdes']';

error = sqrt(sum(abs(errors_vectors).^2,2));

end