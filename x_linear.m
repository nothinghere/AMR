function [ Linj ] = x_linear( ti,tj,Tji,om,t )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Linj1  = om/2*(1/om^2*exp(om*(t-ti)) - ...
        exp(om*(t-tj))*(Tji/om + 1/om^2)).*...
        (1 - stepfun(t,ti));
Linj2 = om/2*(2*(t-ti)/om  -...
         exp(om*(t-tj))*(Tji/om + 1/om^2) + ...
         1*exp(-om*(t-ti))/om^2).*...
         (stepfun(t,ti)-stepfun(t,tj));
Linj3 = om/2*(1/om^2*exp(-om*(t-ti)) + ...
        exp(-om*(t-tj))*(Tji/om - 1/om^2)).*...
        stepfun(t,tj);
Linj = Linj1 + Linj2 + Linj3;


end

