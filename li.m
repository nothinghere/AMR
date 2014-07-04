function [ li ] = li( ti, tj, Tc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
li = [cosh(Tc*(tj-ti)),sinh(Tc*(tj-ti)),-1,0;
    Tc*sinh(Tc*(tj-ti)),Tc*cosh(Tc*(tj-ti)), 0,-Tc];

end

