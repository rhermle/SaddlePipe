clear all
close all

%2D Stokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the L2 error without having the analytical solution
%This will use a numerical solution with small grid spacing
%as the "analytical" solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testPoints = 15;

L2EP = zeros(testPoints,1);
L2EU = zeros(testPoints,1);
L2EV = zeros(testPoints,1);
d = zeros(testPoints,1);

width = 2;
height = 1;
g = 9.8;
mu = 2;
p0 = 200;

%Create the "analytical" solution (small grid spacing)
%function [ P U V X Y numYCells numXCells d] = StokesStaggered(g, numYCells, p0, mu, toGraph, height, width )

[ P U V X Y NUMYCELLS NUMXCELLS DELTA] = StokesStaggered(g, 30, p0, mu, 1, height, width);


% Compute the L2 error
%2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
for i = 4:testPoints
    
    [ p u v x y numYCells numXCells delta] = StokesStaggered(g, i, p0, mu, 0, height, width);
    
    d(i) = delta;
    
    for j=1:numXCells       
        for k = 1:numYCells
            L2EP(i) = L2EP(i) + (p(k,j) - interp2(X,Y,P,x(k,j),y(k,j), 'cubic'))^2;
            L2EU(i) = L2EU(i) + (u(k,j) - interp2(X,Y,U,x(k,j),y(k,j), 'cubic'))^2;
            L2EV(i) = L2EV(i) + (v(k,j) - interp2(X,Y,V,x(k,j),y(k,j), 'cubic'))^2;
        end
    end
    L2EP(i) = sqrt(L2EP(i) / (numXCells * numYCells));
    L2EU(i) = sqrt(L2EU(i) / (numXCells * numYCells));
    L2EV(i) = sqrt(L2EV(i) / (numXCells * numYCells));
end

d = d(4:end);
L2EP = L2EP(4:end);
L2EU = L2EU(4:end);
L2EV = L2EV(4:end);

figure(7)
loglog(d,L2EP,'-',d,d.^2,'--');
title('L2 Error for P (Pressure), analytical, g=0');

figure(8)
loglog(d,L2EU,'-',d,d.^2,'--');
title('L2 Error for U (Horizontal Velocity), analytical, g=0');

figure(9)
loglog(d,L2EV,'-',d,d.^2,'--');
title('L2 Error for V (Vertical Velocity), analytical, g=0');






