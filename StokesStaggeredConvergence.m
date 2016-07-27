clear all
close all

%2D Stokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the L2 error with the actual analytical solution for g = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testPoints = 15;

L2EP = zeros(testPoints,1);
L2EU = zeros(testPoints,1);
L2EV = zeros(testPoints,1);
d = zeros(testPoints,1);

width = 1.0;
height = 1.0;
g = 0;
mu = 2.0;
p0 = 200;


% Compute the L2 error
%2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
for i = 4:testPoints
    
    %function [ P U V X Y numYCells numXCells d] = StokesStaggered(g, numYCells, p0, mu, toGraph, height, width )
    [ p u v x y numYCells numXCells delta] = StokesStaggered(g, i, p0, mu, g, height, width);
    
    d(i) = delta;
    
%     for j=1:numXCells
%         for k = 1:numYCells
%             L2EP(i) = L2EP(i) + (p(k,j) - (x(k,j) * ((100.0 - p0) / width) + p0))^2;
%             L2EU(i) = L2EU(i) + (u(k,j) - ((1/(2*mu)) * ((100.0 - p0) / width) * y(k,j) * (y(k,j) - height)))^2;
%             L2EV(i) = L2EV(i) + (v(k,j) - 0)^2;
%         end
%     end
        

    L2EP(i) = sum(sum((p - (x * ((100.0 - p0) / width) + p0)).^2));
    L2EU(i) = sum(sum((u - ((1/(2*mu)) * ((100.0 - p0) / width) * y) .* (y - height)).^2));
    L2EV(i) = sum(sum((v - 0).^2));
    
    
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
%title('L2 Error for P (Pressure), analytical, g=0');
xlabel('\Delta x');
ylabel('Discrete L2 Error (p)');

figure(8)
loglog(d,L2EU,'-',d,d.^2,'--');
%title('L2 Error for U (Horizontal Velocity), analytical, g=0');
xlabel('\Delta x');
ylabel('Discrete L2 Error (u)');

figure(9)
loglog(d,L2EV,'-',d,d.^2,'--');
%title('L2 Error for V (Vertical Velocity), analytical, g=0');
xlabel('\Delta x');
ylabel('Discrete L2 Error (v)');

[ p u v x y numYCells numXCells delta] = StokesStaggered(g, i, p0, mu, 1, height, width);

figure(10)
surf(x,y,((1/(2*mu)) * ((100 - p0) / width) * y .* (y - height)));
title('u component, analytical');
xlabel('x');
ylabel('y');

figure(11)
surf(x,y,u - ((1/(2*mu)) * ((100 - p0) / width) * y .* (y - height)));
title('u component, computational - analytical');
xlabel('x');
ylabel('y');







