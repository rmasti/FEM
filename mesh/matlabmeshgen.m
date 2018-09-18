% Mesh generator for ring
% Written By Robert Masti
% Output the coordinates for the c++ code
clc, clear, close all


ny = 6;
nx = 8;

% r and theta coordinates

theta = linspace(pi,0, nx+1);
rl = linspace(0.25, 0.5, ny+1);


[r,t] = meshgrid(rl, theta);

x=r.*cos(t);
y=r.*sin(t);

xlower = 0.25*cos(theta);
ylower = 0.25*sin(theta);

figure
plot(x,y,'*r');

xlin=reshape(x,[],1);
ylin=reshape(y,[],1);

plot(xlin,ylin,'r');

fileID=fopen('debugMatlab.msh','w');

fprintf(fileID, '%i\n',(nx+1)*(ny+1));
fprintf(fileID, '%i\n',nx+1);
fprintf(fileID, '%i\n',ny+1);

outvec = vertcat(xlin,ylin);
fprintf(fileID, '%d\n', outvec);
length(outvec)

