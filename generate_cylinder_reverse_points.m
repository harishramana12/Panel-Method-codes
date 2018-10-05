%%
clc;clear;close all;
r = 1.0;
N_panel = 9;
deltheta = (2*pi)/(N_panel-1);
theta= (pi-(deltheta/2)):(deltheta):(2*pi+pi-(deltheta/2));
theta_deg = theta*(180/pi);
%theta(N_panel+1) = theta(1);

x = r*cos(theta);
y = r*sin(theta);

figure(1);hold on;
plot(x,y,'--bx');

%%
fil = 'circular_cylinder.txt';
fid = fopen(fil,'w');

for i = 1:1:length(x)
    
    fprintf(fid,'%d %d %d',x(i),y(i),theta(i));
    fprintf(fid,'\n');
    
end
