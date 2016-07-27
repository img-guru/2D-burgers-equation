close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2D Burger's equation as described by this tutorial:
% http://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lesso
% ns/10_Step_8.ipynb
% Step 8 of the CFD 12-step tutorial by Prof. Lorena Barba
% and that tutorial (but in 1D): 
% http://www.thevisualroom.com/burgers_equation.html
% and that tutorial:
% http://www.thevisualroom.com/2D_burgers_equation.html
% Attempted by Heba A. Shalaby 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx = 41;
ny = 41;
nt = 600;
c = 0.1;
dx = 2/(nx-1);
dy = 2/(ny-1);
sigma = 0.1;
dt = sigma*dx*dy/c; % from the Lorena Barba tutorial
% dt = 0.5 /(nt-1);     % from the visualroom tutorial

x = linspace(0,2,nx);
y = linspace(0,2,ny);

u = ones(ny,nx);
v = ones(ny,nx);
un = ones(ny,nx);
vn = ones(ny,nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u(0.5/dy:1/dy+1 , 0.5/dx:1/dx+1) = 2;     
v(0.5/dy:1/dy+1 , 0.5/dx:1/dx+1) = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(x,y);

%axis([0 2 0 2 0 2 0 2])

s1 = surf(X,Y,v);                       % plot diffusion in x-direction
title('2D Burgers Equation');
colormap(jet);
% view(90,90);
view(80,30);
set(gca,'clim',[1 2]);                  % to have colormap change in the whole jet spectrum between the 1 and 2 axis values 
axis([0,2,0,2,0,2]);                    % limit the axis for viewing

opengl('software')
set(gca,'nextplot','replacechildren','visible','off')
f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,jet(256),'nodither');
im(1,1,1,nt) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:nt
    un = u;
    vn = v;
    
    for i=2:nx-1
        for j=2:ny-1
            
            % find u^n+1_(i,j) in x-direction 
            % forward difference in time & backward different in space
            u(i,j) = un(i,j)-(un(i,j)*dt/dx*(un(i,j)-un(i-1,j)))  ...
                            -(vn(i,j)*dt/dy*(un(i,j)-un(i,j-1)))...
                            + c*dt/(dx^2)*(un(i+1,j)-2*un(i,j)+ un(i-1,j)) ...
                            + c*dt/(dy^2)*(un(i,j+1)-2*un(i,j)+ un(i,j-1));
                        
            % find u^n+1_(i,j) in y-direction 
            % forward difference in time & backward different in space
             v(i,j) = vn(i,j)-(un(i,j)*dt/dx*(vn(i,j)-vn(i-1,j)))  ...
                             -(vn(i,j)*dt/dy*(vn(i,j)-vn(i,j-1))) ...
                             + c*dt/(dx^2)*(vn(i+1,j)-2*vn(i,j)+ vn(i-1,j)) ...
                             + c*dt/(dy^2)*(vn(i,j+1)-2*vn(i,j)+ vn(i,j-1));
                             
        
            u(1,j) = 1;
            u(i,1) = 1;
            u(nx,j) = 1;
            u(i,ny) = 1;
                
            v(1,j) = 1;
            v(i,1) = 1;
            v(nx,j) = 1;
            v(i,ny) = 1;
            
        end
    end
    
    set(s1,'ZData',v);            % update 3D plot animation (x-direction)
    pause(0.001);                 % pause to control animation speed
    drawnow;
    f = getframe(gcf);
    im(:,:,1,t) = rgb2ind(f.cdata,map,'nodither');
   
  
end

imwrite(im,jet(256),'burgers-equation.gif','gif','DelayTime',0,'LoopCount',inf, 'BackgroundColor', 0) %g443800



