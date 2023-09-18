clc; close all; clear all
%%%%%% PARAMETERS %%%%%%%%%%%%%%
Nx=250;                 % No of grids in x direction
dx=10;                  % grid increment
Ny=250;                 % No of grids in y direction
dy=10;                  % grid increment
%%%%%%%SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRCNX = Nx/2;           %source position    
SRCNY = Ny/2;                                      
T = 1000;               % total time
dt=.002;                % time increment
f = 10;                 % frequency of source
t0 = 0.15;              % source term
%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = zeros(Nx,Ny);               % velocity 

xw = 200;


for i = 1:Nx;
    for j = 1:Ny;
        if i < Nx/4;
            c(i,j) = 2500;
        else
            c(i,j) = 2000;
        end
    end
end



r = 2000;               % density
const2 = (1/(dx*r));    % constant used in time updation
%%%%%%%%%%%%%FIELD VARIABLES%%%%%%%%%%%%%%%
p0 = zeros(Nx,Ny);          p2 = zeros(Nx,Ny);  %initializing field variables for pressure 
vx0 = zeros(Nx,Ny);        vx2 = zeros(Nx,Ny);  %initializing field variables for velocity in x direction
vy0 = zeros(Nx,Ny);        vy2 = zeros(Nx,Ny);  %initializing field variables for velocity in y direction
%%%%%%%%%%%%TIME UPDATION%%%%%%%%%%%%%%%%%
for k = 0:T
    t = k*dt; 
     p0(SRCNX,SRCNY) = p0(SRCNX,SRCNY) +(ricker(f,t,t0));
    for i=2:Nx-1
       for j=2:Ny-1 
           const1 = (r*((c(i,j))^2)/dx);  % constant used in time updation

           
       p2(i,j) = ((2.0*dt)/(2.0))*((p0(i,j)*((2.0)/(2.0*dt)))-...
           ((const1)*(vx0(i+1,j)-vx0(i,j)+vy0(i,j+1)-vy0(i,j))));
       vx2(i,j) = ((2.0*dt)/(2.0))*((vx0(i,j)*((2.0)/(2.0*dt)))-...
           ((const2)*(p2(i,j)-p2(i-1,j))));
       vy2(i,j) = ((2.0*dt)/(2.0))*((vy0(i,j)*((2.0)/(2.0*dt)))-...
           ((const2)*(p2(i,j)-p2(i,j-1))));
       end
    end
   
    p0 = p2;
    vx0 = vx2;
    vy0 = vy2;
    %%%%%%%%PLOT%%%%%%%%%%%%%%
    figure(1);imagesc(p2);
    pause(.01)
end




