clc; close all; clear all
%%%%%% PARAMETERS %%%%%%%%%%%%%%
Nx=250;                  % No of grids in x direction
dx=10;                   % grid increment
x = 0:dx:(Nx-1)*dx;      % x axis for ploting
%%%%%%%SOURCE %%%%%%%%%%%%%%%%%%
SRCNX = 1;              %source position    
SRCNY = Nx/2;           %source position
T=1000                  % Total time
dt=.001;                % time increment
f = 10;                 % frequency of source
t0 = 0.1;               % source term
%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%
c = 2000;           %velocity 
r = 2300;           %density
const1 = (r*(c^2)/dx); %constant used in time updation
const2 = (1.0/(dx*r)); %constant used in time updation
%%%%%%%%%%%%%FIELD VARIABLES%%%%%%%%%%%%%%%
p0 = zeros(Nx,1);       p2 = zeros(Nx,1);   %initialize field variables for pressure
v0 = zeros(Nx,1);       v2 = zeros(Nx,1);   %initialize field variables for velocity
%%%%%%%%%%%%%%PML boundary conditions%%%%%% 
d0=150;                                         %maximum value of damping 
w = 15;                                         % width of PML in terms of nodes
a = (1:w);                                      %constant1(Array) 
b = fliplr(a);                                  %constant2(Array) 
sigma = zeros(Nx,1);                            %initialize array for damping
sigma(1:w) = d0*((b/w).^2);                     %Assign damping
sigma((Nx-(w-1)):Nx) = d0*((a/w)).^2;           %Assign damping
%%%%%%%%%%%%TIME UPDATING%%%%%%%%%%%%%%%%%
for k = 0:T
    t = k*dt; 
    
    for i=2:Nx-2
       v2(i,1) = ((2.0*dt)/(2.0+(sigma(i,1)*dt)))*((v0(i,1)*((2.0-(sigma(i,1)*dt))/(2.0*dt)))-...
           ((const2)*(p0(i+1,1)-p0(i,1)))); 
       p2(i,1) = ((2.0*dt)/(2.0+(sigma(i,1)*dt)))*((p0(i,1)*((2.0-(sigma(i,1)*dt))/(2.0*dt)))-...
           ((const1)*(v2(i,1)-v2(i-1,1))));
       
    end
    p2(SRCNY,SRCNX) =  (ricker(f,t,t0));
    p0 = p2;
    v0 = v2;
    
    p0(2,1)=0;
    p2(2,1)=0;
    %%%%%%%%PLOT%%%%%%%%%%%%%%
    plot(x,p2);xlabel('distance(m)');ylabel('pressure (Pa)');   axis([min(x) max(x) -1 1]);
    pause(.01) 
end