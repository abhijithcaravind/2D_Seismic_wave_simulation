%  staggered grid FDM code for elastic wave propagation in 2D fourth order approaximation 
% in spatial differencing first order approaximation in temporal differencing(levander 1988 scheme) 
% where a Fracture is modelled is as an equivalent anisotropic medium (Coates and Schoenberg 1995)   
% Output: 1) 'RecX.out' measuring the Vx component of the seismogram
%         1) 'time.out' corresponfding time measurement 
% 

clc; close all; clear all;

%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 250;                 % No of grids in x direction
Ny = 250;                 % No of grids in y direction
dh = 10;                 % grid increment common for both X, Z directions
%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SRCNX = Nx/2;           % source position    
SRCNY = Ny/2;                                      
T = 1200;               % total time
dt = 0.00009;         % time increment in Seconds
f = 40;               % frequency of source in Hz
t0 = 0.03;            % source term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SS = dt/dh;
X=[];

%%%%%%%%%%%%%%%%%%%%% Medium MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vp = zeros(Nx,Ny);                      % velocity of P wave in m/s
Vp(:,:) = 4200;

Vs = zeros(Nx,Ny);                      % velocity of S wave in m/s
Vs(:,:) = 2700;

r = zeros(Nx,Ny);                       % density in Kg/m3
r(:,:) = 2490;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu     = [];              % Initiate matrix for Lame's parameters
lambda = [];
const0 = [];
b = [];                   % 1/density
alpha = 1;                % percentage of displacement jump in fracture
%%%%%%%%%%%%%%%%%%%%%%%Fracture Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = [];
deln = [];
delt = [];
Zn = (alpha*dh)/(r(1,1)*((Vp(1,1).^2)));  % normal fracture compliance 
Zt = (alpha*dh)/(r(1,1)*(Vs(1,1)^2));     %  tangentioal fracture compliance
L = 10;
%%%%%%%%%%%%%%%%%%%%%%% Absorbing Boundary (sponge (cerjan etal 1985))%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G =[];
for i=1: Nx;
    for j =1: Ny;      
      G(i,j) = 1;          
      if  i <= 25       
          G(i,j) =  exp(-(0.015*(25 - i)).^2); 
      end
      if  i > (Nx-24); 
          i1 = Nx+1-i;
           G(i,j) =  exp(-(0.015*(25 - i1)).^2);
      end
            if  j <= 25    
           G(i,j) =  exp(-(0.015*(25 - j)).^2);
      end
      if  j > (Ny-24); 
          j1 = Ny+1-j;
           G(i,j) =  exp(-(0.015*(25 - j1)).^2);
      end                     
    end
end
%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS USED IN TIME UPDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1: Nx;
    for j =1: Ny;
        mu(i,j)     = r(i,j)*(Vs(i,j)^2);         % Lame's parameters evaluate at each grid oints
        lambda(i,j) = r(i,j)*((Vp(i,j)^2)-(2.0*(Vs(i,j))^2));
        const0(i,j) = lambda(i,j)+(2*mu(i,j));    % constant time updation
        b(i,j) = 1/r(i,j);                        % 1/density
        
        R(i,j) = lambda(i,j)/const0(i,j);
        deln(i,j) = (Zn*const0(i,j))/(L+(Zn*const0(i,j)));
        delt(i,j) = (Zt*mu(i,j))/(L+(Zt*mu(i,j)));         
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%CONSTANTS USED IN TIME UPDATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
const0(1,1);
mu(1,1);
const1 = b*SS;                  % constant used in time updation
const2 = (const0)*SS ;          % constant used in time updation
const3 = (lambda)*SS ;          % constant used in time updation
const4 = mu*SS ;                % constant used in time updation

c1 = 9.0/8;                       % constant used for the fourth order approaximation  
c2 = -1.0/24;                     % constant used for the fourth order approaximation  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% delta matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [16,20,24,28,32,36,40];             % Insert fracture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = length(A);
delnmat = zeros(Nx,Ny);
deltmat = zeros(Nx,Ny);
%delnmat(round(Nx/4)+1,:) = deln(round(Nx/4)+1,:)\

for i = 1:B;
 delnmat(round(Nx/4):round(3*Nx/4),round(SRCNY-A(1,i))) = deln(round(Nx/4):round(3*Nx/4),round(SRCNY-A(1,i)));
 deltmat(round(Nx/4):round(3*Nx/4),round(SRCNY-A(1,i))) = delt(round(Nx/4):round(3*Nx/4),round(SRCNY-A(1,i)));  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STABILITY CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = Vp(1,1)*SS;
if k>.7;
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Rec1=[];
% Rec2=[];
Rec3=[];      % initiate Recoder 
% Rec4=[];
% Rec5=[];
% Rec6=[];
% Rec7=[];
 time=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIELD VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1 = zeros(Nx,Ny);      U2 = zeros(Nx,Ny);       %initializing field variables for velocity in x direction
V1 = zeros(Nx,Ny);      V2 = zeros(Nx,Ny);       %initializing field variables for velocity in z direction
Txx1 = zeros(Nx,Ny);    Txx2 = zeros(Nx,Ny);     %initializing field variables for pressure xx
Tzz1 = zeros(Nx,Ny);    Tzz2 = zeros(Nx,Ny);     %initializing field variables for pressure zz
Txz1 = zeros(Nx,Ny);    Txz2 = zeros(Nx,Ny);     %initializing field variables for pressure xz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TIME UPDATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 0:(2*T)
    t1 = k*dt;
    %%%%%%%%%%%%%%%%%%%%%%%%% Insert source %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Txx1(SRCNX,SRCNY) = Txx1(SRCNX,SRCNY) + (ricker(f,t1,t0));
     Tzz1(SRCNX,SRCNY) = Tzz1(SRCNX,SRCNY) + (ricker(f,t1,t0));
     %Txz1(SRCNX,SRCNY) = Txz1(SRCNX,SRCNY) + (ricker(f,t1,t0));
    for i=5:Nx-4
       for j=5:Ny-4 
           
        k1 =  1.0-delnmat(i,j) ;  
        k2 = 1.0-deltmat(i,j);  
        k3 = (1.0-(R(i,j)*R(i,j)*delnmat(i,j)));
        
        
        U2(i,j) = U1(i,j) + (const1(i,j)*(c1*(Txx1(i+1,j)-Txx1(i-1,j)) + c2*(Txx1(i+3,j)-Txx1(i-3,j)) + ...
        c1*(Txz1(i,j+1)-Txz1(i,j-1))+ c2*(Txz1(i,j+3)-Txz1(i,j-3)))) ; 
    
    
        V2(i+1,j+1) = V1(i+1,j+1) + (const1(i,j)*(c1*(Txz1(i+2,j+1)-Txz1(i,j+1))+ c2*(Txz1(i+4,j+1)-Txz1(i-2,j+1))+ ...
        c1*(Tzz1(i+1,j+2)-Tzz1(i+1,j))+ c2*(Tzz1(i+1,j+4)-Tzz1(i+1,j-2)))) ; 
    
    
        Txx2(i+1,j) = Txx1(i+1,j) + ((const2(i,j)*k1*(c1*(U2(i+2,j)-U2(i,j))) + c2*(U2(i+4,j)-U2(i-2,j))) + ...
        (const3(i,j)*k1*(c1*(V2(i+1,j+1)-V2(i+1,j-1))+c2*(V2(i+1,j+3)-V2(i+1,j-3)))));
        
    
        Txz2(i,j+1) = Txz1(i,j+1)+ (const4(i,j)*k2*(c1*(U2(i,j+2)-U2(i,j)) + c2*(U2(i,j+4)-U2(i,j-2)) + ...
        c1*(V2(i+1,j+1)-V2(i-1,j+1))+ c2*(V2(i+3,j+1)-V2(i-3,j+1))));
    
    
        Tzz2(i+1,j) = Tzz1(i+1,j) + (const2(i,j)*k3*(c1*(V2(i+1,j+1)-V2(i+1,j-1))+ c2*(V2(i+1,j+3)-V2(i+1,j-3))) + ...
        const3(i,j)*k1*(c1*(U2(i+2,j)-U2(i,j)) + c2*(U2(i+4,j)-U2(i-2,j))));        
            
        U2(i,j) = G(i,j)*U2(i,j);
        V2(i,j) = G(i,j)*V2(i,j);
    
    
       end
    end
        U1 = U2;       
        V1 = V2;         
        Txx1 = Txx2;     
        Tzz1 = Tzz2;
        Txz1 = Txz2;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
P2 = V2;
% P3 = U2;
% P2(2:2:end,:)=[]; % Remove even ROWS
% P2(:,1:2:end)=[];% Remove odd COLUMNS        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Recoder Position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                 Rec1(end+1) = (P2(125,40));
%                Rec2(end+1) = (P2(105,42)); 
              Rec3(end+1) = (P2(95,46));
%              Rec4(end+1) = (P2(85,50));
%             Rec5(end+1) = (P2(75,56));
%            Rec6(end+1) = (P2(65,64));
%           Rec7(end+1) = (P2(55,76));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 time(end+1) = t1;
% figure(1);imagesc(P2);
%      colorbar;
%     caxis([-0.00000005  0.00000005]);
    %colormap(gray);
%     pause(.01)
end
 dlmwrite('time.out',time);
%  dlmwrite('Rec1.out',Rec1);
%   dlmwrite('Rec2.out',Rec2);
   dlmwrite('Rec3.out',Rec3);
%     dlmwrite('Rec4.out',Rec4);
%      dlmwrite('Rec5.out',Rec5);
%       dlmwrite('Rec6.out',Rec6);
%        dlmwrite('Rec7.out',Rec7);
