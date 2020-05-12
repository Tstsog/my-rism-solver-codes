function [] = rism1d_solver_charged_spc_h2o_model_hnc
%-----------------------------------------------------------------
% Solves the Reference-interaction site model (RISM) equation for charged
% SPC water three-point model in hypernetted chain approximation
% An interparticle interaction: the Lennard-Jones + Coulomb potentials
% A method: a Picard iteration technique
% Uses: function lf=lsint(f)
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Reference: 
% Date: February 27, 2019
%%%
% Temperature T = 298K & density = 0.03346 A^(-3)
%%%
% ----------------------------------------------------------------
clear;
clc;
format short
Nr = 512; % number of grid point; you may change it
L = 32.0; % grid parameter; you may change it  
dr = L/(Nr+1); % grid parameter
itermax = 4000; tol = 10^(-12); % max number of iteration and tolerance; you may change it
%
alf = 0.750; % convergence accelerating parameter; you may change it
%%%
rho = 0.03346;  % number density in  Ang^(-3); (Ang = angstrom) 
%
%%% SPC/E- water model: Coulomb parameter
gam1 = 100.650;  % plus for H-H atom % gam1 = (1/k_B*T)*(1/4*pi*epsilon_0)*qe_H*qe_H/r
gam2 = 4.*gam1;  % minus for O-O atom
gam3 = -2.*gam1; % plus for H-O
%%%
gam11 = gam1; gam12 = gam3; gam13 = gam1;
gam22 = gam2; gam23 = gam3; gam33 = gam1;
%
aa = 1.000000; % parameter in Coulomb potential 
%
%%% SPC/E- water model: LJ parameter
sigma11 = 0.400;  % Ang for H: sigma_1 = 0.4Ang (Ang = angstrom)
sigma22 = 3.166;  % Ang for O: sigma_3 = 3.166Ang 
sigma33 = 0.400;  % Ang for H: sigma_2 = 0.4Ang
sigma12 = 0.5*(sigma11+sigma22);
sigma13 = 0.5*(sigma11+sigma33);
sigma23 = 0.5*(sigma22+sigma33);
%
ell_12 = 1.000; % ell = 1.0Ang between O-H
ell_13 = 1.633; % ell = 1.633Ang between H-H  
ell_23 = 1.000; % ell = 1.0Ang between O-H 
%%%
T_red_1 = 12.87; % eps_1 = 0.0460 kcal/mol & (0.0019872041/0.0460)*298 = 12.87 for H atom & at T = 298K
T_red_2 = 3.811; % eps_3 = 0.1554 kcal/mol & (0.0019872041/0.1554)*298 = 3.811 for O atom & at T = 298K
T_red_3 = 12.87; % eps_2 = 0.0460 kcal/mol & 12.87 for H atom & at T = 298K 
T_red_12 = sqrt(T_red_1*T_red_2);
T_red_13 = sqrt(T_red_1*T_red_3);
T_red_23 = sqrt(T_red_2*T_red_3);
%%%
r = zeros(Nr,1); 
r1 = zeros(Nr,1); r2 = zeros(Nr,1); r3 = zeros(Nr,1); 
r12 = zeros(Nr,1); r13 = zeros(Nr,1); r23 = zeros(Nr,1);
%
u11_lj = zeros(Nr,1); u22_lj = zeros(Nr,1); u33_lj = zeros(Nr,1); 
u12_lj = zeros(Nr,1); u13_lj = zeros(Nr,1); u23_lj = zeros(Nr,1);
%
u11_lr = zeros(Nr,1); u22_lr = zeros(Nr,1); u33_lr = zeros(Nr,1); 
u12_lr = zeros(Nr,1); u13_lr = zeros(Nr,1); u23_lr = zeros(Nr,1);
%
u11_sr = zeros(Nr,1); u22_sr = zeros(Nr,1); u33_sr = zeros(Nr,1); 
u12_sr = zeros(Nr,1); u13_sr = zeros(Nr,1); u23_sr = zeros(Nr,1);
%
c11_sr = zeros(Nr,1); c22_sr = zeros(Nr,1); c33_sr = zeros(Nr,1);
c12_sr = zeros(Nr,1); c13_sr = zeros(Nr,1); c23_sr = zeros(Nr,1);
%
for i = 1:Nr
    r(i) = i*dr;
    r1(i) = r(i)/sigma11;
    r2(i) = r(i)/sigma22;
    r3(i) = r(i)/sigma33;    
    r12(i) = r(i)/sigma12;   
    r13(i) = r(i)/sigma13;
    r23(i) = r(i)/sigma23;       
%
    u11_lj(i) = (4./T_red_1)*( (1./r1(i))^12 - (1./r1(i))^6) ; % LJ potential for H-H
    u11_lr(i) = (gam11/r(i))*erf(aa*r(i));                     % long-range of Coulomb potential for H-H
    u11_sr(i) = (gam11/r(i))*(1.0 - erf(aa*r(i)));             % short-range of Coulomb potential for H-H
    c11_sr(i) = exp( -u11_lj(i) -u11_sr(i) )  - 1. ;    
%
    u22_lj(i) = (4./T_red_2)*( (1./r2(i))^12 - (1./r2(i))^6) ;
    u22_lr(i) = (gam22/r(i))*erf(aa*r(i)); 
    u22_sr(i) = (gam22/r(i))*(1.0 - erf(aa*r(i)));
    c22_sr(i) = exp( -u22_lj(i) -u22_sr(i) )  - 1. ;
%
    u33_lj(i) = (4./T_red_3)*( (1./r3(i))^12 - (1./r3(i))^6) ;
    u33_lr(i) = (gam33/r(i))*erf(aa*r(i)); 
    u33_sr(i) = (gam33/r(i))*(1.0 - erf(aa*r(i)));
    c33_sr(i) = exp( -u33_lj(i) -u33_sr(i) )  - 1. ;
%
    u12_lj(i) = (4./T_red_12)*( (1./r12(i))^12 - (1./r12(i))^6) ;
    u12_lr(i) = (gam12/r(i))*erf(aa*r(i));  
    u12_sr(i) = (gam12/r(i))*(1.0 - erf(aa*r(i)));
    c12_sr(i) = exp( -u12_lj(i) -u12_sr(i) )  - 1. ; 
%
    u13_lj(i) = (4./T_red_13)*( (1./r13(i))^12 - (1./r13(i))^6) ;
    u13_lr(i) = (gam13/r(i))*erf(aa*r(i));  
    u13_sr(i) = (gam13/r(i))*(1.0 - erf(aa*r(i)));
    c13_sr(i) = exp( -u13_lj(i) -u13_sr(i) )  - 1. ; 
%
    u23_lj(i) = (4./T_red_23)*( (1./r23(i))^12 - (1./r23(i))^6) ;
    u23_lr(i) = (gam23/r(i))*erf(aa*r(i));  
    u23_sr(i) = (gam23/r(i))*(1.0 - erf(aa*r(i)));
    c23_sr(i) = exp( -u23_lj(i) -u23_sr(i) )  - 1. ; 
end
%%%
%%%
h11_k_old = 0.;
h12_k_old = 0.;
h13_k_old = 0.;
h22_k_old = 0.;
h23_k_old = 0.;
h33_k_old = 0.;
%
Nk = Nr; j = (1:1:Nk)'; ii = (1:1:Nr)'; dk = pi/L; k = (j*dk);
%
w12_k = sin(k.*ell_12)./(k.*ell_12); % w(k) for H-O
w13_k = sin(k.*ell_13)./(k.*ell_13); % w(k) for H-H
w23_k = sin(k.*ell_23)./(k.*ell_23); % w(k) for O-H
%%%
c11_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c11_sr.*ii)./j;
c12_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c12_sr.*ii)./j;
c13_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c13_sr.*ii)./j;
c22_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c22_sr.*ii)./j;
c23_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c23_sr.*ii)./j;
c33_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c33_sr.*ii)./j;
%
u11_k_lr = (4.*pi*gam11./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u12_k_lr = (4.*pi*gam12./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u13_k_lr = (4.*pi*gam13./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u22_k_lr = (4.*pi*gam22./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u23_k_lr = (4.*pi*gam23./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u33_k_lr = (4.*pi*gam33./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%
c11_k = c11_k_sr - u11_k_lr;
c12_k = c12_k_sr - u12_k_lr;
c13_k = c13_k_sr - u13_k_lr;
c22_k = c22_k_sr - u22_k_lr;
c23_k = c23_k_sr - u23_k_lr;
c33_k = c33_k_sr - u33_k_lr;
%
%--------linear albegra calculation begins -------------------
%%%
d11 = 1. - c12_k.*rho.*w12_k - c13_k.*rho.*w13_k - c11_k.*rho; % 1st element of [I-rho*w_k*w_k] matrix with size 3x3
d12 = -c12_k.*rho - c22_k.*rho.*w12_k - c23_k.*rho.*w13_k;
d13 = -c13_k.*rho - c23_k.*rho.*w12_k - c33_k.*rho.*w13_k;
%
d21 = -c12_k.*rho - c11_k.*rho.*w12_k - c13_k.*rho.*w23_k;
d22 = 1. - c12_k.*rho.*w12_k - c23_k.*rho.*w23_k - c22_k.*rho;
d23 = -c23_k.*rho - c13_k.*rho.*w12_k - c33_k.*rho.*w23_k;
%
d31 = - c13_k.*rho - c11_k.*rho.*w13_k - c12_k.*rho.*w23_k;
d32 = - c23_k.*rho - c12_k.*rho.*w13_k - c22_k.*rho.*w23_k;
d33 = 1. - c13_k.*rho.*w13_k - c23_k.*rho.*w23_k - c33_k.*rho;
%
Det_k = d11.*d22.*d33 - d11.*d23.*d32 - d12.*d21.*d33 + d12.*d23.*d31 + d13.*d21.*d32 - d13.*d22.*d31; % det[I-rho*w_k*w_k]
%
D11_k = (d22.*d33 - d23.*d32)./Det_k; D12_k = -(d12.*d33 - d13.*d32)./Det_k; D13_k = (d12.*d23 - d13.*d22)./Det_k; % [I-rho*w_k*w_k]^(-1)
D21_k = -(d21.*d33 - d23.*d31)./Det_k; D22_k = (d11.*d33 - d13.*d31)./Det_k; D23_k = -(d11.*d23 - d13.*d21)./Det_k;
D31_k = (d21.*d32 - d22.*d31)./Det_k; D32_k = -(d11.*d32 - d12.*d31)./Det_k; D33_k = (d11.*d22 - d12.*d21)./Det_k;
%%%
h11_k = w12_k.*(c12_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c22_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + ...
        c23_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + w13_k.*(c13_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + ...
        c23_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c33_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + ...
        c11_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c12_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c13_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k); 
h11_k_new = h11_k;
h11_k_n = alf.*h11_k_new + (1.-alf).*h11_k_old;
%
h12_k = w12_k.*(c11_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c12_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + ...
        c13_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + w23_k.*(c13_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + ...
        c23_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c33_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + ...
        c12_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c22_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c23_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k);
h12_k_new = h12_k;
h12_k_n = alf.*h12_k_new + (1.-alf).*h12_k_old;
%
h13_k = w13_k.*(c11_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c12_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + ...
        c13_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + w23_k.*(c12_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + ...
        c22_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c23_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + ...
        c13_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c23_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c33_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k);
h13_k_new = h13_k;
h13_k_n = alf.*h13_k_new + (1.-alf).*h13_k_old;    
%   
h22_k = w12_k.*(c11_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c12_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + ...
        c13_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + w23_k.*(c13_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + ...
        c23_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c33_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + ...
        c12_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c22_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c23_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k);
h22_k_new = h22_k;
h22_k_n = alf.*h22_k_new + (1.-alf).*h22_k_old;    
%
h23_k = w13_k.*(c11_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c12_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + ...
        c13_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + w23_k.*(c12_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + ...
        c22_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c23_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + ...
        c13_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c23_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c33_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k); 
h23_k_new = h23_k;
h23_k_n = alf.*h23_k_new + (1.-alf).*h23_k_old;    
%    
h33_k = w13_k.*(c11_k.*(D31_k + D32_k.*w12_k + D33_k.*w13_k) + c12_k.*(D32_k + D31_k.*w12_k + D33_k.*w23_k) + ...
        c13_k.*(D33_k + D31_k.*w13_k + D32_k.*w23_k)) + w23_k.*(c12_k.*(D31_k + D32_k.*w12_k + D33_k.*w13_k) + ...
        c22_k.*(D32_k + D31_k.*w12_k + D33_k.*w23_k) + c23_k.*(D33_k + D31_k.*w13_k + D32_k.*w23_k)) + ...
        c13_k.*(D31_k + D32_k.*w12_k + D33_k.*w13_k) + c23_k.*(D32_k + D31_k.*w12_k + D33_k.*w23_k) + c33_k.*(D33_k + D31_k.*w13_k + D32_k.*w23_k);
h33_k_new = h33_k;
h33_k_n = alf.*h33_k_new + (1.-alf).*h33_k_old;
%%%
%
%--------linear albegra calculation ends -------------------
%%%
% compute c(r) and h(r) 
h11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h11_k_n.*j);
h12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h12_k_n.*j);
h13_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h13_k_n.*j);
h22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h22_k_n.*j);
h23_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h23_k_n.*j);
h33_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h33_k_n.*j);
%
c11_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c11_k_sr.*j);
c12_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c12_k_sr.*j);
c13_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c13_k_sr.*j);
c22_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c22_k_sr.*j);
c23_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c23_k_sr.*j);
c33_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c33_k_sr.*j);
%%%
G11_sr = h11_r - c11_sr; % indirect correlation functions
G12_sr = h12_r - c12_sr;
G13_sr = h13_r - c13_sr;
G22_sr = h22_r - c22_sr;
G23_sr = h23_r - c23_sr;
G33_sr = h33_r - c33_sr;
%
%Iflag = 0.;  % converfence flag, 1 indicates convergence
%
for iter = 2:itermax
%
    iter
% iteration 2
    G11_sr_old = G11_sr;
    G12_sr_old = G12_sr;
    G13_sr_old = G13_sr;    
    G22_sr_old = G22_sr;
    G23_sr_old = G23_sr;    
    G33_sr_old = G33_sr;    
%%%
    c11_sr_old = exp(-u11_lj -u11_sr + G11_sr_old ) - G11_sr_old - 1. ;    
    c12_sr_old = exp(-u12_lj -u12_sr + G12_sr_old ) - G12_sr_old - 1. ;    
    c13_sr_old = exp(-u13_lj -u13_sr + G13_sr_old ) - G13_sr_old - 1. ;        
    c22_sr_old = exp(-u22_lj -u22_sr + G22_sr_old ) - G22_sr_old - 1. ;  
    c23_sr_old = exp(-u23_lj -u23_sr + G23_sr_old ) - G23_sr_old - 1. ;  
    c33_sr_old = exp(-u33_lj -u33_sr + G33_sr_old ) - G33_sr_old - 1. ;      
%%%
    c11_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c11_sr_old.*ii)./j;
    c12_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c12_sr_old.*ii)./j;
    c13_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c13_sr_old.*ii)./j;
    c22_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c22_sr_old.*ii)./j;
    c23_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c23_sr_old.*ii)./j;
    c33_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c33_sr_old.*ii)./j;
%%%    
    c11_k = c11_k_sr - u11_k_lr;
    c12_k = c12_k_sr - u12_k_lr;
    c13_k = c13_k_sr - u13_k_lr;    
    c22_k = c22_k_sr - u22_k_lr;
    c23_k = c23_k_sr - u23_k_lr;    
    c33_k = c33_k_sr - u33_k_lr;
%%%
%
%--------linear albegra calculation begins -------------------
%
    d11 = 1. - c12_k.*rho.*w12_k - c13_k.*rho.*w13_k - c11_k.*rho; % 1st element of [I-rho*w_k*w_k] matrix
    d12 = -c12_k.*rho - c22_k.*rho.*w12_k - c23_k.*rho.*w13_k;
    d13 = -c13_k.*rho - c23_k.*rho.*w12_k - c33_k.*rho.*w13_k;
%
    d21 = -c12_k.*rho - c11_k.*rho.*w12_k - c13_k.*rho.*w23_k;
    d22 = 1. - c12_k.*rho.*w12_k - c23_k.*rho.*w23_k - c22_k.*rho;
    d23 = -c23_k.*rho - c13_k.*rho.*w12_k - c33_k.*rho.*w23_k;
%
    d31 = - c13_k.*rho - c11_k.*rho.*w13_k - c12_k.*rho.*w23_k;
    d32 = - c23_k.*rho - c12_k.*rho.*w13_k - c22_k.*rho.*w23_k;
    d33 = 1. - c13_k.*rho.*w13_k - c23_k.*rho.*w23_k - c33_k.*rho;
%
    Det_k = d11.*d22.*d33 - d11.*d23.*d32 - d12.*d21.*d33 + d12.*d23.*d31 + d13.*d21.*d32 - d13.*d22.*d31; % det[I-rho*w_k*w_k]
%
    D11_k = (d22.*d33 - d23.*d32)./Det_k; D12_k = -(d12.*d33 - d13.*d32)./Det_k; D13_k = (d12.*d23 - d13.*d22)./Det_k; % [I-rho*w_k*w_k]^(-1)
    D21_k = -(d21.*d33 - d23.*d31)./Det_k; D22_k = (d11.*d33 - d13.*d31)./Det_k; D23_k = -(d11.*d23 - d13.*d21)./Det_k;
    D31_k = (d21.*d32 - d22.*d31)./Det_k; D32_k = -(d11.*d32 - d12.*d31)./Det_k; D33_k = (d11.*d22 - d12.*d21)./Det_k;
%%%
    h11_k = w12_k.*(c12_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c22_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + ...
            c23_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + w13_k.*(c13_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + ...
            c23_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c33_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + ...
            c11_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c12_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c13_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k); 
    h11_k_old = h11_k_n;
    h11_k_new = h11_k;
    h11_k_n = alf.*h11_k_new + (1.-alf).*h11_k_old;   
%
    h12_k = w12_k.*(c11_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c12_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + ...
            c13_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + w23_k.*(c13_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + ...
            c23_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c33_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + ...
            c12_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c22_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c23_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k);
    h12_k_old = h12_k_n;        
    h12_k_new = h12_k;
    h12_k_n = alf.*h12_k_new + (1.-alf).*h12_k_old;
%
    h13_k = w13_k.*(c11_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c12_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + ...
            c13_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + w23_k.*(c12_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + ...
            c22_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c23_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k)) + ...
            c13_k.*(D11_k + D12_k.*w12_k + D13_k.*w13_k) + c23_k.*(D12_k + D11_k.*w12_k + D13_k.*w23_k) + c33_k.*(D13_k + D11_k.*w13_k + D12_k.*w23_k);
    h13_k_old = h13_k_n;        
    h13_k_new = h13_k;
    h13_k_n = alf.*h13_k_new + (1.-alf).*h13_k_old;    
%
    h22_k = w12_k.*(c11_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c12_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + ...
            c13_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + w23_k.*(c13_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + ...
            c23_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c33_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + ...
           c12_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c22_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c23_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k);
    h22_k_old = h22_k_n;       
    h22_k_new = h22_k;
    h22_k_n = alf.*h22_k_new + (1.-alf).*h22_k_old;    
%
    h23_k = w13_k.*(c11_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c12_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + ...
             c13_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + w23_k.*(c12_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + ...
             c22_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c23_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k)) + ...
             c13_k.*(D21_k + D22_k.*w12_k + D23_k.*w13_k) + c23_k.*(D22_k + D21_k.*w12_k + D23_k.*w23_k) + c33_k.*(D23_k + D21_k.*w13_k + D22_k.*w23_k); 
    h23_k_old = h23_k_n;         
    h23_k_new = h23_k;
    h23_k_n = alf.*h23_k_new + (1.-alf).*h23_k_old;    
%    
    h33_k = w13_k.*(c11_k.*(D31_k + D32_k.*w12_k + D33_k.*w13_k) + c12_k.*(D32_k + D31_k.*w12_k + D33_k.*w23_k) + ...
            c13_k.*(D33_k + D31_k.*w13_k + D32_k.*w23_k)) + w23_k.*(c12_k.*(D31_k + D32_k.*w12_k + D33_k.*w13_k) + ...
             c22_k.*(D32_k + D31_k.*w12_k + D33_k.*w23_k) + c23_k.*(D33_k + D31_k.*w13_k + D32_k.*w23_k)) + ...
             c13_k.*(D31_k + D32_k.*w12_k + D33_k.*w13_k) + c23_k.*(D32_k + D31_k.*w12_k + D33_k.*w23_k) + c33_k.*(D33_k + D31_k.*w13_k + D32_k.*w23_k);
    h33_k_old = h33_k_n;         
    h33_k_new = h33_k;
    h33_k_n = alf.*h33_k_new + (1.-alf).*h33_k_old;
%%%
%--------linear albegra calculation ends -------------------
%%%
% compute c(r) and h(r) 
    h11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h11_k_n.*j);
    h12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h12_k_n.*j);
    h13_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h13_k_n.*j);
    h22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h22_k_n.*j);
    h23_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h23_k_n.*j);
    h33_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h33_k_n.*j);
%
    c11_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c11_k_sr.*j);
    c12_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c12_k_sr.*j);
    c13_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c13_k_sr.*j);
    c22_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c22_k_sr.*j);
    c23_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c23_k_sr.*j);
    c33_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c33_k_sr.*j);
%%%
    G11_sr = h11_r - c11_sr;
    G12_sr = h12_r - c12_sr;
    G13_sr = h13_r - c13_sr;
    G22_sr = h22_r - c22_sr;
    G23_sr = h23_r - c23_sr;
    G33_sr = h33_r - c33_sr;
 %%%
    G11_sr_new = G11_sr; % indirect cf
    G12_sr_new = G12_sr;
    G13_sr_new = G13_sr;    
    G22_sr_new = G22_sr;
    G23_sr_new = G23_sr;    
    G33_sr_new = G33_sr;
%%%
    cor = rms(G11_sr_new - G11_sr_old) + rms(G12_sr_new - G12_sr_old) + ...
          rms(G13_sr_new - G13_sr_old) + rms(G22_sr_new - G22_sr_old) + ...
          rms(G23_sr_new - G23_sr_old) + rms(G33_sr_new - G33_sr_old); 
%%%
    if (abs(cor) < tol);
          Iflag = 1.
          break
    end
%
    G11_sr = G11_sr_new;
    G12_sr = G12_sr_new;
    G13_sr = G13_sr_new;    
    G22_sr = G22_sr_new;   
    G23_sr = G23_sr_new;  
    G33_sr = G33_sr_new;       
%  
end
%%%
figure(1) % plots for radial distribution functions for each site-site
hold on
plot(r, h11_r+1, 'b')  % for H-H
plot(r, h12_r+1, 'k')  % for H-O
plot(r, h22_r+1, 'g')  % for O-O 
hold off
axis([0. 16. 0.000 3.00])
set(gca,'FontSize',20)
xlabel('r(Ang)') % ,'fontsize',16
ylabel('g_{ij}(r)','Rotation', 1) %

return
end

%%%
% A function lsint.m is taken from C.T. Kelley, B. M. Pettitt, 
%                      J. Chmp. Phys. \textbf{197}, 491 (2004)
% LSINT
% Fast sine transform with MATLAB's FFT.
%
function lf=lsint(f)
n=length(f);
ft=-fft([0,f']',2*n+2);
lf=imag(ft(2:n+1));
%
return
end

