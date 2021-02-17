%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves excess chemical potential for solvent and
% solvation free energy for solute particle for solute-solvent system 
% in the infinite dilution limit using the Reference interaction site model (RISM) theory
% in the hyppernetted-chain (HNC) approximation.
%%%
%%% --- Reference interaction site model (RISM) theory
% The RISM coupled equations are solved in infinite dilution limit (solute density is zero). 
% 1. solve solvent-solvent equation: 
% solvent_solvent_rism1d_solver.m function solves the solvent-solvent RISM equation.
% 2. solve solute-solvent equation:
% solute_solvent_rism_1d_solver. function solves the solute-solvent RISM equation.
%
% Here:
% solvent is SPC/E water model.
% solute is K+ ion with partial charge +1.
% temparature is T = 298K.
%%%
%
% The Picard iterative method is applied. For solute-solvent case, we first
% solve solute-solvent RISM equation for small gamma parameter (Gamma=40),
% and output from this previous calculation is used for Gamma = 237.50 in next step.
%%%
%
% The excess chemical potential in HNC approximation is
% beta*mu = rho_solvent*int^infinity_0*(0.5*h_vv(r)^{2} - c_vv(r) - 0.5*h_vv(r)*c_vv(r)) * (4*pi*r^2)*dr, 
%                                       where h_vv(r) & c_vv(r) are from solvent-solvent calculation.
%
% For solute-solvent case, in evaluation of beta*mu, one uses h_uv(r) and c_uv from solute-solvent calculation. 
% 
%%%
%%% ---
%
% References: D. Chandler and H. C. Andersen, J. Chem. Phys. v57, p1930 (1972)
%             D. Chandler, J. Chem. Phys. v59, p2742 (1973)
%             L. J. Lowden and D. Chandler, J. Chem. Phys. v59, p6587 (1973)
%             I. S. Joung, T. Luchko, and D. A. Case, J. Chem. Phys. v138, p044103 (2013)
%
% Written by Tsogbayar Tsednee (PhD), Kyungpook National University, South Korea 
% Contact: tsog215@gmail.com
% Date: Feb 16, 2020
%
%
function [] = rism_1D_solute_solvent_equation_hnc_example_1
%
clear all;
format long
clear;
clc;
Nr = 4*256;                      %  number of grid point
L = 32.0; dr = L/(Nr+1);         %  length parameter & step size
itermax = 1000; tol = 10^(-6);   %  number of iteration & tolerance
%
alf = 0.750;                     %  parameter for iteration convergence for solvent-solvent RISM calculation
alf2 = 0.900;                    %  parameter for iteration convergence for solute-solvent RISM calculation
%%%
%
[x_mu_h2o,chi11,chi22,chi33,chi12,chi13,chi23,h11_r,h12_r,h22_r,r] = solvent_solvent_rism1d_solver(Nr,L,dr,itermax,tol,alf);
%
%%% excess chemical potential beta*mu
mu_kcal_over_mol_h2o = x_mu_h2o * 0.5918;                         % in unit of kcal/mol
%
[mu_kap_ion,h14_r,h24_r,h34_r] = solute_solvent_rism_1d_solver(Nr,L,dr,itermax,chi11,chi22,chi33,chi12,chi13,chi23,tol,alf2);
%
%%% solvation free energy (excess chemical potential) beta*mu for solute
%%% particle
mu_kcal_over_mol = mu_kap_ion * 0.5918; % in unit of kcal/mol
%
%%% excess chemical potentias for solvent-solvent and solute-solvent interactions:
[mu_kcal_over_mol_h2o, mu_kcal_over_mol]               % in unit of kcal/mol
%%%%
%
% radial distribution functions for solvent
figure(1)
hold on
plot(r, h11_r+1, 'g')    % for H-H
plot(r, h12_r+1, 'k')    % for H-O
plot(r, h22_r+1, 'b')    % for O-O 
%
% radial distribution functions for solute-solvent
plot(r, h14_r+1, 'g--')    % for H - K+
plot(r, h24_r+1, 'b--')    % for O - K+
%plot(r, h34_r+1, 'ro')     % for H - K+
%
hold off
axis([0. 12. 0.000 4.00])
set(gca,'FontSize',20)
xlabel('r(Ang)') % ,'fontsize',16
ylabel('g_{ij}(r)','Rotation', 1) %   

%
return
end

function [x_mu_h2o,chi11,chi22,chi33,chi12,chi13,chi23,h11_r,h12_r,h22_r,r] = solvent_solvent_rism1d_solver(Nr,L,dr,itermax,tol,alf)
%%%
%%%
rho_red = 0.033460;  % number density, A^(-3)
gam1 = 100.650;  % plus for H-H atom % gam1
gam2 = 4.*gam1;  % minus for O-O atom
gam3 = -2.*gam1; % plus for H-O
%%%
gam11 = gam1;
gam12 = gam3;
gam13 = gam1;
gam22 = gam2;
gam23 = gam3;
gam33 = gam1;
%
aa = 1.000000;
%
%%% SPC/E- water model
sigma11 = 0.400; % Ang for H: sigma_1 = 0.4Ang
sigma22 = 3.166; % Ang for O: sigma_3 = 3.166Ang 
sigma33 = 0.400; % Ang for H: sigma_2 = 0.4Ang
sigma12 = 0.5*(sigma11+sigma22);
sigma13 = 0.5*(sigma11+sigma33);
sigma23 = 0.5*(sigma22+sigma33);
%
ell_12 = 1.000; % ell = 1.0Ang between O-H
ell_13 = 1.633; % ell = 1.633Ang between H-H  
ell_23 = 1.000; % ell = 1.0Ang between O-H 
%%%
%rho_red = 0.10000 - 0.00000;
rho_red_1 = rho_red; % A^(-3) x_1 = 0.5;
%rho_red_2 = rho_red_1; % A^(-3) x_2 = 0.5;
%rho_red_3 = rho_red_1; % A^(-3) x_2 = 0.5;
%%%
T_red_1 = 12.870; % eps_1 = 0.0460 kcal/mol & (0.0019872041/0.0460)*298 = 12.87 for H atom & at T = 298K
T_red_2 = 3.811; % eps_3 = 0.1554 kcal/mol & (0.0019872041/0.1554)*298 = 3.811 for O atom & at T = 298K
T_red_3 = 12.87; % eps_2 = 0.0460 kcal/mol & 12.87 for H atom & at T = 298K 
T_red_12 = sqrt(T_red_1*T_red_2);
T_red_13 = sqrt(T_red_1*T_red_3);
T_red_23 = sqrt(T_red_2*T_red_3);
%%%
%%%
%
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
    u11_lj(i) = (4./T_red_1)*( (1./r1(i))^12 - (1./r1(i))^6) ;
    u11_lr(i) = (gam11/r(i))*erf(aa*r(i));  
    u11_sr(i) = (gam11/r(i))*(1.0 - erf(aa*r(i)));
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
h11_k_old = 0.; h12_k_old = 0.; h13_k_old = 0.;
h22_k_old = 0.; h23_k_old = 0.;
h33_k_old = 0.;
%
Nk = Nr;
j = (1:1:Nk)';
ii = (1:1:Nr)';
dk = pi/L;
k = (j*dk);
%
w12_k = sin(k.*ell_12)./(k.*ell_12); % w(k) for H-O
w13_k = sin(k.*ell_13)./(k.*ell_13); % w(k) for H-H
w23_k = sin(k.*ell_23)./(k.*ell_23); % w(k) for O-H
%%%
c11_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c11_sr.*ii)./j;
c12_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c12_sr.*ii)./j;
c13_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c13_sr.*ii)./j;
%c21_k_sr = c12_k_sr;
c22_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c22_sr.*ii)./j;
c23_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c23_sr.*ii)./j;
%c31_k_sr = c13_k_sr;
%c32_k_sr = c23_k_sr;
c33_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c33_sr.*ii)./j;
%
u11_k_lr = (4.*pi*gam11./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u12_k_lr = (4.*pi*gam12./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u13_k_lr = (4.*pi*gam13./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u22_k_lr = (4.*pi*gam22./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u23_k_lr = (4.*pi*gam23./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u33_k_lr = (4.*pi*gam33./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%%%
%ck = ck_sr - uk_lr;
c11_k = c11_k_sr - u11_k_lr;
c12_k = c12_k_sr - u12_k_lr;
c13_k = c13_k_sr - u13_k_lr;
%c21_k = c12_k; 
c22_k = c22_k_sr - u22_k_lr;
c23_k = c23_k_sr - u23_k_lr;
%c31_k = c13_k;
%c32_k = c23_k;
c33_k = c33_k_sr - u33_k_lr;
%%%%%%%%%%%% 
%%%
%rho = rho_red_1;
%%%
%%% solve rism 3x3 matrix equation
% w*c*w
wd11_k = c11_k + w12_k.*(c12_k + c22_k.*w12_k + c23_k.*w13_k) + w13_k.*(c13_k + c23_k.*w12_k + c33_k.*w13_k) + c12_k.*w12_k + c13_k.*w13_k;
wd12_k = c12_k + w12_k.*(c11_k + c12_k.*w12_k + c13_k.*w13_k) + w23_k.*(c13_k + c23_k.*w12_k + c33_k.*w13_k) + c22_k.*w12_k + c23_k.*w13_k;
wd13_k = c13_k + w13_k.*(c11_k + c12_k.*w12_k + c13_k.*w13_k) + w23_k.*(c12_k + c22_k.*w12_k + c23_k.*w13_k) + c23_k.*w12_k + c33_k.*w13_k;
%
wd22_k = c22_k + w12_k.*(c12_k + c11_k.*w12_k + c13_k.*w23_k) + w23_k.*(c23_k + c13_k.*w12_k + c33_k.*w23_k) + c12_k.*w12_k + c23_k.*w23_k;
wd23_k = c23_k + w13_k.*(c12_k + c11_k.*w12_k + c13_k.*w23_k) + w23_k.*(c22_k + c12_k.*w12_k + c23_k.*w23_k) + c13_k.*w12_k + c33_k.*w23_k;
%
wd33_k = c33_k + w13_k.*(c13_k + c11_k.*w13_k + c12_k.*w23_k) + w23_k.*(c23_k + c12_k.*w13_k + c22_k.*w23_k) + c13_k.*w13_k + c23_k.*w23_k;
%%%
%%% elements of (I-W*C*R)
d11_k = 1. - rho_red_1.*(c11_k + c12_k.*w12_k + c13_k.*w13_k);
d12_k = -rho_red_1.*(c12_k + c22_k.*w12_k + c23_k.*w13_k);
d13_k = -rho_red_1.*(c13_k + c23_k.*w12_k + c33_k.*w13_k);
%
d21_k = -rho_red_1.*(c12_k + c11_k.*w12_k + c13_k.*w23_k);
d22_k = 1. - rho_red_1.*(c22_k + c12_k.*w12_k + c23_k.*w23_k);
d23_k = -rho_red_1.*(c23_k + c13_k.*w12_k + c33_k.*w23_k);
%
d31_k = -rho_red_1.*(c13_k + c11_k.*w13_k + c12_k.*w23_k);
d32_k = -rho_red_1.*(c23_k + c12_k.*w13_k + c22_k.*w23_k);
d33_k = 1. - rho_red_1.*(c33_k + c13_k.*w13_k + c23_k.*w23_k);
%
%
Det_3x3_mat = (d11_k.*d22_k.*d33_k - d11_k.*d23_k.*d32_k - d12_k.*d21_k.*d33_k + ...
               d12_k.*d23_k.*d31_k + d13_k.*d21_k.*d32_k - d13_k.*d22_k.*d31_k);
           
%
id11_k = (d22_k.*d33_k - d23_k.*d32_k)./Det_3x3_mat;
id12_k = -(d12_k.*d33_k - d13_k.*d32_k)./Det_3x3_mat;
id13_k = (d12_k.*d23_k - d13_k.*d22_k)./Det_3x3_mat;
%
id21_k = -(d21_k.*d33_k - d23_k.*d31_k)./Det_3x3_mat;
id22_k = (d11_k.*d33_k - d13_k.*d31_k)./Det_3x3_mat;
id23_k = -(d11_k.*d23_k - d13_k.*d21_k)./Det_3x3_mat;
%
id31_k = (d21_k.*d32_k - d22_k.*d31_k)./Det_3x3_mat;
id32_k = -(d11_k.*d32_k - d12_k.*d31_k)./Det_3x3_mat;
id33_k = (d11_k.*d22_k - d12_k.*d21_k)./Det_3x3_mat;
%
%%%
h11_k = id11_k.*wd11_k + id12_k.*wd12_k + id13_k.*wd13_k;
h11_k_new = h11_k;
h11_k_n = alf.*h11_k_new + (1.-alf).*h11_k_old; 
%
h12_k = id11_k.*wd12_k + id12_k.*wd22_k + id13_k.*wd23_k;
h12_k_new = h12_k;
h12_k_n = alf.*h12_k_new + (1.-alf).*h12_k_old; 
%
h13_k = id11_k.*wd13_k + id12_k.*wd23_k + id13_k.*wd33_k;
h13_k_new = h13_k;
h13_k_n = alf.*h13_k_new + (1.-alf).*h13_k_old; 
%
h22_k = id21_k.*wd12_k + id22_k.*wd22_k + id23_k.*wd23_k;
h22_k_new = h22_k;
h22_k_n = alf.*h22_k_new + (1.-alf).*h22_k_old; 
%
h23_k = id21_k.*wd13_k + id22_k.*wd23_k + id23_k.*wd33_k;
h23_k_new = h23_k;
h23_k_n = alf.*h23_k_new + (1.-alf).*h23_k_old; 
%
h33_k = id31_k.*wd13_k + id32_k.*wd23_k + id33_k.*wd33_k;
h33_k_new = h33_k;
h33_k_n = alf.*h33_k_new + (1.-alf).*h33_k_old; 
%%%
% compute c(r) and h(r) 
c11_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c11_k_sr.*j);
c12_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c12_k_sr.*j);
c13_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c13_k_sr.*j);
%c21_r = c12_r;
c22_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c22_k_sr.*j);
c23_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c23_k_sr.*j);
%
c33_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c33_k_sr.*j);
%
h11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h11_k_n.*j);
h12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h12_k_n.*j);
h13_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h13_k_n.*j);
%h21_r = h12_r;
h22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h22_k_n.*j);
h23_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h23_k_n.*j);
%
h33_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h33_k_n.*j);
%
G11_sr = h11_r - c11_sr; % indirect cf
G12_sr = h12_r - c12_sr;
G13_sr = h13_r - c13_sr;
%G21_r = G12_r;
G22_sr = h22_r - c22_sr;
G23_sr = h23_r - c23_sr;
%
G33_sr = h33_r - c33_sr;
%%%
%Iflag = 0.;  % converfence flag, 1 indicates convergence

%
for iter = 2:itermax
%
%   iter;
% iteration 2
    G11_sr_old = G11_sr; G12_sr_old = G12_sr; G13_sr_old = G13_sr; 
    G22_sr_old = G22_sr; G23_sr_old = G23_sr;
    G33_sr_old = G33_sr; 
%
    c11_sr_old = exp(-u11_lj -u11_sr + G11_sr_old) - G11_sr_old - 1.;
    c12_sr_old = exp(-u12_lj -u12_sr + G12_sr_old) - G12_sr_old - 1.;
    c13_sr_old = exp(-u13_lj -u13_sr + G13_sr_old) - G13_sr_old - 1.;
%
    c22_sr_old = exp(-u22_lj -u22_sr + G22_sr_old) - G22_sr_old - 1.;
    c23_sr_old = exp(-u23_lj -u23_sr + G23_sr_old) - G23_sr_old - 1.;
%
    c33_sr_old = exp(-u33_lj -u33_sr + G33_sr_old) - G33_sr_old - 1.;
%%%
    c11_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c11_sr_old.*ii)./j;
    c12_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c12_sr_old.*ii)./j;
    c13_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c13_sr_old.*ii)./j;
    c22_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c22_sr_old.*ii)./j;
    c23_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c23_sr_old.*ii)./j;
    c33_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c33_sr_old.*ii)./j;
%
    %ck = ck_sr - uk_lr;
    c11_k = c11_k_sr - u11_k_lr;
    c12_k = c12_k_sr - u12_k_lr;
    c13_k = c13_k_sr - u13_k_lr;
    %
    c22_k = c22_k_sr - u22_k_lr;
    c23_k = c23_k_sr - u23_k_lr;
    %
    c33_k = c33_k_sr - u33_k_lr;
    %
    % w*c*w
    wd11_k = c11_k + w12_k.*(c12_k + c22_k.*w12_k + c23_k.*w13_k) + w13_k.*(c13_k + c23_k.*w12_k + c33_k.*w13_k) + c12_k.*w12_k + c13_k.*w13_k;
    wd12_k = c12_k + w12_k.*(c11_k + c12_k.*w12_k + c13_k.*w13_k) + w23_k.*(c13_k + c23_k.*w12_k + c33_k.*w13_k) + c22_k.*w12_k + c23_k.*w13_k;
    wd13_k = c13_k + w13_k.*(c11_k + c12_k.*w12_k + c13_k.*w13_k) + w23_k.*(c12_k + c22_k.*w12_k + c23_k.*w13_k) + c23_k.*w12_k + c33_k.*w13_k;
    %
    wd22_k = c22_k + w12_k.*(c12_k + c11_k.*w12_k + c13_k.*w23_k) + w23_k.*(c23_k + c13_k.*w12_k + c33_k.*w23_k) + c12_k.*w12_k + c23_k.*w23_k;
    wd23_k = c23_k + w13_k.*(c12_k + c11_k.*w12_k + c13_k.*w23_k) + w23_k.*(c22_k + c12_k.*w12_k + c23_k.*w23_k) + c13_k.*w12_k + c33_k.*w23_k;
    %
    wd33_k = c33_k + w13_k.*(c13_k + c11_k.*w13_k + c12_k.*w23_k) + w23_k.*(c23_k + c12_k.*w13_k + c22_k.*w23_k) + c13_k.*w13_k + c23_k.*w23_k;
    %%%
    %%% elements of (I-W*C*R)
    d11_k = 1. - rho_red_1.*(c11_k + c12_k.*w12_k + c13_k.*w13_k);
    d12_k = -rho_red_1.*(c12_k + c22_k.*w12_k + c23_k.*w13_k);
    d13_k = -rho_red_1.*(c13_k + c23_k.*w12_k + c33_k.*w13_k);
    %
    d21_k = -rho_red_1.*(c12_k + c11_k.*w12_k + c13_k.*w23_k);
    d22_k = 1. - rho_red_1.*(c22_k + c12_k.*w12_k + c23_k.*w23_k);
    d23_k = -rho_red_1.*(c23_k + c13_k.*w12_k + c33_k.*w23_k);
    %
    d31_k = -rho_red_1.*(c13_k + c11_k.*w13_k + c12_k.*w23_k);
    d32_k = -rho_red_1.*(c23_k + c12_k.*w13_k + c22_k.*w23_k);
    d33_k = 1. - rho_red_1.*(c33_k + c13_k.*w13_k + c23_k.*w23_k);
    %
    %
    Det_3x3_mat = (d11_k.*d22_k.*d33_k - d11_k.*d23_k.*d32_k - d12_k.*d21_k.*d33_k + ...
                   d12_k.*d23_k.*d31_k + d13_k.*d21_k.*d32_k - d13_k.*d22_k.*d31_k);

    %
    id11_k = (d22_k.*d33_k - d23_k.*d32_k)./Det_3x3_mat;
    id12_k = -(d12_k.*d33_k - d13_k.*d32_k)./Det_3x3_mat;
    id13_k = (d12_k.*d23_k - d13_k.*d22_k)./Det_3x3_mat;
    %
    id21_k = -(d21_k.*d33_k - d23_k.*d31_k)./Det_3x3_mat;
    id22_k = (d11_k.*d33_k - d13_k.*d31_k)./Det_3x3_mat;
    id23_k = -(d11_k.*d23_k - d13_k.*d21_k)./Det_3x3_mat;
    %
    id31_k = (d21_k.*d32_k - d22_k.*d31_k)./Det_3x3_mat;
    id32_k = -(d11_k.*d32_k - d12_k.*d31_k)./Det_3x3_mat;
    id33_k = (d11_k.*d22_k - d12_k.*d21_k)./Det_3x3_mat;
    %
    %%%
    h11_k = id11_k.*wd11_k + id12_k.*wd12_k + id13_k.*wd13_k;
    h11_k_old = h11_k_n;
    h11_k_new = h11_k;
    h11_k_n = alf.*h11_k_new + (1.-alf).*h11_k_old; 
    %
    h12_k = id11_k.*wd12_k + id12_k.*wd22_k + id13_k.*wd23_k;
    h12_k_old = h12_k_n;
    h12_k_new = h12_k;
    h12_k_n = alf.*h12_k_new + (1.-alf).*h12_k_old; 
    %
    h13_k = id11_k.*wd13_k + id12_k.*wd23_k + id13_k.*wd33_k;
    h13_k_old = h13_k_n;
    h13_k_new = h13_k;
    h13_k_n = alf.*h13_k_new + (1.-alf).*h13_k_old; 
    %
    h22_k = id21_k.*wd12_k + id22_k.*wd22_k + id23_k.*wd23_k;
    h22_k_old = h22_k_n;
    h22_k_new = h22_k;
    h22_k_n = alf.*h22_k_new + (1.-alf).*h22_k_old; 
    %
    h23_k = id21_k.*wd13_k + id22_k.*wd23_k + id23_k.*wd33_k;
    h23_k_old = h23_k_n;
    h23_k_new = h23_k;
    h23_k_n = alf.*h23_k_new + (1.-alf).*h23_k_old; 
    %
    h33_k = id31_k.*wd13_k + id32_k.*wd23_k + id33_k.*wd33_k;
    h33_k_old = h33_k_n;
    h33_k_new = h33_k;
    h33_k_n = alf.*h33_k_new + (1.-alf).*h33_k_old; 
    %%%
    % compute c(r) and h(r) 
    c11_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c11_k_sr.*j);
    c12_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c12_k_sr.*j);
    c13_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c13_k_sr.*j);
    %c21_r = c12_r;
    c22_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c22_k_sr.*j);
    c23_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c23_k_sr.*j);
    %
    c33_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c33_k_sr.*j);
    %
    h11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h11_k_n.*j);
    h12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h12_k_n.*j);
    h13_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h13_k_n.*j);
    %h21_r = h12_r;
    h22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h22_k_n.*j);
    h23_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h23_k_n.*j);
    %
    h33_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h33_k_n.*j);
    %
    G11_sr_new = h11_r - c11_sr; % indirect cf
    G12_sr_new = h12_r - c12_sr;
    G13_sr_new = h13_r - c13_sr;
    %G21_r = G12_r;
    G22_sr_new = h22_r - c22_sr;
    G23_sr_new = h23_r - c23_sr;
    %
    G33_sr_new = h33_r - c33_sr;
    
    %%%
    %
    cor = 1.*rms(G11_sr_new - G11_sr_old) + 1.*rms(G12_sr_new - G12_sr_old) + 1.*rms(G13_sr_new - G13_sr_old) + ...
          1.*rms(G22_sr_new - G22_sr_old) + 1.*rms(G23_sr_new - G23_sr_old) + ...
          1.*rms(G33_sr_new - G33_sr_old);
%
    G11_sr = G11_sr_new;
    G12_sr = G12_sr_new;
    G13_sr = G13_sr_new;   
    G22_sr = G22_sr_new; 
    G23_sr = G23_sr_new;    
    G33_sr = G33_sr_new; 
%
    if (cor < tol)        
        Iflag = 1.
    break
    end
end
%%%
c11_r = c11_sr; c12_r = c12_sr; c13_r = c13_sr;
c22_r = c22_sr; c23_r = c23_sr;
c33_r = c33_sr;
%
c11_r_full = c11_sr - u11_lr ;
c12_r_full = c12_sr - u12_lr ;
c13_r_full = c13_sr - u13_lr ;
c22_r_full = c22_sr - u22_lr ;
c23_r_full = c23_sr - u23_lr ;
c33_r_full = c33_sr - u33_lr ;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculation of beta*mu begins
sm4_h2o = 0.; % beta*mu
hr = zeros(Nr,1); cr = zeros(Nr,1); hrcr = zeros(Nr,1);
for i = 1:Nr
%%%
    hr(i) = h11_r(i)^2 + 2.*h12_r(i)^2 + 2.*h13_r(i)^2 + h22_r(i)^2 + 2.*h23_r(i)^2 + h33_r(i)^2;
%
    cr(i) = c11_r(i) + 2.*c12_r(i) + 2.*c13_r(i) + c22_r(i) + 2.*c23_r(i) + c33_r(i);
%
%    hrcr(i) = h11_r(i)*c11_r(i) + 2.*h12_r(i)*c12_r(i) + 2.*h13_r(i)*c13_r(i) + h22_r(i)*c22_r(i) + 2.*h23_r(i)*c23_r(i) + h33_r(i)*c33_r(i);
    hrcr(i) = h11_r(i)*c11_r_full(i) + 2.*h12_r(i)*c12_r_full(i) + 2.*h13_r(i)*c13_r_full(i) + h22_r(i)*c22_r_full(i) + 2.*h23_r(i)*c23_r_full(i) + h33_r(i)*c33_r_full(i) ;
%                            
    sm4_h2o = sm4_h2o + (4*pi)*rho_red_1*(0.5*hr(i) - cr(i) - 0.5*hrcr(i))*r(i)^2*dr;    %  
%%%
end
%%%
x_mu_h2o = sm4_h2o; %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculation of beta*mu ends
%%%
%%% calculation of susceptibility for water-model begins ---------------------------------
w11_k = ones(Nr,1);
chi11 = w11_k + rho_red.*h11_k_n;
chi22 = w11_k + rho_red.*h22_k_n;
chi33 = w11_k + rho_red.*h33_k_n;
%
chi12 = w12_k + rho_red.*h12_k_n;
chi13 = w13_k + rho_red.*h13_k_n;
chi23 = w23_k + rho_red.*h23_k_n;
%%% calculation of susceptibility for water-model ends ---------------------------------

%%%
return
end

%%% solute-solvent rism_1d calculation begins -----------------

function [mu_kap_ion,h14_r,h24_r,h34_r] = solute_solvent_rism_1d_solver(Nr,L,dr,itermax,chi11,chi22,chi33,chi12,chi13,chi23,tol,alf2)
%%%
%%% SPC/E- water model
sigma11 = 0.400; % Ang for H: sigma_1 = 0.4Ang
sigma22 = 3.166; % Ang for O: sigma_3 = 3.166Ang 
sigma33 = 0.400; % Ang for H: sigma_2 = 0.4Ang
%
T_red_1 = 12.870; % eps_1 = 0.0460 kcal/mol & (0.0019872041/0.0460)*298 = 12.87 for H atom & at T = 298K
T_red_2 = 3.811; % eps_3 = 0.1554 kcal/mol & (0.0019872041/0.1554)*298 = 3.811 for O atom & at T = 298K
T_red_3 = 12.87; % eps_2 = 0.0460 kcal/mol & 12.87 for H atom & at T = 298K 
%
rho_red = 0.033460;  % 
%%% SPC/E- water model
%
aa = 1.000000;
%%%
sigma44 = 3.332; % Ang for K+ ion: sigma_4 = 3.332Ang 
%%%
sigma14 = 0.5*(sigma11+sigma44);
sigma24 = 0.5*(sigma22+sigma44);
sigma34 = 0.5*(sigma33+sigma44);
%
T_red_4 = 5.922; % eps_2 = 0.1000 kcal/mol & (0.0019872041/0.1000)*298 = 5.921868217999999 for Ka+ ion & at T = 298K 
%
T_red_14 = sqrt(T_red_1*T_red_4);
T_red_24 = sqrt(T_red_2*T_red_4);
T_red_34 = sqrt(T_red_3*T_red_4);
%%%
%%% Ka+ ion part
gam14 = 30.;%237.50;%-237.50; % for H - Ka+ ion
gam24 = -2.*gam14; % for O - Ka+ ion
gam34 = gam14; % for H - Ka+ ion
%gam44 = 0;%560.40; % for Ka+ and Ka+
%%%
%
r = zeros(Nr,1); 
r4 = zeros(Nr,1); 
r14 = zeros(Nr,1); r24 = zeros(Nr,1); r34 = zeros(Nr,1);
%
%u44_lj = zeros(Nr,1); 
u14_lj = zeros(Nr,1); u24_lj = zeros(Nr,1); u34_lj = zeros(Nr,1);
%
%u44_lr = zeros(Nr,1); 
u14_lr = zeros(Nr,1); u24_lr = zeros(Nr,1); u34_lr = zeros(Nr,1);
%
%u44_sr = zeros(Nr,1); 
u14_sr = zeros(Nr,1); u24_sr = zeros(Nr,1); u34_sr = zeros(Nr,1);
%
%c44_sr = zeros(Nr,1);
c14_sr = zeros(Nr,1); c24_sr = zeros(Nr,1); c34_sr = zeros(Nr,1);
%
for i = 1:Nr
    r(i) = i*dr;
    r4(i) = r(i)/sigma44;        
    r14(i) = r(i)/sigma14;    
    r24(i) = r(i)/sigma24;    
    r34(i) = r(i)/sigma34;    
%
%    u44_lj(i) = (4./T_red_4)*( (1./r4(i))^12 - (1./r4(i))^6) ; 
%    u44_lr(i) = (gam44/r(i))*erf(aa*r(i)); 
%    u44_sr(i) = (gam44/r(i))*(1.0 - erf(aa*r(i)));
%    c44_sr(i) = exp( -u44_lj(i) -u44_sr(i) ) - 1. ;    
%    
    u14_lj(i) = (4./T_red_14)*( (1./r14(i))^12 - (1./r14(i))^6) ;
    u14_lr(i) = (gam14/r(i))*erf(aa*r(i));  
    u14_sr(i) = (gam14/r(i))*(1.0 - erf(aa*r(i)));
    c14_sr(i) = exp( -u14_lj(i) -u14_sr(i) )  - 1. ;         
%    
    u24_lj(i) = (4./T_red_24)*( (1./r24(i))^12 - (1./r24(i))^6) ;
    u24_lr(i) = (gam24/r(i))*erf(aa*r(i));  
    u24_sr(i) = (gam24/r(i))*(1.0 - erf(aa*r(i)));
    c24_sr(i) = exp( -u24_lj(i) -u24_sr(i) ) - 1. ; 
%    
    u34_lj(i) = (4./T_red_34)*( (1./r34(i))^12 - (1./r34(i))^6) ;
    u34_lr(i) = (gam34/r(i))*erf(aa*r(i));  
    u34_sr(i) = (gam34/r(i))*(1.0 - erf(aa*r(i)));
    c34_sr(i) = exp( -u34_lj(i) -u34_sr(i) )  - 1. ;     
%    
end
%%%
h14_k_old = 0.; h24_k_old = 0.;  h34_k_old = 0.;
%
Nk = Nr;
j = (1:1:Nk)';
ii = (1:1:Nr)';
dk = pi/L;
k = (j*dk);

%
c14_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c14_sr.*ii)./j;
c24_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c24_sr.*ii)./j;
c34_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c34_sr.*ii)./j;
%c44_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c44_sr.*ii)./j;
%
u14_k_lr = (4.*pi*gam14./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u24_k_lr = (4.*pi*gam24./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u34_k_lr = (4.*pi*gam34./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%u44_k_lr = (4.*pi*gam44./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%
%ck = ck_sr - uk_lr;
c14_k = c14_k_sr - u14_k_lr;
c24_k = c24_k_sr - u24_k_lr;
c34_k = c34_k_sr - u34_k_lr;
%c44_k = c44_k_sr - u44_k_lr;
%ck = ck_sr - uk_lr;
%
h14_k = c14_k.*chi11;
h14_k_new = h14_k;
h14_k_n = alf2.*h14_k_new + (1.-alf2).*h14_k_old;
%
h24_k = c24_k.*chi22;
h24_k_new = h24_k;
h24_k_n = alf2.*h24_k_new + (1.-alf2).*h24_k_old;
%
h34_k = c34_k.*chi33;
h34_k_new = h34_k;
h34_k_n = alf2.*h34_k_new + (1.-alf2).*h34_k_old;
%%%
% compute c(r) and h(r) 
c14_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c14_k_sr.*j);
c24_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c24_k_sr.*j);
c34_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c34_k_sr.*j);
%c44_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c44_k_sr.*j);
%
h14_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h14_k_n.*j);
h24_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h24_k_n.*j);
h34_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h34_k_n.*j);
%h44_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h44_k.*j);
%
G14_sr = h14_r - c14_sr;
G24_sr = h24_r - c24_sr;
G34_sr = h34_r - c34_sr;
%G44_sr = h44_r - c44_sr;
%%%
%Iflag = 0.;  % converfence flag, 1 indicates convergence
%
for iter = 2:itermax
%
%   iter
% iteration 2
    G14_sr_old = G14_sr;
    G24_sr_old = G24_sr;
    G34_sr_old = G34_sr; 
%    G44_sr_old = G44_sr;
%
    c14_sr_old = exp(-u14_lj -u14_sr + G14_sr_old) - G14_sr_old - 1.;    
    c24_sr_old = exp(-u24_lj -u24_sr + G24_sr_old) - G24_sr_old - 1.;    
    c34_sr_old = exp(-u34_lj -u34_sr + G34_sr_old) - G34_sr_old - 1.;    
%    c44_sr_old = exp(-u44_lj + G44_sr_old) - G44_sr_old - 1.;
    %
    c14_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c14_sr_old.*ii)./j;
    c24_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c24_sr_old.*ii)./j;
    c34_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c34_sr_old.*ii)./j;
%%%
%ck = ck_sr - uk_lr;
    %ck = ck_sr - uk_lr;
    c14_k = c14_k_sr - u14_k_lr;
    c24_k = c24_k_sr - u24_k_lr;
    c34_k = c34_k_sr - u34_k_lr;
%    c44_k = c44_k_sr - u44_k_lr;
    %
    h14_k = c14_k.*chi11 + c24_k.*chi12 + c34_k.*chi13;
    h14_k_old = h14_k_n; 
    h14_k_new = h14_k;
    h14_k_n = alf2.*h14_k_new + (1.-alf2).*h14_k_old;
    %
    h24_k = c14_k.*chi12 + c24_k.*chi22 + c34_k.*chi23;
    h24_k_old = h24_k_n; 
    h24_k_new = h24_k;
    h24_k_n = alf2.*h24_k_new + (1.-alf2).*h24_k_old;
    %
    h34_k = c14_k.*chi13 + c24_k.*chi23 + c34_k.*chi33;
    h34_k_old = h34_k_n; 
    h34_k_new = h34_k;
    h34_k_n = alf2.*h34_k_new + (1.-alf2).*h34_k_old;    
%%%

%%%
    % compute c(r) and h(r) 
% compute c(r) and h(r) 
c14_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c14_k_sr.*j);
c24_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c24_k_sr.*j);
c34_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c34_k_sr.*j);
%c44_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c44_k_sr.*j);
%
h14_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h14_k_n.*j);
h24_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h24_k_n.*j);
h34_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h34_k_n.*j);
%h44_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h44_k.*j);
%
%
% compute c(r) and h(r) 
%c14_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c14_k.*j);
%c24_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c24_k.*j);
%c34_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c34_k.*j);
%c44_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c44_k.*j);

%
G14_sr_new = h14_r - c14_sr;
G24_sr_new = h24_r - c24_sr;
G34_sr_new = h34_r - c34_sr;
    %
    cor = 1.*rms(G14_sr_new - G14_sr_old) + ...
          1.*rms(G24_sr_new - G24_sr_old) + ...
          1.*rms(G34_sr_new - G34_sr_old);    
%
    G14_sr = G14_sr_new;       
    G24_sr = G24_sr_new;
    G34_sr = G34_sr_new;    
%    G44_sr = G44_sr_new;    
%
    if (cor < tol)        
        Iflag = 1.
    break
    end


end
%   

%%% 2nd -----------------------------------------------------
gam14 = 237.50;%-237.50; % for H - Ka+ ion
gam24 = -2.*gam14; % for O - Ka+ ion
gam34 = gam14; % for H - Ka+ ion
gam44 = 0;%560.40; % for Ka+ and Ka+

%
r = zeros(Nr,1); 
r4 = zeros(Nr,1); 
r14 = zeros(Nr,1); r24 = zeros(Nr,1); r34 = zeros(Nr,1);
%
u44_lj = zeros(Nr,1); 
u14_lj = zeros(Nr,1); u24_lj = zeros(Nr,1); u34_lj = zeros(Nr,1);
%
u44_lr = zeros(Nr,1); 
u14_lr = zeros(Nr,1); u24_lr = zeros(Nr,1); u34_lr = zeros(Nr,1);
%
u44_sr = zeros(Nr,1); 
u14_sr = zeros(Nr,1); u24_sr = zeros(Nr,1); u34_sr = zeros(Nr,1);
%
%e44_sr = zeros(Nr,1);
%e14_sr = zeros(Nr,1); e24_sr = zeros(Nr,1); e34_sr = zeros(Nr,1);
%
for i = 1:Nr
    r(i) = i*dr;
    r4(i) = r(i)/sigma44;        
    r14(i) = r(i)/sigma14;    
    r24(i) = r(i)/sigma24;    
    r34(i) = r(i)/sigma34;    
%
    u44_lj(i) = (4./T_red_4)*( (1./r4(i))^12 - (1./r4(i))^6) ; 
    u44_lr(i) = (gam44/r(i))*erf(aa*r(i)); 
    u44_sr(i) = (gam44/r(i))*(1.0 - erf(aa*r(i)));
%    
    u14_lj(i) = (4./T_red_14)*( (1./r14(i))^12 - (1./r14(i))^6) ;
    u14_lr(i) = (gam14/r(i))*erf(aa*r(i));  
    u14_sr(i) = (gam14/r(i))*(1.0 - erf(aa*r(i)));
%    
    u24_lj(i) = (4./T_red_24)*( (1./r24(i))^12 - (1./r24(i))^6) ;
    u24_lr(i) = (gam24/r(i))*erf(aa*r(i));  
    u24_sr(i) = (gam24/r(i))*(1.0 - erf(aa*r(i)));
%    
    u34_lj(i) = (4./T_red_34)*( (1./r34(i))^12 - (1./r34(i))^6) ;
    u34_lr(i) = (gam34/r(i))*erf(aa*r(i));  
    u34_sr(i) = (gam34/r(i))*(1.0 - erf(aa*r(i)));
%    
end
%%%
%%%
%
u14_k_lr = (4.*pi*gam14./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u24_k_lr = (4.*pi*gam24./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
u34_k_lr = (4.*pi*gam34./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%u44_k_lr = (4.*pi*gam44./k.^2).*exp(-k.^2./(4.*aa^2)); % for (gamma/r)*(1.-exp(-aa*r))
%
%Iflag = 0.;  % converfence flag, 1 indicates convergence
%
for iter = 2:itermax
%
%   iter
% iteration 2
    G14_sr_old = G14_sr;
    G24_sr_old = G24_sr;
    G34_sr_old = G34_sr; 
%    G44_sr_old = G44_sr;
%
    c14_sr_old = exp(-u14_lj -u14_sr + G14_sr_old) - G14_sr_old - 1.;    
    c24_sr_old = exp(-u24_lj -u24_sr + G24_sr_old) - G24_sr_old - 1.;    
    c34_sr_old = exp(-u34_lj -u34_sr + G34_sr_old) - G34_sr_old - 1.;    
    %
    c14_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c14_sr_old.*ii)./j;
    c24_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c24_sr_old.*ii)./j;
    c34_k_sr = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c34_sr_old.*ii)./j;
%%%
%ck = ck_sr - uk_lr;
    %ck = ck_sr - uk_lr;
    c14_k = c14_k_sr - u14_k_lr;
    c24_k = c24_k_sr - u24_k_lr;
    c34_k = c34_k_sr - u34_k_lr;
    %
    h14_k = c14_k.*chi11 + c24_k.*chi12 + c34_k.*chi13;
    h14_k_old = h14_k_n; 
    h14_k_new = h14_k;
    h14_k_n = alf2.*h14_k_new + (1.-alf2).*h14_k_old;
    %
    h24_k = c14_k.*chi12 + c24_k.*chi22 + c34_k.*chi23;
    h24_k_old = h24_k_n; 
    h24_k_new = h24_k;
    h24_k_n = alf2.*h24_k_new + (1.-alf2).*h24_k_old;
    %
    h34_k = c14_k.*chi13 + c24_k.*chi23 + c34_k.*chi33;
    h34_k_old = h34_k_n; 
    h34_k_new = h34_k;
    h34_k_n = alf2.*h34_k_new + (1.-alf2).*h34_k_old;    
%%%

%%%
    % compute c(r) and h(r) 
% compute c(r) and h(r) 
c14_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c14_k_sr.*j);
c24_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c24_k_sr.*j);
c34_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c34_k_sr.*j);
%c44_sr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c44_k_sr.*j);
%
h14_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h14_k_n.*j);
h24_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h24_k_n.*j);
h34_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h34_k_n.*j);
%h44_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h44_k.*j);
%
%
% compute c(r) and h(r) 
%c14_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c14_k.*j);
%c24_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c24_k.*j);
%c34_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c34_k.*j);
%c44_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c44_k.*j);
%
G14_sr_new = h14_r - c14_sr;
G24_sr_new = h24_r - c24_sr;
G34_sr_new = h34_r - c34_sr;
    %
    cor = 1.*rms(G14_sr_new - G14_sr_old) + ...
          1.*rms(G24_sr_new - G24_sr_old) + ...
          1.*rms(G34_sr_new - G34_sr_old);    
%
    G14_sr = G14_sr_new;       
    G24_sr = G24_sr_new;
    G34_sr = G34_sr_new;    
%    G44_sr = G44_sr_new;    
%
    if (cor < tol)        
        Iflag = 1.
    break
    end


end
%   
%%%
%rho_red_4 = 0.; % dilute limit
h41_r = h14_r; h42_r = h24_r; h43_r = h34_r;
%c41_r = c14_r; c42_r = c24_r; c43_r = c34_r;
c41_r = c14_sr - u14_lr ; c42_r = c24_sr  - u24_lr ; c43_r = c34_sr  - u34_lr ;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculation of beta*mu for Ka+ ion begins
mu_kap_ion = 0.;
for i = 1:Nr
    mu_kap_ion = mu_kap_ion + rho_red *(0.5*h41_r(i)^2 - c41_r(i) - 0.5*h41_r(i)*c41_r(i)) * 4.*pi*r(i)*r(i)*dr + ...
                              rho_red *(0.5*h42_r(i)^2 - c42_r(i) - 0.5*h42_r(i)*c42_r(i)) * 4.*pi*r(i)*r(i)*dr + ...
                              rho_red *(0.5*h43_r(i)^2 - c43_r(i) - 0.5*h43_r(i)*c43_r(i)) * 4.*pi*r(i)*r(i)*dr ;
%                              rho_red_4 *(0.5*h44_r(i)^2 - c44_r(i) - 0.5*h44_r(i)*c44_r(i)) * 4.*pi*r(i)*r(i)*dr;
end
%mu_kap_ion;   %                
%%% calculation of excess chemical potential for Ka+ ion ends
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




return
end