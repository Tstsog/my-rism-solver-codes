function [] = rism1d_solver_heteronuclear_diatomic_liquid_hnc
%-----------------------------------------------------------------
% Solves the Reference-interaction site model (RISM) equation for uncharged
% heteronuclear diatomic liquid in hypernetted chain approximation
% An interparticle interaction: the Lennard-Jones potential
% A method: a Picard iteration technique
% Uses: function lf=lsint(f)
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Reference: F. Hirata, B. M. Pettitt and P. J. Rossky, J. Chem. Phys. 77, 509 (1982);
% Date: November 7, 2018
%%%
% Temperature T = 210K & density  = 0.018 A^(-3)
%%%
% ----------------------------------------------------------------
clear;
clc;
format short
Nr = 256; % number of grid point; you may change it
L = 16.0; % grid parameter; you may change it  
dr = L/(Nr+1); % grid parameter
itermax = 4000; tol = 10^(-5); % max number of iteration and tolerance; you may change it
%
alf = 0.9500 ; % convergence acceleration parameter; you may change it
%%%
% LJ parameters for HCl-like model
sigma11 = 0.400; % Ang for H: sigma_1 = 0.4Ang (Ang = angstrom)
sigma22 = 3.353; % Ang for Cl: sigma_2 = 3.353Ang
sigma12 = 0.5*(sigma11+sigma22);
%
ell = 1.30; % distance between two sites, in angstrom
%%%
rho_red_1 = 0.018; % Ang^(-3) x_1 = 0.5;
rho_red_2 = 0.018; % Ang^(-3) x_2 = 0.5;
%%%
T_red_1 = 10.500; % eps_1/k_B = 20K for H & T = 210K
T_red_2 = 0.8108; % eps_1/k_B = 259K for Cl & T = 210K
T_red_12 = sqrt(T_red_1*T_red_2);
%%%
r = zeros(Nr,1); r2 = zeros(Nr,1); r1 = zeros(Nr,1); r12 = zeros(Nr,1);
U11 = zeros(Nr,1); c11 = zeros(Nr,1);
U12 = zeros(Nr,1); c12 = zeros(Nr,1);
U22 = zeros(Nr,1); c22 = zeros(Nr,1);
for i = 1:Nr
    r(i) = i*dr;
    r1(i) = r(i)/sigma11;
    r2(i) = r(i)/sigma22;    
    r12(i) = r(i)/sigma12;        
%
    U11(i) = (4./T_red_1)*( (1./r1(i))^12 - (1./r1(i))^6) ;
    c11(i) = exp(-U11(i)) - 1.;
%
    U22(i) = (4./T_red_2)*( (1./r2(i))^12 - (1./r2(i))^6) ;
    c22(i) = exp(-U22(i)) - 1.;
%
    U12(i) = (4./T_red_12)*( (1./r12(i))^12 - (1./r12(i))^6) ;
    c12(i) = exp(-U12(i)) - 1.;  
end
%%%
h_k_mat_old = 0.;
%
Nk = Nr; j = (1:1:Nk)'; ii = (1:1:Nr)'; dk = pi/L; k = (j*dk);
%
w_k = sin(k.*ell)./(k.*ell); % w(k), the Fourier transforms of the intramolecular correlation functions
%
w11 = 1.*diag(ones(Nr,1)); w12 = diag(w_k); w21 = w12; w22 = w11;
omega_k = [w11, w12; w21, w22]; % matrix form
%%%
rho_mat = [(rho_red_1).*diag(ones(Nr,1)), diag(zeros(Nr,1)); diag(zeros(Nr,1)), (rho_red_2).*diag(ones(Nr,1))];
%
c11_k = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c11.*ii)./j;
c12_k = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c12.*ii)./j;
c21_k = c12_k;
c22_k = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c22.*ii)./j;
c_k_mat = [diag(c11_k), diag(c12_k); diag(c21_k), diag(c22_k)]; 
%
% solve SSOZ: H(k) = (I - W(k)*C(k)*Rho)^(-1)*W(k)*C(k)*W(k), 
% (SSOZ = site-site Ornstein-Zernike)
% M_denom (= I - W(k)*C(k)*Rho ) in k-space
M_denom = sparse(diag(ones(2.*Nr,1)) - omega_k*c_k_mat*rho_mat); % rho is equal
h_k_mat = M_denom\sparse(omega_k*c_k_mat*omega_k); % omega_k.*c_k_mat.*inv(M_denom).*omega_k
%
h_k_mat_new = h_k_mat;
h_k_mat_n = alf.*h_k_mat_new + (1.-alf).*h_k_mat_old;
%
h11_k = diag(h_k_mat_n(1:Nr,1:Nr)); 
h12_k = diag(h_k_mat_n(1:Nr,Nr+1:2*Nr)); 
%h21_k = h12_k;
h22_k = diag(h_k_mat_n(Nr+1:2*Nr,Nr+1:2*Nr) );
%%%
% compute c(r) and h(r) 
c11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c11_k.*j);
c12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c12_k.*j);
%c21_r = c12_r;
c22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c22_k.*j);
%
h11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h11_k.*j);
h12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h12_k.*j);
%h21_r = h12_r;
h22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h22_k.*j);
%
G11_r = h11_r - c11_r; % indirect cf
G12_r = h12_r - c12_r;
%G21_r = G12_r;
G22_r = h22_r - c22_r;
%%%
%Iflag = 0.;  % converfence flag, 1 indicates convergence
%
for iter = 2:itermax
%
    iter
% iteration 2
    G11_r_old = G11_r;
    G12_r_old = G12_r;
    G22_r_old = G22_r;
%
    c11_r_old = exp(-U11 + G11_r_old ) - G11_r_old - 1.;
    c12_r_old = exp(-U12 + G12_r_old ) - G12_r_old - 1.;
    c22_r_old = exp(-U22 + G22_r_old ) - G22_r_old - 1.;
%
    c11_k = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c11_r_old.*ii)./j;
    c12_k = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c12_r_old.*ii)./j;
    c21_k = c12_k;
    c22_k = 4*pi*(dr^3*(Nr+1)./pi).*lsint(c22_r_old.*ii)./j;
    c_k_mat = [diag(c11_k), diag(c12_k); diag(c21_k), diag(c22_k)]; % rho = equal
%%%
%   solving SSOZ: H(k) = (I - W(k)*C(k)*Rho)^(-1)*W(k)*C(k)*W(k)
%   M_denom (= I - W(k)*C(k)*Rho ) in k-space
    M_denom = sparse(diag(ones(2.*Nr,1)) - omega_k*c_k_mat*rho_mat); % rho is equal
    h_k_mat = M_denom\sparse(omega_k*c_k_mat*omega_k); % omega_k.*c_k_mat.*inv(M_denom).*omega_k
%
    h_k_mat_old = h_k_mat_n;
    h_k_mat_new = h_k_mat;
    h_k_mat_n = alf*h_k_mat_new + (1.-alf)*h_k_mat_old;
    %
    h11_k = diag(h_k_mat_n(1:Nr,1:Nr)); 
    h12_k = diag(h_k_mat_n(1:Nr,Nr+1:2*Nr)); 
    h22_k = diag(h_k_mat_n(Nr+1:2*Nr,Nr+1:2*Nr) );
    %%%
    % compute c(r) and h(r) 
    c11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c11_k.*j);
    c12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c12_k.*j);
    c22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(c22_k.*j);
    %
    h11_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h11_k.*j);
    h12_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h12_k.*j);
    %h21_r = h12_r;
    h22_r = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(h22_k.*j);
    %
    G11_r = h11_r - c11_r; % indirect cf
    G12_r = h12_r - c12_r;
    G22_r = h22_r - c22_r;
    %%%
    G11_r_new = G11_r; % indirect cf
    G12_r_new = G12_r;
    G22_r_new = G22_r;
    %
    cor = rms(G11_r_new - G11_r_old) + rms(G12_r_new - G12_r_old) + ...
          rms(G22_r_new - G22_r_old);    
%
    if (rms(cor) < tol);
         Iflag = 1.
         break
    end
%
    G11_r = G11_r_new;
    G12_r = G12_r_new;
    G22_r = G22_r_new;    
%
end
%%%
figure(1)
hold on
plot(r, h11_r+1, 'b') % g11 is for (H-H) 
plot(r, h12_r+1, 'k') % g12 is for (H-Cl)
plot(r, h22_r+1, 'g') % g22 is for (Cl-Cl)
hold off
axis([0. 16. 0. 3.00])
set(gca,'FontSize',18)
xlabel('r(Ang)') % ,'fontsize',16
ylabel('g_{ij}(r)','Rotation', 1) %
%
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