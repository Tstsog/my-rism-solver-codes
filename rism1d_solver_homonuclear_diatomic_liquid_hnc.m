function [] = rism1d_solver_homonuclear_diatomic_liquid_hnc
%-----------------------------------------------------------------
% Solves the Reference-interaction site model (RISM) equation for uncharged
% diatomic liquid in hypernetted chain approximation
% An interparticle interaction: the Lennard-Jones potential
% A method: a Picard iteration technique
% Uses: function lf=lsint(f)
%
% Written by Tsogbayar Tsednee (PhD), California State University Northridge 
% Contact: tsog215@gmail.com
% Reference: E. Johnson and R. P. Hazoume, J. Chem. Phys. 70, 1599 (1979);
% Date: October 21, 2018
%%%
% Temperature T = 328K & density  = 0.00633 A^(-3)
%%%
% ----------------------------------------------------------------
clear;
clc;
format short
Nr = 512; % number of grid point; you may change it
L = 32.0; % grid parameter; you may change it  
dr = L/(Nr+1); % grid parameter
itermax = 4000; tol = 10^(-14); % max number of iteration and tolerance; you may change it
%
alf = 0.850 ; % convergence acceleration parameter; you may change it
%%%
% LJ parameters for benzene type model
ell = 1.756; % bond length, l = 1.756 angstrom; (Ang = angstrom)
sigma = 3.5; % in angstrom
%
rho_red = 0.00633; % in A^(-3)
T_red = 328/77; % eps/k_B = 77K & T = 328/77 = 4.260 & T = 328K
%
r = zeros(Nr,1); u_lj = zeros(Nr,1); c = zeros(Nr,1); 
for i = 1:Nr
    r(i) = i*dr;
    u_lj(i) = (4./T_red)*( (sigma/r(i))^12 - (sigma/r(i))^6 );
%    
    c(i) = exp( -u_lj(i) ) - 1. ;
end
%%%
Nk = Nr; j = (1:1:Nk)'; ii = (1:1:Nr)'; dk = pi/L; k = (j*dk);
%
w_k = sin(k.*ell)./(k.*ell); % w(k), the Fourier transforms of the intramolecular correlation functions
%
hk_old = 0.;
%
ck = 4.*pi*(dr^3*(Nr+1)./pi).*lsint(c.*ii)./j;
%%%
hk = (ck.*(1.0+w_k).^2./(1.-2.*rho_red.*(1.0+w_k).*ck)); % solving RISM or scalar SSOZ equation in k-space; (SSOZ = site-site Ornstein-Zernike)
%
hk_new = hk;
hk_n = alf*hk_new + (1.-alf)*hk_old;
% compute c(r) and gam(r) 
cr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(ck.*j);
hr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(hk_n.*j);
%%%
gam = hr - cr; % indirect correlation function, gamma(r) = h(r) - c(r)
%%%
%
%Iflag = 0.;  % converfence flag, 1 indicates convergence
%
for iter = 2:itermax
%    
    iter
    %%% iter 2
    gam_old = gam;
%
    cr_old = exp(-u_lj + gam_old) - gam_old - 1. ;    
    ck = 4.*pi.*(dr^3*(Nr+1)./pi).*lsint(cr_old.*ii)./j;
%%%
    hk = (ck.*(1.0+w_k).^2./(1.-2.*rho_red.*(1.0+w_k).*ck)); % solving SSOZ equation in k-space
%    
    hk_old = hk_n;
    hk_new = hk;
    hk_n = alf*hk_new + (1.-alf)*hk_old;  
% compute c(r) and gam(r) 
    cr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(ck.*j);
    hr = (1./(2*pi)^(3))*4*pi.*(pi^2./(ii.*(Nr+1)^2.*dr^3)).*lsint(hk_n.*j);
%%%
    gam = hr - cr; % indirect correlation function
%%%
    gam_new = gam;
%    
    if (rms(gam_new - gam_old) < tol);
         Iflag = 1.
         break
    end
%
    gam = gam_new;
%%%
end
%
gr = exp(-u_lj + gam); %    radial correlation function
%
%figure(1)
plot(r,gr, '-b')
axis([0. 12. 0.000 1.30])
set(gca,'FontSize',20)
xlabel('r(Ang)') % ,'fontsize',16
ylabel('g(r)','Rotation', 1) %

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

