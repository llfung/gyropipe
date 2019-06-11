function [Drr,Drz,Dzz,er,ez]=gtd2d_libinter_cfunvec(Grr_norm,Gzr_norm,Grz_norm,Gzz_norm) %#codegen
%% GTD Library Interpolation
% Taking G in (r,theta,z) coordinate form and get diffusivity using
% tranformed interpolation method. Method is, first tranform to S=1, e1=0,
% then diagonalise, then recreate Diffusivity in (1,2,3), then rotate back.
%% Load Database
persistent loadmat G12_mesh G11_mesh G21_mesh G22_mesh
    if isempty(loadmat)
        loadmat=coder.load('GTD_beta_22_GTD_libv2_32_101_2D.mat');
        [G12_mesh,G21_mesh,G11_mesh,G22_mesh]=ndgrid(loadmat.G12_loop,loadmat.G21_loop,loadmat.G11_loop,loadmat.G22_loop);
    end
    
%% Get omega in (r,theta,z)
Drr=zeros(size(Grr_norm,1),1);
Drz=zeros(size(Grr_norm,1),1);
Dzz=zeros(size(Grr_norm,1),1);
er=zeros(size(Grr_norm,1),1);
ez=zeros(size(Grr_norm,1),1);

Gzr_norm=-Gzr_norm;
Grz_norm=-Grz_norm;

omega_e3=Grz_norm-Gzr_norm;

%% Interpolate Database
Drr(:,1)=interpn(G12_mesh,G21_mesh,G11_mesh,G22_mesh,loadmat.res_array(:,:,:,:,1),Grz_norm,Gzr_norm,Grr_norm,Gzz_norm,'linear');
Drz(:,1)=-interpn(G12_mesh,G21_mesh,G11_mesh,G22_mesh,loadmat.res_array(:,:,:,:,2),Grz_norm,Gzr_norm,Grr_norm,Gzz_norm,'linear');
Dzz(:,1)=interpn(G12_mesh,G21_mesh,G11_mesh,G22_mesh,loadmat.res_array(:,:,:,:,3),Grz_norm,Gzr_norm,Grr_norm,Gzz_norm,'linear');

%% Interpolation eavg
er(:,1)=interp1(loadmat.omega_e3_col,loadmat.e1_col,omega_e3,'spline');
ez(:,1)=-interp1(loadmat.omega_e3_col,loadmat.e2_col,omega_e3,'spline');

end

