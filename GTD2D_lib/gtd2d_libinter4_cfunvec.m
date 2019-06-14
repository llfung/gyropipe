function [Drr,Drz,Dzz,er,ez]=gtd2d_libinter4_cfunvec(Grr_norm,Gzr_norm,Grz_norm,Gzz_norm) %#codegen
%% GTD Library Interpolation
% Taking G in (r,theta,z) coordinate form and get diffusivity using
% tranformed interpolation method. Method is, first tranform to S=1, e1=0,
% then diagonalise, then recreate Diffusivity in (1,2,3), then rotate back.
%% Load Database
persistent loadmat G12_mesh G11_mesh G21_mesh G22_leng G22_step
    if isempty(loadmat)
        loadmat=coder.load('GTD_beta_22_GTD_libv2_32_101_2D.mat');
        [G21_mesh,G12_mesh,G11_mesh]=meshgrid(loadmat.G21_loop,loadmat.G12_loop,loadmat.G11_loop);
        G22_leng=length(loadmat.G22_loop);
        G22_step=loadmat.G22_loop(floor(G22_leng/2))-loadmat.G22_loop(floor(G22_leng/2)-1);
    end
    
%% Get omega in (r,theta,z)
inp_siz=size(Grr_norm,1);
Drr=zeros(inp_siz,1);
Drz=zeros(inp_siz,1);
Dzz=zeros(inp_siz,1);
er=zeros(inp_siz,1);
ez=zeros(inp_siz,1);

Gzz_maxi=min(ceil((max(Gzz_norm)-loadmat.G22_loop(1))/G22_step)+3,G22_leng);
Gzz_mini=max(floor((min(Gzz_norm)-loadmat.G22_loop(1))/G22_step)-3,1);
Drr_int=zeros(inp_siz,G22_leng);
Drz_int=zeros(inp_siz,G22_leng);
Dzz_int=zeros(inp_siz,G22_leng);


Gzr_norm=-Gzr_norm;
Grz_norm=-Grz_norm;

omega_e3=Grz_norm-Gzr_norm;

%% Interpolate Database
for i=Gzz_mini:Gzz_maxi
    Drr_int(:,i)= interp3(G21_mesh,G12_mesh,G11_mesh,loadmat.res_array(:,:,:,i,1),Gzr_norm,Grz_norm,Grr_norm,'spline');
    Drz_int(:,i)=-interp3(G21_mesh,G12_mesh,G11_mesh,loadmat.res_array(:,:,:,i,2),Gzr_norm,Grz_norm,Grr_norm,'spline');
    Dzz_int(:,i)= interp3(G21_mesh,G12_mesh,G11_mesh,loadmat.res_array(:,:,:,i,3),Gzr_norm,Grz_norm,Grr_norm,'spline');
end
for i=1:inp_siz
    Drr(i)=interp1(loadmat.G22_loop(Gzz_mini:Gzz_maxi),Drr_int(i,Gzz_mini:Gzz_maxi),Gzz_norm(i),'spline');
    Drz(i)=interp1(loadmat.G22_loop(Gzz_mini:Gzz_maxi),Drz_int(i,Gzz_mini:Gzz_maxi),Gzz_norm(i),'spline');
    Dzz(i)=interp1(loadmat.G22_loop(Gzz_mini:Gzz_maxi),Dzz_int(i,Gzz_mini:Gzz_maxi),Gzz_norm(i),'spline');
end
%% Interpolation eavg
er(:,1)=interp1(loadmat.omega_e3_col,loadmat.e1_col,omega_e3,'spline');
ez(:,1)=-interp1(loadmat.omega_e3_col,loadmat.e2_col,omega_e3,'spline');
end

