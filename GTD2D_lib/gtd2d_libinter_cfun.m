function Diffrz=gtd2d_libinter_cfun(Grz_flat) %#codegen
%% GTD Library Interpolation
% Taking G in (r,theta,z) coordinate form and get diffusivity using
% tranformed interpolation method. Method is, first tranform to S=1, e1=0,
% then diagonalise, then recreate Diffusivity in (1,2,3), then rotate back.
%% Load Database
persistent loadmat G21_mesh2 G12_mesh2 e1_array2 e2_array2 G12_mesh G11_mesh G21_mesh
    if isempty(loadmat)
        loadmat=coder.load('GTD_beta_22_GTD_lib_32_101_2D.mat');
        [G12_mesh,G11_mesh,G21_mesh]=meshgrid(loadmat.G12_loop',loadmat.G11_loop,loadmat.G21_loop);
        [G21_mesh2,G12_mesh2]=meshgrid(loadmat.G21_loop',loadmat.G12_loop);
        e1_array2=reshape(loadmat.e1_array(1,:,:),numel(loadmat.G12_loop),numel(loadmat.G21_loop));
        e2_array2=reshape(loadmat.e2_array(1,:,:),numel(loadmat.G12_loop),numel(loadmat.G21_loop));
    end

%% Get omega in (r,theta,z)
Diffrz=zeros(5,1);
% Grz=reshape([Grz_flat -Grz_flat(1)],2,2);

% G11=Grz(1,1);
% G12=-Grz(1,2);
% G21=-Grz(2,1);
G11=Grz_flat(1);
G21=-Grz_flat(2);
G12=-Grz_flat(3);
%% Interpolate Database
for i=1:loadmat.NBvar
     Diffrz(i)=interp3(G12_mesh,G11_mesh,G21_mesh,loadmat.res_array(:,:,:,i),G12,G11,G21,'spline');
end

%% Interpolation eavg
Diffrz(4)=interp2(G21_mesh2,G12_mesh2,e1_array2,G21,G12,'spline');
Diffrz(5)=interp2(G21_mesh2,G12_mesh2,e2_array2,G21,G12,'spline');

Diffrz(5)=-Diffrz(5);
Diffrz(2)=-Diffrz(2);
end
