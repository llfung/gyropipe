%% Script to generator state.cdf.in

%% Initiazing Variables and Attributes
i_N=120;
i_K=32;
i_M=1;
i_Mp=1;
ReIm=2;

t     = 0.0;
Re    = 0.401070456591576;
Ri    = 4.322797117392170;
dr    = 0.668211770763543;
Vs    = 0.314159265358979;
Pe    = 1/Vs/Vs*dr;
alpha = 0.15;

dt=1e-8;
dtcor=1e-8;
dtcfl=1e-8;
nint=0.5;


%% Creating Variables
h=zeros(i_N,i_K*i_M,ReIm);
ur=zeros(i_N,i_K*i_M,ReIm);
ut=zeros(i_N,i_K*i_M,ReIm);
uz=zeros(i_N,i_K*i_M,ReIm);
r=zeros(i_N,1);

for i=1:i_N
    r(i)=(1+cos(pi*(i_N-i)/i_N))/2;
end
%% Modifying initial condition
epsilon=1e-11;
h(:,2:end,:)=(2*rand(i_N,i_K*i_M-1,ReIm)-1)*epsilon;


%% Creating NETCDF File
nfid=netcdf.create('state.cdf.in','CLOBBER');

%% Define Variables and Dimensions
% Attributes
global_id=netcdf.getConstant('NC_GLOBAL');

netcdf.putAtt(nfid,global_id,'t',t);
netcdf.putAtt(nfid,global_id,'Re',Re);
netcdf.putAtt(nfid,global_id,'Ri',Ri);
netcdf.putAtt(nfid,global_id,'dr',dr);
netcdf.putAtt(nfid,global_id,'Vs',Vs);
netcdf.putAtt(nfid,global_id,'Pe',Pe);
netcdf.putAtt(nfid,global_id,'alpha',alpha);

% Dimensions
r_id   =netcdf.defDim(nfid,'r' ,i_N);
H_id   =netcdf.defDim(nfid,'H' ,i_K*i_M);
ReIm_id=netcdf.defDim(nfid,'ReIm' ,ReIm);

r_var_id    =netcdf.defVar(nfid,    'r','NC_DOUBLE',r_id);
dt_var_id   =netcdf.defVar(nfid,   'dt','NC_DOUBLE',[]);
dtcor_var_id=netcdf.defVar(nfid,'dtcor','NC_DOUBLE',[]);
dtcfl_var_id=netcdf.defVar(nfid,'dtcfl','NC_DOUBLE',[]);
nint_var_id =netcdf.defVar(nfid, 'nint','NC_DOUBLE',[]);

dim=[r_id,H_id,ReIm_id];
Ur_id       =netcdf.defVar(nfid,'Ur','NC_DOUBLE',dim);
Ut_id       =netcdf.defVar(nfid,'Ut','NC_DOUBLE',dim);
Uz_id       =netcdf.defVar(nfid,'Uz','NC_DOUBLE',dim);
 T_id       =netcdf.defVar(nfid, 'T','NC_DOUBLE',dim);

netcdf.putAtt(nfid,Ur_id,'K',i_K);
netcdf.putAtt(nfid,Ur_id,'M',i_M);
netcdf.putAtt(nfid,Ur_id,'Mp',i_Mp);
netcdf.putAtt(nfid,Ut_id,'K',i_K);
netcdf.putAtt(nfid,Ut_id,'M',i_M);
netcdf.putAtt(nfid,Ut_id,'Mp',i_Mp);
netcdf.putAtt(nfid,Uz_id,'K',i_K);
netcdf.putAtt(nfid,Uz_id,'M',i_M);
netcdf.putAtt(nfid,Uz_id,'Mp',i_Mp);
netcdf.putAtt(nfid, T_id,'K',i_K);
netcdf.putAtt(nfid, T_id,'M',i_M);
netcdf.putAtt(nfid, T_id,'Mp',i_Mp);

%% Write Variables
netcdf.endDef(nfid);

netcdf.putVar(nfid,r_var_id, r);
netcdf.putVar(nfid,   Ur_id,ur);
netcdf.putVar(nfid,   Ut_id,ut);
netcdf.putVar(nfid,   Uz_id,uz);
netcdf.putVar(nfid,    T_id, h);

netcdf.putVar(nfid,   dt_var_id,   dt);
netcdf.putVar(nfid,dtcor_var_id,dtcor);
netcdf.putVar(nfid,dtcfl_var_id,dtcfl);
netcdf.putVar(nfid, nint_var_id, nint);

% Close netcdf file
netcdf.close(nfid);
