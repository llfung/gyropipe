function Diff_rtz=gtd_direct_cfun(Grtz_flat) %#codegen

%Input reshape
Grtz=reshape([Grtz_flat -Grtz_flat(1)-Grtz_flat(5)],3,3);

%% Define discretization parameters
settings=struct( ...
   'n_phi',16,...% Has to be even for FFT. 2^N recommended
   'n_theta',31,...% Had better to be 1+(multiples of 5,4,3 or 2) for Newton-Cote
   'N_total',0,...
   'tn_phi',0,...
   'tn_theta',0,...
   'dtheta',0.0,...
   'dphi',0.0,...
   'theta',zeros(31,1),...
   'phi',zeros(1,16),...
   'kphi',zeros(1,9),...
   'integrad',zeros(1,31),...
   'beta',2.2,'S',1.0,...
   'G11',0.0,'G12',0.0,'G13',0.0,...
   'G21',0.0,'G22',0.0,'G23',0.0,...
   'G31',0.0,'G32',0.0,'G33',0.0,...
   'omega_e1',0.0,'omega_e2',0.0,'omega_e3',0.0...
   );

settings.N_total=settings.n_phi*settings.n_theta;
settings.tn_phi=settings.n_phi/2+1;
settings.tn_theta=(floor(settings.n_theta/2)+1);

settings.dtheta=(pi/(settings.n_theta-1));
settings.dphi=2*pi/(settings.n_phi);

settings.theta=[0:settings.dtheta:pi]';
settings.phi=[0:settings.dphi:(2*pi-settings.dphi)];
settings.kphi=[0:(settings.tn_phi-1)];

e_all_field=zeros(settings.n_theta,settings.n_phi,3);
e_all_field(:,:,1)=sin(settings.theta)*cos(settings.phi);
e_all_field(:,:,2)=cos(settings.theta)*ones(size(settings.phi));
e_all_field(:,:,3)=sin(settings.theta)*sin(settings.phi);

b_endings=[1 settings.n_theta settings.n_theta+1 2*settings.n_theta 2*settings.n_theta+1 3*settings.n_theta];

%% Newton-Cote Integrand
    if mod(settings.n_theta,5)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[19 repmat([75 50 50 75 38],1,(settings.n_theta-6)/5)  75 50 50 75 19]*5/288*settings.dtheta*2*pi;
    elseif mod(settings.n_theta,4)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[7 repmat([32 12 32 14],1,(settings.n_theta-5)/4)  32 12 32 7]*2/45*settings.dtheta*2*pi;
    elseif mod(settings.n_theta,3)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[1 repmat([3 3 2],1,(settings.n_theta-4)/3) 3 3 1]*3/8*settings.dtheta*2*pi;
    elseif mod(settings.n_theta,2)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[1 repmat([4 2],1,(N-3)/2) 4 1]/3*settings.dtheta*2*pi;
    else
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0]*settings.dtheta*pi*2; %% Trapezoid Rule
    end
    
%% Define input parameters
settings.beta=2.2;

% Rotate
settings.G11=Grtz(1,1);
settings.G12=-Grtz(1,3);
settings.G13=Grtz(1,2);
settings.G21=-Grtz(3,1);
settings.G22=Grtz(3,3);
settings.G23=-Grtz(3,2);
settings.G31=Grtz(2,1);
settings.G32=-Grtz(2,3);
settings.G33=-settings.G11-settings.G22;
settings.omega_e1=settings.G23-settings.G32;
settings.omega_e2=settings.G31-settings.G13;
settings.omega_e3=settings.G12-settings.G21;

% Magnitude
settings.S=1;
    
    %% Solving for f(e)
    LHSL=l_fd_lhs_coder(settings);
    LHS_F=l_fd_lhs_bc_coder(LHSL,settings);
    LHS_F(settings.N_total/2,:)=kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi);
    RHS_F=zeros(settings.n_theta*settings.n_phi,1);
    RHS_F(settings.N_total/2,1)=1;
    f_sol=LHS_F\RHS_F;
    f0=transpose(reshape(f_sol,settings.n_phi,settings.n_theta));
    
    %% Finding e_avg and FK Diffusitivity
    ebracket=eavg(f0,e_all_field,settings.integrad);   
    
    %% Solving for b(e)
%     LHS=b_fd_lhs_coder(settings);
        %% bG term
        G=[settings.G11 settings.G21 settings.G31;settings.G12 settings.G22 settings.G32;settings.G13 settings.G23 settings.G33];
        bG=settings.S*kronsp(sparse(G),speye(settings.n_theta*settings.n_phi));

        %% Combined
        LHS=-kronsp(speye(3),LHSL)-bG;
        LHS=b_fd_lhs_bc_coder(LHS,settings);
%     f0_col=reshape(transpose(f0),(settings.n_theta)*settings.n_phi,1);
    RHS_mat=[f0.*(e_all_field(:,:,1)-ebracket.e1);f0.*(e_all_field(:,:,2)-ebracket.e2);f0.*(e_all_field(:,:,3)-ebracket.e3)];
    RHS_mat(b_endings,:)=0;
    RHS=reshape(transpose(RHS_mat),3*(settings.n_theta)*settings.n_phi,1);

    % Normalisation (Legacy)
%        LHS=[LHS; kron(speye(3),kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi))];
%        RHS=[RHS;0;0;0];
    % Normalisation (Sq. Matrix)
    LHS(settings.N_total/2,:)=[kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi) zeros(1,2*settings.N_total)];
    LHS(settings.N_total/2*3,:)=[zeros(1,settings.N_total) kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi) zeros(1,settings.N_total)];
    LHS(settings.N_total/2*5,:)=[zeros(1,2*settings.N_total) kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi)];
    RHS([settings.N_total/2 settings.N_total/2*3 settings.N_total/2*5],1)=0;
    RHS=sparse(LHS)\RHS;
    b=transpose(reshape(RHS,settings.n_phi,settings.n_theta*3));
    
    b1=b(b_endings(1):b_endings(2),:);
    b2=b(b_endings(3):b_endings(4),:);
    b3=b(b_endings(5):b_endings(6),:);

    %% Finding Diffusitivity

    b1G=b1*settings.G11+b2*settings.G21+b3*settings.G31;
    b2G=b1*settings.G12+b2*settings.G22+b3*settings.G32;
    b3G=b1*settings.G13+b2*settings.G23+b3*settings.G33;
    
    D11t=b1.*e_all_field(:,:,1)+settings.S*b1.*b1G./f0; %#ok<*PFBNS>
    D12t=b1.*e_all_field(:,:,2)+settings.S*b1.*b2G./f0;
    D13t=b1.*e_all_field(:,:,3)+settings.S*b1.*b3G./f0;
    D21t=b2.*e_all_field(:,:,1)+settings.S*b2.*b1G./f0;
    D22t=b2.*e_all_field(:,:,2)+settings.S*b2.*b2G./f0;
    D23t=b2.*e_all_field(:,:,3)+settings.S*b2.*b3G./f0;
    D31t=b3.*e_all_field(:,:,1)+settings.S*b3.*b1G./f0;
    D32t=b3.*e_all_field(:,:,2)+settings.S*b3.*b2G./f0;
    D33t=b3.*e_all_field(:,:,3)+settings.S*b3.*b3G./f0;
    
    %Symmetric Part
    D12t=(D12t+D21t)/2;
    D23t=(D23t+D32t)/2;
    D13t=(D13t+D31t)/2;
    
    Diff=zeros(3,3);
    Diff(1,1)=settings.integrad*mean(D11t,2);
    Diff(1,2)=settings.integrad*mean(D12t,2);Diff(2,1)=settings.integrad*mean(D12t,2);
    Diff(1,3)=settings.integrad*mean(D13t,2);Diff(3,1)=settings.integrad*mean(D13t,2);
    Diff(2,2)=settings.integrad*mean(D22t,2);
    Diff(2,3)=settings.integrad*mean(D23t,2);Diff(3,2)=settings.integrad*mean(D23t,2);
    Diff(3,3)=settings.integrad*mean(D33t,2);   

%% Rotate Diffusivity and avg swimming back to (r,theta,z)
Diff_rtz=zeros(1,9);
Diff_rtz(1,1)=Diff(1,1);
Diff_rtz(1,2)=Diff(1,3);
Diff_rtz(1,3)=-Diff(1,2);
Diff_rtz(1,4)=Diff(3,3); %Dtt
Diff_rtz(1,5)=-Diff(3,2); %Dtz
Diff_rtz(1,6)=Diff(2,2); %Dzz
Diff_rtz(1,7)=ebracket.e1; %er
Diff_rtz(1,8)=ebracket.e3; %et
Diff_rtz(1,9)=-ebracket.e2; %ez
end