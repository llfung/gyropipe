function Diff_rtz=gtd_eig_cfun(Grtz_flat) %#codegen
% persistent Grtz settings e_all_field b_endings LHSL LHS_F RHS_F f_sol f0 ebracket G W D G_flag G_dia_flag ...
%     ebracket_array RHS_mat RHS LHS B ifind G_small bG2 LHS2 RHS_temp1 RHS_temp2 RHS_mat2 RHS2 b1 b2 b3 bG ...
%     LHS3 RHS_mat3 RHS3 b Vinv Bc RHSc RHS_matc LHSc LHS_Fc b1_col b2_col b3_col ...
%     b1G b2G b3G D11t D12t D13t D21t D23t D22t D31t D32t D33t Diff

%Input reshape
Grtz=reshape([Grtz_flat -Grtz_flat(1)-Grtz_flat(5)],3,3);

%% Define discretization parameters
settings=struct( ...
   'n_phi',16,...% Has to be even for FFT. 2^N recommended
   'n_theta',31,...% Had better to be 1+(multiples of 5,4,3 or 2) for Newton-Cote
   'N_total',16*31,...
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

% settings.N_total=settings.n_phi*settings.n_theta;
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
%     if mod(settings.n_theta,5)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[19 repmat([75 50 50 75 38],1,(settings.n_theta-6)/5)  75 50 50 75 19]*5/288*settings.dtheta*2*pi;
%     elseif mod(settings.n_theta,4)==1
%         settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[7 repmat([32 12 32 14],1,(settings.n_theta-5)/4)  32 12 32 7]*2/45*settings.dtheta*2*pi;
%     elseif mod(settings.n_theta,3)==1
%         settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[1 repmat([3 3 2],1,(settings.n_theta-4)/3) 3 3 1]*3/8*settings.dtheta*2*pi;
%     elseif mod(settings.n_theta,2)==1
%         settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[1 repmat([4 2],1,(N-3)/2) 4 1]/3*settings.dtheta*2*pi;
%     else
%         settings.integrad=[0 sin(settings.theta(2:end-1)') 0]*settings.dtheta*pi*2; %% Trapezoid Rule
%     end
    
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
    %% eig G
    G=[settings.G11 settings.G12 settings.G13;settings.G21 settings.G22 settings.G23;settings.G31 settings.G32 settings.G33];
    [W,D]=eig(G);
    if cond(W)>1e3
        %% Zero out
        G=G';
        G_flag=abs(G/norm(G))<1e-6;
        G(G_flag)=0;
        G_dia_flag=[G_flag(1,2) && G_flag(1,3);G_flag(2,1) && G_flag(2,3);G_flag(3,1) && G_flag(3,2)];
        if any(G_dia_flag)
            B=zeros(settings.N_total,3);
            ebracket_array=zeros(1,3);
            ebracket_array(1)=ebracket.e1;
            ebracket_array(2)=ebracket.e2;
            ebracket_array(3)=ebracket.e3;
            for i=1:3
                if G_dia_flag(i)
                    RHS_mat=e_all_field(:,:,i)-ebracket_array(i);
                    RHS_mat=RHS_mat.*f0;
                    RHS_mat([1,settings.n_theta],:)=0;
                    RHS=reshape(transpose(RHS_mat),(settings.n_theta)*settings.n_phi,1);

                    LHS=-LHSL-speye(settings.N_total)*G(i,i)*settings.S;
                    LHS_F=l_fd_lhs_bc_coder(LHS,settings);
                    LHS_F(settings.N_total/2,:)=kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi);
                    RHS(settings.N_total/2,1)=0;
                    B(:,i)=LHS_F\RHS;
                end
            end
            if sum(G_dia_flag)~=3
            if sum(G_dia_flag)==2
                ifind=find(~G_dia_flag,1);
                i=ifind(1);
                RHS_mat=e_all_field(:,:,i)-ebracket_array(i);
                RHS_mat=RHS_mat.*f0;
                for j=1:3
                    if G_dia_flag(j)
                        RHS_mat=RHS_mat+G(i,j)*transpose(reshape(B(:,j),settings.n_phi,settings.n_theta));
                    end
                end
                RHS_mat([1,settings.n_theta],:)=0;
                RHS=reshape(transpose(RHS_mat),(settings.n_theta)*settings.n_phi,1);

                LHS=-LHSL-speye(settings.N_total)*G(i,i)*settings.S;
                LHS_F=l_fd_lhs_bc_coder(LHS,settings);
                LHS_F(settings.N_total/2,:)=kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi);
                RHS(settings.N_total/2,1)=0;
                B(:,i)=LHS_F\RHS;
            else
                i_=find(~G_dia_flag,2);
                j_f=find(G_dia_flag,1);
                j_=j_f(1);
                G_small=G(i_,i_);
                %% 2x2 G
                bG2=settings.S*kronsp(sparse(G_small),speye(settings.n_theta*settings.n_phi));

                %% Combined
                LHS2=-kronsp(speye(2),LHSL)-bG2;
                LHS2=btwo_fd_lhs_bc_coder(LHS2,settings);
                RHS_temp1=f0.*(e_all_field(:,:,i_(1))-ebracket_array(i_(1)))+G(i_(1),j_)*transpose(reshape(B(:,j_),settings.n_phi,settings.n_theta));
                RHS_temp2=f0.*(e_all_field(:,:,i_(2))-ebracket_array(i_(2)))+G(i_(2),j_)*transpose(reshape(B(:,j_),settings.n_phi,settings.n_theta));
                RHS_mat2=[RHS_temp1;RHS_temp2];
                RHS_mat2(b_endings(1:4),:)=0;
                RHS2=reshape(transpose(RHS_mat2),2*(settings.n_theta)*settings.n_phi,1);

                % Normalisation (Sq. Matrix)
                LHS2(settings.N_total/2,:)=[kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi) zeros(1,settings.N_total)];
                LHS2(settings.N_total/2*3,:)=[zeros(1,settings.N_total) kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi)];
                RHS2([settings.N_total/2 settings.N_total/2*3],1)=0;
                RHS2=LHS2\RHS2;
                B(:,i_(1))=RHS2(1:settings.N_total,1);
                B(:,i_(2))=RHS2(settings.N_total+1:2*settings.N_total,1);
            end
            end
            b1=transpose(reshape(B(:,1),settings.n_phi,settings.n_theta));
            b2=transpose(reshape(B(:,2),settings.n_phi,settings.n_theta));
            b3=transpose(reshape(B(:,3),settings.n_phi,settings.n_theta));
        else
            %% Slowest Method
            bG=settings.S*kronsp(sparse(G),speye(settings.n_theta*settings.n_phi));

            %% Combined
            LHS3=-kronsp(speye(3),LHSL)-bG;
            LHS3=b_fd_lhs_bc_coder(LHS3,settings);
            RHS_mat3=[f0.*(e_all_field(:,:,1)-ebracket.e1);f0.*(e_all_field(:,:,2)-ebracket.e2);f0.*(e_all_field(:,:,3)-ebracket.e3)];
            RHS_mat3(b_endings,:)=0;
            RHS3=reshape(transpose(RHS_mat3),3*(settings.n_theta)*settings.n_phi,1);

            % Normalisation (Sq. Matrix)
            LHS3(settings.N_total/2,:)=[kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi) zeros(1,2*settings.N_total)];
            LHS3(settings.N_total/2*3,:)=[zeros(1,settings.N_total) kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi) zeros(1,settings.N_total)];
            LHS3(settings.N_total/2*5,:)=[zeros(1,2*settings.N_total) kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi)];
            RHS3([settings.N_total/2 settings.N_total/2*3 settings.N_total/2*5],1)=0;
            RHS3=LHS3\RHS3;
            b=transpose(reshape(RHS3,settings.n_phi,settings.n_theta*3));

            b1=b(b_endings(1):b_endings(2),:);
            b2=b(b_endings(3):b_endings(4),:);
            b3=b(b_endings(5):b_endings(6),:);
        end
    else
        Vinv=transpose(W)\eye(3);
        Bc=complex(zeros(settings.N_total,3));
        for i=1:3
            RHS_matc=W(1,i)*(e_all_field(:,:,1)-ebracket.e1)+W(2,i)*(e_all_field(:,:,2)-ebracket.e2)+W(3,i)*(e_all_field(:,:,3)-ebracket.e3);
            RHS_matc=RHS_matc.*f0;
            RHS_matc([1,settings.n_theta],:)=0;
            RHSc=reshape(transpose(RHS_matc),(settings.n_theta)*settings.n_phi,1);
            
            LHSc=-LHSL-speye(settings.N_total)*D(i,i)*settings.S;
            LHS_Fc=l_fd_lhs_bc_coder(LHSc,settings);
            LHS_Fc(settings.N_total/2,:)=kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi);
            RHSc(settings.N_total/2,1)=0;
            Bc(:,i)=LHS_Fc\RHSc;
        end
        
        b1_col=real(Vinv(1,1)*Bc(:,1)+Vinv(1,2)*Bc(:,2)+Vinv(1,3)*Bc(:,3));
        b2_col=real(Vinv(2,1)*Bc(:,1)+Vinv(2,2)*Bc(:,2)+Vinv(2,3)*Bc(:,3));
        b3_col=real(Vinv(3,1)*Bc(:,1)+Vinv(3,2)*Bc(:,2)+Vinv(3,3)*Bc(:,3));
        
        b1=transpose(reshape(b1_col,settings.n_phi,settings.n_theta));
        b2=transpose(reshape(b2_col,settings.n_phi,settings.n_theta));
        b3=transpose(reshape(b3_col,settings.n_phi,settings.n_theta));
    end
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