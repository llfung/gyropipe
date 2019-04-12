function [LHS]=l_fd_lhs_coder(settings)  %#codegen
% LHS Matrix for Operator (-L) for solving f(e)
%     persistent n_theta n_phi cot_term d1phicos_wide d1phicos d1phi_wide d1phi d1phisin_wide d1phisin MM sin2_term d2phi_wide d2phi
    n_theta=31;
    n_phi=16;
%     N_total=16*31;
    
    %% Vorticity Term
    cot_term=[0;cot(settings.theta(2:(n_theta-1)));0];
    
    d1phicos_wide=spdiags(cos(settings.phi')*[-1/60 3/20 -3/4 0 3/4 -3/20 1/60],0:6,n_phi,n_phi+6); %Central Difference 6th order in phi
    d1phicos=d1phicos_wide(:,3+(1:n_phi));
    d1phicos(1:3,end-2:end)=d1phicos(1:3,end-2:end)+d1phicos_wide(1:3,1:3);
    d1phicos(end-2:end,1:3)=d1phicos(end-2:end,1:3)+d1phicos_wide(end-2:end,end-2:end);
    LHS=-settings.omega_e1*settings.S/2/settings.dphi*...
        kronsp(spdiags(cot_term,0,n_theta,n_theta)...
        ,d1phicos);
    
    d1phi_wide=spdiags(ones(n_phi,1)*[-1/60 3/20 -3/4 0 3/4 -3/20 1/60],0:6,n_phi,n_phi+6); %Central Difference 6th order in phi
    d1phi=d1phi_wide(:,3+(1:n_phi));
    d1phi(1:3,end-2:end)=d1phi(1:3,end-2:end)+d1phi_wide(1:3,1:3);
    d1phi(end-2:end,1:3)=d1phi(end-2:end,1:3)+d1phi_wide(end-2:end,end-2:end);
    LHS=LHS-settings.omega_e2*settings.S/2/settings.dphi*...
        kronsp(speye(n_theta,n_theta)...
        ,d1phi);
    
    d1phisin_wide=spdiags(sin(settings.phi')*[-1/60 3/20 -3/4 0 3/4 -3/20 1/60],0:6,n_phi,n_phi+6); %Central Difference 6th order in phi
    d1phisin=d1phisin_wide(:,3+(1:n_phi));
    d1phisin(1:3,end-2:end)=d1phisin(1:3,end-2:end)+d1phisin_wide(1:3,1:3);
    d1phisin(end-2:end,1:3)=d1phisin(end-2:end,1:3)+d1phisin_wide(end-2:end,end-2:end);
    LHS=LHS-settings.omega_e3*settings.S/2/settings.dphi*...
        kronsp(spdiags(cot_term,0,n_theta,n_theta)...
        ,d1phisin);   
    
    LHS=LHS-settings.omega_e1*settings.S/4/settings.dtheta*...
        kronsp(spdiags(ones(n_theta,1)*[-1 1],[-1 1],n_theta,n_theta)...
        ,spdiags(sin(settings.phi'),0,n_phi,n_phi));
    LHS=LHS+settings.omega_e3*settings.S/4/settings.dtheta*...
        kronsp(spdiags(ones(n_theta,1)*[-1 1],[-1 1],n_theta,n_theta)...
        ,spdiags(cos(settings.phi'),0,n_phi,n_phi));    
    %% Gyrotactic Term
    MM=transpose(spdiags([-sin(settings.theta)/(2*settings.dtheta)...
        -2*cos(settings.theta) sin(settings.theta)/(2*settings.dtheta)],...
        -1:1,n_theta,n_theta));
    LHS=LHS-settings.beta*kronsp(MM,speye(n_phi,n_phi));
    %% Diffusive (Laplacian) Term
%     pre_ll=sparse(n_theta,n_theta,n_theta*3-2);
    MM=transpose(...
        spdiags(...
        [cot_term/...
        (2*settings.dtheta)+settings.dtheta^(-2) ...
        -2*settings.dtheta^(-2)*ones(n_theta,1) ...
        -cot_term/(2*settings.dtheta)+settings.dtheta^(-2)]...
        ,-1:1,n_theta,n_theta));
    
    sin2_term=[0;sin(settings.theta(2:(n_theta-1))).^(-2);0];
    d2phi_wide=spdiags(ones(n_phi,1)*[1/90 -3/20 3/2 -49/18 3/2 -3/20 1/90],0:6,n_phi,n_phi+6);
    d2phi=d2phi_wide(:,3+(1:n_phi));
    d2phi(1:3,end-2:end)=d2phi(1:3,end-2:end)+d2phi_wide(1:3,1:3);
    d2phi(end-2:end,1:3)=d2phi(end-2:end,1:3)+d2phi_wide(end-2:end,end-2:end);
    
    LHS=LHS+kronsp(MM,speye(n_phi,n_phi)) ...
        +kronsp(spdiags(sin2_term,0,n_theta,n_theta),...
        d2phi)*settings.dphi^(-2);

end
    