function LHS=l_fd_lhs_bc_coder(LHS,settings)
% Applying Boundary Condition for LHS Matrix for Operator (-L) for solving f(e)
%   By assuming df\dtheta=0 at theta=0,pi, in 6th order,
%   using average over phi 
% persistent BC1 BC2
    n_theta=31;
    n_phi=16;
    N_total=16*31;


    BC1=kronsp(sparse([-49/20 zeros(1,n_theta-1);...
        zeros(n_theta-2,n_theta);...
        zeros(1,n_theta-1) 49/20]),speye(n_phi))*(1/settings.dtheta);
    BC2=kronsp(sparse([0 6 -15/2 20/3 -15/4 6/5 -1/6 zeros(1,n_theta-7);...
        zeros(n_theta-2,n_theta);...
        zeros(1,n_theta-7) 1/6 -6/5 15/4 -20/3 15/2 -6 0]),sparse(ones(n_phi,n_phi)))*(1/settings.dtheta/n_phi);

    LHS(1:n_phi,:)=0;
    LHS(n_phi*(n_theta-1)+1:N_total,:)=0;
    
    LHS=LHS+BC1+BC2;
end