% Grtz=rand(1,8)*2-1;
%  Grtz=[0.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0];
% Grtz=[1.0 -1.0 -0.5 4.0 -3.0 -1.5 2.0 1.0];
% for i=1:32
% [Diff_rtz_d]=gtd_direct_cfun(Grtz);
%  [Diff_rtz_e]=gtd_eig_cfun(Grtz);
% end

Grtz=[1.0 -1.0 -0.5 4.0 -3.0 -1.5 2.0 1.0];
for i=1:32*18
[Diff_rtz_e]=gtd_eig_cfun(Grtz);
end

% Grtz=[0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0];
%  [Diff_rtz_e]=gtd_eig_cfun(Grtz);
% max((Diff_rtz_e-Diff_rtz_d)./norm(Diff_rtz_d))
