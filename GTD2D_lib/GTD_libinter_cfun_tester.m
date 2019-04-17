
Grtz=[rand(1,1)*2-1 0 rand(1,1)*2-1 0 0 0 rand(1,1)*2-1 0];

Grz=[Grtz(1) Grtz(3) Grtz(7)];

Diff_rz=gtd2d_libinter_cfun(Grz);

