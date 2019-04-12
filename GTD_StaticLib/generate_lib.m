cfg=coder.config('lib');
cfg.BuildConfiguration='Faster Runs';
cfg.FilePartitionMethod='SingleFile';

GTD_direct_cfun_tester;
inp=coder.typeof(Grtz);

codegen -config cfg gtd_eig_cfun -args {inp}
