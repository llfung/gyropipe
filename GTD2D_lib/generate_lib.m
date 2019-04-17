cfg=coder.config('lib');
cfg.BuildConfiguration='Faster Runs';
cfg.FilePartitionMethod='SingleFile';

GTD_libinter_cfun_tester;
inp=coder.typeof(Grz);

codegen -config cfg gtd2d_libinter_cfun -args {inp}
