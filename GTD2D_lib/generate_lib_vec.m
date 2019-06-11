cfg=coder.config('lib');
cfg.BuildConfiguration='Faster Runs';
cfg.FilePartitionMethod='SingleFile';

maxsize=1080;
codegen -config cfg gtd2d_libinter_cfunvec -args {coder.typeof(0,[maxsize 1],[1 0]),coder.typeof(0,[maxsize 1],[1 0]),coder.typeof(0,[maxsize 1],[1 0]),coder.typeof(0,[maxsize 1],[1 0])}
