# gyropipe
Gyrotactic Suspension in a vertical Pipe DNS, based on OpenPipeFlow
## Compiling `gtd_eig_cfun.a`
0. Make sure you use > MATLAB R2019a
1. Open MATLAB and `cd` to `$(project_folder)/GTD_StaticLib`
2. Then execute the following:
```
generate_lib
```
3. Copy the `gtd_eig_cfun.a` from `$(project_folder)/GTD_StaticLib/codegen/lib/gtd_eig_cfun` to `$(project_folder)`
4. Now you can `make runall` as normal.