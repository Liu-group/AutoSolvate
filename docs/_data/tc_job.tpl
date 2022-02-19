basis        lacvps_ecp
method       ub3lyp
dispersion   yes
guess        ./scr/ca0 ./scr/cb0
scf          diis+a
threall      1e-13
convthre     2e-5
xtol         1-4
dftd         no
maxit        200
dftgrid      1
charge       6
spinmult     6
scrdir       ./scr
keep_scr     yes
end
