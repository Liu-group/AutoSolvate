sander -O -i qmmmmin.in -o qmmmmin.out -p water_solvated.prmtop -c mm.ncrst -r qmmm.ncrst  -inf qmmmmin.info -x water_solvated-qmmmmin.netcdf
sed -i '3 a guess        ./scr/ca0 ./scr/cb0' tc_job.tpl
sander -O -i qmmmmin.in -o qmmmmin.out -p water_solvated.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmmin.info -x water_solvated-qmmmmin.netcdf
sander -O -i qmmmheat.in -o qmmmheat.out -p water_solvated.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmheat.info -x water_solvated-qmmmheat.netcdf
sander -O -i qmmmnve.in -o qmmmnve.out -p water_solvated.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmnve.info -x water_solvated-qmmmnve.netcdf
sander -O -i qmmmnvt.in -o qmmmnvt.out -p water_solvated.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmnvt.info -x water_solvated-qmmmnvt.netcdf
