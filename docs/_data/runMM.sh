sander -O -i mmmin.in -o mmmin.out -p water_solvated.prmtop -c water_solvated.inpcrd -r mm.ncrst -inf mmmin.info
pmemd.cuda -O -i mmheat.in -o mmheat.out -p water_solvated.prmtop -c mm.ncrst -r mm.ncrst -x water_solvated-heat.netcdf -inf mmheat.info
pmemd.cuda -O -i mmnve.in -o mmnve.out -p water_solvated.prmtop -c mm.ncrst -r mm.ncrst -x water_solvated-mmnve.netcdf -inf mmnve.info
pmemd.cuda -O -i mmnpt.in -o mmnpt.out -p water_solvated.prmtop -c mm.ncrst -r mm.ncrst -x water_solvated-mmnpt.netcdf -inf mmnpt.info
