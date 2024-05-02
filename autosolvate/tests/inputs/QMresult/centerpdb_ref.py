import sys
import subprocess
number = sys.argv[1]
ifile = open(number + 'ref.in','w')
ifile.write('parm '+ number + '_nptlast.pdb\n')
ifile.write('trajin ' + number + '_nptlast.pdb\n')
ifile.write('center :1 mass\n')
ifile.write('image center familiar\n')
ifile.write('trajout  ' + number + '_nptlast-center.pdb pdb\n')
ifile.write('go\n')
ifile.close()


