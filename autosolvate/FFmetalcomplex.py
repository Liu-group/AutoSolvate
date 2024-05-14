import os
import sys
from . import autoMCPB
from . import QM_gen
from . import boxgen_metal
import getopt
from . import gen_esp
import subprocess

support_method = ['B3LYP']
support_basisset = ['6-31G','6-31G*','LANL2DZ','6-31GS']


class genFF():
    def __init__(self, filename,metal_charge, chargefile, solvent,slv_count,folder,
                 mode,spinmult,software,basisset,method,cubesize,closeness,
                 totalcharge,nprocs,QMexe,amberhome,solvent_off,solvent_frcmod,opt):
        self.metal_charge = metal_charge
        self.filename = filename
        self.xyzfile = filename + '.xyz'
        self.chargefile = chargefile
        self.mode = mode
        self.spinmult = spinmult
        self.software = software
        self.basisset = basisset
        self.method = method
        self.totalcharge = totalcharge
        self.nprocs = nprocs
        self.QMexe = QMexe
        self.amberhome = amberhome
        self.solvent_off=solvent_off
        self.solvent_frcmod=solvent_frcmod
        self.opt = opt
        self.solvent = solvent
        self.cubesize = cubesize
        self.closeness = closeness
        self.slv_count = slv_count
        self.folder = folder

    
    def inputCheck(self):
        if os.path.exists(self.xyzfile):
            print('Find ', self.xyzfile)
        else:
            raise TypeError('Error:','can not find ',self.xyzfile,'as initial input')
            sys.exit()

        if self.mode in ['M','m']:
            print('Start to use manual mode for charge assignment')
            if os.path.exists(self.chargefile):
                print('Read the chargefile of ', self.chargefile)
            else:
                raise TypeError('Error: can not find', self.chargefile, 'as chargefile for manual mode')
                sys.exit()
        
        elif self.mode in ['A','a']:
            print('Start to use auto mode for charge assignment')
        
        else:
            raise TypeError('Error: Invalid inputs for --mode or -m')
            sys.exit()

        try:
            int(self.spinmult)
            print('The spin multiplicity is assigned as', self.spinmult)
        except ValueError:
            print('Error: Invalid inputs spin multiplicity')
            sys.exit()
        
        if self.method.upper() in support_method:
            print('Use ', self.method,'to carry out QM calculation')
        else:
            raise TypeError('Erorr',self.method,'is not supported')
        

        if self.basisset.upper() in support_basisset:
            print('Use ', self.basisset,'to carry out QM calculation')
        else:
            print('Erorr',self.basisset,'is not supported')
        
        try:
            int(self.metal_charge)
            print('The charge of metal is assigned as', self.metal_charge)
        except ValueError:
            print('Error: Invalid input for the metal charge')
            sys.exit()
        
        try:
            int(self.nprocs)
            print('The number of procs is assigned as',self.nprocs)
        except ValueError:
            print('Error: Invalid input for the number of procs')
            sys.exit()
        
        if self.totalcharge.uppper() in ['DEFAULT']:
            print('The charge of total system is set as default, will be calculated after the charge assignments for each ligand and metal')
        
        else:
            try:
                int(self.totalcharge)
                print('The totalcharge is assigned as',self.totalcharge)
            except ValueError:
                print('Error: Invalid input for the number of totalcharge')
                sys.exit()
        
        if self.opt not in  ["Y","N"]:
            raise TypeError('Error: Invalid input for opt option')
        
        if self.software not in ['orca','gau','g03','g09','gms']:
            raise TypeError('Error: Invalid software of QM, please use : orca, gau, g03, g09')

    def check_freq_result(self):
        finish = 'notconverged'
        if self.software == 'orca':
            input = self.filename + '_small_fc.orcaout'
            hess = self.filename + '_small_fc.orca.hess'
            if os.path.exists(input):
                if os.path.exists(hess):
                    with open(input,'r') as f:
                        for line in f:
                            if 'ORCA TERMINATED NORMALLY' in line:
                                finish = 'converged'
        if self.software == 'gms':
            input = self.filename + '_small_fc.log'
            if os.path.exists(input):
                    with open(input,'r') as f:
                        for line in f:
                            if 'GAMESS TERMINATED NORMALLY' in line:
                                finish = 'converged'

        if self.software in ['gau','g03','g09']:
            input = self.filename + '_small_fc.log'
            if os.path.exists(input):
                    with open(input,'r') as f:
                        for line in f:
                            if ' Normal termination of Gaussian' in line:
                                finish = 'converged'            
        
        return finish
    
    def check_opt_result(self):
        finish = 'notconverged'
        if self.software == 'gms':
            optout = self.filename + '_small_opt.log'
            with open(optout,'r') as f:
                for line in f:
                    if  'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
                  #      print('QM opt is terminated normally!')
                        finish = 'converged'
        if self.software == 'orca':
            optout  = self.filename + '_small_opt.orcaout'
            with open(optout,'r') as f:
                for line in f:
                    if 'ORCA TERMINATED NORMALLY' in line:
                  #      print('QM opt is terminated normally!')
                        finish = 'converged'
        if self.software in ['gau','g03','g09']:
            optout = self.filename + '_small_opt.log'
            with open(optout,'r') as f:
                for line in f:
                    if 'Normal termination of Gaussian' in line:
                        finish = 'converged'
        return finish

    def check_mk_result(self):
        finish = 'notconverged'
        if self.software in ['gms']:
            optout = self.filename + '_large_mk.log'
            if os.path.exists(optout):
                with open(optout,'r') as f:
                    for line in f:
                        if  'EXECUTION OF GAMESS TERMINATED NORMALLY' in line:
                        #      print('QM opt is terminated normally!')
                            finish = 'converged'
        elif self.software in ['orca']:
            optout  = self.filename + '_large_mk.orcaout'
            if  os.path.exists(optout):
                with open(optout,'r') as f:
                    for line in f:
                        if 'ORCA TERMINATED NORMALLY' in line:
                            finish = 'converged'
        
        elif self.software in ['gau','g03','g09']:
            optout  = self.filename + '_large_mk.log'
            if  os.path.exists(optout):
                with open(optout,'r') as f:
                    for line in f:
                        if 'Normal termination of Gaussian' in line:
                            finish = 'converged'
        return finish    
            
    def runQM_opt(self):
        if self.software == 'gms':
            optinput = self.filename + '_small_opt.inp'
            optout = self.filename + '_small_opt.log'
            if os.path.exists(optinput):
                pass
           #     print('Find ', optinput)
            else:
                raise TypeError('Erorr, can not find', optinput )
             #   sys.exit()
            
            if os.path.exists(optout):
                check = self.check_opt_result()
                print(check)
                if check == 'notconverged':
                    cmd = self.QMexe + ' ' + optinput +' > ' + optout
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                    stdout, stderr = proc.communicate()
                    proc.wait()
                elif check == 'converged':
                    print('Next to submit Freq calculation')
                
            else:
                cmd = self.QMexe + ' ' + optinput +' > ' + optout
                print('------------------start to submit opt job of ', optinput)
                print(cmd)
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                stdout, stderr = proc.communicate()
                proc.wait()
        if self.software == 'orca':
            optinput = self.filename + '_small_opt.orca'
            optout = self.filename + '_small_opt.orcaout'
            if os.path.exists(optinput):
                pass
           #     print('Find ', optinput)
            else:
                raise TypeError('Erorr, can not find', optinput )
              #  sys.exit()
            
            if os.path.exists(optout):
                check = self.check_opt_result()
                if check == 'notconverged':
                    cmd = self.QMexe + ' ' + optinput +' > ' + optout
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                    stdout, stderr = proc.communicate()
                    print(stdout)
                    proc.wait()
                elif check == 'converged':
                     print('Next to submit Freq calculation')

            else:
                cmd = self.QMexe + ' ' + optinput +' > ' + optout
                print('------------------start to submit opt job of ', optinput)
                print(cmd)
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                stdout, stderr = proc.communicate()
                proc.wait()
        if self.software in ['gau','g09','g03']:
            optinput = self.filename + '_small_opt.com'
            if os.path.exists(optinput):
                pass
            else:
                raise TypeError('Erorr, can not find', optinput )
               # sys.exit()
            optout = self.filename + '_small_opt.log'
            cmd = self.QMexe + ' < ' +  optinput  + ' > ' + optout
            if os.path.exists(optout):
                check = self.check_opt_result()
                if check == 'notconverged':
                    cmd = self.QMexe + ' ' + optinput +' > ' + optout
                    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                    stdout, stderr = proc.communicate()
                    print(stdout)
                    proc.wait()
                elif check == 'converged':
                     print('Next to submit Freq calculation')
            else:
                proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                stdout, stderr = proc.communicate()
                proc.wait()

    
    def runQM_freq(self):
        if self.software == 'orca':
            freqinp = self.filename + '_small_fc.orca'
            freqout = self.filename + '_small_fc.orcaout'
            optout = self.filename + '_small_opt.orcaout'
            check = self.check_freq_result()
            if check == 'notconverged':
                if os.path.exists(freqinp):
                    print('Find',freqinp)
                else:
                    print('Erorr, can not find', freqinp )
                    sys.exit()
                if os.path.exists(optout):
                    print('Start to check ',optout)
                    check = self.check_opt_result()
                    if check == 'converged':
                    #    print('Geometry opt is converged, start to calculate freq')
                        cmd = self.QMexe + ' ' + freqinp + ' > ' + freqout
                        print(cmd)
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    if check == 'notconverged':
                        print('Erorr: Geometry opt is converged! please restart QM opt calculation')
                      #  sys.exit()        
                else:
                    raise TypeError('Error, can not find',optout, 'to get optimized coordinates for Freq calculation')
            elif check == 'converged':
                print('Freq calculation is finished, start to run QM charge calculation')
        
        if self.software == 'gms':
            freqinp = self.filename + '_small_fc.inp'
            freqout = self.filename + '_small_fc.log'
            optout = self.filename + '_small_opt.log'
            check = self.check_freq_result()
            if check == 'notconverged':
                if os.path.exists(freqinp):
                    print('Find',freqinp)
                else:
                    print('Erorr, can not find', freqinp )
                    sys.exit()
                if os.path.exists(optout):
               #    print('Start to check!!!! ',optout)
                    check = self.check_opt_result()
                    if check == 'converged':
                        print('Geometry opt is converged, start to calculate freq')
                        cmd = self.QMexe + ' ' + freqinp + ' > ' + freqout
                        print(cmd)
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    if check == 'notconverged':
                        print('Erorr: Geometry opt is converged! please restart QM opt calculation')
                        sys.exit()        
                else:
                    raise TypeError('Error, can not find',optout, 'to get optimized coordinates for Freq calculation')
                  #  sys.exit()
            elif check == 'converged':
                print('Freq calculation is finished, start to run QM charge calculation')
        
        if self.software in ['g03','g09','gau']:
            freqinp = self.filename + '_small_fc.com'
            freqout = self.filename + '_small_fc.log'
            optout = self.filename + '_small_opt.log'
            check = self.check_freq_result()
            if check == 'notconverged':
                if os.path.exists(freqinp):
                    print('Find',freqinp)
                else:
                    print('Erorr, can not find', freqinp )
                    sys.exit()
                if os.path.exists(optout):
               #    print('Start to check!!!! ',optout)
                    check = self.check_opt_result()
                    if check == 'converged':
                        print('Geometry opt is converged, start to calculate freq')
                        cmd = self.QMexe + ' ' + freqinp + ' > ' + freqout
                        print(cmd)
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)   
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    if check == 'notconverged':
                        print('Erorr: Geometry opt is converged! please restart QM opt calculation')
                        sys.exit()        
                else:
                    raise TypeError('Error, can not find',optout, 'to get optimized coordinates for Freq calculation')
                    sys.exit()
            elif check == 'converged':
                print('Freq calculation is finished, start to run QM charge calculation')

    def runQM_mk(self):
        print('runingg QM mk')
        if self.software in[ 'gms' ]:
            mkinp = self.filename + '_large_mk.inp'
            mkout = self.filename + '_large_mk.log'
            check_mk = self.check_mk_result()
            if check_mk == 'notconverged':
                if self.opt == 'Y':
                    check = self.check_opt_result()
                    if check == 'converged':
                        cmd = self.QMexe + ' ' + mkinp + ' > ' + mkout
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    if check == 'notconverged':
                        print('Erorr: Geometry opt is converged! please restart QM opt calculation')
                        sys.exit()
                elif self.opt == 'N':
                    if os.path.exists(mkinp):
                        cmd = self.QMexe + ' ' + mkinp + ' > ' + mkout
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    else:
                        raise TypeError('Erorr: cant find',mkinp)
            elif check_mk == 'converged':
                raise TypeError('charge calculation is converged!')
        
        elif self.software in ['orca','ORCA']:
            mkinp = self.filename + '_large_mk.orca'
            mkout = self.filename + '_large_mk.orcaout'
            check_mk = self.check_mk_result()
            if check_mk == 'notconverged':
                if self.opt == 'Y':
                    check = self.check_opt_result()
                    if check == 'converged':
                        cmd = self.QMexe + ' ' + mkinp + ' > ' + mkout
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    if check == 'notconverged':
                        print('Erorr: Geometry opt is converged! please restart QM opt calculation')
                        sys.exit()
                elif self.opt == 'N':
                    if os.path.exists(mkinp):
                        cmd = self.QMexe + ' ' + mkinp + ' > ' + mkout
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    else:
                        raise TypeError('Erorr: cant find',mkinp)
        
        elif self.software in ['g03','g09','gau']:
            mkinp = self.filename + '_large_mk.com'
            mkout = self.filename + '_large_mk.log'
            check_mk = self.check_mk_result()

            if check_mk == 'notconverged':
          #      print('yyyyyyyyyyyy')
                if self.opt == 'Y':
                    check = self.check_opt_result()
                    if check == 'converged':
                        cmd = self.QMexe + ' < ' + mkinp + ' > ' + mkout
                    #    print(cmd)
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = proc.communicate()
                        print(stderr)
                        proc.wait()
                    if check == 'notconverged':
                        raise TypeError('Erorr: Geometry opt is converged! please restart QM opt calculation')
                        sys.exit()
                elif self.opt == 'N':
                    if os.path.exists(mkinp):
                        cmd = self.QMexe + ' < ' + mkinp + ' > ' + mkout
                        proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = proc.communicate()
                        proc.wait()
                    else:
                        print('Erorr: cant find',mkinp)
                
    def build(self):
        print (self.folder)
        os.makedirs(self.folder, exist_ok=True)
        os.system('cp ' + self.xyzfile + ' ' + self.folder)
        cwd = os.getcwd()
        print(cwd)
        os.chdir(self.folder)

        print('******************** start to generate inputs for MCPB.py -s 1 ********************')
        step1 = autoMCPB.AutoMCPB(filename=self.filename,metal_charge=self.metal_charge, spinmult=self.spinmult,amberhome=self.amberhome,
                       mode=self.mode,chargefile=self.chargefile,round='1',software=self.software)
        step1.build()
        print('******************** Finish generating inputs for MCPB.py -s 1 ********************')
        print('******************** start to QM calculations for', self.software + '_small_opt', self.software + '_small_fc',self.software + '_large_mk','********************')
        step2 = QM_gen.QM_inputs_gen(filename=self.filename,software=self.software,caltype='opt',multi=self.spinmult,method=self.method,
                            basisset=self.basisset,metal='default',totalcharge=self.totalcharge,nprocs=self.nprocs,opt=self.opt)
        step2.build()
        self.runQM_opt()
        step3 = QM_gen.QM_inputs_gen(filename=self.filename,software=self.software,caltype='fc',multi=self.spinmult,method=self.method,
                            basisset=self.basisset,metal='default',totalcharge=self.totalcharge,nprocs=self.nprocs,opt=self.opt)
        step3.build()
        self.runQM_freq()

        step4 = QM_gen.QM_inputs_gen(filename=self.filename,software=self.software,caltype='mk',multi=self.spinmult,method=self.method,
                            basisset=self.basisset,metal='default',totalcharge=self.totalcharge,nprocs=self.nprocs,opt=self.opt)
        step4.build()
        #print('#######################')
        self.runQM_mk()
        checkfreq = self.check_freq_result()
        print('******************** start to generate inputs for MCPB.py -s 2 ********************')

        if self.software == 'gau':
            if os.path.exists(self.filename + '_small_opt.chk'):
                cmd = self.QMexe[:-3] + 'formchk ' + self.filename + '_small_opt.chk' + ' ' + self.filename + '_small_opt.fchk'
                subprocess.call(cmd,shell=True)
            else:
                raise TypeError('can not find ' + self.filename + '_small_opt.chk')

        if checkfreq == 'converged':
            step5=autoMCPB.AutoMCPB(filename=self.filename,metal_charge=self.metal_charge, spinmult=self.spinmult,amberhome=self.amberhome,
                        mode=self.mode,chargefile=self.chargefile,round='2',software=self.software)
            step5.build()
        #    print('nnnnnnnnnnnnnn')
        elif checkfreq == 'notconverged':
            raise TypeError('Erorr: Freq calculation is not finished!')
            sys.exit()
        if os.path.exists(self.filename + '_mcpbpy.frcmod'):
            print(self.filename + '_mcpbpy.frcmod', 'is generated by MCPB.py -s')
        else:
            raise TypeError('Error:',self.filename + '_mcpbpy.frcmod is not generated')
            sys.exit()
        print('******************** Finish generating frcmod files for metal bonds ********************')
        print('********************    start to generate inputs for MCPB.py -s 3   ********************')
        checkmk = self.check_mk_result()
        if checkmk == 'converged':
            if self.software != 'orca':
                step6 = autoMCPB.AutoMCPB(filename=self.filename,metal_charge=self.metal_charge, spinmult=self.spinmult,
                mode=self.mode,chargefile=self.chargefile,round='3',software=self.software,amberhome=self.amberhome)
                step6.build()
            else:
                if self.opt == 'Y':
                    orcaout = self.filename + '_small_opt.orcaout'
                    trjxyz = self.filename +'_small_opt.orca_trj.xyz'
                    crds = QM_gen.get_crds_from_orcaopt(orcaout=orcaout,trjxyz=trjxyz)
                    elements = []
                    newcrds = []
                    for line in crds:
                        newcrds.append([float(line[1]),float(line[2]),float(line[3])])
                        elements.append(line[0].capitalize())
                elif self.opt == 'N':
                    originalxyz = self.filename + '.xyz'
                    newcrds, elements = gen_esp.read_xyz(originalxyz)
                
                orcapath  = self.QMexe + '_vpot'
                gbw = self.filename +'_large_mk.orca.gbw'
                density = self.filename +'_large_mk.orca.scfp'
                out = self.filename +'_large_mk.orcaespout'
                espf = self.filename +'_large_mk.orca.esp'
                gen_esp.gen_grids(crds=newcrds,elements=elements,orcapath=orcapath,gbw=gbw,denisty=density,out=out)
                gen_esp.convertoesp(espin=out,espout=espf,crds=newcrds)
                ### using MCPB_orca to generate resp1.in resp2.in ####
                cmd = 'MCPB_orca ' + self.filename + '_MCPB.in 3'
                with open('respinputgen.log', 'w') as output_file:
                    subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
                
                """
                cmd1 = self.amberhome + "resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg \
                        -e %s -s resp1_calc.esp" %espf
                cmd2 = self.amberhome + "resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg \
                        -t resp2.chg -e %s -s resp2_calc.esp" %espf
                with open('resp1.log',"w") as output_file:
                    subprocess.call(cmd1,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
                with open('resp2.log',"w") as output_file:
                    subprocess.call(cmd2,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
                """ 

        elif checkmk == 'notconverged':
            print('Erorr: charge calculation is not finished!')
            sys.exit()
        if os.path.exists('resp2.chg'):
            step7 = autoMCPB.AutoMCPB(filename=self.filename,metal_charge=self.metal_charge, spinmult=self.spinmult,
                        mode=self.mode,chargefile=self.chargefile,round='4',software=self.software,amberhome=self.amberhome)
            step7.build()
        else:
            print('Erorr: can not find resp2.chg')
            sys.exit()
        if os.path.exists(self.filename + '_dry.prmtop'):
            step8 = autoMCPB.AutoMCPB(filename=self.filename,metal_charge=self.metal_charge, spinmult=self.spinmult,
                        mode=self.mode,chargefile=self.chargefile,round='5',software=self.software,amberhome=self.amberhome)
            step8.build()
        else:
            print('Erorr: can not find ' + self.filename + '_dry.prmtop')
        
        step9 = boxgen_metal.solventBoxBuilderMetal(pdb_prefix=self.filename,totalcharge=self.totalcharge,solvent=self.solvent,
                                              solvent_frcmod=self.solvent_frcmod, solvent_off=self.solvent_off,closeness=self.closeness,
                                             outputFile='',amberhome=self.amberhome,cube_size=self.cubesize,slv_count=self.slv_count)
        step9.build()
        os.chdir(cwd)
"""
        pdb_prefix = self.filename
        solvent = self.solvent
        output = ''
        charge = self.totalcharge
        cube_size = self.cubesize
        amberhome = self.amberhome
        closeness = self.closeness
        solvent_off = self.solvent_off
        solvent_frcmod = self.solvent_frcmod
        cmd = [
        "boxgen_metal",
        "-m", str(pdb_prefix),
        "-s", str(solvent),
        "-o", str(output),
        "-c", str(charge),
        "-b", str(cube_size),
        "-a", str(amberhome),
        "-t", str(closeness),
        "-l", str(solvent_off),
        "-p", str(solvent_frcmod)
               ]
        try:
           result = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            cmd = [
                "python -m autosolvate.boxgen_metal",
                "-m", str(pdb_prefix),
                "-s", str(solvent),
                "-o", str(output),
                "-c", str(charge),
                "-b", str(cube_size),
                "-a", str(amberhome),
                "-t", str(closeness),
                "-l", str(solvent_off),
                "-p", str(solvent_frcmod)
                    ]
            result = subprocess.run(cmd, capture_output=True, text=True)
        if result.stderr:
            print( result.stderr)
"""

def startFFgen(argumentList):
    options = 'hm:c:u:v:f:x:k:r:G:n:l:p:s:A:Q:e:b:t:z:'
    long_options = ["help",'filename=','metal_charge=','spin=','mode=', 'folder=',
                    'chargefile=','totalcharge=','software=','nprocs=','method=','cubesize=','closeness='
                    'qmexe=','basisset=','solventoff=','solventfrcmod=','amberhome=','opt=','solvent=']
    arguments, values = getopt.getopt(argumentList,options,long_options)
    filename = None
    metal_charge = None
    mode='A'
    spinmult=None
    chargefile = None
    software = 'gms'
    spinmult = '1'
    totalcharge= 'default'
    nprocs='8'
    method='B3LYP'
    QMexe = None
    solvent_off=''
    solvent_frcmod = ''
    basisset='6-31G*'
    amberhome='$AMBERHOME/bin/'
    opt = 'Y'
    gamessexe = None
    solvent = 'water'
    closeness = 'automated'
    cubesize = 54
    slv_count = 210*8
    folder = './'
    

    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            message = '''Usage: autoMCPB.py [OPTIONS]
            Options:
                -h, --help              Display this help message.
                -m, --main              metalcomplex xyz file including extension.
                -c, --metal_charge      Set the metal charge.
                -u, --spin              Set the spin default: 1
                -v, --mode              Set the mode (A/M).
                -f, --chargefile        Specify the charge file if --mode is M.
                                        chargefile example:
                                        LG1 -1 #ligand name charge 
                -x, --software          gau,g09,g16,gms,orca, default:gms
                -k, --totalcharge       total charge of the whole system, the default is caluated after charge assignment
                -r, --nprocs             procs to run orca QM calculation, if -x orca
                -G  --qmexe             path to QM exe e.g. /opt/orca/5.0.2/orca
                -n  --method            method of QM default:B3LYP
                -l  --solventoff        path to the custom solvent .off library file
                -p  --solventfrcmod     path to the custom solvent .frcmod file
                -s  --basisset          basisset used in QM calculation, default: 6-31G*
                -A  --amberhome         path of AmberTools bin default: $AMBERHOME/bin/
                -Q  --opt               use(Y) or not use(N) the QM optimized structure for charge calculation default: Y
                -e  --solvent           name of solvent (water, methanol, chloroform, nma, acetonitrile)
                -b  --cubesize          size of solvent cube in angstroms default: 56 
                -t  --closeness         Solute-solvent closeness setting default value is used if the option is not specified
                -z  --folder            the path of outputfiles, default: the current folder './'


                '''
            print(message)
            sys.exit()
        elif currentArgument in ("-m","--main"):
            filename = os.path.splitext(str(currentValue))[0]
        elif currentArgument in ("-c",'--metal_charge'):
            metal_charge = int(currentValue)
        elif currentArgument in ('-u','--spin'):
            spinmult=int(currentValue)
        elif currentArgument in ('-v','--mode'):
            mode = str(currentValue)
        elif currentArgument in ('-f','--chargefile'):
            chargefile = str(currentValue)
        elif currentArgument in ('-x','--software'):
            software = str(currentValue)
        elif currentArgument in ('-k','--totalcharge'):
            totalcharge = str(currentValue)
        elif currentArgument in ('-r','--nprocs'):
            nprocs = str(currentValue)
        elif currentArgument in ('-G','--qmexe'):
            QMexe = str(currentValue)
        elif currentArgument in ('-n','--method'):
            method = str(currentValue)
        elif currentArgument in ('-l','--solventoff'):
            solvent_off = str(currentValue)
        elif currentArgument in ('-p','--solventfrcmod'):
            solvent_frcmod  = str(currentValue)
        elif currentArgument in ('-s','--basisset'):
            basisset  = str(currentValue)
        elif currentArgument in ('-A','--amberhome'):
            amberhome = str(currentValue)
        elif currentArgument in ('-Q','-opt'):
            opt = str(currentValue)
        elif currentArgument in ('-e','-solvent'):
            solvent = str(currentValue)
        elif currentArgument in ('-b','-cubesize'):
            cubesize = str(currentValue)
        elif currentArgument in ('-t','-closeness'):
            closeness = str(currentValue)
        elif currentArgument in ('-z','-folder'):
            folder = str(currentValue)

    builder = genFF(filename=filename,metal_charge=metal_charge, spinmult=spinmult,basisset=basisset,closeness=closeness,folder = folder,
                    mode=mode,chargefile=chargefile,software=software,solvent_frcmod=solvent_frcmod,amberhome=amberhome,cubesize=cubesize,
                    totalcharge=totalcharge,nprocs=nprocs,QMexe=QMexe,method=method,solvent_off=solvent_off,opt=opt,solvent=solvent,slv_count=slv_count)
 #   print(arguments)

    builder.build()
    
if __name__ == '__main__':
    argumentList = sys.argv[1:]
   # print(argumentList)
    startFFgen(argumentList)
                    


    
    


    



    