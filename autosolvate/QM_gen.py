    
import subprocess
import numpy as np
from glob import glob
import sys
import getopt
import collections
from . import GAMESS_input_gen

metalbasisset = ['LANL2DZ','SBKJC']
vdw = {"H":1.1,"Li":2.2,"Be":1.9,"B":1.8,"C":1.7,"N":1.55,"O":1.52,"F":1.47,"Na":2.4,"Mg":2.2,"Al":2.1,"Si":2.1,
        "P":1.8,"S":1.8,"Cl":1.75,"K":2.8,"Ca":2.31,"Sc":2.3,"Ti":2.15,"V":2.05,
        "Cr":2.05,"Mn":2.05,"Fe":2.05,"Co":2.0,"Ni":2.0,"Cu":2.0,"Zn":2.1,"Ga":2.1,
        "As":2.05,"Se":1.9,"Br":1.83,"Rb":2.9,"Sr":2.55,"Y":2.4,"Zr":2.3,"Nb":2.15,
        "Mo":2.1,"Tc":2.05,"Ru":2.05,"Rh":2.0,"Pd":2.05,"Ag":2.1,"Cd":2.2,"In":2.2,"I":1.98,"W":2.1,"Os":2.1,"Ir":2.1}

def get_crds_from_gmsopt(inp):
    cor_data = []
    start = False
    start_line = "EQUILIBRIUM GEOMETRY LOCATED"
    end_line = "INTERNUCLEAR DISTANCES (ANGS.)"
    cor_data = []
    with open(inp, 'r') as f:
        for line in f:
            if start_line in line:
                start = True
                continue  
            if start:
                if end_line in line:
                    break  
                if "ATOM" in line or "------" in line:
                    continue  
                if len(line.split()) == 5: 
                    cor_data.append(line)
    return cor_data

def get_crds_from_gmsinp(inp):
    crds = []
    with open(inp,'r') as f:
        b = 0
        e = 0
        data = f.readlines()
        for sort, line in enumerate(data):
            line = line.strip()
            if '$DATA' in line:
                b = sort + 3
            elif len(line.split()) == 1:
                if '$END' in line:
                    e = sort
        gms_cord = data[b:e]
        for line in gms_cord:
            crds.append(line)
    return crds


def get_crds_from_xyz(inp):
    crds = []
    with open(inp,'r') as f:
        data = f.readlines()
        for line in data[2:]:
            line = line.strip()
            atom = line.split()[0].capitalize()
            x = line.split()[1]
            y = line.split()[2]
            z = line.split()[3]
            crds.append([atom,x,y,z])
    return crds
         
def get_crds_from_orcaopt(orcaout,trjxyz):
    check = False
    atomnumber = None
    crds = []
    #### check opt finish #####
    try:
        with open(orcaout) as f:
            for line in f:
                if " ****ORCA TERMINATED NORMALLY****" in line:
                    check = True
                    print("ORCA opt terminated normally")

    except FileNotFoundError:
        print("Erorr:",orcaout, "file not found")

    if check:
        try:
            with open(trjxyz, 'r') as f:
                firstline = f.readline()
                atomnumber = list(firstline.strip())
                if len(atomnumber[0].split()) == 1:
                    atomnumber = int(firstline)
                    queue = collections.deque(f,maxlen=atomnumber)
                    queue = list(queue)
                    if len(queue) != atomnumber:
                        print('Erorr: Reading orca trj file has an error, please check the space and newline character.')
                        sys.exit()
                    else:
                        for crd in queue:
                            crd = crd.strip()
                            if len(crd.split()) != 4:
                                print('Erorr: wrong format of orca trj xyzfile')
                                sys.exit()
                            else:
                                crds.append(crd.split()[0:4])
    
        except FileNotFoundError:
            print("Erorr:",trjxyz, "file not found")
        
        except ValueError:
            print("Erorr: ","wrong format of orca trj xyzfile")
        
    else:
        print("Erorr: ORCA opt is not terminated normally")
   # print(crds)    
    return crds



class QM_inputs_gen():
    def __init__(self, filename,software,caltype,multi,method,basisset,totalcharge,metal,nprocs,opt):
        self.software = software
        self.caltype = caltype
        self.multi = multi
        self.method = method
        self.basisset = basisset
        self.filename = filename
        self.info = filename + '.info'
        self.totalcharge = totalcharge
        self.metal = metal
        self.nprocs = nprocs
        self.opt = opt

        if self.totalcharge in ['default','Default'] or self.metal in ['default','Default']:
            try:
                with open(self.info, 'r') as f:
                    data = f.readlines()

                begin, end = 0, len(data) 
                for sort, line in enumerate(data):
                    if '@charges' in line:
                        begin = sort + 1
                    elif '@connection_infos' in line:
                        end = sort
                newtotalcharge = 0  

                for line in data[begin:end]:
                    if 'LG' not in line:
                        parts = line.split()
                        if len(parts) == 2:
                            new_metal = parts[0].capitalize()
                            self.metal = new_metal
                            metalcharge = int(parts[1])
                            newtotalcharge += metalcharge
                    else:
                        parts = line.split()
                        if len(parts) == 2:
                            ligand = parts[0].capitalize()
                            ligandchg = int(parts[1])
                            newtotalcharge += ligandchg
                self.totalcharge = newtotalcharge

            except FileNotFoundError:
                print('Error, cant find',self.info, ' or the info file is not correct, please assign the metal name and total charge manually')
                sys.exit()
            except IOError:
                print('Error, cant find',self.info, ' or the info file is not correct, please assign the metal name and total charge manually')
                sys.exit()
            except IndexError:
                print('Error, cant find',self.info, ' or the info file is not correct, please assign the metal name and total charge manually')
                sys.exit()
            except ValueError:
                print('Error, cant find',self.info, ' or the info file is not correct, please assign the metal name and total charge manually')
                sys.exit()
    
    def write_gaussian_opt(self,crds,out,multi,totalcharge,nprocs,basisset,method):
        chkfile = self.filename + '_small_opt.chk'
        ofile = open(out,'w')
        ofile.write('%Chk=' + chkfile+'\n')
        ofile.write('%Mem=3000MB' + '\n')
        ofile.write('%NProcShared=' + str(nprocs) + '\n')
        ofile.write('# ' + method + '/'  + basisset + ' Geom=PrintInputOrient Integral=(Grid=UltraFine) Opt\n\n')
        ofile.write('CLR\n\n')
        ofile.write(str(totalcharge) + ' ' + str(multi)+ '\n')
        for line in crds:
            print( "%-6s   %8.3f %8.3f %8.3f" %(line[0], float(line[1]), float(line[2]),float(line[3])), file=ofile)
        ofile.write('\n\n\n\n')
        ofile.close()

    def write_gaussian_mk_from_opt(self,out,nprocs,basisset,method):
        if self.metal.capitalize() in vdw:
            vdwradii = vdw[self.metal.capitalize()]
        else:
            raise TypeError('lack the vdw-radii of ' + self.metal.capitalize())
            sys.exit()
        chkfile = self.filename + '_small_opt.chk'
        ofile = open(out,'w')
        ofile.write('%Chk=' + chkfile+'\n')
        ofile.write('%Mem=3000MB' + '\n')
        ofile.write('%NProcShared=' + str(nprocs) + '\n')
        ofile.write('# ' + method + '/'  + basisset + ' Integral=(Grid=UltraFine) Pop(MK,ReadRadii) Geom=AllCheckpoint\nIOp(6/33=2,6/42=6)\n')
        ofile.write('\n')
        ofile.write(self.metal.capitalize() + ' ' +  str(vdwradii) + '\n\n\n\n')
        ofile.close()
    
    def write_gaussian_mk(self,crds,out,multi,totalcharge,nprocs,basisset,method):
        if self.metal.capitalize() in vdw:
            vdwradii = vdw[self.metal.capitalize()]
        else:
            raise TypeError('lack the vdw-radii of ' + self.metal.capitalize())
            sys.exit()
        chkfile = self.filename + '_small_opt.chk'
        ofile = open(out,'w')
        ofile.write('%Chk=' + chkfile+'\n')
        ofile.write('%Mem=3000MB' + '\n')
        ofile.write('%NProcShared=' + str(nprocs) + '\n')
        ofile.write('# ' + method + '/'  + basisset + ' Integral=(Grid=UltraFine) Pop(MK,ReadRadii)  Geom=AllCheckpoint\nIOp(6/33=2,6/42=6)\n\n')
        ofile.write('CLR\n\n')
        ofile.write(str(totalcharge) + ' ' + str(multi)+ '\n')
        for line in crds:
            print("%-6s   %8.3f %8.3f %8.3f" %(line[0], float(line[1]), float(line[2]), float(line[3])),file=ofile)
        ofile.write('\n')
        ofile.write(self.metal.capitalize() + ' ' +   str(vdwradii))
        ofile.write('\n\n\n\n')
        ofile.close() 

    def write_gaussian_fc(self,out,nprocs,basisset,method):
        chkfile = self.filename + '_small_opt.chk'
        ofile = open(out,'w')
        ofile.write('%Chk=' + chkfile+'\n')
        ofile.write('%Mem=3000MB' + '\n')
        ofile.write('%NProcShared=' + str(nprocs) + '\n')
        ofile.write('# ' + method + '/'  + basisset + ' Freq Geom=AllCheckpoint Guess=Read\nIntegral=(Grid=UltraFine) IOp(7/33=1)\n\n\n')
        ofile.close()

    def write_orca_opt(self,crds,out,multi,totalcharge,nprocs,basisset,method):
        ofile = open(out,'w')
        ofile.write('! ' + method + ' ' + basisset + ' OPT\n')
        ofile.write('\n')
        ofile.write("%pal\nnprocs " + str(nprocs) + "\nend\n\n")
        ofile.write('* xyz ' + str(totalcharge) + ' ' + str(multi) + '\n')
        for line in crds:
            atomname = line[0]
            x  = line[1]
            y = line[2]
            z = line[3]
            ofile.write(atomname + '  ' + x + '  ' + y  + '  ' + z + '\n')
        ofile.write('*')
        ofile.close()
    
    def write_orca_fc(self,crds,out,multi,totalcharge,nprocs,basisset,method):
        ofile = open(out,'w')
        ofile.write('! ' + method + ' ' + basisset + ' Freq\n\n')
        ofile.write("%pal\nnprocs " + str(nprocs) + "\nend\n\n")
        ofile.write('* xyz ' + str(totalcharge) + ' ' + str(multi) + '\n')
        for line in crds:
            atomname = line[0]
            x  = line[1]
            y = line[2]
            z = line[3]
            ofile.write(atomname + '  ' + x + '  ' + y  + '  ' + z + '\n')
        ofile.write('*')
        ofile.close()
    
    def write_orca_mk(self,crds,out,multi,totalcharge,basisset,nprocs, method):
        ofile = open(out,'w')
        ofile.write('! ' + method + ' ' + basisset + ' keepdens\n\n')
        ofile.write("%pal\nnprocs " + str(nprocs) + "\nend\n\n")
        ofile.write('* xyz ' + str(totalcharge) + ' ' + str(multi) + '\n')
        for line in crds:
            atomname = line[0]
            x  = line[1]
            y = line[2]
            z = line[3]
            ofile.write(atomname + '  ' + x + '  ' + y  + '  ' + z + '\n')
        ofile.write('*')
        ofile.close()  
  
    def write_orca_mk_old(self,inp,crds,out,multi,totalcharge,basiset,method): ##### not going to use this function
        ####get ### elements
        elements = []
        with open(inp,'r') as f:
            b = 0
            e = 0
            data = f.readlines()
            for sort, line in enumerate(data):
                line = line.strip()
                if '$DATA' in line:
                    b = sort + 3
                elif len(line.split()) == 1:
                    if '$END' in line:
                        e = sort
            gms_cord = data[b:e]
            for line in gms_cord:
              #  print(line)
                ele = line.split()[1]
                elements.append(ele)

        if basiset.upper() in ['6-31G','6-31G*','6-31GS']:
            ofile = open(out,'w')
            ofile.write(' $SYSTEM MEMDDI=400 MWORDS=200 $END\n')
            line = ' $CONTRL DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=1'  + ' $END\n'
            if int(multi) > 1:
                line = ' $CONTRL SCFTYP=UHF DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       MAXIT=200 $END\n'
            ofile.write(line)
            ofile.write(' $ELPOT IEPOT=1 WHERE=PDC $END\n')
            ofile.write(' $PDC PTSEL=CONNOLLY CONSTR=NONE $END\n')
            if basiset in ['6-31G*','6-31g*','6-31gd','6-31GS']:
                ofile.write( ' $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n')
            elif basiset in ['6-31G','6-31g']:
                ofile.write( ' $BASIS GBASIS=N31 NGAUSS=6 $END\n')
            ofile.write(' $DATA\n')
            ofile.write('Cluster/' + basiset +'\n')
            ofile.write('C1\n')
            for i, line in enumerate(crds):
                atom = line[0]
                number = elements[i].capitalize()
                x = line[1]
                y = line[2]
                z = line[3]
                ofile.write(atom.ljust(6) +number.rjust(4)+ x.rjust(12)+ y.rjust(12)+z.rjust(12) + '\n')
            ofile.write(' $END')
            ofile.close()
        
        elif basiset.upper() in metalbasisset:
            ECPinput = GAMESS_input_gen.GAMESS_ECP_input_gen(filename=self.filename,software='orca',caltype='mk',
                                                             multi=self.multi,method=self.method,
                                                             basisset=self.basisset,totalcharge=self.totalcharge,
                                                             metal=self.metal,nprocs=self.nprocs)
            ECPinput.write_orca_mk(crds=crds,multi=multi,
                                   totalcharge=totalcharge,basiset=basiset,method=method,out=out)
    
    def write_gms_opt(self,inp,out,multi,totalcharge,basiset,method):  ## just for 6-31g*
        if basiset in ['6-31gs','6-31GS','6-31G*','6-31g*']:
            ofile = open(out,'w')
            with open(inp,'r') as f:
                data = f.readlines()
                data_b = 0
                data_e = 0
                buer =False
                for sort, line in enumerate(data):
                    if 'DATA' in line:
                        if len(line.split()) == 1:
                            buer = True
                            data_b = sort + 1
                         #   print(data_b)
                    if buer == True:
                        if 'END' in line:
                            if len(line.split()) == 1:
                                data_e = sort
                                break
                cor_data = data[data_b:data_e]

                elements = []
                for atom in cor_data[2:]:
                    elements.append(atom.split()[0])
                
                ofile.write( " $SYSTEM MEMDDI=400 MWORDS=200 $END\n")
                if int(multi) > 1:
                    ofile.write(" $CONTRL SCFTYP=UHF DFTTYP=" + method + " RUNTYP=OPTIMIZE\n      ICHARG=" + str(totalcharge) + " MULT=" + str(multi)+" MAXIT=200 $END\n")
                else:
                    ofile.write(" $CONTRL DFTTYP=" + method + " RUNTYP=OPTIMIZE\n        ICHARG=" + str(totalcharge) + " MULT=" + str(multi)+" MAXIT=200 $END\n")
                ofile.write(' $STATPT NSTEP=1000 $END\n')
            
                ofile.write(' $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n')
                ofile.write(' $DATA\n')
                for line in cor_data[:2]:
                    ofile.write(line)
                for line in cor_data[2:]:
                    ofile.write(' ' + line)
                ofile.write(' $END\n')
        elif basiset.upper() in metalbasisset:
            print('begin to generat ECP')
            ECPinput = GAMESS_input_gen.GAMESS_ECP_input_gen(filename=self.filename,software='gms',caltype='opt',
                                                    multi=self.multi,method=self.method,
                                                    basisset=self.basisset,totalcharge=self.totalcharge,
                                                    metal=self.metal,nprocs=self.nprocs)
            ECPinput.write_gms_opt(inp,out,multi,totalcharge,basiset,method)

    
    def write_gms_fc(self, inp, out, multi, totalcharge, basiset, method):
        if basiset in ['6-31G*','6-31gs','6-31g*','6-31GS']:
            start = False
            start_line = "EQUILIBRIUM GEOMETRY LOCATED"
            end_line = "INTERNUCLEAR DISTANCES (ANGS.)"
            cor_data = []
            with open(inp, 'r') as f:
                for line in f:
                    if start_line in line:
                        start = True
                        continue  
                    if start:
                        if end_line in line:
                            break  
                        if "ATOM" in line or "------" in line:
                            continue  
                        if len(line.split()) == 5: 
                            cor_data.append(line)  
            
            ofile = open(out,'w')
            ofile.write(' $SYSTEM MEMDDI=400 MWORDS=200 $END\n')
            if int(self.multi) > 1:
                line = ' $CONTRL SCFTYP=UHF DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       NOSYM=1 RUNTYP=HESSIAN MAXIT=200$END\n'
                ofile.write(line)
            else:
                line = ' $CONTRL DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       NOSYM=1 RUNTYP=HESSIAN MAXIT=200 $END\n'
                ofile.write(line)
            ofile.write(' $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n')
            ofile.write(' $DATA\n')
            ofile.write('Cluster/6-31G\n')
            ofile.write('C1\n')
            for line in cor_data:
                ofile.write(' ' + line)
            ofile.write(' $END\n')
            ofile.close()
        elif basiset in metalbasisset:
            ECPinput = GAMESS_input_gen.GAMESS_ECP_input_gen(filename=self.filename,software='gms',caltype='fc',
                                                    multi=self.multi,method=self.method,
                                                    basisset=self.basisset,totalcharge=self.totalcharge,
                                                    metal=self.metal,nprocs=self.nprocs)
            ECPinput.write_gms_fc(inp,out,multi,totalcharge,basiset,method)


    def write_gms_mk(self,crds,out,multi,totalcharge,basiset,method):
        cor_data = crds
        if basiset.upper() in ['6-31G','6-31G*','6-31GS']:
            ofile = open(out,'w')
            ofile.write(' $SYSTEM MEMDDI=400 MWORDS=200 $END\n')
            line = ' $CONTRL DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=1'  + ' $END\n'
            if int(multi) > 1:
                line = ' $CONTRL SCFTYP=UHF DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       MAXIT=200 $END\n'
            ofile.write(line)
            ofile.write(' $ELPOT IEPOT=1 WHERE=PDC $END\n')
            ofile.write(' $PDC PTSEL=CONNOLLY CONSTR=NONE $END\n')
            if basiset in ['6-31G*','6-31g*','6-31gd','6-31GS']:
                ofile.write( ' $BASIS GBASIS=N31 NGAUSS=6 NDFUNC=1 $END\n')
            elif basiset in ['6-31G','6-31g']:
                ofile.write( ' $BASIS GBASIS=N31 NGAUSS=6 $END\n')
            ofile.write(' $DATA\n')
            ofile.write('Cluster/' + basiset +'\n')
            ofile.write('C1\n')
            for i in cor_data:
                ofile.write(i)
            ofile.write(' $END')
        
        elif basiset.upper() in metalbasisset:
            ECPinput = GAMESS_input_gen.GAMESS_ECP_input_gen(filename=self.filename,software='gms',caltype='mk',
                                                    multi=self.multi,method=self.method,
                                                    basisset=self.basisset,totalcharge=self.totalcharge,
                                                    metal=self.metal,nprocs=self.nprocs)
            ECPinput.write_gms_mk(out,crds=crds,multi=self.multi,totalcharge=self.totalcharge,basiset=self.basisset,method=self.method)

    def build(self):
        if self.software in ['orca','ORCA','Orca']:
            if self.caltype in ['opt','OPT']:
                crds = get_crds_from_xyz(self.filename + '_final.xyz')
                out = self.filename + '_small_opt.orca'
                self.write_orca_opt(crds=crds,out=out,multi=self.multi,totalcharge=self.totalcharge,
                                    nprocs=self.nprocs,basisset=self.basisset,method=self.method)
            elif self.caltype in ['fc','FC','Fc']:
                orcainput = self.filename + '_small_opt.orca'
                orcaoutput = self.filename + '_small_opt.orcaout'
                out = self.filename + '_small_fc.orca'
                trjxyz = orcainput + '_trj.xyz'
                crds = get_crds_from_orcaopt(orcaout=orcaoutput,trjxyz=trjxyz)
                self.write_orca_fc(crds=crds,out=out,multi=self.multi,totalcharge=self.totalcharge,
                                   nprocs=self.nprocs,basisset=self.basisset,method=self.method)
            elif self.caltype in ['MK','mk','Mk']:
                if self.opt == 'Y':
                    orcainput = self.filename + '_small_opt.orca'
                    orcaoutput = self.filename + '_small_opt.orcaout'
                    out = self.filename + '_large_mk.orca'
                    trjxyz = orcainput + '_trj.xyz'
                    inp = self.filename + '_small_opt.inp'
                    crds = get_crds_from_orcaopt(orcaout=orcaoutput,trjxyz=trjxyz)
                    self.write_orca_mk(out=out, crds=crds,multi=self.multi, nprocs=self.nprocs,
                                    totalcharge=self.totalcharge,basisset=self.basisset,method=self.method)
                elif self.opt == 'N':
                    inp = self.filename + '_small_opt.inp'
                    out = self.filename + '_large_mk.orca'
                    crds_origin = get_crds_from_gmsinp(inp)
                    new_crds = []
                    for line in crds_origin:
                        new_crds.append([line.split()[0],line.split()[2],line.split()[3],line.split()[4]])
                    self.write_orca_mk( out=out, crds=new_crds,multi=self.multi,  nprocs=self.nprocs,
                                    totalcharge=self.totalcharge,basiset=self.basisset,method=self.method)
                    
        elif self.software in ['gms','GAMESS']:
            if self.caltype in ['MK','mk','Mk']:
                if self.opt == 'Y':
                    out = self.filename + '_large_mk.inp'
                    inp = self.filename + '_small_opt.log'
                    crds = get_crds_from_gmsopt(inp)
                    #  print(crds)
                    if self.filename+'_small_opt.logback' not in glob('*'):
                        cmd = 'cp ' + inp + ' ' + self.filename+'_small_opt.logback'
                        subprocess.call(cmd,shell=True)
                    self.write_gms_mk(out=out,multi=self.multi,crds=crds,
                                totalcharge=str(self.totalcharge),basiset=self.basisset,method=self.method )
                elif self.opt == 'N':
                    out = self.filename + '_large_mk.inp'
                    inp = self.filename + '_small_opt.inp'
                    crds = get_crds_from_gmsinp(inp)
                    self.write_gms_mk(out=out,multi=self.multi,crds=crds,
                                totalcharge=str(self.totalcharge),basiset=self.basisset,method=self.method)
            
            elif self.caltype in ['opt','OPT']:
                inp=self.filename+'_small_opt.inpback'
                if inp not in glob('*'):
                    cmd = 'cp ' + self.filename+'_small_opt.inp ' + inp
                    subprocess.call(cmd,shell=True)
                out = self.filename+'_small_opt.inp'
                self.write_gms_opt(inp=inp,out=out,multi=self.multi, totalcharge=self.totalcharge,
                                    basiset=self.basisset,method=self.method)

            elif self.caltype in ['FC','fc']:
                inp=self.filename+'_small_opt.log'
                if self.filename+'_small_opt.logback' not in glob('*'):
                    cmd = 'cp ' + inp + ' ' + self.filename+'_small_opt.logback'
                    subprocess.call(cmd,shell=True)
                out = self.filename+'_small_fc.inp'
                self.write_gms_fc(inp=inp,out=out,multi=self.multi, totalcharge=self.totalcharge,
                                    basiset=self.basisset,method=self.method)
        
        elif self.software in ['gau','g09','g03']:
            if self.caltype in ['opt','OPT']:
                inp = self.filename + '_small_opt.inp'
                if self.filename + '_small_opt.inpback' not in glob('*'):
                    cmd = 'cp ' + inp + ' ' + self.filename + '_small_opt.inpback'
                    subprocess.call(cmd,shell=True)
                crds = get_crds_from_gmsinp(inp + 'back')
                new_crds = []
                out = self.filename + '_small_opt.com'
                for line in crds:
                    if len(line.split()) == 5:
                        atomname = line.split()[0]
                        x = float(line.split()[2])
                        y = float(line.split()[3])
                        z = float(line.split()[4])
                        new_crds.append([atomname,x,y,z])
                self.write_gaussian_opt(out=out,crds=new_crds,multi=self.multi,totalcharge=self.totalcharge,
                                        nprocs=self.nprocs,basisset=self.basisset,method=self.method)
            
            elif self.caltype in ['fc','Fc','FC']:
                out = self.filename + '_small_fc.com'
                self.write_gaussian_fc(out=out,nprocs=self.nprocs,basisset=self.basisset,method=self.method)
            

            elif self.caltype in ['mk','MK','Mk']:
                out =  self.filename + '_large_mk.com'
                if self.opt == 'Y':
                    self.write_gaussian_mk_from_opt(out=out,nprocs=self.nprocs,basisset=self.basisset,method=self.method)
                elif self.opt == 'N':
                    inp = self.filename + '_small_opt.inp'
                    if self.filename + '_small_opt.inpback' not in glob('*'):
                        cmd = 'cp ' + inp + ' ' + self.filename + '_small_opt.logback'
                        subprocess.call(cmd,shell=True)
                    crds = get_crds_from_gmsinp(inp + 'back')
                    new_crds = []
                    out = self.filename + '_small_opt.com'
                    for line in crds:
                        if len(line.split()) == 5:
                            atomname = line.split()[0]
                            x = float(line.split()[2])
                            y = float(line.split()[3])
                            z = float(line.split()[4])
                            new_crds.append(atomname,x,y,z)
                    self.write_gaussian_mk(crds=new_crds,out=out,multi=self.multi,totalcharge=self.totalcharge,nprocs=self.nprocs,
                                           basisset=self.basisset,method=self.method)



                


def startautoQM(argumentList):
    options = "hn:c:u:x:t:m:b:p:l:Q:"
    long_options = ["help",'filename=','totalcharge=','multi=','software=','caltype=','metal=','basisset=','nprocs=','orcapath=','opt=']
    arguments, values = getopt.getopt(argumentList,options,long_options)
    software = 'gms'
    filename = None
    totalcharge = str(0)
    multi = '1'
    caltype = None
    metal = 'default'
    method = 'B3LYP'
    basisset = 'ECP'
    nprocs = '16'
    opt = 'Y'

    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            message = '''Usage: QM_inputs_gen.py [OPTIONS]
            Options:
                -h, --help              Display this help message.
                -n, --filename          Specify the prefix of metal complex.xyz file.
                -t, --caltype           opt,mk or fc
                -u, --multi              spin multiplicity default: 1
                -x, --software          g09,g16 or gms
                -m, --metal             name of the metal, default: Default reading filename.info file
                -c, --totalcharge       totalcharge of the system
                -b, --basisset          ECP,allECP(only used in opt process),6-31G,the basisset, default: ECP, (ECP is only applied to metal
                -p, --nprocs            number of procs (only when orca))
                -Q  --opt               use(Y) or not use(N) the QM optimized structure for charge calculation default: Y

                '''
            print(message)
            sys.exit()
        elif currentArgument in ("-n","--filename"):
            print ("metal_complex filename", currentValue)
            filename = str(currentValue)
            print ('system name is ' + filename)
        elif currentArgument in ("-u",'--spin'):
            multi = str(currentValue)
            print('spin multiplicity = ' + multi)
        elif currentArgument in ("-x","--software"):
            software = str(currentValue)
            print('the input file is for ' + software)
        elif currentValue in ("-m","--metal"):
            metal = str(currentValue).capitalize()
            if metal in ['default','Default']:
                print('Reading the metal name from info files generated by autoMCPB')
            else:
                print('the metal is ' + metal)
        elif currentArgument in ("-c","--totalcharge"):
            totalcharge = str(currentValue)
            if totalcharge in ['default','Default']:
                print('Reading the charge from info files generated by autoMCPB')
            else:
                print('the total charge is set as', totalcharge)
        elif currentArgument in ("-t","--caltype"):
            caltype = str(currentValue)
            print('the QM input type is for opt')
        elif currentArgument in ('-b',"--basisset"):
            basisset = str(currentValue)
        elif currentArgument in ('-p',"--nprocs"):
            nprocs = str(currentValue)
        elif currentArgument in ('-Q','--opt'):
            opt = str(currentValue)

    builder = QM_inputs_gen(filename=filename,software=software,caltype=caltype,multi=multi,method=method,
                            basisset=basisset,metal=metal,totalcharge=totalcharge,nprocs=nprocs,opt=opt)
  #  print(builder.caltype,'!!!!!!!!!!!!!!!!!!!!!!')
    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startautoQM(argumentList)
  #  print(argumentList)
                    

                

                 
            









      