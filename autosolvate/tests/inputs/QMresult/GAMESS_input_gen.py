import basis_set_exchange as bse
import sys
from mendeleev import element
import collections

def get_basis_set(atom,basisset):
    basissetout = ' $' + atom + '\n'
    ECPout = ''
    atom_number = element(atom.capitalize()).atomic_number
    basis_str = bse.get_basis(basisset,elements=[atom_number],fmt='GAMESS_US').split('\n')
    ECP = False
    basisset_b = ''
    basisset_e = ''
    ecp_b = ''
    ecp_e = ''
    for line in basis_str:
        if 'ECP' in line:
            ECP = True
    if ECP:
        for sort, line in enumerate(basis_str):
            if '$DATA' in line:
                basisset_b = sort + 3
            elif '$ECP' in line:
                basisset_e = sort - 2
                ecp_b = sort + 1
                ecp_e = len(basis_str) - 1
        basissetpart = basis_str[basisset_b:basisset_e]
        ECPpart = basis_str[ecp_b:ecp_e]
        for line in basissetpart:
            if line != '':
                if 'END' not in line:
                  basissetout = basissetout + ' ' + line + '\n'
        basissetout = basissetout + '\n'
        for line in ECPpart:
            if line != '':
                if 'END' not in line:
                    ECPout = ECPout + ' ' + line + '\n'
    
    else:
        for sort, line in enumerate(basis_str):
            if '$DATA' in line:
                basisset_b = sort + 3
        basisset_e = len(basis_str) -1
        basissetpart = basis_str[basisset_b:basisset_e]
        for line in basissetpart:
            if line != '':
                if 'END' not in line:
                    basissetout = basissetout + ' ' + line + '\n'
        ECPout = atom.upper() + '-ECP none'
    
    return basissetout, ECPout

class GAMESS_ECP_input_gen():
    def __init__(self, filename,software,caltype,multi,method,basisset,totalcharge,metal,nprocs):
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
        
    def write_orca_mk(self,crds,out,multi,totalcharge,basiset,method):
    #    print('now start to generate GAMES ECP inputs')
        elements = []
        elements_all = []
        elements_number = []
        for line in crds:
            atomname = line[0].capitalize()
            atomnumber = str(element(atomname).atomic_number)
            elements.append(atomname)
            elements_number.append(atomnumber)
            elements_all.append(atomname)
        ofile = open(out,'w')
        ofile.write(' $SYSTEM MEMDDI=400 MWORDS=200 $END\n')
        line = ' $CONTRL DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=1'  + ' $END\n'
        if int(multi) > 1:
            line = ' $CONTRL SCFTYP=UHF DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       PP=READ MAXIT=200 $END\n'
        ofile.write(line)
        ofile.write(' $ELPOT IEPOT=1 WHERE=PDC $END\n')
        ofile.write(' $PDC PTSEL=CONNOLLY CONSTR=NONE $END\n')
        ofile.write(' $BASIS BASNAM(1)=' + elements[0] + ',')
        basnam_all = []
        for i in range(0, len(elements[1:]),10): 
            basnam = ''                     
            for  j in elements[1:][i:i+10]:
                basnam = basnam  + j + ','
            basnam_all.append(basnam)

            if len(basnam_all) == 0:
                basnam = ''
                print('Error: no atom data found in QM inputs file!')
            elif len(basnam_all) == 1 :
                basnam = ''
                basnam = basnam + basnam_all[0] + ' $END\n'
            elif len(basnam_all) == 2:
                basnam = ''
                basnam = basnam + basnam_all[0] + '\n     ' + basnam_all[1][:-1] + ' $END\n'
            elif len(basnam_all) > 2:
                basnam = ''
                basnam = basnam + basnam_all[0] 
                for i in range(len(basnam_all))[1:-1]:
                    basnam = basnam + '    \n    ' + basnam_all[i]
                basnam = basnam + '     \n    ' + basnam_all[-1][:-1] + ' $END\n'
        ofile.write(basnam)
        ofile.write(' $DATA\n')
        ofile.write('Cluster/' + basiset +'\n')
        ofile.write('C1\n')
        for i, line in enumerate(crds):
            atom = line[0]
            number = elements_number[i].capitalize()
            x = line[1]
            y = line[2]
            z = line[3]
            ofile.write(atom.ljust(6) +number.rjust(4)+ x.rjust(12)+ y.rjust(12)+z.rjust(12) + '\n')
        ofile.write(' $END\n')
        ECP_all = []
     #   print(elements_all)
        for line in elements_all:
            atom= line
            basis, ecp = get_basis_set(atom, self.basisset)
            ofile.write(basis)
            ECP_all.append(ecp)
            ofile.write('\n $END\n')
        ofile.write(' $ECP\n')
        for line in ECP_all:
            ofile.write(' ' +line + '\n')
        ofile.write(' $END\n')
        ofile.close()
    
    def write_gms_mk(self,out,multi,totalcharge,basiset,method,crds):
        elements = []
        print(crds)
        for line in crds:
            ele = line.split()[0]
            elements.append(ele)
        ofile = open(out,'w')
        ofile.write(' $SYSTEM MEMDDI=400 MWORDS=200 $END\n')
        line = ' $CONTRL DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=1'  + ' $END\n'
        if int(multi) > 1:
            line = ' $CONTRL SCFTYP=UHF DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       MAXIT=200 PP=READ $END\n'
        ofile.write(line)
        ofile.write(' $ELPOT IEPOT=1 WHERE=PDC $END\n')
        ofile.write(' $PDC PTSEL=CONNOLLY CONSTR=NONE $END\n')
        ofile.write(' $BASIS BASNAM(1)=' + elements[0] + ',')
        basnam_all = []
        for i in range(0, len(elements[1:]),10): 
            basnam = ''                     
            for  j in elements[1:][i:i+10]:
                basnam = basnam  + j + ','
            basnam_all.append(basnam)

            if len(basnam_all) == 0:
                basnam = ''
                print('Error: no atom data found in QM inputs file!')
            elif len(basnam_all) == 1 :
                basnam = ''
                basnam = basnam + basnam_all[0] + ' $END\n'
            elif len(basnam_all) == 2:
                basnam = ''
                basnam = basnam + basnam_all[0] + '\n     ' + basnam_all[1][:-1] + ' $END\n'
            elif len(basnam_all) > 2:
                basnam = ''
                basnam = basnam + basnam_all[0] 
                for i in range(len(basnam_all))[1:-1]:
                    basnam = basnam + '    \n    ' + basnam_all[i]
                basnam = basnam + '     \n    ' + basnam_all[-1][:-1] + ' $END\n'
        ofile.write(basnam)
        ofile.write(' $DATA\n')
        ofile.write('Cluster/' + basiset +'\n')
        ofile.write('C1\n')
        for line in crds:
            ofile.write(line)
        ofile.write(' $END\n')
        ECP_all = []
        print(elements)
        for line in elements:
            atom= line
            basis, ecp = get_basis_set(atom, self.basisset)
            ofile.write(basis)
            ECP_all.append(ecp)
            ofile.write('\n $END\n')
        ofile.write(' $ECP\n')
        for line in ECP_all:
            ofile.write(' ' +line + '\n')
        ofile.write(' $END\n')
        ofile.close()
    
    def write_gms_opt(self,inp,out,multi,totalcharge,basiset,method):
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
                      #  print(data_b)
                if buer == True:
                    if 'END' in line:
                        if len(line.split()) == 1:
                            data_e = sort
                            break
            cor_data = data[data_b:data_e]
        #    print(cor_data)

            elements = []
            for atom in cor_data[2:]:
                #print(atom)
                elements.append(atom.split()[0])
       #     print(elements,'########')
            
            ofile.write( " $SYSTEM MEMDDI=400 MWORDS=200 $END\n")
            if int(multi) > 1:
                ofile.write(" $CONTRL SCFTYP=UHF DFTTYP=" + method + " RUNTYP=OPTIMIZE\n      ICHARG=" + str(totalcharge) + " MULT=" + str(multi)+" MAXIT=200 PP=READ $END\n")
            else:
                ofile.write(" $CONTRL DFTTYP=" + method + " RUNTYP=OPTIMIZE\n        ICHARG=" + str(totalcharge) + " MULT=" + str(multi)+" MAXIT=200 $END\n")
            ofile.write(' $STATPT NSTEP=1000 $END\n')
        ofile.write(' $BASIS BASNAM(1)=' + elements[0] + ',')
        basnam_all = []
        for i in range(0, len(elements[1:]),10): 
            basnam = ''                     
            for  j in elements[1:][i:i+10]:
                basnam = basnam  + j + ','
            basnam_all.append(basnam)

            if len(basnam_all) == 0:
                basnam = ''
                print('Error: no atom data found in QM inputs file!')
            elif len(basnam_all) == 1 :
                basnam = ''
                basnam = basnam + basnam_all[0] + ' $END\n'
            elif len(basnam_all) == 2:
                basnam = ''
                basnam = basnam + basnam_all[0] + '\n     ' + basnam_all[1][:-1] + ' $END\n'
            elif len(basnam_all) > 2:
                basnam = ''
                basnam = basnam + basnam_all[0] 
                for i in range(len(basnam_all))[1:-1]:
                    basnam = basnam + '    \n    ' + basnam_all[i]
                basnam = basnam + '     \n    ' + basnam_all[-1][:-1] + ' $END\n'
        ofile.write(basnam)
        ofile.write(' $DATA\n')
        ofile.write('Cluster/' + basiset +'\n')
        ofile.write('C1\n')
        for line in cor_data[2:]:
            ofile.write(line)
        ofile.write(' $END\n')
        ECP_all = []
        for line in elements:
            atom= line
            basis, ecp = get_basis_set(atom, self.basisset)
            ofile.write(basis)
        #    print(basis)
            ECP_all.append(ecp)
            ofile.write('\n $END\n')
        ofile.write(' $ECP\n')
        for line in ECP_all:
            ofile.write(' ' +line + '\n')
        ofile.write(' $END\n')
        ofile.close()
    
    def write_gms_fc(self,inp,out,multi,totalcharge,basiset,method):
        start = False
        start_line = "EQUILIBRIUM GEOMETRY LOCATED"
        end_line = "INTERNUCLEAR DISTANCES (ANGS.)"
        cor_data = []
        elements = []
        with open(inp, 'r') as f:
            for line in f:
             #   print(line)
                if start_line in line:
                    start = True
                    continue  
                if start:
                    if end_line in line:
                        break  
                    if "ATOM" in line or "------" in line:
                        continue  
                    if len(line.split()) == 5:
                        ele = line.split()[0] 
                        cor_data.append(line)
                        elements.append(ele)
        ofile = open(out,'w')
        ofile.write(' $SYSTEM MEMDDI=400 MWORDS=200 $END\n')
        line = ' $CONTRL DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       NOSYM=1 RUNTYP=HESSIAN MAXIT=200 PP=READ $END\n'
        if int(multi) > 1:
            line = ' $CONTRL SCFTYP=UHF DFTTYP=' + method.upper() + ' ICHARG=' + str(totalcharge) + ' MULT=' + str(multi) + '\n       NOSYM=1 RUNTYP=HESSIAN MAXIT=200 PP=READ $END\n'
        ofile.write(line)
        ofile.write( ' $SCF DIRSCF=.TRUE. $END\n $CPHF CPHF=AO $END\n')
        ofile.write(' $BASIS BASNAM(1)=' + elements[0] + ',')
        basnam_all = []
        for i in range(0, len(elements[1:]),10): 
            basnam = ''                     
            for  j in elements[1:][i:i+10]:
                basnam = basnam  + j + ','
            basnam_all.append(basnam)

            if len(basnam_all) == 0:
                basnam = ''
                print('Error: no atom data found in QM inputs file!')
            elif len(basnam_all) == 1 :
                basnam = ''
                basnam = basnam + basnam_all[0] + ' $END\n'
            elif len(basnam_all) == 2:
                basnam = ''
                basnam = basnam + basnam_all[0] + '\n     ' + basnam_all[1][:-1] + ' $END\n'
            elif len(basnam_all) > 2:
                basnam = ''
                basnam = basnam + basnam_all[0] 
                for i in range(len(basnam_all))[1:-1]:
                    basnam = basnam + '    \n    ' + basnam_all[i]
                basnam = basnam + '     \n    ' + basnam_all[-1][:-1] + ' $END\n'
        ofile.write(basnam)
        ofile.write(' $DATA\n')
        ofile.write('Cluster/' + basiset +'\n')
        ofile.write('C1\n')
        for line in cor_data:
            ofile.write(line)
        ofile.write(' $END\n')
        ECP_all = []
        for line in elements:
            atom= line
            basis, ecp = get_basis_set(atom, self.basisset)
            ofile.write(basis)
            ECP_all.append(ecp)
            ofile.write('\n $END\n')
        ofile.write(' $ECP\n')
        for line in ECP_all:
            ofile.write(' ' +line + '\n')
        ofile.write(' $END\n')
        ofile.close()
    
    def write_gms_opt(self,inp,out,multi,totalcharge,basiset,method):
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
                      #  print(data_b)
                if buer == True:
                    if 'END' in line:
                        if len(line.split()) == 1:
                            data_e = sort
                            break
            cor_data = data[data_b:data_e]
        #    print(cor_data)

            elements = []
            for atom in cor_data[2:]:
                #print(atom)
                elements.append(atom.split()[0])
        #    print(elements,'########')
            
            ofile.write( " $SYSTEM MEMDDI=400 MWORDS=200 $END\n")
            if int(multi) > 1:
                ofile.write(" $CONTRL SCFTYP=UHF DFTTYP=" + method + " RUNTYP=OPTIMIZE\n      ICHARG=" + str(totalcharge) + " MULT=" + str(multi)+" MAXIT=200 PP=READ $END\n")
            else:
                ofile.write(" $CONTRL DFTTYP=" + method + " RUNTYP=OPTIMIZE\n        ICHARG=" + str(totalcharge) + " MULT=" + str(multi)+" MAXIT=200 $END\n")
            ofile.write(' $STATPT NSTEP=1000 $END\n')
        ofile.write(' $BASIS BASNAM(1)=' + elements[0] + ',')
        basnam_all = []
        for i in range(0, len(elements[1:]),10): 
            basnam = ''                     
            for  j in elements[1:][i:i+10]:
                basnam = basnam  + j + ','
            basnam_all.append(basnam)

            if len(basnam_all) == 0:
                basnam = ''
                print('Error: no atom data found in QM inputs file!')
            elif len(basnam_all) == 1 :
                basnam = ''
                basnam = basnam + basnam_all[0] + ' $END\n'
            elif len(basnam_all) == 2:
                basnam = ''
                basnam = basnam + basnam_all[0] + '\n     ' + basnam_all[1][:-1] + ' $END\n'
            elif len(basnam_all) > 2:
                basnam = ''
                basnam = basnam + basnam_all[0] 
                for i in range(len(basnam_all))[1:-1]:
                    basnam = basnam + '    \n    ' + basnam_all[i]
                basnam = basnam + '     \n    ' + basnam_all[-1][:-1] + ' $END\n'
        ofile.write(basnam)
        ofile.write(' $DATA\n')
        ofile.write('Cluster/' + basiset +'\n')
        ofile.write('C1\n')
        for line in cor_data[2:]:
            ofile.write(line)
        ofile.write(' $END\n')
        ECP_all = []
        for line in elements:
            atom= line
            basis, ecp = get_basis_set(atom, self.basisset)
            ofile.write(basis)
        #    print(basis)
            ECP_all.append(ecp)
            ofile.write('\n $END\n')
        ofile.write(' $ECP\n')
        for line in ECP_all:
            ofile.write(' ' +line + '\n')
        ofile.write(' $END\n')
        ofile.close()
    
    



