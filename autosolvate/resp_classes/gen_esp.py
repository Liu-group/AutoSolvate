import numpy as np
import subprocess
from pyvdwsurface import vdwsurface

AtoB = 1/0.52917721092

def read_xyz(xyz):
    crds = []
    elements = []
    ifile = open(xyz,'r')
    data = ifile.readlines()
    ifile.close()
    for line in data[2:]:
        ele = str(line.split()[0])
        x = float(line.split()[1])
        y = float(line.split()[2])
        z = float(line.split()[3])
        crds.append([x,y,z])
        elements.append(ele.capitalize())
  #  print(crds, elements)
    return crds, elements

def gen_grids(crds,elements,orcapath,gbw,denisty,out):
    gridall = []
    atoms = np.array(crds,dtype=float)
    elements_bytes = [element.encode('utf-8') for element in elements]
    points_1 = vdwsurface(atoms, elements=elements_bytes, density=1,scale_factor = 2.0)
    new_point_1 =  [[x *AtoB for x in sublist] for sublist in points_1]
    gridall.extend(new_point_1)

    points_2 = vdwsurface(atoms, elements=elements_bytes, density=1,scale_factor = 1.4)
    new_point_2=  [[x *AtoB for x in sublist] for sublist in points_2]
    gridall.extend(new_point_2)

    points_3 = vdwsurface(atoms, elements=elements_bytes, density=1,scale_factor = 1.6)
    new_point_3=  [[x *AtoB for x in sublist] for sublist in points_3]
    gridall.extend(new_point_3)

    points_4 = vdwsurface(atoms, elements=elements_bytes, density=1,scale_factor = 1.8)
    new_point_4=  [[x *AtoB for x in sublist] for sublist in points_4]
    gridall.extend(new_point_4)

   # print(gridall)

    with open('esp.xyz','w') as f:
        f.write(str(len(gridall))+'\n')
        for line in gridall:
            for crd in line:
                f.write(str(crd) + '  ')
            f.write('\n')
    cmd = orcapath + ' ' + gbw  + ' ' + denisty + ' esp.xyz ' + out
    print(cmd)
    with open('esp_gen.log','w') as f:
        subprocess.call(cmd,shell=True,stdout=f, stderr=subprocess.STDOUT)

def convertoesp(espin,espout,crds):
    ifile = open(espin,'r')
    data = ifile.readlines()
    ifile.close()
    length_crds = len(crds)
    length_esp = int(data[0])
    with open(espout,'w') as f:
        print("%5d%5d%5d" %(length_crds, length_esp, 0), file=f)
        for i in crds:
            print("%16s %15.7E %15.7E %15.7E" %(' ', float(i[0])*AtoB, float(i[1])*AtoB, float(i[2])*AtoB), file=f)
        for i in data[1:]:
            line = i.split()
            print("%16.7E %15.7E %15.7E %15.7E" %(float(line[3]), float(line[0]), float(line[1]), float(line[2])), file=f)
