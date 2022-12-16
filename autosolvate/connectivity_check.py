import getopt, sys, os
import subprocess
from openbabel import pybel
from openbabel import openbabel as ob




def read_xyz_multiple(fname):
    f = open(fname, "r")
    lines = f.readlines()
    f.close()
    curlineidx = 0
    elements, datas = [], []
    while curlineidx < len(lines) and lines[curlineidx].strip():
        natom = int(lines[curlineidx])
        datastr = lines[curlineidx+2:curlineidx+2+natom]
        data = [[0., 0., 0.] for i in range(natom)]
        element = []
        for i, line in enumerate(datastr):
            temp = line.split()
            if len(temp) != 4:
                continue
            a, x, y, z = temp
            element.append(a)
            data[i][0] = float(x)
            data[i][1] = float(y)
            data[i][2] = float(z)
        elements.append(element)
        datas.append(data)
        curlineidx = curlineidx+2+natom
    return elements, datas

def write_xyz(fname, element, data):
    f = open(fname, "w")
    f.write(str(len(element)) + "\n")
    f.write("\n")
    for i in range(len(element)):
        f.write(f"{element[i]}    {data[i][0]:>12.8f}    {data[i][1]:>12.8f}    {data[i][2]:>12.8f}\n")
    f.close()

def getResLabel(i:int):
    i = i % 26
    A2Z = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    return "SL"+A2Z[i]

def getFragments(filename:str):
    """check the connectivity of a pdb file. output each fragment as separate pdb file"""
    basename, ext = os.path.splitext(filename)
    pdbname = basename + "_temp.pdb"
    xyzname = basename + "_temp.xyz"
    ext = ext[1:]
    subprocess.run(f"obabel -i {ext} {filename} -o xyz > {xyzname}", shell = True)
    subprocess.run(f"obabel -i {ext} {filename} -o pdb > {pdbname}", shell = True)

    element, data = read_xyz_multiple(f"{basename}.xyz")
    element, data = element[0], data[0]
    with open(pdbname, "r") as f:
        connectdata = [line for line in f if line.startswith("CONECT")]
    connectdata = [[int(dat)-1 for dat in cline.split()[1:]] for cline in connectdata]
    connectdata = {cline[0]:cline[1:] for cline in connectdata}
    natom = len(element)

    visited = [0 for i in range(len(element))]
    fragments = []
    def dfs(i, residxs):
        if visited[i]:
            return
        residxs.append(i)
        visited[i] = 1
        for j in connectdata[i]:
            dfs(j, residxs)
    for i in range(natom):
        if visited[i]:
            continue
        residxs = []
        dfs(i, residxs)
        fragments.append(sorted(residxs))

    nfragment = len(fragments)
    allelem, alldata = [], []
    for i, fragment in enumerate(nfragment):
        sepelem, sepdata = [], []
        for idx in fragment:
            sepelem.append(element[idx])
            sepdata.append(data[idx])
        fragxyzname = basename + f"-fragment{i}.xyz"
        write_xyz(fragxyzname, sepelem, sepdata)
        allelem += sepelem
        alldata += sepdata
    write_xyz(xyzname, allelem, alldata)

    inst = ob.OBConversion()
    inst.SetInAndOutFormats("xyz", "pdb")
    allmol = pybel.readfile(xyzname).__next__()
    for i in range(nfragment):
        mol = pybel.readfile(basename + f"-fragment{i}.xyz").__next__()
        mol:pybel.Molecule
        res = mol.OBMol.GetResidue(0)
        res.SetName(getResLabel(i))
        inst.WriteFile(mol, basename + f"-fragment{i}.pdb")
        res = 
    
    
    

        
    return fragments

if __name__ == "__main__":
    result = getFragments("01.xyz")
    print(result)