from tkinter import *
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import messagebox
from tkinter.ttk import *
import nglview
import glob
import os
import subprocess
import shutil
import imolecule

this_dir, this_filename = os.path.split(__file__)

def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol


def view_w_ngl(xyz):
    shutil.copy(xyz, '_autosolvate_input.xyz')
    SCRIPT_PATH = os.path.join(this_dir, "nb_scripts", "view_autosolvate_input.ipynb")
    shutil.copy(SCRIPT_PATH, './')

    subprocess.call('jupyter notebook view_autosolvate_input.ipynb &', shell=True)


def view_w_imol(xyz):
    imolecule.draw(open(xyz).read(), format='xyz')


def cleanUp():
    if os.path.exists('_autosolvate_input.xyz'):
        os.remove('_autosolvate_input.xyz')
    if os.path.exists('view_autosolvate_input.ipynb'):
        os.remove('view_autosolvate_input.ipynb')


colwidth = [25, 25, 7, 7, 50]

class boxgenGUI():
    def __init__(self,master):
        self.master = master
        self.master.title("Automated solvated box structure and MD parameter generator")
        self.master.geometry('820x600')
        self.full_enumeration = False
        self.cube_size_max = 100

        self.lbl00 = Label(self.master, text="Enter solute xyz file path", width=colwidth[0])
        self.lbl00.grid(column=0, row=0)
        
        self.txt01 = Entry(self.master, width=colwidth[3])
        self.txt01.grid(column=1, row=0, columnspan=3, sticky=W+E)
        
        self.lbl10 = Label(self.master, text="Current solute xyz file path:", width=colwidth[0])
        self.lbl10.grid(column=0, row=1)
        
        self.lbl11 = Label(self.master, text="Waiting...", width=colwidth[4])
        self.lbl11.grid(column=1, row=1, columnspan=6, sticky=W+E)
        
        self.xyzfile = StringVar()
        self.txt_torsion = []

        
        def set_xyz():
            my_filetypes = [('xyz files', '.xyz')]
            mypath = self.txt01.get()
            if mypath !="" and os.path.exists(mypath):
                    self.xyzfile.set(self.txt01.get())
            else:
                answer = filedialog.askopenfilename(parent=self.master,
                                    initialdir=os.getcwd(),
                                    title="No file path entered.\n Please select a file:",
                                    filetypes=my_filetypes)
                self.xyzfile.set(answer)
                self.txt01.delete(0,END)
                self.txt01.insert(0,answer)
            res = self.xyzfile.get()
            self.lbl11.configure(text=res)
        
        
        self.btn02 = Button(self.master, text="Set solute xyz", command=set_xyz, width=colwidth[3])
        self.btn02.grid(column=4, row=0,columnspan=3,sticky=W+E)
        
        # Second set of buttons: enter select the way to view the molecule
        self.lbl20 = Label(self.master, text="Select visualization method", width=colwidth[0])
        self.lbl20.grid(column=0, row=2)
        
        # Adding combobox drop down list 
        self.n = StringVar()
        self.molecule_chosen = Combobox(self.master, textvariable=self.n, width=colwidth[3])
        self.molecule_chosen['values'] = ('imolecule',
                                     'nglview',)
        
        self.molecule_chosen.current(0)
        self.molecule_chosen.grid(column=1, row=2,columnspan=3,sticky=W+E)
        
        
        def view_xyz():
            if self.molecule_chosen.get() == 'imolecule':
                view_w_imol(self.xyzfile.get())
            else:
                view_w_ngl(self.xyzfile.get())
        
        
        self.btn22 = Button(self.master, text="View", command=view_xyz, width=colwidth[3])
        self.btn22.grid(column=4, row=2,columnspan=3,sticky=W+E)
        
        # Add the radio button to control whether to run autosolvate with srun
        # TODO (will adjust the location of this option)
        self.srunuse = BooleanVar()
        
        self.lbl03 = Label(self.master, text="Use srun:", width=colwidth[0])
        self.lbl03.grid(column=0, row=3)
        
        self.rad30 = Radiobutton(self.master, text='Yes', value=True, variable=self.srunuse, width=colwidth[3])
        self.rad30.grid(column=1, row=3)
        
        self.rad31 = Radiobutton(self.master, text='No', value=False, variable=self.srunuse)
        self.rad31.grid(column=2, row=3)
        
        self.srunuse.set(False)
        
        # Verbose output
        self.verbose = BooleanVar()
        
        self.lbl04 = Label(self.master, text="Verbose output", width=colwidth[0])
        self.lbl04.grid(column=0, row=4)
        
        self.rad14 = Radiobutton(self.master, text='Yes', value=True, variable=self.verbose, width=colwidth[3])
        self.rad14.grid(column=1, row=4)
        
        self.rad24 = Radiobutton(self.master, text='No', value=False, variable=self.verbose)
        self.rad24.grid(column=2, row=4)
        
        self.verbose.set(False)
        
        # Adding combobox drop down list for selecting solvent 
        self.lbl050 = Label(self.master, text="Select solvent", width=colwidth[0])
        self.lbl050.grid(column=0, row=5)
        self.n5 = StringVar()
        self.solvent = Combobox(self.master, textvariable=self.n5, width=colwidth[3])
        self.solvent['values'] = ('water',
                                     'methanol',
                                     'chloroform',
                                     'nma',
                                     'acetonitrile')
        
        self.solvent.current(0)
        self.solvent.grid(column=1, row=5,columnspan=3,sticky=W+E)
        
        # Adding combobox drop down list for selecting charge method
        self.lbl060 = Label(self.master, text="Select charge method for force field fitting", width=colwidth[0])
        self.lbl060.grid(column=0, row=6)
        self.n6 = StringVar()
        self.charge_method = Combobox(self.master, textvariable=self.n6, width=colwidth[3])
        self.charge_method['values'] = ('bcc',
                                     'gaussian')
        
        self.charge_method.current(0)
        self.charge_method.grid(column=1, row=6,columnspan=3,sticky=W+E)

        # Output File path
        self.lbl08 = Label(self.master, text="Output directory", width=colwidth[0])
        self.lbl08.grid(column=0, row=8)
        
        self.txt18 = Entry(self.master)
        self.txt18.grid(column=1, row=8, columnspan=3, sticky=W+E)
        
        self.output_path = StringVar()
        
        def set_output_path():
            mypath = self.txt18.get()
            if mypath !="" and os.path.exists(mypath):
                    self.output_path.set(self.txt18.get())
            else:
                answer = filedialog.askdirectory(title="No file path entered.\n Please select a file:")
                self.output_path.set(answer)
                self.txt18.delete(0,END)
                self.txt18.insert(0,answer)
            res = self.output_path.get()
        
        self.btn28 = Button(self.master, text="Set", command=set_output_path, width=colwidth[3])
        self.btn28.grid(column=4, row=8,columnspan=3,sticky=W+E)
        
        # Terachem input file template
        self.lbl09 = Label(self.master, text="Output file name prefix", width=colwidth[0])
        self.lbl09.grid(column=0, row=9)
        
        
        self.txt19 = Entry(self.master)
        self.txt19.grid(column=1, row=9, columnspan=3, sticky=W+E)
        
        self.output_prefix = StringVar()
        
        def set_prefix():
            if self.txt19.get()!= "":
                    self.output_prefix.set(self.txt19.get())
            else:
                default_output_prefix = self.solvent.get()+"_solvated"
                self.output_prefix.set(default_output_prefix)
                self.txt19.insert(0,default_output_prefix)
        
        self.btn29 = Button(self.master, text="Set", command=set_prefix, width=colwidth[3])
        self.btn29.grid(column=4, row=9,columnspan=3,sticky=W+E)
        
        # Solute charge 
        self.lbl010 = Label(self.master, text="solute charge", width=colwidth[0])
        self.lbl010.grid(column=0, row=10)
        
        
        self.txt110 = Entry(self.master)
        self.txt110.grid(column=1, row=10, columnspan=3, sticky=W+E)
        
        self.charge = IntVar()
        
        def set_charge():
            answer = 0
            try: 
                answer = int(self.txt110.get())
            except ValueError:
                answer = simpledialog.askinteger(parent=self.master,
                                                 title="Dialog",
                                                 prompt="Charge must be an integer! Please re-enter.") 
                self.txt110.delete(0,END)
                self.txt110.insert(0,answer)
            self.charge.set(answer)
            res = self.charge.get()
        
        self.btn210 = Button(self.master, text="Set", command=set_charge, width=colwidth[3])
        self.btn210.grid(column=4, row=10,columnspan=3,sticky=W+E)

        # Solute spin multiplicity
        self.lbl011 = Label(self.master, text="solute spin multiplicity", width=colwidth[0])
        self.lbl011.grid(column=0, row=11)
        
        
        self.txt111 = Entry(self.master)
        self.txt111.grid(column=1, row=11, columnspan=3, sticky=W+E)
        
        self.spin_multiplicity = IntVar()
        
        def set_spin_multiplicity():
            answer = 0
            try: 
                answer = int(self.txt111.get())
            except ValueError:
                answer = simpledialog.askinteger(parent=self.master,
                                     title="Dialog",
                                     prompt="Spin multiplicity must be a positive integer! Please re-enter.")
                self.txt111.delete(0,END)
                self.txt111.insert(0,answer)
            self.spin_multiplicity.set(answer)
            res = self.spin_multiplicity.get()
        
        self.btn211 = Button(self.master, text="Set", command=set_spin_multiplicity, width=colwidth[3])
        self.btn211.grid(column=4, row=11,columnspan=3,sticky=W+E)

        # Solvent cube size
        self.lbl012 = Label(self.master, text="Solvent cube size (Angstrom)", width=colwidth[0])
        self.lbl012.grid(column=0, row=12)
        
        
        self.txt112 = Entry(self.master)
        self.txt112.grid(column=1, row=12, columnspan=3, sticky=W+E)
        
        self.cube_size = DoubleVar()
        
        def set_cubesize():
            answer = 0
            try: 
                answer = int(self.txt112.get())
            except ValueError:
                answer = simpledialog.askfloat(parent=self.master,
                                               title="Dialog",
                                               prompt="Solvent cube size must be a float! Please re-enter.")
                self.txt112.delete(0,END)
                self.txt112.insert(0,answer)
            self.cube_size.set(answer)
            res = self.cube_size.get()
        
        self.btn212 = Button(self.master, text="Set", command=set_cubesize, width=colwidth[3])
        self.btn212.grid(column=4, row=12,columnspan=3,sticky=W+E)
        
        # Sanity check to make sure that all required buttons are set
        def GUI_input_sanity_check():
            boxgen_error =0 
            if self.xyzfile.get() =="":
                  print("Solute molecule xyz file must be provided!")
                  boxgen_error = 1
            else:
                  print("Solute molecule xyz file: ", self.xyzfile.get())
            # Check charge method based on spin multiplicity
            if self.spin_multiplicity.get() > 1 and  self.charge_method_chose.get()!='gaussian':
                print("{:s} charge method cannot ".format(self.charge_method_chose.get()) +
                    + "handle open-shell system with spin multiplicity"
                    + "{:d}".format(self.spin_multiplicity.get()))
                boxgen_error = 2
            if self.cube_size.get() < 0 or self.cube_size.get() > self.cube_size_max :
                print("Solvent box size (Angstrom) must be a positive"
                    + "number no bigger than {:d}".format(self.cube_size_max))
                boxgen_error = 3 
            # TODO: add check for number of electrons and spin multiplicity
            
            return boxgen_error

        def write_boxgen_input():
            cmd = "python " + os.path.join(this_dir,"../set_up_solvated_AmberMD_ob3.py")
        
            if self.xyzfile.get() != "":
               cmd += " -m " + self.xyzfile.get()
            if self.solvent.get() != "":
                cmd += " -s " + self.solvent.get()
            if self.verbose.get() == 1:
                cmd += " -v "
            if self.charge_method.get() != "":
                cmd += " -g " + self.charge_method.get()
            if self.charge.get() != "":
                cmd += " -c {:d}".format(self.charge.get())
            if self.spin_multiplicity.get() != "":
                cmd += " -u {:d}".format(self.spin_multiplicity.get())
            if self.cube_size.get() != "":
                cmd += " -b {:f}".format(self.cube_size.get())
            if self.output_prefix !="":
                cmd += " -o {:s}".format(self.output_prefix.get())
            return cmd
        
        
        # Execute  python command to run autosolvate
        def execute():
            boxgen_error = GUI_input_sanity_check()
            if boxgen_error == 0:
                cmd = write_boxgen_input()
                res = "Congratulations! AutoSolvate command generated: \n" + cmd
                messagebox.showinfo(title="Confirmation", message=res)
                question = "Do you want to continue to generate the "\
                         + "solvent box structure and force field parameters?"
                answer = messagebox.askyesno(title="Confirmation", message=question)
                if answer == True:
                    subprocess.call(cmd, shell=True)
        
            else:
                res = "Error found in input.\n"
                if boxgen_error == 1:
                    res += "Please enter xyz file path"
                elif boxgen_error == 2:
                    res += "Please modify charge method. Only gaussian RESP charge\n"\
                          + "fitting can handle open-shell system with spin multiplicity > 1"
                elif boxgen_error == 3:
                    res += "Please re-define solvent box size. Solvent box size (Angstrom)"\
                         + "must be a positive number no bigger than {:d}".format(self.cube_size_max)
                messagebox.showerror(message=res)
        
        # Generate solvated box and MD parameter files
        self.btn132 = Button(self.master, text="Generate Solvent box and MD parameters! ", command=execute)
        self.btn132.grid(column=0, row=21, columnspan=3, sticky=W+E, padx=10, pady=5)

### START cluster extraction window ###
class clusterGUI():
    def __init__(self,master):
        self.master = master
        master.title("Microsolvated cluster extraction")
        master.geometry('820x800')
     #TODO: link to scripts that post process MD trajectories

### END Cluster extraction window ###

## The master GUI of AutoSolvate where we select what task to do ##
class autosolvateGUI():
    def __init__(self,master):
        self.master = master
        self.padx = 30
        master.title("Welcome to AutoSolvate!")
        master.geometry('400x80')
        self.lbl00 = Label(master, text="Please select the task",width=20)
        self.lbl00.grid(column=0, row=0,  sticky=W+E, padx=self.padx)
        self.lbl00.configure(anchor="center")

        # Adding combobox drop down list 
        self.task = StringVar()
        self.task_chosen = Combobox(master, textvariable=self.task, width=30)
        self.task_chosen['values'] = ('Solvated box and MD parameter generation',
                                     'Microsolvated cluster extraction') 
        self.task_chosen.current(0)
        self.task_chosen.grid(column=0, row=1,  sticky=W+E, padx=self.padx)

        # Create new window to do the task
        def create_task_window():
            if self.task_chosen.get() == 'Solvated box and MD parameter generation':
                self.master2 = Toplevel(self.master)
                my_gui = boxgenGUI(self.master2)
            else:
                self.master3 = Toplevel(self.master)
                my_gui = clusterGUI(self.master3)
        
        self.btn02 = Button(master, text="Go!", command=create_task_window, width=20)
        self.btn02.grid(column=0, row=2,  sticky=W+E, padx=self.padx)

# This part does not need to be modified
if __name__ == '__main__':
    window = Tk()
    my_gui = autosolvateGUI(window)
    window.mainloop()

    cleanUp()
