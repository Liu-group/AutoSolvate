import tkinter as tk
from tkinter import *
from tkinter import filedialog
from tkinter import simpledialog
from tkinter import messagebox
from tkinter.ttk import *
from PIL import ImageTk, Image
import nglview
import glob
import os
import subprocess
import shutil
import imolecule
import pkg_resources

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

class baseGUI():
    def __init__(self,master):
        self.master = master
        self.padx = 5
        self.irow = 0

    def display_logo(self):
        path = pkg_resources.resource_filename('autosolvate', 'GUI/images/logo.png')

        #Creates a Tkinter-compatible photo image, which can be used everywhere Tkinter expects an image object. Scale the image to fix into the current window
        self.master.update()
        #win_width = self.master.winfo_width() - self.padx*2
        win_width = 360
        logo_padx = (self.master.winfo_width() - win_width)/2.0
        img = Image.open(path)
        zoom = win_width/img.size[0]
        #multiple image size by zoom
        pixels_x, pixels_y = tuple([int(zoom * x)  for x in img.size])
        scaled_img = ImageTk.PhotoImage(img.resize((pixels_x, pixels_y)))

        #The Label widget is a standard Tkinter widget used to display a text or image on the screen.
        self.logo = Label(self.master, image = scaled_img)
        self.logo.image = scaled_img
        self.logo.grid(column=0, row=self.irow, columnspan=6, sticky=W+E, padx=logo_padx)
        self.irow += 1


class boxgenGUI(baseGUI):
    def __init__(self,master):
        super().__init__(master)
        self.master.title("Automated solvated box structure and MD parameter generator")
        self.master.geometry('820x500')
        self.cube_size_max = 100
        
        self.display_logo()

        ### Add a solute path
        self.lbl00 = Label(self.master, text="Solute xyz file path", width=colwidth[0])
        self.lbl00.grid(column=0, row=self.irow)
        
        self.txt01 = Entry(self.master, width=colwidth[3])
        self.txt01.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.xyzfile = StringVar()
        self.txt_torsion = []
        
        def set_xyz():
            my_filetypes = [('xyz files', '.xyz'),
                            ('all files', '.*')]
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
        self.btn02.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### Display solute path
        self.lbl10 = Label(self.master, text="Current solute xyz file path:", width=colwidth[0])
        self.lbl10.grid(column=0, row=self.irow)
        
        self.lbl11 = Label(self.master, text="Waiting...", width=colwidth[4])
        self.lbl11.grid(column=1, row=self.irow, columnspan=6, sticky=W+E)
        self.irow += 1
        
        ### Visualize the solute molecule
        self.lbl20 = Label(self.master, text="Select visualization method", width=colwidth[0])
        self.lbl20.grid(column=0, row=self.irow)
        
        # Adding combobox drop down list 
        self.vtxt = StringVar()
        self.view_chosen = Combobox(self.master, textvariable=self.vtxt, width=colwidth[3])
        self.view_chosen['values'] = ('imolecule',
                                     'nglview')
        
        self.view_chosen.current(0)
        self.view_chosen.grid(column=1, row=self.irow,columnspan=3,sticky=W+E)
        
        
        def view_xyz():
            if self.view_chosen.get() == 'imolecule':
                view_w_imol(self.xyzfile.get())
            else:
                view_w_ngl(self.xyzfile.get())
        
        
        self.btn22 = Button(self.master, text="View", command=view_xyz, width=colwidth[3])
        self.btn22.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1
        
        ### Verbose output or not
        self.verbose = BooleanVar()
        
        self.lbl04 = Label(self.master, text="Verbose output", width=colwidth[0])
        self.lbl04.grid(column=0, row=self.irow)
        
        self.rad14 = Radiobutton(self.master, text='Yes', value=True, variable=self.verbose, width=colwidth[3])
        self.rad14.grid(column=1, row=self.irow)
        
        self.rad24 = Radiobutton(self.master, text='No', value=False, variable=self.verbose)
        self.rad24.grid(column=2, row=self.irow)
        
        self.verbose.set(False)
        self.irow += 1
        
        ### Adding combobox drop down list for selecting solvent 
        self.lbl050 = Label(self.master, text="Select solvent", width=colwidth[0])
        self.lbl050.grid(column=0, row=self.irow)
        self.stxt = StringVar()
        self.solvent = Combobox(self.master, textvariable=self.stxt, width=colwidth[3])
        self.solvent['values'] = ('water',
                                     'methanol',
                                     'chloroform',
                                     'nma',
                                     'acetonitrile')
        
        self.solvent.current(0)
        self.solvent.grid(column=1, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1
        
        ### Adding combobox drop down list for selecting charge method
        self.lbl060 = Label(self.master, text="Select charge method for\n force field fitting", width=colwidth[0])
        self.lbl060.grid(column=0, row=self.irow)
        self.chargetxt = StringVar()
        self.charge_method = Combobox(self.master, textvariable=self.chargetxt, width=colwidth[3])
        self.charge_method['values'] = ('bcc',
                                     'resp')
        
        self.charge_method.current(0)
        self.charge_method.grid(column=1, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Output File path
        self.lbl08 = Label(self.master, text="Output directory", width=colwidth[0])
        self.lbl08.grid(column=0, row=self.irow)
        
        self.txt18 = Entry(self.master)
        self.txt18.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.output_path = StringVar()
        
        def set_output_path():
            mypath = self.txt18.get()
            if mypath !="" and os.path.exists(mypath):
                    self.output_path.set(self.txt18.get())
            else:
                answer = filedialog.askdirectory(title="No file path entered.\n Please select a path:")
                self.output_path.set(answer)
                self.txt18.delete(0,END)
                self.txt18.insert(0,answer)
            res = self.output_path.get()
        
        self.btn28 = Button(self.master, text="Set", command=set_output_path, width=colwidth[3])
        self.btn28.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1
        
        ### Terachem input file template
        self.lbl09 = Label(self.master, text="Output file name prefix", width=colwidth[0])
        self.lbl09.grid(column=0, row=self.irow)
        
        
        self.txt19 = Entry(self.master)
        default_output_prefix = self.solvent.get()+"_solvated"
        self.txt19.insert(0,default_output_prefix)
        self.txt19.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.output_prefix = StringVar(value=default_output_prefix)
        
        def set_prefix():
            if self.txt19.get()!= "":
                    self.output_prefix.set(self.txt19.get())
            else:
                default_output_prefix = self.solvent.get()+"_solvated"
                self.output_prefix.set(default_output_prefix)
                self.txt19.insert(0,default_output_prefix)
        
        self.btn29 = Button(self.master, text="Set", command=set_prefix, width=colwidth[3])
        self.btn29.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Solute charge 
        self.lbl010 = Label(self.master, text="solute charge", width=colwidth[0])
        self.lbl010.grid(column=0, row=self.irow)
        
        
        self.txt110 = Entry(self.master)
        defaultCharge = 0
        self.txt110.insert(0, str(defaultCharge))
        self.txt110.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.charge = IntVar(value=defaultCharge)
        
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
        
        self.btn210 = Button(self.master, text="Set", command=set_charge, width=colwidth[3])
        self.btn210.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Solute spin multiplicity
        self.lbl011 = Label(self.master, text="solute spin multiplicity", width=colwidth[0])
        self.lbl011.grid(column=0, row=self.irow)
        
        
        self.txt111 = Entry(self.master)
        defaultSpin = 1
        self.txt111.insert(0, str(defaultSpin))
        self.txt111.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.spin_multiplicity = IntVar(value=defaultSpin)
        
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
        self.btn211.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)

        self.irow += 1
        
        ### Solvent cube size
        self.lbl012 = Label(self.master, text="Solvent cube size (Angstrom)", width=colwidth[0])
        self.lbl012.grid(column=0, row=self.irow)
        
        
        self.txt112 = Entry(self.master)
        defaultCubesize = 54
        self.txt112.insert(0, str(defaultCubesize))
        self.txt112.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.cube_size = DoubleVar(value=defaultCubesize)
        
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
        self.btn212.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Reminder that the following options are not required
        #self.lblMsg1 = Label(self.master, text="The following entry is optional", font='Helvetica 18 bold')
        #self.lblMsg1.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        #self.irow += 1
        #msg2 = "AmberTools have been installed automatically when installing autosolvate\n"
        #msg2 += "You will want to specify another AMBERHOME directory only if you have another version\n"
        #msg2 += "of Amber installed and you prefer to use that. Otherwise, simply skip this entry."
        #self.lblMsg2 = Label(self.master, text=msg2)
        #self.lblMsg2.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        #self.irow += 1
        
        ### AMBERHOME path
        self.lbl013 = Label(self.master, text="AMBERHOME directory (optional)", width=colwidth[0])
        self.lbl013.grid(column=0, row=self.irow)
        
        self.txt113 = Entry(self.master)
        self.txt113.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.amberhome = StringVar()
        
        def set_amberhome():
            mypath = self.txt113.get()
            if mypath !="" and os.path.exists(mypath):
                    self.amberhome.set(self.txt113.get())
            else:
                answer = filedialog.askdirectory(title="No file path entered.\n Please select a path:")
                self.amberhome.set(answer)
                self.txt113.delete(0,END)
                self.txt113.insert(0,answer)
            res = self.amberhome.get()
        
        self.btn213 = Button(self.master, text="Set", command=set_amberhome, width=colwidth[3])
        self.btn213.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Reminder that the following options are required if using RESP
        #self.lblMsg3 = Label(self.master, text="The following entries are required if using RESP for charge fitting", font='Helvetica 18 bold')
        #self.lblMsg3.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        #self.irow += 1
        #msg4 = "If you skip the specifying the following entries, default path for Gaussian 16 will be set.\n"
        #msg4 = "That may not work because the default path may not exist.\n"
        #msg4 = "If failed, you should re-run AutoSolvate and enter the Gaussian path.\n"
        #self.lblMsg4 = Label(self.master, text=msg4)
        #self.lblMsg4.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        #self.irow += 1
        ### Gaussian EXE
        self.lbl014 = Label(self.master, text="Select gaussian exe (optional)", width=colwidth[0])
        self.lbl014.grid(column=0, row=self.irow)
        self.gstxt = StringVar()
        self.gaussianexe = Combobox(self.master, textvariable=self.gstxt, width=colwidth[3])
        self.gaussianexe['values'] = ('None',
                                      'g09',
                                      'g16')
        
        self.gaussianexe.current(0)
        self.gaussianexe.grid(column=1, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Gaussian path
        self.lbl015 = Label(self.master, text="Gaussian EXE directory (optional)", width=colwidth[0])
        self.lbl015.grid(column=0, row=self.irow)
        
        self.txt115 = Entry(self.master)
        self.txt115.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.gaussiandir = StringVar()
        
        def set_gaussiandir():
            mypath = self.txt115.get()
            if mypath !="" and os.path.exists(mypath):
                    self.gaussiandir.set(self.txt115.get())
            else:
                answer = filedialog.askdirectory(title="No file path entered.\n Please select a path:")
                self.gaussiandir.set(answer)
                self.txt115.delete(0,END)
                self.txt115.insert(0,answer)
            res = self.gaussiandir.get()
        
        self.btn215 = Button(self.master, text="Set", command=set_gaussiandir, width=colwidth[3])
        self.btn215.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Sanity check to make sure that all required buttons are set
        def GUI_input_sanity_check():
            boxgen_error =0 
            if self.xyzfile.get() =="":
                  print("Solute molecule xyz file must be provided!")
                  boxgen_error = 1
            else:
                  print("Solute molecule xyz file: ", self.xyzfile.get())
            # Check charge method based on spin multiplicity
            if self.spin_multiplicity.get() > 1 and  self.charge_method.get()!='resp':
                print("{:s} charge method cannot ".format(self.charge_method.get()) 
                    + "handle open-shell system with spin multiplicity"
                    + "{:d}".format(self.spin_multiplicity.get()))
                boxgen_error = 2
            if self.cube_size.get() < 0 or self.cube_size.get() > self.cube_size_max :
                print("Solvent box size (Angstrom) must be a positive"
                    + "number no bigger than {:d}".format(self.cube_size_max))
                boxgen_error = 3 
            # TODO: add check for number of electrons and spin multiplicity
            if self.amberhome.get()=="":
                print("WARNING: AMBERHOME not provided from GUI.")
            if self.charge_method.get()=='resp':
                if self.gaussianexe.get()=='None':
                    print("WARNING: Gaussian exe file name not specified for RESP charge fitting.")
                    print("WARNING: Will use default value with the risk to fail later.")
                if self.gaussiandir.get()=="":
                    print("WARNING: Gaussian exe directory not specified for RESP charge fitting.")
                    print("WARNING: Will use default value with the risk to fail later.")
            
            return boxgen_error

        def write_boxgen_input():
            cmd = "autosolvate boxgen"
        
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
            if self.output_prefix.get() !="":
                cmd += " -o {:s}".format(self.output_prefix.get())
            if self.gaussianexe.get() !="None":
                cmd += " -e {:s}".format(self.gaussianexe.get())
            if len(self.gaussiandir.get()) > 0:
                cmd += " -d {:s}".format(self.gaussiandir.get())
            if len(self.amberhome.get()) > 0:
                cmd += " -a {:s}".format(self.amberhome.get())
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
        self.btn132.grid(column=0, row=self.irow, columnspan=3, sticky=W+E, padx=self.padx, pady=5)
        self.irow += 1

### START MD automation window ###
class mdGUI(baseGUI):
    def __init__(self, master):
        super().__init__(master)
        self.master.title("MD simulation automation")
        self.master.geometry('820x730')
        self.display_logo()
        self.padx = 10
        self.pady = 5

        ### classical MD control block starts 
        self.lblMain = Label(self.master, text="Essential control options", font='Helvetica 18 bold')
        self.lblMain.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        self.irow += 1

        ### Enter .parmtop and .inpcrd filename prefix 
        self.lbl00 = Label(self.master, text="File name prefix for existing\n.inpcrd and .parmtop files", width=colwidth[0])
        self.lbl00.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txt01 = Entry(self.master)
        self.txt01.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.prefix = StringVar()
        
        def set_prefix():
            mypath = self.txt01.get()
            inpcrd = "{:s}.inpcrd".format(mypath)
            prmtop = "{:s}.prmtop".format(mypath)
            if mypath !="" and os.path.exists(inpcrd) and os.path.exists(prmtop):
                    self.prefix.set(self.txt01.get())
            else:
                msg = "Invalid prefix provided!\n"
                answer = ""
                answerValid = False
                if mypath !="":
                    if not os.path.exists(inpcrd):
                        msg+= "File " + inpcrd + " does not exist!\n"
                    if not os.path.exists(prmtop):
                        msg+= "File " + prmtop + " does not exist!\n"
                    msg+= "Please re-try\n"
                else:
                    msg = "Empty file prefix. Please enter.\n"
                while (not answerValid):
                    answer = simpledialog.askstring(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    inpcrd = "{:s}.inpcrd".format(answer)
                    prmtop = "{:s}.prmtop".format(answer)
                    msg = "Invalid prefix provided!\n"
                    answerValid = True
                    if not os.path.exists(inpcrd):
                        msg+= "File " + inpcrd + " does not exist!\n"
                        answerValid = False
                    if not os.path.exists(prmtop):
                        msg+= "File " + prmtop + " does not exist!\n"
                        answerValid = False
                    msg+= "Please re-try\n"

                self.prefix.set(answer)
                self.txt01.delete(0,END)
                self.txt01.insert(0,answer)
        
        self.btn02 = Button(self.master, text="Set file prefix", command=set_prefix, width=colwidth[3])
        self.btn02.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### set Temperature
        self.lblTemp = Label(self.master, text="Temperature (K)", width=colwidth[0])
        self.lblTemp.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtTemp = Entry(self.master)
        defaultTemp = 298
        self.txtTemp.insert(0,str(defaultTemp))
        self.txtTemp.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.Temp = DoubleVar(value=defaultTemp)
        
        def set_Temp():
            answer = defaultTemp
            try: 
                answer = float(self.txtTemp.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "Temperature must be a positive float!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askfloat(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = float(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False

            self.Temp.set(answer)
            self.txtTemp.delete(0,END)
            self.txtTemp.insert(0,answer)
        
        self.btnTemp = Button(self.master, text="Set", command=set_Temp, width=colwidth[3])
        self.btnTemp.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### set Pressure
        self.lblPressure = Label(self.master, text="Pressure (bar)", width=colwidth[0])
        self.lblPressure.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtPressure = Entry(self.master)
        defaultPressure = 1
        self.txtPressure.insert(0,str(defaultPressure))
        self.txtPressure.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.Pressure = DoubleVar(value=defaultPressure)
        
        def set_Pressure():
            answer = defaultPressure
            try: 
                answer = float(self.txtPressure.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "Pressure must be a positive float (unit: bar)!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askfloat(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = float(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False

            self.Pressure.set(answer)
            self.txtPressure.delete(0,END)
            self.txtPressure.insert(0,answer)
        
        self.btnPressure = Button(self.master, text="Set", command=set_Pressure, width=colwidth[3])
        self.btnPressure.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### classical MD control block starts 
        self.lblMD1 = Label(self.master, text="Classical MD control options", font='Helvetica 18 bold')
        self.lblMD1.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        self.irow += 1

        ### set number of steps for MM minimization
        self.lblMMMinSteps = Label(self.master, text="MM minimization steps", width=colwidth[0])
        self.lblMMMinSteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtMMMinSteps = Entry(self.master)
        defaultMMMinSteps = 1000
        self.txtMMMinSteps.insert(0,str(defaultMMMinSteps))
        self.txtMMMinSteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.MMMinSteps = IntVar(value=defaultMMMinSteps)
        
        def set_MMMinSteps():
            answer = defaultMMMinSteps
            try: 
                answer = int(self.txtMMMinSteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "MM minimization steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False

            self.MMMinSteps.set(answer)
            self.txtMMMinSteps.delete(0,END)
            self.txtMMMinSteps.insert(0,answer)
        
        self.btnMMMinSteps = Button(self.master, text="Set", command=set_MMMinSteps, width=colwidth[3])
        self.btnMMMinSteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1


        ### set number of steps for MM heat up
        self.lblMMHeatSteps = Label(self.master, text="MM heat up steps", width=colwidth[0])
        self.lblMMHeatSteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtMMHeatSteps = Entry(self.master)
        defaultMMHeatSteps = 1000
        self.txtMMHeatSteps.insert(0,str(defaultMMHeatSteps))
        self.txtMMHeatSteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.MMHeatSteps = IntVar(value=defaultMMHeatSteps)
        
        def set_MMHeatSteps():
            answer = defaultMMHeatSteps
            try: 
                answer = int(self.txtMMHeatSteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "MM heat up steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False

            self.MMHeatSteps.set(answer)
            self.txtMMHeatSteps.delete(0,END)
            self.txtMMHeatSteps.insert(0,answer)

        self.btnMMHeatSteps = Button(self.master, text="Set", command=set_MMHeatSteps, width=colwidth[3])
        self.btnMMHeatSteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### set number of steps for MM NPT pressure equilibration 
        self.lblMMNPTSteps = Label(self.master, text="MM NPT pressure equilibration steps", width=colwidth[0])
        self.lblMMNPTSteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtMMNPTSteps = Entry(self.master)
        defaultMMNPTSteps = 1000
        self.txtMMNPTSteps.insert(0,str(defaultMMNPTSteps))
        self.txtMMNPTSteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.MMNPTSteps = IntVar(value=defaultMMNPTSteps)
        
        def set_MMNPTSteps():
            answer = defaultMMNPTSteps
            try: 
                answer = int(self.txtMMNPTSteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "MM NPT equilibration steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False

            self.MMNPTSteps.set(answer)
            self.txtMMNPTSteps.delete(0,END)
            self.txtMMNPTSteps.insert(0,answer)
        
        self.btnMMNPTSteps = Button(self.master, text="Set", command=set_MMNPTSteps, width=colwidth[3])
        self.btnMMNPTSteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1
        
        ### set number of steps for MM NVE production run steps 
        self.lblMMNVESteps = Label(self.master, text="MM NVE production run steps", width=colwidth[0],foreground='blue')
        self.lblMMNVESteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtMMNVESteps = Entry(self.master)
        defaultMMNVESteps = 0
        self.txtMMNVESteps.insert(0,str(defaultMMNVESteps))
        self.txtMMNVESteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.MMNVESteps = IntVar(value=defaultMMNVESteps)
        
        def set_MMNVESteps():
            answer = defaultMMNVESteps
            try: 
                answer = int(self.txtMMNVESteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "MM NVE steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
        
            self.MMNVESteps.set(answer)
            self.txtMMNVESteps.delete(0,END)
            self.txtMMNVESteps.insert(0,answer)
        
        style = Style()
        self.btnMMNVESteps = Button(self.master, text="Set", command=set_MMNVESteps, width=colwidth[3])
        self.btnMMNVESteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### Do QM/MM or not
        self.doQMMM = BooleanVar(value=False)
        
        self.lbldoQMMM = Label(self.master, text="Do QMM/MM ?", width=colwidth[0])
        self.lbldoQMMM.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.radQMMMyes = Radiobutton(self.master, text='Yes', value=True, variable=self.doQMMM, width=colwidth[3])
        self.radQMMMyes.grid(column=1, row=self.irow)
        
        self.radQMMMno = Radiobutton(self.master, text='No', value=False, variable=self.doQMMM)
        self.radQMMMno.grid(column=2, row=self.irow)
        
        self.irow += 1

        self.lblQMMMreminder = Label(self.master, text="If \"Yes\", answer questions in the next section", font='Helvetica 14 bold', foreground='red')
        self.lblQMMMreminder.grid(column=0, row=self.irow, columnspan=3, sticky=W, padx=self.padx, pady=self.pady)
        self.irow += 1
        

        ### QM/MM MD control block starts
        self.lblQMMM1 = Label(self.master, text="QM/MM control options", font='Helvetica 18 bold')
        self.lblQMMM1.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        self.irow += 1

        ### Solute charge 
        self.lblcharge = Label(self.master, text="solute charge", width=colwidth[0])
        self.lblcharge.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        
        self.txtcharge = Entry(self.master)
        defaultCharge = 0
        self.txtcharge.insert(0, str(defaultCharge))
        self.txtcharge.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.charge = IntVar(value=defaultCharge)
        
        def set_charge():
            answer = 0
            try: 
                answer = int(self.txtcharge.get())
            except ValueError:
                answer = simpledialog.askinteger(parent=self.master,
                                                 title="Dialog",
                                                 prompt="Charge must be an integer! Please re-enter.") 
                self.txtcharge.delete(0,END)
                self.txtcharge.insert(0,answer)
            self.charge.set(answer)
        
        self.btncharge = Button(self.master, text="Set", command=set_charge, width=colwidth[3])
        self.btncharge.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Solute spin multiplicity
        self.lblspin = Label(self.master, text="solute spin multiplicity", width=colwidth[0])
        self.lblspin.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        
        self.txtspin = Entry(self.master)
        defaultSpin = 1
        self.txtspin.insert(0, str(defaultSpin))
        self.txtspin.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        
        self.spin_multiplicity = IntVar(value=defaultSpin)
        
        def set_spin_multiplicity():
            answer = 0
            try: 
                answer = int(self.txtspin.get())
            except ValueError:
                answer = simpledialog.askinteger(parent=self.master,
                                     title="Dialog",
                                     prompt="Spin multiplicity must be a positive integer! Please re-enter.")
                self.txtspin.delete(0,END)
                self.txtspin.insert(0,answer)
            self.spin_multiplicity.set(answer)
            res = self.spin_multiplicity.get()
        
        self.btnspin = Button(self.master, text="Set", command=set_spin_multiplicity, width=colwidth[3])
        self.btnspin.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)

        self.irow += 1

        # DFT functional for QMMM 
        self.lblfunc = Label(self.master, text="Select QM method", width=colwidth[0])
        self.lblfunc.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        self.functxt = StringVar()
        self.func_chosen = Combobox(self.master, textvariable=self.functxt, width=colwidth[3])
        self.func_chosen['values'] = ("b3lyp",
                                      "hf",
                                      "case",
                                      "dftb",
                                      "svwn",
                                      "blyp",
                                      "bhandhlyp",
                                      "pbe",
                                      "revpbe",
                                      "pbe0",
                                      "revpbe0",
                                      "wpbe",
                                      "wpbeh",
                                      "bop",
                                      "mubop",
                                      "camb3lyp",
                                      "b97",
                                      "wb97",
                                      "wb97x",
                                      "wb97xd3")
        
        self.func_chosen.current(0)
        self.func_chosen.grid(column=1, row=self.irow,columnspan=3,sticky=W+E)
        self.irow += 1

        ### Set number of QM/MM minimization steps
        self.lblQMMMminSteps = Label(self.master, text="QM/MM minimization steps", width=colwidth[0],foreground='blue')
        self.lblQMMMminSteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtQMMMminSteps = Entry(self.master)
        defaultQMMMminSteps = 1000
        self.txtQMMMminSteps.insert(0,str(defaultQMMMminSteps))
        self.txtQMMMminSteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.QMMMminSteps = IntVar(value=defaultQMMMminSteps)
        
        def set_QMMMminSteps():
            answer = defaultQMMMminSteps
            try: 
                answer = int(self.txtQMMMminSteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "QM/MM minimization steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
        
            self.QMMMminSteps.set(answer)
            self.txtQMMMminSteps.delete(0,END)
            self.txtQMMMminSteps.insert(0,answer)
        
        style = Style()
        self.btnQMMMminSteps = Button(self.master, text="Set", command=set_QMMMminSteps, width=colwidth[3])
        self.btnQMMMminSteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### Set number of QM/MM heatup steps
        self.lblQMMMheatSteps = Label(self.master, text="QM/MM heatup steps", width=colwidth[0])
        self.lblQMMMheatSteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtQMMMheatSteps = Entry(self.master)
        defaultQMMMheatSteps = 1000
        self.txtQMMMheatSteps.insert(0,str(defaultQMMMheatSteps))
        self.txtQMMMheatSteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.QMMMheatSteps = IntVar(value=defaultQMMMheatSteps)
        
        def set_QMMMheatSteps():
            answer = defaultQMMMheatSteps
            try: 
                answer = int(self.txtQMMMheatSteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "QM/MM minimization steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
        
            self.QMMMheatSteps.set(answer)
            self.txtQMMMheatSteps.delete(0,END)
            self.txtQMMMheatSteps.insert(0,answer)
        
        style = Style()
        self.btnQMMMheatSteps = Button(self.master, text="Set", command=set_QMMMheatSteps, width=colwidth[3])
        self.btnQMMMheatSteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### Set number of QM/MM NVT equlibration steps
        self.lblQMMMeqNVTSteps = Label(self.master, text="QM/MM NVT run steps", width=colwidth[0])
        self.lblQMMMeqNVTSteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtQMMMeqNVTSteps = Entry(self.master)
        defaultQMMMeqNVTSteps = 1000
        self.txtQMMMeqNVTSteps.insert(0,str(defaultQMMMeqNVTSteps))
        self.txtQMMMeqNVTSteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.QMMMeqNVTSteps = IntVar(value=defaultQMMMeqNVTSteps)
        
        def set_QMMMeqNVTSteps():
            answer = defaultQMMMeqNVTSteps
            try: 
                answer = int(self.txtQMMMeqNVTSteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "QM/MM minimization steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
        
            self.QMMMeqNVTSteps.set(answer)
            self.txtQMMMeqNVTSteps.delete(0,END)
            self.txtQMMMeqNVTSteps.insert(0,answer)
        
        self.btnQMMMeqNVTSteps = Button(self.master, text="Set", command=set_QMMMeqNVTSteps, width=colwidth[3])
        self.btnQMMMeqNVTSteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### set number of steps for QMMM NVE production run steps 


        self.lblQMMMNVESteps = Label(self.master, text="QMMM NVE steps", width=colwidth[0],foreground='blue')
        self.lblQMMMNVESteps.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtQMMMNVESteps = Entry(self.master)
        defaultQMMMNVESteps = 0
        self.txtQMMMNVESteps.insert(0,str(defaultQMMMNVESteps))
        self.txtQMMMNVESteps.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.QMMMNVESteps = IntVar(value=defaultQMMMNVESteps)
        
        def set_QMMMNVESteps():
            answer = defaultQMMMNVESteps
            try: 
                answer = int(self.txtQMMMNVESteps.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "MM NVE steps must be a positive integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askinteger(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
        
            self.QMMMNVESteps.set(answer)
            self.txtQMMMNVESteps.delete(0,END)
            self.txtQMMMNVESteps.insert(0,answer)
        
        self.btnQMMMNVESteps = Button(self.master, text="Set", command=set_QMMMNVESteps, width=colwidth[3])
        self.btnQMMMNVESteps.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### Job control block starts
        self.lblJobControl = Label(self.master, text="Simulation job control options", font='Helvetica 18 bold')
        self.lblJobControl.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        self.irow += 1

        ### Use srun or not
        self.srunuse = BooleanVar(value=False)
        
        self.lblsrunuse = Label(self.master, text="Use srun to execute commands ?", width=colwidth[4])
        self.lblsrunuse.grid(column=0, row=self.irow,  columnspan=2, sticky=W, padx=self.padx)
        
        self.radsrunyes = Radiobutton(self.master, text='Yes', value=True, variable=self.srunuse, width=colwidth[3])
        self.radsrunyes.grid(column=2, row=self.irow)
        
        self.radsrunno = Radiobutton(self.master, text='No', value=False, variable=self.srunuse)
        self.radsrunno.grid(column=3, row=self.irow)
        
        self.irow += 1


        ### Use GPU accelerated MM or not
        self.pmemduse = BooleanVar(value=False)
        
        self.lblpmemduse = Label(self.master, text="Use GPU accelerated Amber MM (PMEMD instead of Sander) ?", width=colwidth[4])
        self.lblpmemduse.grid(column=0, row=self.irow, columnspan=2, sticky=W, padx=self.padx)
        
        self.radpmemdyes = Radiobutton(self.master, text='Yes', value=True, variable=self.pmemduse, width=colwidth[3])
        self.radpmemdyes.grid(column=2, row=self.irow)
        
        self.radpmemdno = Radiobutton(self.master, text='No', value=False, variable=self.pmemduse)
        self.radpmemdno.grid(column=3, row=self.irow)
        
        self.irow += 1

        ### Dry run mode or not
        self.dryrun = BooleanVar(value=True)
        
        dryrunmsg = "Dry run mode?\n(Only generate the commands "
        dryrunmsg += "to run MD programs without executing.)"
        self.lbldryrun = Label(self.master, text=dryrunmsg, width=colwidth[4])
        self.lbldryrun.grid(column=0, row=self.irow,  columnspan=2,  sticky=W, padx=self.padx)
        
        self.raddryrunyes = Radiobutton(self.master, text='Yes', value=True, variable=self.dryrun, width=colwidth[3])
        self.raddryrunyes.grid(column=2, row=self.irow)
        
        self.raddryrunno = Radiobutton(self.master, text='No', value=False, variable=self.dryrun)
        self.raddryrunno.grid(column=3, row=self.irow)
        
        self.irow += 1

        def GUI_input_sanity_check():
            mdrun_error =0 
            print("### Start MDrun input option fact check ###")
            if self.prefix.get() =="":
                  print("Input file (.prmtop) prefix must be provided!")
                  mdrun_error = 1
            else:
                  print("Input file prefix: ", self.prefix.get())
            # Check temperature
            if self.Temp.get() <=0:
                print("Temperature specified: {:.2f} (K) is not valid\n".format(self.Temp.get() )) 
                mdrun_error = 2
            if self.Pressure.get() <=0:
                print("Pressure specified: {:.4f} (bar) is not valid\n".format(self.Pressure.get()))
                mdrun_error = 3 
            return mdrun_error

        def write_mdrun_input():
            cmd = "autosolvate mdrun"
        
            cmd += " -f " + self.prefix.get()
            cmd += " -t {:.2f}".format(self.Temp.get())
            cmd += " -p {:.4f}".format(self.Pressure.get())
            #cmd += " -m {:d} ".format(self.MMMinSteps.get())
            cmd += " -m {:d} ".format(self.MMHeatSteps.get())
            cmd += " -n {:d} ".format(self.MMNPTSteps.get())
            #cmd += " -e {:d} ".format(self.MMNVESteps.get())
            if self.doQMMM.get() == True:
                 cmd += " -l {:d} ".format(self.QMMMminSteps.get())
                 cmd += " -o {:d} ".format(self.QMMMheatSteps.get())
                 cmd += " -s {:d} ".format(self.QMMMeqNVTSteps.get())
                 #cmd += " -s {:d} ".format(self.QMMMNVESteps.get())
                 cmd += " -q {:d} ".format(self.charge.get())
                 cmd += " -u {:d} ".format(self.spin_multiplicity.get())
                 cmd += " -k {:s} ".format(self.func_chosen.get())
            if self.srunuse.get() == True:
                cmd += "-r"
            if self.pmemduse.get() == True:
                cmd += "-x"
            if self.dryrun.get() == True:
                cmd += "-d"
            return cmd

        ### Execute  python command to generate MD input files and jobscripts
        def execute():
            mdrun_error = GUI_input_sanity_check()
            if mdrun_error == 0:
                cmd = write_mdrun_input()
                res = "Congratulations! MDrun command line generated: \n" + cmd
                messagebox.showinfo(title="Confirmation", message=res, icon='info')
                question = "Do you want to continue to generate the "\
                         + "MD input and job scripts?"
                answer = messagebox.askyesno(title="Confirmation", message=question, icon='question')
                if answer == True:
                    subprocess.call(cmd, shell=True)
                    if self.dryrun.get() == True:
                       res = "MD input and job script generation finished! \n"
                       res += "Please check the job scripts in runMM.sh and runQMMM.sh\n"
                    else:
                       res = "MD input generation and simulation finished! \n"
                    res += "Thanks for using AutoSolvate. Have a nice day!\n"
                    messagebox.showinfo(title="Success", message=res, icon='info')

        
            else:
                res = "Error found in input.\n"
                res += "Please check the input setting and retry\n"
                messagebox.showerror(message=res)
        
        #### Generate MD simulation input files and job scripts
        self.btnGo = Button(self.master, text="Generate MD input files and job scripts! ", command=execute)
        self.btnGo.grid(column=0, row=self.irow, columnspan=3, sticky=W+E, padx=self.padx, pady=5)
        self.irow += 1
### END MD automation window ###

### START cluster extraction window ###
class clusterGUI(baseGUI):
    def __init__(self,master):
        super().__init__(master)
        master.title("Microsolvated cluster extraction")
        master.geometry('820x280')
        self.display_logo()
        ### Enter .prmtop filename 
        self.lblParm = Label(self.master, text="Path for .prmtop file", width=colwidth[0])
        self.lblParm.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtParm = Entry(self.master)
        self.txtParm.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.prmtop = StringVar()
        
        def set_prmtop():
            mypath = self.txtParm.get()
            my_filetypes = [('Amber prmtop', '.prmtop'),
                            ('all files', '.*')]
            if mypath !="" and os.path.exists(mypath) and 'prmtop' in mypath:
                    self.prmtop.set(self.txtParm.get())
            else:
                msg = "Invalid prmtop file path provided!\nPlease select a path:"
                answer = filedialog.askopenfilename(parent=self.master,
                                                    initialdir=os.getcwd(),
                                                    title="Request", 
                                                    message=msg,
                                                    filetypes=my_filetypes)
                self.prmtop.set(answer)
                self.txtParm.delete(0,END)
                self.txtParm.insert(0,answer)
        
        self.btnParm = Button(self.master, text="Set", command=set_prmtop, width=colwidth[3])
        self.btnParm.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1
      
        ### Enter .netcdf filename 
        self.lblMdcrd = Label(self.master, text="Path for .netcdf file", width=colwidth[0])
        self.lblMdcrd.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtMdcrd = Entry(self.master)
        self.txtMdcrd.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.netcdf = StringVar()
        
        def set_netcdf():
            mypath = self.txtMdcrd.get()
            my_filetypes = [('Amber .netcdf trajectory file', '.netcdf'),
                            ('all files', '.*')]
            if mypath !="" and os.path.exists(mypath) and 'netcdf' in mypath:
                    self.netcdf.set(self.txtMdcrd.get())
            else:
                msg = "Invalid netcdf file path provided!\nPlease select a path:"
                answer = filedialog.askopenfilename(parent=self.master,
                                                    initialdir=os.getcwd(),
                                                    title="Request", 
                                                    message=msg,
                                                    filetypes=my_filetypes)
                self.netcdf.set(answer)
                self.txtMdcrd.delete(0,END)
                self.txtMdcrd.insert(0,answer)
        
        self.btnMdcrd = Button(self.master, text="Set", command=set_netcdf, width=colwidth[3])
        self.btnMdcrd.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1
        ### set first frame to extract 
        self.lblStartFrame = Label(self.master, text="ID of first frame to extract (ID starts from 0)", width=colwidth[0])
        self.lblStartFrame.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtStartFrame = Entry(self.master)
        defaultStartFrame = 0
        self.txtStartFrame.insert(0,str(defaultStartFrame))
        self.txtStartFrame.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.StartFrame = IntVar(value=defaultStartFrame)
        
        def set_StartFrame():
            answer = defaultStartFrame
            try: 
                answer = int(self.txtStartFrame.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "StartFrame must be a non-negative integer!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askint(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
      
            self.StartFrame.set(answer)
            self.txtStartFrame.delete(0,END)
            self.txtStartFrame.insert(0,answer)
        
        self.btnStartFrame = Button(self.master, text="Set", command=set_StartFrame, width=colwidth[3])
        self.btnStartFrame.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1
        ### set Extraction Interval
        self.lblInterval = Label(self.master, text="Cluster extraction interval (steps)", width=colwidth[0])
        self.lblInterval.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtInterval = Entry(self.master)
        defaultInterval = 1
        self.txtInterval.insert(0,str(defaultInterval))
        self.txtInterval.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.Interval = IntVar(value=defaultInterval)
        
        def set_Interval():
            answer = defaultInterval
            try: 
                answer = int(self.txtInterval.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "Interval must be a positive int!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askint(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = int(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False
      
            self.Interval.set(answer)
            self.txtInterval.delete(0,END)
            self.txtInterval.insert(0,answer)
        
        self.btnInterval = Button(self.master, text="Set", command=set_Interval, width=colwidth[3])
        self.btnInterval.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### set shell size
        self.lblShellSize = Label(self.master, text="Shell thickness (Angstrom)", width=colwidth[0])
        self.lblShellSize.grid(column=0, row=self.irow, sticky=W, padx=self.padx)
        
        self.txtShellSize = Entry(self.master)
        defaultShellSize = 4.0
        self.txtShellSize.insert(0,str(defaultShellSize))
        self.txtShellSize.grid(column=1, row=self.irow, columnspan=3, sticky=W+E)
        self.ShellSize = DoubleVar(value=defaultShellSize)
        
        def set_ShellSize():
            answer = defaultShellSize
            try: 
                answer = float(self.txtShellSize.get())
            except ValueError:
                anserValid = False
                while (not answerValid):
                    msg = "ShellSize (in Angstrom) must be a positive float!\n"
                    msg += "Please re-enter.\n"
                    answer = simpledialog.askfloat(parent=self.master,
                                               title="Dialog",
                                               prompt=msg)
                    answerValid = True
                    try: 
                        answer = float(answer)
                    except ValueError:
                        answerValid = False
                    if answerValid:
                        if answer <= 0 :
                            answerValid = False

            self.ShellSize.set(answer)
            self.txtShellSize.delete(0,END)
            self.txtShellSize.insert(0,answer)
        
        self.btnShellSize = Button(self.master, text="Set", command=set_ShellSize, width=colwidth[3])
        self.btnShellSize.grid(column=4, row=self.irow,columnspan=3,sticky=W+E)
        
        self.irow += 1

        ### Use srun or not
        self.srunuse = BooleanVar(value=False)
        
        self.lblsrunuse = Label(self.master, text="Use srun to execute commands ?", width=colwidth[4])
        self.lblsrunuse.grid(column=0, row=self.irow,  columnspan=2, sticky=W, padx=self.padx)
        
        self.radsrunyes = Radiobutton(self.master, text='Yes', value=True, variable=self.srunuse, width=colwidth[3])
        self.radsrunyes.grid(column=2, row=self.irow)
        
        self.radsrunno = Radiobutton(self.master, text='No', value=False, variable=self.srunuse)
        self.radsrunno.grid(column=3, row=self.irow)
        
        self.irow += 1

        def GUI_input_sanity_check():
            clusterrun_error =0 
            print("### Start ClusterGen input option fact check ###")
            if self.prmtop.get() =="":
                  print("Input prmtop file must be provided!")
                  clusterrun_error = 1
            else:
                  print("Input prmtop: ", self.prmtop.get())
            if self.netcdf.get() =="":
                  print("Input netcdf file must be provided!")
                  clusterrun_error = 2
            else:
                  print("Input netcdf: ", self.netcdf.get())
            if self.ShellSize.get() <=0:
                print("Solvent shell thickness specified: {:.2f} (Angstrom) is not valid\n".format(self.ShellSize.get() )) 
                clusterrun_error = 3
            return clusterrun_error

        def write_clusterrun_input():
            cmd = "autosolvate clusterrun"
        
            cmd += " -f " + self.prmtop.get()
            cmd += " -t " + self.netcdf.get()
            cmd += " -a {:d}".format(self.StartFrame.get())
            cmd += " -i {:d}".format(self.Interval.get())
            cmd += " -s {:.4f}".format(self.ShellSize.get())
            if self.srunuse.get() == True:
                cmd += "-r"
            return cmd

        def execute():
            clusterrun_error = GUI_input_sanity_check()
            if clusterrun_error == 0:
                cmd = write_clusterrun_input()
                res = "Congratulations! ClusterGen command line generated: \n" + cmd
                messagebox.showinfo(title="Confirmation", message=res)
                question = "Do you want to continue to generate the "\
                         + "solvated cluster trajectory file?"
                answer = messagebox.askyesno(title="Confirmation", message=question)
                if answer == True:
                    subprocess.call(cmd, shell=True)
                    res = "Cluster extraction finished!\n"
                    res += "Please check the xyz file(s) generated in the current folder.\n"
                    res += "Thanks for using AutoSolvate. Have a nice day!\n"
                    messagebox.showinfo(title="Success", message=res, icon='info')
        
            else:
                res = "Error found in input.\n"
                res += "Please check the input setting and retry\n"
                messagebox.showerror(message=res)

        #### Do cluster extraction
        self.btnGo = Button(self.master, text="Generate solvated cluster! ", command=execute)
        self.btnGo.grid(column=0, row=self.irow, columnspan=3, sticky=W+E, padx=self.padx, pady=5)
        self.irow += 1
### END Cluster extraction window ###

## The master GUI of AutoSolvate where we select what task to do ##
class autosolvateGUI(baseGUI):
    def __init__(self,master):
        super().__init__(master)
        self.master.title("Welcome to AutoSolvate!")
        self.master.geometry('360x180')

        # Display logo
        self.display_logo()

        ### select the task to do
        self.lbl00 = Label(master, text="Please select the task",font='Helvetica 20 bold',width=20)
        self.lbl00.grid(column=0, row=self.irow,  sticky=W+E, padx=self.padx)
        self.lbl00.configure(anchor="center")
        self.irow += 1
        # Adding combobox drop down list 
        self.task = StringVar()
        self.task_chosen = Combobox(master, textvariable=self.task, width=30)
        self.task_chosen['values'] = ('Solvated box and MD parameter generation',
                                      'MD automation',
                                     'Microsolvated cluster extraction') 
        self.task_chosen.current(0)
        self.task_chosen.grid(column=0, row=self.irow,  sticky=W+E, padx=self.padx)
        self.irow += 1

        # Create new window to do the task
        def create_task_window():
            if self.task_chosen.get() == 'Solvated box and MD parameter generation':
                self.master2 = Toplevel(self.master)
                my_gui = boxgenGUI(self.master2)
            elif self.task_chosen.get() == 'MD automation':
                self.master3 = Toplevel(self.master)
                my_gui = mdGUI(self.master3)
            else:
                self.master4 = Toplevel(self.master)
                my_gui = clusterGUI(self.master4)
        
        self.btn02 = Button(master, text="Go!", command=create_task_window, width=20)
        self.btn02.grid(column=0, row=self.irow,  sticky=W+E, padx=self.padx)
        self.irow += 1

# This part does not need to be modified
if __name__ == '__main__':
    window = Tk()
    my_gui = autosolvateGUI(window)
    window.mainloop()

    cleanUp()
