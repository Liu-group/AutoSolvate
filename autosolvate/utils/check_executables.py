import os
import subprocess
import logging
from logging import DEBUG, INFO, WARN, WARNING, CRITICAL
logging.basicConfig(level = WARNING, force = True, handlers=[])

logger              = logging.getLogger(name = "CheckExecutable")
# output_handler      = logging.FileHandler(filename = "log.txt", mode = "a", encoding="utf-8")
output_handler      = logging.StreamHandler()
output_formater     = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
output_handler.setFormatter(output_formater)
if len(logger.handlers) == 0:
    logger.addHandler(output_handler)

ambertools = [
    "antechamber",
    "parmchk2",
    "tleap",
    "pdb4amber",
    "MCPB.py",
    "ForceBalance",
    "mdgx",
]

def get_packmol(packmol_exe = ""):
    if packmol_exe and os.path.exists(packmol_exe):
        return packmol_exe
    packmol_exe = "packmol"
    res = subprocess.getoutput("which {}".format(packmol_exe))
    if not res.find("no {} in".format(packmol_exe)):
        return packmol_exe
    logger.error("ERROR: cannot find packmol. The custom solvation functionality will be disabled!")
    return packmol_exe

def get_amberhome(amberhome = ""):
    if amberhome and os.path.exists(amberhome):
        return amberhome
    cmd = ["echo", "$AMBERHOME"]
    proc=subprocess.Popen(cmd,
                    universal_newlines=True,
                    stdout=subprocess.PIPE)
    output=proc.communicate()[0]
    if len(output)<2:
        logger.error("ERROR: AMBERHOME not defined")
        logger.error("ERROR: Please set up AMBERHOME environment or specify through command line option -a")
        logger.error("ERROR: Exiting...")
        raise ModuleNotFoundError("ERROR: AMBERHOME not defined. {}".format(output))
    else:
        logger.warning("WARNING: AMBERHOME detected: ", output)
        amberhome = output
    return amberhome

def get_gaussian(gaussian_exe = ""):
    if gaussian_exe and os.path.exists(gaussian_exe):
        return gaussian_exe
    logger.warning("WARNING: Gaussian executable name is not specified for RESP charge fitting!")
    logger.warning("WARNING: Using g16 by default. If failed later, please rerun with the option -e specified!")
    gaussian_exe = 'g16'
    return gaussian_exe

def get_gaussian_dir(gaussian_dir = ""):
    """
    Question:
    --------------
    What is gaussian dir for? Does it work in antechamber? Will this variable be discarded in future GaussianDocker?
    """
    if gaussian_dir and os.path.exists(gaussian_dir):
        return gaussian_dir
    logger.warning("WARNING: Gaussian executable directory is not specified for RESP charge fitting!")
    logger.warning("WARNING: Setting to default path: /opt/packages/gaussian/g16RevC.01/g16/")
    logger.warning("WARNING: If failed later, please rerun with the option -d specified!")
    return gaussian_dir

def get_terachem_module(terachem_package = ""):
    terachem_package = "TeraChem"
    traversed_path = [terachem_exe, ]
    if terachem_package and os.path.exists(terachem_package):
        return terachem_package
    
    terachem_exe = "terachem"
    traversed_path.append(terachem_exe)
    result = subprocess.getoutput("which terachem")
    if result.find("no terachem in") == -1:
        return terachem_package
    
    res = subprocess.getoutput("module load TeraChem")
    terachem_package = "TeraChem"
    if res.find("ERROR") == -1:
        return terachem_package
    
    tcmodules = subprocess.getoutput("module avail TeraChem")
    tcmodules = [s.replace("(default)", "", 1) for s in tcmodules.split() if s.find("TeraChem") != -1]
    for tcmodule in tcmodules:
        res = subprocess.getoutput("module load {tcmodule}".format(tcmodule = tcmodules[0]))
        terachem_package = tcmodule
        traversed_path.append(tcmodule)
        if res.find("ERROR") == -1:
            return terachem_package
    else:
        logger.critical("No terachem in: " + " ".join(traversed_path))
        raise ModuleNotFoundError("No terachem in: " + " ".join(traversed_path))

def get_terachem(terachem_exe = ""):
    traversed_path = [terachem_exe, ]
    if terachem_exe and os.path.exists(terachem_exe):
        return terachem_exe
    
    terachem_exe = "terachem"
    traversed_path.append(terachem_exe)
    result = subprocess.getoutput("which terachem")
    if result.find("no terachem in") == -1:
        return terachem_exe
    
    res = subprocess.getoutput("module load TeraChem")
    if res.find("ERROR") == -1:
        return terachem_exe
    
    tcmodules = subprocess.getoutput("module avail TeraChem")
    tcmodules = [s.replace("(default)", "", 1) for s in tcmodules.split() if s.find("TeraChem") != -1]
    for tcmodule in tcmodules:
        res = subprocess.getoutput("module load {tcmodule}".format(tcmodule = tcmodules[0]))
        traversed_path.append(tcmodule)
        if res.find("ERROR") == -1:
            return terachem_exe
    else:
        logger.critical("No terachem in: " + " ".join(traversed_path))
        raise ModuleNotFoundError("No terachem in: " + " ".join(traversed_path))

