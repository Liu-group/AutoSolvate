import autosolvate 
import os 

autosolvate.startboxgen(["-m", "naphthalene_neutral.xyz",
                            "-s", "acetonitrile",
                            "-c", 0,
                            "-u", 1,
                            "-g", "bcc",
                            "-o", "nap_neutral_MeCN"])

