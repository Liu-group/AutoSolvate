######################################################################################################
# test example 4: Naphthalene Radical in Chloroform: box gen: test pdb, prmtop 
# author: Patrick Li (2022-10-18) 
# path: autosolvate/tests/test_case4.py 
######################################################################################################
# @Date   : May 9th 2024 
# @Author : Patrick Li
# @Note   : Disabled because this test case requires the gaussian to be installed 
######################################################################################################

# import os
# import pytest
# import numpy as np

# import autosolvate
# from . import helper_functions as hp


# def test_example_4(tmpdir):
#     testName = "test_example_4"
    
#     autosolvate.startboxgen([
#         "-m", hp.get_input_dir("naphthalene_radical.xyz"),
#         "-s", "chloroform",
#         "-c", "1", 
#         "-u", "2",
#         "-g", "resp", 
#         "-o", "nap_radical_chcl3",
#         ])

#     out = "nap_radical_chcl3" 
#     ref = hp.get_reference_dir("nap_radical_chcl3")

    
#     for suffix in [".pdb", ".prmtop", ".inpcrd"]:
#         assert os.path.exists(out + suffix)
    
#     assert hp.compare_pdb(out, ref)
#     assert hp.compare_inpcrd_prmtop(out, ref)

    
