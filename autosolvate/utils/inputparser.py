import os 
import re
import json
import shutil 



def parse_json_input(jsonpath:str):
    with open(jsonpath, "r") as f:
        data = json.load(f)
    