from config import useQGTagging
import os

lines = ""
folder = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(folder, "Category.cc")) as f:
    f = iter(f)
    for line in f:
        if "float WP" in line:
            break
    for a in 1, 2:
        for line in f:
            lines += line
            if "}" in line:
                break
lines = lines.replace("{", ":").replace("}", "").replace("  if", "if").replace("  else", "else")

def WP_VBF2j(ZZMass=125):
    exec lines
    return WP_VBF2j
def WP_VBF1j(ZZMass=125):
    exec lines
    return WP_VBF1j
def WP_ZHh(ZZMass=125):
    exec lines
    return WP_ZHh
def WP_WHh(ZZMass=125):
    exec lines
    return WP_WHh
