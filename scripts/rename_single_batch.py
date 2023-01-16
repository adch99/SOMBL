import os

filelist = os.glob.glob("data/mbl_40x40_*_TU1_TD1_N100_BS100_B1_grgrstar_*.dat)
for file in filelist:
    parts = file.split("_")
    newparts = []
    for part in parts:
        if part == "BS100" or part == "B1":
            continue
        else:
            newparts.append(part)
    newname = "_".join(newparts)
    print "Old:", file,"\nNew:", newname, "\n\n"
