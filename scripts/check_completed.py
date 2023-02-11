import os
import subprocess
import sys

FLOAT_SIZE = 781256
DOUBLE_SIZE = 2*FLOAT_SIZE

def get_last_run(logfilename):
    cmd = f'grep -E "Run [0-9]{{1,3}} done" {logfilename} | tail -n 1 | cut -d" " -f 2'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    last_run = int(result.stdout[:-1])
    return last_run

def get_status(logfilename):
    cmd = f'grep -i -e "kill" -e "var" -e "PBS" -e "cannot" {logfilename}'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    output = result.stdout[:-1]
    if output == "":
        return "success"
    else:
        return "killed"

def get_filetype(datafilename):
    cmd = f'file -b {datafilename}'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    output = result.stdout[:-1]
    if "ASCII" in output:
        return "text"
    elif output == "data":
        return "binary"
    else:
        print(f"File {datafilename} has unknown type {output}", file=sys.stderr)
        # return output
        # Assume it's just binary
        return "binary"

def get_filesize(datafilename):
    cmd = f'du {datafilename} | cut -f 1'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    size = int(result.stdout[:-1])
    return size



def get_logfilename(c, w, batchsize, batch, numruns):
    # Ex: logs/mbl_100x100_W9.0_C1.6_TU1.0_TD1.0_N100_BS10_B9.log
    name = f"logs/mbl_100x100_W{w:.1f}_C{c:.1f}_TU1.0_TD1.0_N{numruns}_BS{batchsize}_B{batch}.log"
    return name

def get_datafilename(c, w, batchsize, batch, numruns):
    # Ex: data/mbl_100x100_W9_C1.6_TU1_TD1_N100_BS10_B9_grgrstar_a0_b0_full.dat
    name = f"data/mbl_100x100_W{w}_C{c:g}_TU1_TD1_N{numruns}_BS{batchsize}_B{batch}_grgrstar_a0_b0_full.dat"
    return name

def check_completed(logfilename, batch, batchsize):
    stop = get_last_run(logfilename) # get the last run
    if stop == batch*batchsize:
        return True
    else:
        return False

def check_type(datafilename):
    filetype = get_filetype(datafilename)
    if filetype == "binary":
        filesize = get_filesize(datafilename)
        if filesize == FLOAT_SIZE:
            return ("binary", "float")
        elif filesize == DOUBLE_SIZE:
            return ("binary", "double")
        else:
            print(f"Can't figure out type for {datafilename}", file=sys.stderr)
            print(f"Size: {filesize}, Type: {filetype}", file=sys.stderr)
            return ("binary", "")
    elif filetype == "text":
        return ("text", "double")

    else:
        print(f"Can't figure out type for {datafilename}", file=sys.stderr)
        print(f"Type: {filetype}", file=sys.stderr)
        return ("", "")

def print_csv(datarow):
    varlist = [datarow["c"], datarow["w"], datarow["batch"],
            datarow["batchsize"], datarow["numruns"],
            datarow["status"], datarow["filetype"],
            datarow["size"]]
    row = ",".join([str(var) for var in varlist])
    print(row)

def run(c, w, batchsize, batch, numruns):
    datafilename = get_datafilename(c, w, batchsize, batch, numruns)
    logfilename = get_logfilename(c, w, batchsize, batch, numruns)
    # print(datafilename)
    # print(logfilename)

    datarow = {
        "c": c,
        "w": w,
        "batch": batch,
        "batchsize": batchsize,
        "numruns": numruns
    }
    # Check if the files exist first
    logexists = os.path.isfile(logfilename)
    dataexists = os.path.isfile(datafilename)

    if dataexists and logexists:
        complete = check_completed(logfilename, batch, batchsize)
        ftype, dtype = check_type(datafilename)
        datarow["status"] = str(complete)
        datarow["filetype"] = ftype
        datarow["size"] = get_filesize(datafilename)
        # datarow["dtype"] = dtype
        print_csv(datarow)
    
    elif dataexists and (not logexists):
        ftype, dtype = check_type(datafilename)
        datarow["status"] = "na"
        datarow["filetype"] = ftype
        # data["dtype"] = dtype
        datarow["size"] = get_filesize(datafilename)
        print_csv(datarow)

    elif (not dataexists) and (not logexists):
        datarow["status"] = "na"
        datarow["filetype"] = "na"
        # datarow["dtype"] = "na"
        datarow["size"] = "na"
        
    else:
        print("Something is wrong: datafile doesn't exist but logfile does.", file=sys.stderr)
        print(f"Datafile: {datafilename}", file=sys.stderr)
        print(f"Logfile: {logfilename}", file=sys.stderr)




def main():
    # c = 1.6
    # w = 9
    batchsize = 10
    batch = 9
    numruns = 100
    crange = [x*0.1 for x in range(20)]
    for c in crange:
        for w in range(8, 18):
            for batch in range(1, batchsize+1):
                run(c, w, batchsize, batch, numruns)


if __name__ == "__main__":
    main()