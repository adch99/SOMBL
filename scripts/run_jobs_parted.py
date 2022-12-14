import os

size = 100
runs = 90
batchsize = 15
numbatches = 6
hopup = 1.0
hopdn = 1.0


coupLow = 0
coupHigh = 2
coupStep = 0.1

disLow = 8
disHigh = 18
disStep = 1

numCoups = 21
numDis = 11
numRuns = numCoups * numDis

start = 0
stop = 80

for index in range(start, stop):
    coupIdx = index // numDis
    disIdx = index % numDis
    coupling = coupLow + coupStep*coupIdx
    disorder = disLow + disStep*disIdx

    for batchnum in range(1, numbatches+1):
        jobfilename = "jobs/mbl_%dx%d_W%.1f_C%.1f_TU%.1f_TD%.1f_N%d_BS%d_B%d.pbs" % (size,size,disorder,coupling,hopup,hopdn,runs,batchsize,batchnum)
        os.system("echo %s" % jobfilename)
        # os.system("qsub %s" % jobfilename)
