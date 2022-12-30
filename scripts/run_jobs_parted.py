import os

size = 60
runs = 100
batchsize = 10
numbatches = 10
hopup = 1.0
hopdn = 1.0


coupLow = 0.3
coupHigh = 0.5
coupStep = 0.1

disLow = 13
disHigh = 15
disStep = 1

numCoups = 3
numDis = 3
numRuns = numCoups * numDis

start = 0
stop = 9

basename = "mbl_average"

for index in range(start, stop):
    coupIdx = index // numDis
    disIdx = index % numDis
    coupling = coupLow + coupStep*coupIdx
    disorder = disLow + disStep*disIdx

    for batchnum in range(1, numbatches+1):
        subs = (basename, size, size,
                disorder, coupling,
                hopup, hopdn, runs,
                batchsize, batchnum)
        jobfilename = "jobs/%s_%dx%d_W%.1f_C%.1f_TU%.1f_TD%.1f_N%d_BS%d_B%d.pbs" % subs
        os.system("echo %s" % jobfilename)
        # os.system("qsub %s" % jobfilename)
