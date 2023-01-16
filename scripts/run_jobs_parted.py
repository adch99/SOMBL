import os

size = 100
runs = 100
batchsize = 10
numbatches = 10
hopup = 1.0
hopdn = 1.0


coupLow = 0.0
coupHigh = 2.0
coupStep = 0.1

disLow = 8
disHigh = 18
disStep = 1

numCoups = 21
numDis = 11
numRuns = numCoups * numDis

start = 0
stop = 22

basename = "mbl"

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

#     subs = (basename, size, size,
#             disorder, coupling,
#             hopup, hopdn, runs)
#     jobfilename = "jobs/%s_%dx%d_W%.1f_C%.1f_TU%.1f_TD%.1f_N%d.pbs" % subs
#     os.system("echo %s" % jobfilename)
    # os.system("qsub %s" % jobfilename)
