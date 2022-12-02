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
stop = 480

for index in range(start, stop):
    coupIdx = index // numDis
    disIdx = index % numDis
    coupling = coupLow + coupStep*coupIdx
    disorder = disLow + disStep*disIdx

    for batchnum in range(1, numbatches+1):
        jobfilename = f"jobs/mbl_{size}x{size}_W{disorder:.1f}_C{coupling:.1f}_TU{hopup:.1f}_TD{hopdn:.1f}_N{runs}_BS{batchsize}_B{batchnum}.pbs"
        os.system(f"qsub {jobfilename}")
