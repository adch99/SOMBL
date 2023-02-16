import os


param_vals = [
    (0.1, 9, [2, 10]),
    (0.1, 10, range(1, 11)),
    (0.1, 11, range(1, 11)),
    (0.1, 12, [1, 3, 4, 5, 6, 7, 8, 9, 10]),
    (0.1, 13, [1, 2, 3, 4]),
    (0.4, 9, [6]),
    (0.4, 10, [6, 7]),
    (0.4, 11, [1, 3]),
    (0.4, 13, [1, 9]),
    (0.4, 15, [2, 8, 9]),
    (0.5, 13, [8]),
    (1.8, 11, range(1, 11)),
    (1.8, 12, [1, 2, 3, 4, 5]),
    (1.8, 13, range(1, 11)),
    (1.8, 14, [1, 2, 3, 4, 5, 8, 9, 10]),
    (1.8, 15, range(1, 11)),
    (1.8, 16, [1, 2, 6, 7, 8, 9, 10]),
    (1.8, 17, range(1, 11)),
    (1.8, 18, [2, 3, 4, 5, 6, 7, 8, 9, 10]),
    (1.9, 8, [1, 2, 3, 4, 5, 6, 9, 10])
]


size = 100
runs = 100
batchsize = 10
numbatches = 10
hopup = 1.0
hopdn = 1.0


basename = "mbl"

for index in range(len(param_vals)):
    coupling = param_vals[index][0]
    disorder = param_vals[index][1]

    for batchnum in param_vals[index][2]:
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
