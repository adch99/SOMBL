import numpy as np


def alt_up_empty(length):
    xUp = np.arange(0, length, 2)
    yUp = np.arange(0, length, 2)
    upList = xUp + length * yUp[:, np.newaxis]
    downList = None
    return upList.flatten(), downList


def alt_up_down():
    pass


def alt_updown_empty():
    pass


def full_updown():
    pass


def outputPattern(upList, downList, filename):
    if upList is None:
        upNum = 0
        upListStr = "\n"
    else:
        upNum = len(upList)
        upListStr = ",".join([str(x) for x in sorted(upList)]) + "\n"

    if downList is None:
        downNum = 0
        downListStr = "\n"
    else:
        downNum = len(downList)
        downListStr = ",".join([str(x) for x in sorted(downList)]) + "\n"

    with open(filename, "w") as ofile:
        ofile.write("#NUMUP\n")
        ofile.write(str(upNum) + "\n")
        ofile.write("#NUMDN\n")
        ofile.write(str(downNum) + "\n")
        ofile.write("#SPINSUP\n")
        ofile.write(upListStr)
        ofile.write("#SPINSDN")
        ofile.write(downListStr)


def main():
    length = 3
    filename = f"data/alt_up_empty_L{length}.dat"
    upList, downList = alt_up_empty(length)
    outputPattern(upList, downList, filename)


if __name__ == "__main__":
    main()
