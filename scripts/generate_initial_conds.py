import numpy as np


def get_index(x, y, length):
    return x + y * length


def alt_up_empty(length):
    xUp = np.arange(0, length, 2)
    yUp = np.arange(0, length, 2)
    upList = xUp + length * yUp[:, np.newaxis]
    downList = None
    return upList.flatten(), downList


def alt_up_down(length):
    upList = []
    downList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
            else:
                downList.append(get_index(x, y, length))
    return np.array(upList), np.array(downList)


def alt_updown_empty(length):
    upList = []
    downList = []
    for x in range(0, length, 2):
        for y in range(0, length, 2):
            upList.append(get_index(x, y, length))
            downList.append(get_index(x, y, length))
    return np.array(upList), np.array(downList)


def altn_up_empty(length):
    upList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
    return upList, None

def altn_updown_empty(length):
    upList = []
    downList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
                downList.append(get_index(x, y, length))
    return np.array(upList), np.array(downList)

def altn_altupdown_updown(length):
    upList = []
    downList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
                downList.append(get_index(x, y, length))
            elif ((x + y) % 4 == 1):
                upList.append(get_index(x, y, length))
            else:
                downList.append(get_index(x, y, length))

    return np.array(upList), np.array(downList)


def altn_altupdown_empty(length):
    upList = []
    downList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                continue
            elif ((x + y) % 4 == 1):
                upList.append(get_index(x, y, length))
            else:
                downList.append(get_index(x, y, length))

    return np.array(upList), np.array(downList)



def altn_random_updown(length):
    upList = []
    downList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
                downList.append(get_index(x, y, length))
            else:
                if np.random.randint(2) == 0:
                    upList.append(get_index(x, y, length))
                else:
                    downList.append(get_index(x, y, length))

    return np.array(upList), np.array(downList)


def altn_randomequal_updown(length):
    upList = []
    downList = []
    rng = np.random.default_rng()
    subLatSize = length*length//2
    choices = rng.choice(np.arange(subLatSize), subLatSize//2, replace=False)
    isUp = np.full(subLatSize, False, dtype=bool)
    isUp[choices] = True
    ctr = 0
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
                downList.append(get_index(x, y, length))
            else:
                if isUp[ctr]:
                    upList.append(get_index(x, y, length))
                else:
                    downList.append(get_index(x, y, length))
                ctr += 1

    return np.array(upList), np.array(downList)



def altn_random_empty(length):
    upList = []
    downList = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                continue
            else:
                if np.random.randint(2) == 0:
                    upList.append(get_index(x, y, length))
                else:
                    downList.append(get_index(x, y, length))

    return np.array(upList), np.array(downList)

def altn_randomequal_empty(length):
    upList = []
    downList = []
    rng = np.random.default_rng()
    subLatSize = length*length//2
    choices = rng.choice(np.arange(subLatSize), subLatSize//2, replace=False)
    isUp = np.full(subLatSize, False, dtype=bool)
    isUp[choices] = True
    ctr = 0
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                continue
            else:
                if isUp[ctr]:
                    upList.append(get_index(x, y, length))
                else:
                    downList.append(get_index(x, y, length))
                ctr += 1

    return np.array(upList), np.array(downList)



def pnjunction(length):
    upList = []
    downList = []
    for x in range(length//2):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
                downList.append(get_index(x, y, length))
            else:
                upList.append(get_index(x, y, length))

    for x in range(length//2, length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                upList.append(get_index(x, y, length))
                downList.append(get_index(x, y, length))
            else:
                downList.append(get_index(x, y, length))

    return np.array(upList), np.array(downList)





def full_updown():
    pass


def single_up(length):
    index = get_index(length//2, length//2, length)
    return np.array([index]), []


def single_updown(length):
    index = get_index(length//2, length//2, length)
    return np.array([index]), np.array([index])


def adj_up(length):
    index1 = get_index(length//2, length//2, length)
    index2 = get_index(length//2 + 1, length//2, length)
    return np.array([index1, index2]), []


def adj_updown(length):
    index1 = get_index(length//2, length//2, length)
    index2 = get_index(length//2 + 1, length//2, length)
    return np.array([index1, index2]), np.array([index1, index2])


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
        ofile.write("#SPINSDN\n")
        ofile.write(downListStr)



def main():
    length = 60
    # filename = f"data/pnjunction_L{length}.dat"
    # filename = f"data/testing_sample_initial_cond_L{length}.dat"

    filename = f"data/alt_up_down_L{length}.dat"
    upList, downList = alt_up_down(length)
    # upList = [2, 3, 5, 9, 11, 12, 15]
    # downList = [1, 2, 4, 5, 7, 8, 11, 14, 15]
    print(f"Up: {upList} Down: {downList}")
    outputPattern(upList, downList, filename)

    # for n in range(1, 2):
    #     filename = f"data/altn_randomequal_updown_n{n}_L{length}.dat"
    #     upList, downList = altn_randomequal_updown(length)
    #     outputPattern(upList, downList, filename)


if __name__ == "__main__":
    main()
