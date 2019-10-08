sW = 0.23119 ** 0.5
cW = (1 - sW**2) ** 0.5
Lambda1 = 10000
mZ = 91.1876
e = 0.313328534329

def ghw2(ghz1, ghz2, ghz4, ghz1prime2):
    return cW**2 * ghz2

def ghw4(ghz1, ghz2, ghz4, ghz1prime2):
    return cW**2 * ghz4

def ghw1prime2(ghz1, ghz2, ghz4, ghz1prime2):
    return Lambda1**2 / (cW**2 - sW**2) * (ghz1prime2/Lambda1**2 + 2 * sW**2 * (-ghz2 / mZ**2))

def ghzgs1prime2(ghz1, ghz2, ghz4, ghz1prime2):
    return Lambda1**2 / (cW**2 - sW**2) * 2 * sW * cW * (ghz1prime2/Lambda1**2 - ghz2 / mZ**2)

def deltacz(ghz1):
    return 0.5 * ghz1 - 1

def czz(ghz2):
    return -2 * sW**2 * cW**2 / e**2 * ghz2

def czbox(ghz1prime2):
    return mZ**2 * sW**2 / e**2 * ghz1prime2 / Lambda1**2

def czztilde(ghz4):
    return -2 * sW**2 * cW**2 / e**2 * ghz4
