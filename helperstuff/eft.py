sW = 0.23119
cW = (1 - sW**2) ** 0.5
Lambda1 = 10000
mZ = 91.1876

def ghw2(ghz1, ghz2, ghz4, ghz1prime2):
    return cW**2 * ghz2

def ghw4(ghz1, ghz2, ghz4, ghz1prime2):
    return cW**2 * ghz4

def ghw1prime2(ghz1, ghz2, ghz4, ghz1prime2):
    return Lambda1**2 / (cW**2 - sW**2) * (ghz1prime2/Lambda1**2 + 2 * sW**2 * (-ghz2 / mZ**2))

def ghzgs1prime2(ghz1, ghz2, ghz4, ghz1prime2):
    return Lambda1**2 / (cW**2 - sW**2) * 2 * sW * cW * (ghz1prime2/Lambda1**2 - ghz2 / mZ**2)
