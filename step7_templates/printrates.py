from helperstuff.combinehelpers import getrates

if __name__ == "__main__":
    for analysis in "fa3", "fa2", "fL1":
        print analysis
        for flavor in "2e2mu", "4e", "4mu":
            print flavor
            print
            print getrates(flavor, analysis)
            print
        print
