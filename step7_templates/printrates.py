from helperstuff.combinehelpers import getrates

if __name__ == "__main__":
    for flavor in "2e2mu", "4e", "4mu":
        print flavor
        print
        print getrates(flavor, "fa3")
        print
