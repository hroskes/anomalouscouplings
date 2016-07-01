from helperstuff.combinehelpers import getrate

if __name__ == "__main__":
  for i in range(5):
    for analysis in "fa3", "fa2", "fL1":
        print analysis
        for flavor in "2e2mu", "4e", "4mu":
            print flavor
            print
            #print getrates(flavor)
            print getrate(flavor, "ggZZ")
            print
        print
