from collections import OrderedDict
import json

#Configuration options
files = ["ggH0+_0.root", "ggH0+_1.root", "ggH0+_2.root", "ggH0+_3.root", "ggH0+_4.root", "ggH0-_0.root", "ggH0-_1.root", "ggH0-_2.root", "ggH0-_3.root", "ggHL1_0.root", "ggHL1_1.root", "ggHL1_2.root", "ggHL1_3.root", "ggHL1_4.root", "ggHa2_0.root", "ggHa2_1.root", "ggHa2_2.root", "ggHa2_3.root", "ggHa2_4.root", "ggHfL10.5_0.root", "ggHfL10.5_1.root", "ggHfL10.5_2.root", "ggHfL10.5_3.root", "ggHfL10.5_4.root", "ggHfa20.5_0.root", "ggHfa20.5_1.root", "ggHfa20.5_2.root", "ggHfa20.5_3.root", "ggHfa20.5_4.root", "ggHfa30.5_0.root", "ggHfa30.5_1.root", "ggHfa30.5_2.root", "ggHfa30.5_3.root"]

def addtemplates(jsondict, name, weight, flavorselection, domirror):
    postprocessing = [
      OrderedDict([("type", "smooth"), ("kernel", "adaptive"), ("entriesperbin", 50)]),
      OrderedDict([("type", "reweight"), ("axes",[0,1,2])])
    ]
    if not domirror:
      postprocessing.append(
        OrderedDict([("type", "floor")])
      )

    jsondict["templates"].append(OrderedDict([
      ("name", name),
      ("files", files),
      ("tree", "SelectedTree"),
      ("variables", ["D_0minus_decay","D_CP_decay","D_bkg_0plus"]),
      ("weight", weight),
      ("selection", "ZZMass>105 && ZZMass<140 && "+flavorselection),
      ("binning",OrderedDict([
        ("type", "fixed"),
        ("bins",[50,0.,1.,50,-0.5,0.5,50,0.,1.])
      ])),
      ("conserveSumOfWeights", True),
      ("postprocessing",postprocessing)
    ]))

    if domirror:
      jsondict["templates"].append(OrderedDict([
        ("name", name+"Mirror"),
        ("templatesum",[
          OrderedDict([("name", name),("factor", 1.)])
        ]),
        ("postprocessing",[
         OrderedDict([("type", "mirror"), ("axis", 1)]),
          OrderedDict([("type", "floor")])
        ])
      ]))


def maketemplates(flavor):

  jsondict = OrderedDict([
    ("inputDirectory", "/work-zfs/lhc/heshy/anomalouscouplings/Summer2015_VBF/maketemplates/step4_withdiscriminants/"),
    ("outputFile", "templates_{}.root".format(flavor)),
    #template definitions
    ("templates", [])
  ])

  if flavor == "2e2mu":
    flavorproduct = 13*13*11*11
  elif flavor == "4mu":
    flavorproduct = 13*13*13*13
  elif flavor == "4e":
    flavorproduct = 11*11*11*11
  else:
    raise ValueError("Bad flavor {}".format(flavor))
  flavorselection = "Z1Flav*Z2Flav == {}".format(flavorproduct)

  addtemplates(jsondict, "template0PlusAdapSmooth", "MC_weight_ggH_g1", flavorselection, domirror=True)
  addtemplates(jsondict, "template0MinusAdapSmooth", "MC_weight_ggH_g4", flavorselection, domirror=True)
  addtemplates(jsondict, "templateMixAdapSmooth", "MC_weight_ggH_g1g4", flavorselection, domirror=False)
  #template interference (use non-mirrored inputs), anti-mirror done afterwards
  jsondict["templates"].append(OrderedDict([
    ("name", "templateIntAdapSmoothMirror"),
    ("templatesum",[
      OrderedDict([("name", "templateMixAdapSmooth"),("factor", 1.)]),
      OrderedDict([("name", "template0PlusAdapSmooth"),("factor", -1.)]),
      OrderedDict([("name", "template0MinusAdapSmooth"),("factor", -1.)])
    ]),
    ("postprocessing",[
      OrderedDict([("type", "mirror"), ("antisymmetric", True), ("axis", 1)])
    ])
  ]))

  with open("templates_{}.json".format(flavor), "w") as f:
    f.write(json.dumps(jsondict, indent=4, separators=(',', ': ')))


if __name__ == "__main__":
  for flavor in "2e2mu", "4e", "4mu":
    maketemplates(flavor)
