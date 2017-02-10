import json
from numbers import Number
import os

import config
from enums import Analysis, Category, Channel, EnumItem, MultiEnum, MyEnum, ProductionMode
from utilities import getnesteddictvalue, setnesteddictvalue

class YieldSystematic(MyEnum):
    enumname = "yieldsystematic"
    enumitems = (
                 EnumItem("BTag"),
                 EnumItem("JEC"),
                )

class YieldValue(MultiEnum):
    enumname = "yieldvalue"
    enums = (Analysis, Category, Channel, ProductionMode)

    yieldsfile = os.path.join(config.repositorydir, "data", "yields.json")

    @classmethod
    def getyieldsdict(cls, trycache=True):
      import globals
      if globals.yieldsdict_cache is None or not trycache:
        try:
          with open(cls.yieldsfile) as f:
            jsonstring = f.read()
        except IOError:
          try:
            os.makedirs(os.path.dirname(cls.yieldsfile))
          except OSError:
            pass
          with open(cls.yieldsfile, "w") as f:
            f.write("{}\n")
            jsonstring = "{}"
        globals.yieldsdict_cache = json.loads(jsonstring)
      return globals.yieldsdict_cache

    @classmethod
    def writeyieldsdict(cls):
      dct = cls.getyieldsdict()
      jsonstring = json.dumps(dct, sort_keys=True, indent=4, separators=(',', ': '))
      with open(cls.yieldsfile, "w") as f:
        f.write(jsonstring)

    @property
    def keys(self):
        return (
                str(self.analysis),
                str(self.category),
                str(self.channel),
                str(self.productionmode),
               )

    @property
    def value(self):
        return getnesteddictvalue(self.getyieldsdict(), *self.keys, default=None)

    @value.setter
    def value(self, value):
        setnesteddictvalue(self.getyieldsdict(), *self.keys, value=value)
        assert self.value == value

    def __float__(self):
        return self.__value

class YieldSystematicValue(MultiEnum):
    enumname = "yieldsystematicvalue"
    enums = (YieldSystematic, Analysis, Category, Channel, ProductionMode)

    yieldsystematicsfile = os.path.join(config.repositorydir, "data", "categorysystematics.json")

    @classmethod
    def getyieldsystematicsdict(cls, trycache=True):
      import globals
      if globals.yieldsystematicsdict_cache is None or not trycache:
        try:
          with open(cls.yieldsystematicsfile) as f:
            jsonstring = f.read()
        except IOError:
          try:
            os.makedirs(os.path.dirname(cls.yieldsystematicsfile))
          except OSError:
            pass
          with open(cls.yieldsystematicsfile, "w") as f:
            f.write("{}\n")
            jsonstring = "{}"
        globals.yieldsystematicsdict_cache = json.loads(jsonstring)
      return globals.yieldsystematicsdict_cache

    @classmethod
    def writeyieldsystematicsdict(cls):
      dct = cls.getyieldsystematicsdict()
      jsonstring = json.dumps(dct, sort_keys=True, indent=4, separators=(',', ': '))
      with open(cls.yieldsystematicsfile, "w") as f:
        f.write(jsonstring)

    @property
    def keys(self):
        return (
                str(self.yieldsystematic),
                str(self.analysis),
                str(self.category),
                str(self.channel),
                str(self.productionmode),
               )

    @property
    def value(self):
        return getnesteddictvalue(self.getyieldsystematicsdict(), *self.keys, default=None)

    @value.setter
    def value(self, value):
        setnesteddictvalue(self.getyieldsystematicsdict(), *self.keys, value=value)
        assert self.value == value

    def __str__(self):
        if self.value is None or self.value == 0:
            return "-"
        if isinstance(self.value, Number):
            if not -1 < self.value < 1: raise ValueError("{!r} value {} should be between 0 and 1!".format(self, self.value))
            return str(self.value)
        if not (hasattr(self.value, __len__) and len(self.value) == 2):
            raise ValueError("{!r} value '{!r}' should be None, a number, or a list (tuple, etc.) of length 2".format(self, self.value))
        if not all(-1 < _ < 1 for _ in self.value):
            raise ValueError("Both elements of {!r} value {} should be between -1 and 1!".format(self, self.value))
        return "{}/{}".format(1+self.value[0], 1+self.value[1])
