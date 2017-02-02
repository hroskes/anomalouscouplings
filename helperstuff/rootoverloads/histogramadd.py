import ROOT

def __add__(self, other):
    result = self.Clone(self.GetName()+"+"+other.GetName())
    result.Add(other)
    return result

def __iadd__(self, other):
    self.Add(other)
    return self

def __sub__(self, other):
    result = self.Clone(self.GetName()+"+"+other.GetName())
    result.Add(other, -1)
    return result

def __isub__(self, other):
    self.Add(other, -1)
    return self

def __mul__(self, other):
    result = self.Clone(self.GetName()+"*"+str(other))
    result.Scale(other)
    return result

__rmul__ = __mul__

def __imul__(self, other):
    self.Scale(other)
    return self

def __div__(self, other):
    result = self.Clone(self.GetName()+"/"+str(other))
    result.Scale(1.0/other)
    return result

def __idiv__(self, other):
    self.Scale(1.0/other)
    return self


ROOT.TH1.__add__ = __add__
ROOT.TH1.__iadd__ = __iadd__
ROOT.TH1.__sub__ = __sub__
ROOT.TH1.__isub__ = __isub__
ROOT.TH1.__mul__ = __mul__
ROOT.TH1.__rmul__ = __rmul__
ROOT.TH1.__imul__ = __imul__
ROOT.TH1.__div__ = __div__
ROOT.TH1.__idiv__ = __idiv__
