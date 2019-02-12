import collections
import array
import ROOT

class ExtendedCounter(collections.Counter):
    """For one reason or another I can't get Counter.__add__ and __sub__
       to work here"""

    def __add__(self, other):
        result = ExtendedCounter(self)
        for item in other:
            if item not in result:
                result[item] = 0
            result[item] += other[item]
        return result

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + -other

    def __mul__(self, other):
        """multiply by a scalar"""
        result = ExtendedCounter(self)
        for item in result:
            result[item] *= other
        return result

    def _rmul__(self, other):
        """multiply by a scalar"""
        return self * other

    def __div__(self, other):
        """divide by a scalar"""
        result = ExtendedCounter(self)
        for item in result:
            result[item] /= other
        return result

    def zero(self):
        minvalue = min(self.values())
        for key in self:
            self[key] -= minvalue

    def TGraph(self):
        """make a TGraph with the data in self.  Keys and values have to be numbers."""
        items = self.items()
        items.sort(key = lambda x: x[0])
        keysvalues = zip(*items)
        keys = keysvalues[0]
        values = keysvalues[1]
        x = array.array("d", keys)
        y = array.array("d", values)
        g = ROOT.TGraph(len(self), x, y)
        g.SetTitle("")
        return g

    def TGraph2D(self):
        items = self.items()
        items.sort(key = lambda x: x[0])
        xx, yy, zz = zip(*((_[0][0], _[0][1], _[1]) for _ in items))
        x = array.array("d", xx)
        y = array.array("d", yy)
        z = array.array("d", zz)
        g = ROOT.TGraph2D(len(self), x, y, z)
        g.SetTitle("")

        return g
