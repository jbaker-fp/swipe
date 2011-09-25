# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.4
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.



from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_swipe', [dirname(__file__)])
        except ImportError:
            import _swipe
            return _swipe
        if fp is not None:
            try:
                _mod = imp.load_module('_swipe', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _swipe = swig_import_helper()
    del swig_import_helper
else:
    import _swipe
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def new_doublea(*args):
  return _swipe.new_doublea(*args)
new_doublea = _swipe.new_doublea

def delete_doublea(*args):
  return _swipe.delete_doublea(*args)
delete_doublea = _swipe.delete_doublea

def doublea_getitem(*args):
  return _swipe.doublea_getitem(*args)
doublea_getitem = _swipe.doublea_getitem

def doublea_setitem(*args):
  return _swipe.doublea_setitem(*args)
doublea_setitem = _swipe.doublea_setitem
import numpy as NP
from math import log, fsum, isnan
from bisect import bisect_left

## helper functions

def _mean(x):
    """ 
    Compute the mean of x
    """
    return fsum(x) / len(x)

def _var(x):
    """
    Compute the variance of x 
    """
    my_mean = mean(x)
    s = 0.
    for i in x:
        s += (i - my_mean) ** 2
    return s / len(x) - 1

def _regress(x, y):
    """
    Compute the intercept and slope for y ~ x
    """
    solution = NP.linalg.lstsq(NP.vstack((NP.ones(len(x)), x)).T, y)
    return solution[0]

## the class itself

class Swipe(object):
    """
    Wrapper class representing a SWIPE' p extraction
    """

    def __init__(self, path, pmin=100., pmax=600., st=.3, dt=0.001, mel=False):
        """
        Class constructor:

        path = either a file object pointing to a wav file, or a string path
        pmin = minimum frequency in Hz
        pmax = maximum frequency in Hz
        st = frequency cutoff (must be between [0.0, 1.0]
        dt = samplerate in seconds
        show_nan = if True, voiceless intervals are returned, marked as nan.
        """
        # Get Python path, just in case someone passed a file object
        f = path if isinstance(path, str) else path.name
        # Obtain the vector itself
        P = pyswipe(f, pmin, pmax, st, dt)
        # get function
        conv = None
        if mel: conv = lambda hz: 1127.01048 * log(1. + hz / 700.)
        else: conv = lambda hz: hz
        # generate
        tt = 0.
        self.t = []
        self.p = []
        if P.x < 1: 
            raise ValueError('Failed to read audio')
        for i in range(P.x):
            val = doublea_getitem(P.v, i)
            if not isnan(val):
                self.t.append(tt)
                self.p.append(conv(doublea_getitem(P.v, i)))
            tt += dt

    def __str__(self):
        return '<Swipe pitch track with %d points>' % len(self.t)

    def __len__(self):
        return len(self.t)

    def __iter__(self):
        return iter(zip(self.t, self.p))

    def __getitem__(self, t):
        """ 
        Takes a  argument and gives the nearest sample 
        """
        if self.t[0] <= 0.:
            raise ValueError, 'Time less than 0'
        i = bisect(self.t, t)
        if self.t[i] - t > t - self.t[i - 1]:
            return self.p[i - 1]
        else:
            return self.p[i]

    def _bisect(self, tmin=None, tmax=None):
        """ 
        Helper for bisection
        """
        if not tmin:
            if not tmax:
                raise ValueError, 'At least one of tmin, tmax must be defined'
            else:
                return (0, bisect(self.t, tmax))
        elif not tmax:
            return (bisect(self.t, tmin), len(self.t))
        else:
            return (bisect(self.t, tmin), bisect(self.t, tmax))

    def slice(self, tmin=None, tmax=None):
        """ 
        Slice out samples outside of s [tmin, tmax] inline 
        """
        if tmin or tmax:
            (i, j) = self._bisect(tmin, tmax)
            self.t = self.t[i:j]
            self.p = self.p[i:j]
        else:
            raise ValueError, 'At least one of tmin, tmax must be defined'

    def mean(self, tmin=None, tmax=None):
        """ 
        Return pitch mean 
        """
        if tmin or tmax:
            (i, j) = self._bisect(tmin, tmax)
            return mean(self.p[i:j])
        else:
            return mean(self.p)

    def var(self, tmin=None, tmax=None):
        """ 
        Return pitch variance 
        """
        if tmin or tmax:
            (i, j) = self._bisect(tmin, tmax)
            return var(self.p[i:j])
        else:
            return var(self.p)

    def sd(self, tmin=None, tmax=None):
        """ 
        Return pitch standard deviation 
        """
        return sqrt(self.var(tmin, tmax))

    def regress(self, tmin=None, tmax=None):
        """ 
        Return the linear regression intercept and slope for pitch ~ time. I
        wouldn't advise using this on raw p, but rather the Mel frequency 
        option: e.g., call Swipe(yourfilename, min, max, mel=True). The reason 
        for this is that Mel frequency is log-proportional to p in Hertz, 
        and I find log pitch is much closer to satisfying the normality 
        assumption.
        """
        if tmin or tmax:
            (i, j) = self._bisect(tmin, tmax)
            return _regress(self.t[i:j], self.p[i:j])
        else:
            return _regress(self.t, self.p)

class vector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, vector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, vector, name)
    __repr__ = _swig_repr
    __swig_setmethods__["x"] = _swipe.vector_x_set
    __swig_getmethods__["x"] = _swipe.vector_x_get
    if _newclass:x = _swig_property(_swipe.vector_x_get, _swipe.vector_x_set)
    __swig_setmethods__["v"] = _swipe.vector_v_set
    __swig_getmethods__["v"] = _swipe.vector_v_get
    if _newclass:v = _swig_property(_swipe.vector_v_get, _swipe.vector_v_set)
    def __init__(self): 
        this = _swipe.new_vector()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _swipe.delete_vector
    __del__ = lambda self : None;
vector_swigregister = _swipe.vector_swigregister
vector_swigregister(vector)


def pyswipe(*args):
  return _swipe.pyswipe(*args)
pyswipe = _swipe.pyswipe
# This file is compatible with both classic and new-style classes.


