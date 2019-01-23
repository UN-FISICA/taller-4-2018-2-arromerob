# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 18:29:52 2018

@author: FAMILIA
"""

class Derivada:
    def __init__(self, f, metodo ="adelante", dx= 0.001):
        self.f=f
        if metodo in {"adelante","central","extrapolada","segunda"}: 
            self.metodo=metodo
        else: raise TypeError(str(metodo)+" no es un nombre de método")
        self.dx=dx
    def calc(self,x):
        if self.metodo=="adelante": return (self.f(x+self.dx)-self.f(x))/self.dx
        elif self.metodo=="central": return (self.f(x+self.dx/2)-self.f(x-self.dx/2))/self.dx
        elif self.metodo=="extrapolada": 
            f4=2*(self.f(x+self.dx/4)-self.f(x-self.dx/4))/self.dx
            f2=(self.f(x+self.dx/2)-self.f(x-self.dx/2))/self.dx
            return (4*f4-f2)/3 
        elif self.metodo=="segunda": 
            f1=(self.f(x+self.dx)-self.f(x))/self.dx
            f2=(self.f(x-self.dx)-self.f(x))/self.dx
            return (f1+f2)/self.dx
class Zeros:
    def __init__(self, f, metodo, error=1e-4, max_iter=100):
        self.f=f
        if metodo in {"newton","bisectriz","interpolacion","newton-sp","fsolve-sp","brentq-sp"}: 
            self.metodo=metodo
        else: raise TypeError(str(metodo)+" no es un nombre de método")
        self.error=error
        self.max_iter=max_iter
    def zero(self,vi):
        if self.metodo=="newton": 
            i=0; ei=abs(self.f(vi)); x=vi
            while i<self.max_iter and ei>self.error:
                fprima=(self.f(x+self.error)-ei)/self.error
                x-= ei/fprima
                ei=abs(self.f(x))
                i+=1
            return x
        elif self.metodo=="bisectriz": 
            i=0; a,b=vi; fa=self.f(a); fb=self.f(b); 
            ei=min(abs(self.f((a+b)/2)),abs(fa),abs(fb)); 
            if (self.f(a)>0 and self.f(b)<0) or (self.f(b)>0 and self.f(a)<0) or ei<=self.error:
                while i<self.max_iter and ei>self.error:
                    i+=1; m=(a+b)/2; fm=self.f(m); print(ei)
                    if (fa>0 and fm<0) or (fm>0 and fa<0): 
                        b,fb=m,self.f(m);
                    else: 
                        a,fa=m,self.f(m); 
                    ei=min(abs(fa),abs(fb))
                ei=min(abs(fa),abs(fb),abs(self.f((a+b)/2))); 
                if abs(fa)==ei: return a
                else: 
                    if abs(fb)==ei: return b
                    else: return (a+b)/2
            else: raise TypeError("Con tus valores no se puede hacer bisectriz")
        elif self.metodo=="interpolacion": 
            i=0; a,b=vi; fa=self.f(a); fb=self.f(b); 
            ei=min(abs(fa),abs(fb))
            if (self.f(a)>0 and self.f(b)<0) or (self.f(b)>0 and self.f(a)<0) or ei<=self.error:
                while i<self.max_iter and ei>self.error:
                    n=(b-a)/(fb-fa); i+=1
                    m=a-fa*n; fm=self.f(m)
                    if (fa>0 and fm<0) or (fm>0 and fa<0): b,fb=m,self.f(m)
                    else: a,fa=m,self.f(m)
                    ei=abs(fm)
                if abs(fa)==ei:return a
                else: return b
            else: raise TypeError("Con tus valores no se puede hacer bisectriz")
        from scipy import optimize as sp
        if self.metodo=="newton-sp": 
            return sp.newton(self.f,vi)
        elif self.metodo=="fsolve-sp": 
            return sp.fsolve(self.f,vi)
        elif self.metodo=="brentq-sp":
            return sp.brentq(self.f,vi[0],vi[1])
def pol(x):
    return x*x-2
x=Zeros(pol,"newton-sp"); print(x.zero(2))