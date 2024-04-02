import numpy as np

def rep_seq(x, rep=10):
    L = len(x) * rep
    res = np.zeros(L, dtype=x.dtype)
    idx = np.arange(len(x)) * rep
    for i in np.arange(rep):
         res[idx + i] = x
    return res

a=np.linspace(0,9,10)
print(a)
new= rep_seq(a,5)
print("new",new)