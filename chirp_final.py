import numpy as np

def repet_val(array,rep):
    res = []
    for valor in array:
        res.extend([valor] * int(rep))
    return res

def valor_n(arreglo, cada_n):
    nuevo_arreglo = []
    contador = 0
    for valor in arreglo:
        if contador % cada_n == 0:
            nuevo_arreglo.append(valor)
        contador += 1
    return nuevo_arreglo



def rep_seq(x, rep=10):
     L = len(x) * rep
     res = np.zeros(L, dtype=x.dtype)
     idx = np.arange(len(x)) * rep
     for i in np.arange(rep):
         res[idx + i] = x
     return(res)
