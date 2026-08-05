import matplotlib.pyplot as plt
import numpy as np
import os

def grid(teff=range(10000,25000,5000),logg=range(750,900,50)):

    for t in teff:
        for g in logg:
            ts = str(t)
            gs = str(g)
            gv = float(g)/100.
            gw = str(gv)

            mo = 't' + ts + 'g' + gs
            print(mo)

            a = "sed 's/20000 8.0/" + ts + "  " + gw + "/' zl.5 >" + mo + "l.5"
            os.system(a)
            print(a)
            b = './Tl ' + mo + 'l'
            os.system(b)
            print(b)

            
