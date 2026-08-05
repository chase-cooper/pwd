import numpy as np
import os,shutil,subprocess,sys
from consts import *

def newSpeciesDict():
    """
    Return a species dictionary with every species' mode set to 0 except hydrogen.
    """
    return {
        'h':    {'mode': 2,'abn': 0.},       # Only hydrogen is on by default
        'he':   {'mode': 0,'abn': 0.},
        'li':   {'mode': 0,'abn': 0.},
        'be':   {'mode': 0,'abn': 0.},
        'b':    {'mode': 0,'abn': 0.},
        'c':    {'mode': 0,'abn': 0.},
        'n':    {'mode': 0,'abn': 0.},
        'o':    {'mode': 0,'abn': 0.},
        'fl':   {'mode': 0,'abn': 0.},
        'ne':   {'mode': 0,'abn': 0.},
        'na':   {'mode': 0,'abn': 0.},
        'mg':   {'mode': 0,'abn': 0.},
        'al':   {'mode': 0,'abn': 0.},
        'si':   {'mode': 0,'abn': 0.},
        'p':    {'mode': 0,'abn': 0.},
        's':    {'mode': 0,'abn': 0.},
        'cl':   {'mode': 0,'abn': 0.},
        'ar':   {'mode': 0,'abn': 0.},
        'k':    {'mode': 0,'abn': 0.},
        'ca':   {'mode': 0,'abn': 0.},
        'sc':   {'mode': 0,'abn': 0.},
        'ti':   {'mode': 0,'abn': 0.},
        'v':    {'mode': 0,'abn': 0.},
        'cr':   {'mode': 0,'abn': 0.},
        'mn':   {'mode': 0,'abn': 0.},
        'fe':   {'mode': 0,'abn': 0.},
        'co':   {'mode': 0,'abn': 0.},
        'ni':   {'mode': 0,'abn': 0.}
    }

def writeUnit5(path:str,teff:float,log_g:float,species_dict:dict=newSpeciesDict(),is_lte:bool=True,is_ltgray:bool=True):
    """
    Helper function that writes the unit 5 files for TLUSTY model runs.

    Inputs:
    - path (str):           the location of the new unit 5 file
    - teff (float):         the temperature value to be written
    - log_g (float):        the value of log(g) to be written
    - species_dict (dict):  the dictionary containing all elements, their modes of treatment, and
                            their abundances.
    - is_lte (bool):        whether or not the model to be run is an LTE model
    - is_ltgray (bool):     whether or not the model to be run is an LTGray model. If false, a 
                            previous model output must be provided as fort.8.
    """

    # First block: temp, gravity, LTE/LTGRAY flags, flag file, frequency points
    l = 'T' if is_lte else 'F'
    g = 'T' if is_ltgray else 'F'

    f = open(path+'/fort.5','w+')
    f.write(str(int(teff))+'\t'+np.format_float_positional(log_g,precision=3)+'\t! T_eff, log(g)\n')
    f.write(f'{l}\t{g}       ! LTE, LTGRAY\n')
    f.write('\'cwd.flag\' ! Add\'l params\n')
    f.write('200        ! frequency points(?)\n')

    # Second block: species, modes of computation, and abundances
    num_exp = len(species_dict.keys())
    f.write(str(num_exp)+'\t\t! atoms\n')
    f.write('* mode abn modpf\n')
    for key in species_dict.keys():
        # if key=='h':continue
        f.write(f'  {species_dict[key]['mode']}   {species_dict[key]['abn']}   0   ! {key}\n')
    
    # Third block: Explicit ions
    f.write('*\n* explicit ions\n*\n')
    f.write('* iat iz nlevs ilast ilvlin nonstd typion file\n*\n')

    # Include files for relevant species. Do not set mode=2 to any species but these!!!
    if species_dict['h']['mode']==2:
        f.write("""  1   0   16   0   0     0   ' H 1' './data/h1s16.dat'
  1   1   1    1   0     0   ' H 2' ' '\n""")
    if species_dict['he']['mode']==2:
        f.write("""  2   0   24   0   0     0   'He 1' './data/he1.dat'
  2   1   20   0   0     0   'He 2' './data/he2.dat'
  2   2   1    1   0     0   'He 3' ' '\n""")
    if species_dict['c']['mode']==2:
        f.write("""  6   0   40   0   0     0   ' C 1' './data/c1.dat'
  6   1   22   0   0     0   ' C 2' './data/c2.dat'
  6   2   46   0   0     0   ' C 3' './data/c3_34+12lev.dat'
  6   3   1    1   0     0   ' C 4' ' '\n""")
#   6   3   25   0   0     0   ' C 4' './data/c4.dat'
#   6   4   1    1   0     0   ' C 5' ' '\n""")
    if species_dict['n']['mode']==2:
        f.write("""  7   0   34   0   0     0   ' N 1' './data/n1.dat'
  7   1   42   0   0     0   ' N 2' './data/n2_32+10lev.dat'
  7   2   1    1   0     0   ' N 3' ' '\n""")
#   7   2   32   0   0     0   ' N 3' './data/n3.dat'
#   7   3   48   0   0     0   ' N 4' './data/n4_34+14lev.dat'
#   7   4   16   0   0     0   ' N 5' './data/n5.dat'
#   7   5   1    1   0     0   ' N 6' ' '\n""")
    if species_dict['o']['mode']==2:
        f.write("""  8   0   33   0   0     0   ' O 1' './data/o1_23+10lev.dat'
  8   1   48   0   0     0   ' O 2' './data/o2_36+12lev.dat'
  8   2   41   0   0     0   ' O 3' './data/o3_28+13lev.dat'
  8   3   1    1   0     0   ' O 4' ' '\n""")
#   8   3   39   0   0     0   ' O 4' './data/o4.dat'
#   8   4   6    0   0     0   ' O 5' './data/o5.dat'
#   8   5   1    1   0     0   ' O 6' ' '\n""")
    if species_dict['na']['mode']==2:
        f.write("""  11  0   32   0   0     0   'Na 1' './data/na1.t'
  11  1   8    1   0     0   'Na 2' ' '\n""")
    if species_dict['mg']['mode']==2:
        f.write("""  12  1   25   0   0     0   'Mg 2' './data/mg2.dat'
  12  2   1    1   0     0   'Mg 3' ' '\n""")
#     if species_dict['al']['mode']==2:
#         f.write("""  13  1   29   0   0     0   'Al 2' './data/al2_20+9lev.dat'
#   13  2   23   0   0     0   'Al 3' './data/al3_19+4lev.dat'
#   13  3   1    1   0     0   'Al 4' ' '\n""")
    if species_dict['si']['mode']==2:
        f.write("""  14  1   40   0   0     0   'Si 2' './data/si2_36+4lev.dat'
  14  2   30   0   0     0   'Si 3' './data/si3.dat'
  14  3   23   0   0     0   'Si 4' './data/si4.dat'
  14  4   1    1   0     0   'Si 5' ' '\n""")
    if species_dict['s']['mode']==2:
        f.write("""  16  1   33   0   0     0   ' S 2' 'data/s2_23+10lev.dat'
  16  2   41   0   0     0   ' S 3' './data/s3_29+12lev.dat'
  16  3   38   0   0     0   ' S 4' './data/s4_33+5lev.dat'
  16  4   25   0   0     0   ' S 5' './data/s5_20+5lev.dat'
  16  5   1    1   0     0   ' S 6' ' '\n""")
#     if species_dict['ca']['mode']==2:
#         f.write("""  20  1   
# """)
    if species_dict['fe']['mode']==2:
        f.write("""  26  1   36   0   0     -1  'Fe 2' './data/fe2va.dat'
   0  0                             './data/gf2601.gam'
                                    './data/gf2601.lin'
                                    './data/fe2p_14+11lev.rap'
  26  2   50   0   0     -1  'Fe 3' './data/fe3va.dat'
   0  0                             './data/gf2602.gam'
                                    './data/gf2602.lin'
                                    './data/fe3p_22+7lev.rap'
  26  3   1    1   0     0   'Fe 4' ' '\n""")
#   26  3   43   0   0     -1  'Fe 4' './data/fe4va.dat'
#    0  0                             './data/gf2603.gam'
#                                     './data/gf2603.lin'
#                                     './data/fe4p_21+11lev.rap'
#   26  4   42   0   0     -1  'Fe 5' './data/fe5va.dat'
#    0  0                             './data/gf2604.gam'
#                                     './data/gf2604.lin'
#                                     './data/fe5p_19+11lev.rap'
#   26  5    1   1   0     0   'Fe 6' ' '\n""")

    f.write("  0   0   0   -1   0     0   '    ' ' '")
    f.close()

