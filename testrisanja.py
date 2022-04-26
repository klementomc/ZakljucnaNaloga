import numpy as np
import roboticstoolbox as rp

panda = rp.models.DH.Panda()

q = [0,0,0,0,0,0,0]

panda.plot(q)