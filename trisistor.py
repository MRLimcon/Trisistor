import trisistor as tr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
plt.style.use('ggplot')

n = 100
polos = tr.trisis.geracao_polos(n)
campo_base = tr.trisis.geracao_campo(polos,n)
elecbur = tr.trisis.geracao_bur_elec(n)

for i in range(200):
    if i == 50:
        t = 1
    else:
        t = 0
    #print(i,"teste")
    campo = tr.trisis.atualizacao_campo(elecbur,campo_base,n)
    #print(i)
    elecbur = tr.trisis.atualizacao_elecbur(elecbur,campo,t,n)

densidade = tr.trisis.c_densidade(elecbur,polos,n)
#print(densidade)

plt.quiver(campo_base[:,:,0],campo_base[:,:,1])
plt.show()
plt.imshow(densidade, interpolation='none')
plt.show()
#plt.imshow(polos, interpolation='none')
#plt.show()
plt.imshow(elecbur, interpolation='none')
plt.show()