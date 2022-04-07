import ternary_modules.ternary_calculations as calculations
import ternary_modules.utils as utils
import sys
import matplotlib.pyplot as plt

plt.style.use('ggplot')

if len(sys.argv) <= 2:
    print("Coloque os argumentos para a ativação")
    print("1 para ativado ou 0 para desativado e depois a quantidade de iterações")
    print("Exemplo: python ./trisistor.py 1 500")
    exit()
is_activated = int(sys.argv[1])
time = int(sys.argv[2])

n = 100
polos = utils.utils.generate_doping(n, is_activated)

campo_base = utils.utils.generate_field(polos, n)

# plt.quiver(campo_base[:, :, 0], campo_base[:, :, 1])
# plt.title("Campo elétrico inicial")
# plt.show()
elecbur = utils.utils.generate_electrons(n)

for i in range(time):
    if i % 20 == 0:
        t = 1
    else:
        t = 0

    campo = calculations.trisistor.update_field(elecbur, campo_base, n)
    elecbur = calculations.trisistor.update_positions(elecbur, campo, t, n)

densidade = utils.utils.calculate_density(elecbur, polos, n)

# Plots
plt.quiver(campo[:, :, 0], campo[:, :, 1])
plt.title("Campo elétrico final")
plt.show()

plt.imshow(densidade, interpolation='none', vmin=-10, vmax=10)
plt.title("Densidade eletrônica em cada ponto")
plt.show()

plt.imshow(elecbur, interpolation='none')
plt.title("Localização dos elétrons e dos buracos")
plt.show()
