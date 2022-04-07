import os
import pip

pip.main(["install", "numpy"])
pip.main(["install", "matplotlib"])

os.system("cd ternary_modules && f2py -c \
    --opt='-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -fopenmp' \
    -lgomp \
    ternary_calculations.f95 -m ternary_calculations")

os.system("cd ternary_modules && f2py -c \
    --opt='-O3 -ftree-vectorize -march=native -fno-range-check \
        -floop-nest-optimize -fPIC -fopenmp' \
    -lgomp \
    utils.f95 -m utils")
