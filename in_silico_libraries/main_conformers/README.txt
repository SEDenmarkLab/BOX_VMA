PART 2: Conformer Generation

This makes conformers for all of the ligands.

1. Run generate_conformers.py in rdkit2023 environment (see main_library for list of dependencies).

takes awhile so tune to your resource

2. Run conformer_post_processing.py in rdkit2023 environment.

this makes a nice log output file and makes sure nothing went wrong.

3. Run mode_stereochem_checker.py in rdkit2023 environment.

checks stereochem and moves the good files

4. Run package_conformers_delete_linker.py in mollienv

this will delete the dibromomethylene linker from the conformers and package them as a molli mlib object for ASO calculation


DEPENDENCIES:

# packages in environment at /home/colen2/anaconda3/envs/mollienv:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
appdirs                   1.4.4                    pypi_0    pypi
attrs                     22.1.0                   pypi_0    pypi
bidict                    0.21.4                   pypi_0    pypi
bokeh                     2.4.3              pyhd8ed1ab_3    conda-forge
brotli                    1.0.9                h166bdaf_8    conda-forge
brotli-bin                1.0.9                h166bdaf_8    conda-forge
brotlipy                  0.7.0                    pypi_0    pypi
bzip2                     1.0.8                h7f98852_4    conda-forge
ca-certificates           2022.12.7            ha878542_0    conda-forge
cairo                     1.16.0            ha61ee94_1014    conda-forge
certifi                   2022.12.7                pypi_0    pypi
cffi                      1.15.1                   pypi_0    pypi
charset-normalizer        3.0.1                    pypi_0    pypi
click                     8.1.3           unix_pyhd8ed1ab_2    conda-forge
cloudpickle               2.2.1              pyhd8ed1ab_0    conda-forge
colorama                  0.4.6                    pypi_0    pypi
contourpy                 1.0.7                    pypi_0    pypi
cryptography              39.0.0                   pypi_0    pypi
cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
cytoolz                   0.12.0                   pypi_0    pypi
dask                      2023.1.0           pyhd8ed1ab_0    conda-forge
dask-core                 2023.1.0           pyhd8ed1ab_0    conda-forge
distributed               2023.1.0           pyhd8ed1ab_0    conda-forge
expat                     2.5.0                h27087fc_0    conda-forge
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
fontconfig                2.14.2               h14ed4e7_0    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
fonttools                 4.39.2                   pypi_0    pypi
freetype                  2.12.1               hca18f0e_1    conda-forge
fsspec                    2023.1.0           pyhd8ed1ab_0    conda-forge
gettext                   0.21.1               h27087fc_0    conda-forge
h5py                      3.7.0                    pypi_0    pypi
heapdict                  1.0.1                      py_0    conda-forge
icu                       70.1                 h27087fc_0    conda-forge
idna                      3.4                pyhd8ed1ab_0    conda-forge
imageio                   2.25.0                   pypi_0    pypi
jinja2                    3.1.2              pyhd8ed1ab_1    conda-forge
joblib                    1.2.0              pyhd8ed1ab_0    conda-forge
jpeg                      9e                   h166bdaf_2    conda-forge
kiwisolver                1.4.4                    pypi_0    pypi
lcms2                     2.14                 hfd0df8a_1    conda-forge
ld_impl_linux-64          2.39                 hcc3a1bd_1    conda-forge
lerc                      4.0.0                h27087fc_0    conda-forge
libblas                   3.9.0           16_linux64_openblas    conda-forge
libbrotlicommon           1.0.9                h166bdaf_8    conda-forge
libbrotlidec              1.0.9                h166bdaf_8    conda-forge
libbrotlienc              1.0.9                h166bdaf_8    conda-forge
libcblas                  3.9.0           16_linux64_openblas    conda-forge
libdeflate                1.17                 h0b41bf4_0    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-ng                 12.2.0              h65d4601_19    conda-forge
libgfortran-ng            12.2.0              h69a702a_19    conda-forge
libgfortran5              12.2.0              h337968e_19    conda-forge
libglib                   2.74.1               h606061b_1    conda-forge
libgomp                   12.2.0              h65d4601_19    conda-forge
libiconv                  1.17                 h166bdaf_0    conda-forge
libjpeg-turbo             2.1.4                h166bdaf_0    conda-forge
liblapack                 3.9.0           16_linux64_openblas    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libopenblas               0.3.21          pthreads_h78a6416_3    conda-forge
libpng                    1.6.39               h753d276_0    conda-forge
libsqlite                 3.40.0               h753d276_0    conda-forge
libstdcxx-ng              12.2.0              h46fd767_19    conda-forge
libtiff                   4.5.0                h6adf6a1_2    conda-forge
libuuid                   2.32.1            h7f98852_1000    conda-forge
libwebp-base              1.2.4                h166bdaf_0    conda-forge
libxcb                    1.13              h7f98852_1004    conda-forge
libxml2                   2.10.3               hca2bb57_3    conda-forge
libzlib                   1.2.13               h166bdaf_4    conda-forge
locket                    1.0.0              pyhd8ed1ab_0    conda-forge
lz4                       4.2.0                    pypi_0    pypi
lz4-c                     1.9.3                h9c3ff4c_1    conda-forge
markupsafe                2.1.2                    pypi_0    pypi
matplotlib                3.6.3                    pypi_0    pypi
matplotlib-base           3.7.1           py310he60537e_0    conda-forge
molli                     1.0.0a6                  pypi_0    pypi
msgpack                   1.0.4                    pypi_0    pypi
msgpack-python            1.0.4           py310hbf28c38_1    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
ncurses                   6.3                  h27087fc_1    conda-forge
networkx                  3.0                      pypi_0    pypi
numpy                     1.24.1                   pypi_0    pypi
openbabel                 3.1.1           py310heaf86c6_5    conda-forge
openjpeg                  2.5.0                hfec8fc6_2    conda-forge
openssl                   3.1.0                h0b41bf4_0    conda-forge
packaging                 22.0                     pypi_0    pypi
pandas                    1.5.3                    pypi_0    pypi
partd                     1.3.0              pyhd8ed1ab_0    conda-forge
patsy                     0.5.3              pyhd8ed1ab_0    conda-forge
pcre2                     10.40                hc3806b6_0    conda-forge
pillow                    9.4.0                    pypi_0    pypi
pip                       22.3.1             pyhd8ed1ab_0    conda-forge
pixman                    0.40.0               h36c2ea0_0    conda-forge
pooch                     1.6.0                    pypi_0    pypi
psutil                    5.9.4                    pypi_0    pypi
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
pycparser                 2.21               pyhd8ed1ab_0    conda-forge
pyopenssl                 23.0.0             pyhd8ed1ab_0    conda-forge
pyparsing                 3.0.9              pyhd8ed1ab_0    conda-forge
pysocks                   1.7.1              pyha2e5f31_6    conda-forge
python                    3.10.8          h4a9ceb5_0_cpython    conda-forge
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python_abi                3.10                    3_cp310    conda-forge
pytz                      2022.7.1           pyhd8ed1ab_0    conda-forge
pyvista                   0.37.0                   pypi_0    pypi
pyyaml                    5.4.1                    pypi_0    pypi
readline                  8.1.2                h0f457ee_0    conda-forge
requests                  2.28.2                   pypi_0    pypi
scikit-learn              1.2.2                    pypi_0    pypi
scipy                     1.10.0                   pypi_0    pypi
scooby                    0.7.0                    pypi_0    pypi
seaborn                   0.12.2               hd8ed1ab_0    conda-forge
seaborn-base              0.12.2             pyhd8ed1ab_0    conda-forge
setuptools                66.1.1             pyhd8ed1ab_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
sortedcontainers          2.4.0              pyhd8ed1ab_0    conda-forge
statsmodels               0.13.5                   pypi_0    pypi
tblib                     1.7.0              pyhd8ed1ab_0    conda-forge
threadpoolctl             3.1.0              pyh8a188c0_0    conda-forge
tk                        8.6.12               h27826a3_0    conda-forge
toolz                     0.12.0             pyhd8ed1ab_0    conda-forge
tornado                   6.2                      pypi_0    pypi
tqdm                      4.64.1                   pypi_0    pypi
typing_extensions         4.4.0              pyha770c72_0    conda-forge
tzdata                    2022g                h191b570_0    conda-forge
unicodedata2              15.0.0                   pypi_0    pypi
urllib3                   1.26.14            pyhd8ed1ab_0    conda-forge
vtk                       9.2.5                    pypi_0    pypi
wheel                     0.38.4             pyhd8ed1ab_0    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.0.10               h7f98852_0    conda-forge
xorg-libsm                1.2.3             hd9c2040_1000    conda-forge
xorg-libx11               1.8.4                h0b41bf4_0    conda-forge
xorg-libxau               1.0.9                h7f98852_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h0b41bf4_2    conda-forge
xorg-libxrender           0.9.10            h7f98852_1003    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h0b41bf4_1003    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
yaml                      0.2.5                h7f98852_2    conda-forge
zict                      2.2.0              pyhd8ed1ab_0    conda-forge
zlib                      1.2.13               h166bdaf_4    conda-forge
zstd                      1.5.2                h3eb15da_6    conda-forge