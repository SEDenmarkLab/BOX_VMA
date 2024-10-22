PART 1: In Silico Library Generation and Preprocessing for Conformer Generation

The main library contains everything except for the focused analogues.

1. Run python generate_library.py in cdxmlparser environment

Parses chemdraws, makes combinatorial library, and minimizes the library with obabel UFF force field. Output has CBr2 linker. Saved raw molli output in in_silico_library folder, and UFF minimized structures in in_silico_library_uff_min. There are 96160 in silico ligands.

2. Run python move_stereochem_checker.py in rdkit2023 environment

Checks files for the correct chirality and move good files to in_silico_library_uff_min_checked1 folder. Some files fail RDKit's CIPlabelling for runtime error. These structures are huge and RDKit struggles with CIP assignments for the adamantyl substituents. The failed files (47 files) were excluded, and are found in in_silico_library_uff_min_checked1 folder. Only 3 files actually had problematic stereochemistry, which were not moved from in_silico_library_uff_min folder.

The total size of the checked library is 96110.


DEPENDENCIES:

# packages in environment at /home/colen2/anaconda3/envs/cdxmlparser:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
asttokens                 2.2.1              pyhd8ed1ab_0    conda-forge
attrs                     22.1.0                   pypi_0    pypi
backcall                  0.2.0              pyh9f0ad1d_0    conda-forge
backports                 1.0                pyhd8ed1ab_3    conda-forge
backports.functools_lru_cache 1.6.4              pyhd8ed1ab_0    conda-forge
bidict                    0.21.4                   pypi_0    pypi
bzip2                     1.0.8                h7f98852_4    conda-forge
ca-certificates           2022.12.7            ha878542_0    conda-forge
certifi                   2022.12.7                pypi_0    pypi
charset-normalizer        3.1.0                    pypi_0    pypi
colorama                  0.4.6                    pypi_0    pypi
contourpy                 1.0.7                    pypi_0    pypi
cycler                    0.11.0                   pypi_0    pypi
debugpy                   1.6.7                    pypi_0    pypi
decorator                 5.1.1              pyhd8ed1ab_0    conda-forge
executing                 1.2.0              pyhd8ed1ab_0    conda-forge
fonttools                 4.39.3                   pypi_0    pypi
h5py                      3.7.0                    pypi_0    pypi
idna                      3.4                      pypi_0    pypi
imageio                   2.28.0                   pypi_0    pypi
importlib-metadata        6.6.0              pyha770c72_0    conda-forge
importlib_metadata        6.6.0                hd8ed1ab_0    conda-forge
ipykernel                 6.14.0                   pypi_0    pypi
ipython                   8.4.0                    pypi_0    pypi
jedi                      0.18.2             pyhd8ed1ab_0    conda-forge
jupyter-core              5.3.0                    pypi_0    pypi
jupyter_client            8.2.0              pyhd8ed1ab_0    conda-forge
jupyter_core              5.3.0           py310hff52083_0    conda-forge
kiwisolver                1.4.4                    pypi_0    pypi
ld_impl_linux-64          2.40                 h41732ed_0    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-ng                 12.2.0              h65d4601_19    conda-forge
libgomp                   12.2.0              h65d4601_19    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libsodium                 1.0.18               h36c2ea0_1    conda-forge
libsqlite                 3.40.0               h753d276_1    conda-forge
libstdcxx-ng              12.2.0              h46fd767_19    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libzlib                   1.2.13               h166bdaf_4    conda-forge
matplotlib                3.7.1                    pypi_0    pypi
matplotlib-inline         0.1.6              pyhd8ed1ab_0    conda-forge
molli                     1.0.0a7                  pypi_0    pypi
msgpack                   1.0.5                    pypi_0    pypi
ncurses                   6.3                  h27087fc_1    conda-forge
nest-asyncio              1.5.6              pyhd8ed1ab_0    conda-forge
networkx                  3.1                      pypi_0    pypi
numpy                     1.24.3                   pypi_0    pypi
openssl                   3.1.0                hd590300_2    conda-forge
packaging                 22.0                     pypi_0    pypi
parso                     0.8.3              pyhd8ed1ab_0    conda-forge
pexpect                   4.8.0              pyh1a96a4e_2    conda-forge
pickleshare               0.7.5                   py_1003    conda-forge
pillow                    9.5.0                    pypi_0    pypi
pip                       23.1.2             pyhd8ed1ab_0    conda-forge
platformdirs              3.5.0              pyhd8ed1ab_0    conda-forge
pooch                     1.7.0                    pypi_0    pypi
prompt-toolkit            3.0.38             pyha770c72_0    conda-forge
psutil                    5.9.5                    pypi_0    pypi
ptyprocess                0.7.0              pyhd3deb0d_0    conda-forge
pure_eval                 0.2.2              pyhd8ed1ab_0    conda-forge
py3dmol                   2.0.1.post1              pypi_0    pypi
pygments                  2.15.1             pyhd8ed1ab_0    conda-forge
pyparsing                 3.0.9                    pypi_0    pypi
python                    3.10.10         he550d4f_0_cpython    conda-forge
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python_abi                3.10                    3_cp310    conda-forge
pyvista                   0.38.5                   pypi_0    pypi
pyyaml                    5.4.1                    pypi_0    pypi
pyzmq                     25.0.2                   pypi_0    pypi
readline                  8.2                  h8228510_1    conda-forge
requests                  2.29.0                   pypi_0    pypi
scipy                     1.10.1                   pypi_0    pypi
scooby                    0.7.1                    pypi_0    pypi
setuptools                67.7.2             pyhd8ed1ab_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
stack_data                0.6.2              pyhd8ed1ab_0    conda-forge
tk                        8.6.12               h27826a3_0    conda-forge
tornado                   6.3                      pypi_0    pypi
tqdm                      4.64.1                   pypi_0    pypi
traitlets                 5.9.0              pyhd8ed1ab_0    conda-forge
typing-extensions         4.5.0                hd8ed1ab_0    conda-forge
typing_extensions         4.5.0              pyha770c72_0    conda-forge
tzdata                    2023c                h71feb2d_0    conda-forge
urllib3                   1.26.15                  pypi_0    pypi
vtk                       9.2.6                    pypi_0    pypi
wcwidth                   0.2.6              pyhd8ed1ab_0    conda-forge
wheel                     0.40.0             pyhd8ed1ab_0    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
zeromq                    4.3.4                h9c3ff4c_1    conda-forge
zipp                      3.15.0             pyhd8ed1ab_0    conda-forge


# packages in environment at /home/colen2/anaconda3/envs/rdkit2023:
#
# Name                    Version                   Build  Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
at-spi2-atk               2.38.0               h0630a04_3    conda-forge
at-spi2-core              2.40.3               h0630a04_0    conda-forge
atk-1.0                   2.38.0               hd4edc92_1    conda-forge
boost                     1.78.0          py311h59ea3da_4    conda-forge
boost-cpp                 1.78.0               h75c5d50_1    conda-forge
brotli                    1.0.9                h166bdaf_8    conda-forge
brotli-bin                1.0.9                h166bdaf_8    conda-forge
brotlipy                  0.7.0                    pypi_0    pypi
bzip2                     1.0.8                h7f98852_4    conda-forge
ca-certificates           2022.12.7            ha878542_0    conda-forge
cairo                     1.16.0            ha61ee94_1014    conda-forge
cairomm-1.0               1.14.4               h09cb3b9_0    conda-forge
certifi                   2022.12.7          pyhd8ed1ab_0    conda-forge
cffi                      1.15.1                   pypi_0    pypi
charset-normalizer        2.1.1              pyhd8ed1ab_0    conda-forge
colorama                  0.4.6              pyhd8ed1ab_0    conda-forge
contourpy                 1.0.7                    pypi_0    pypi
cryptography              39.0.2                   pypi_0    pypi
cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
dbus                      1.13.6               h5008d03_3    conda-forge
epoxy                     1.5.10               h166bdaf_1    conda-forge
expat                     2.5.0                h27087fc_0    conda-forge
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
fontconfig                2.14.2               h14ed4e7_0    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
fonttools                 4.38.0                   pypi_0    pypi
freetype                  2.12.1               hca18f0e_1    conda-forge
fribidi                   1.0.10               h36c2ea0_0    conda-forge
gdk-pixbuf                2.42.10              h05c8ddd_0    conda-forge
gettext                   0.21.1               h27087fc_0    conda-forge
glib-tools                2.74.1               h6239696_1    conda-forge
gmp                       6.2.1                h58526e2_0    conda-forge
graph-tool                2.46            py311h35b2f40_0    conda-forge
graph-tool-base           2.46            py311hab43870_0    conda-forge
graphite2                 1.3.13            h58526e2_1001    conda-forge
greenlet                  2.0.1                    pypi_0    pypi
gtk3                      3.24.36              h1e7d460_0    conda-forge
harfbuzz                  6.0.0                h8e241bc_0    conda-forge
hicolor-icon-theme        0.17                 ha770c72_2    conda-forge
icu                       70.1                 h27087fc_0    conda-forge
idna                      3.4                pyhd8ed1ab_0    conda-forge
jpeg                      9e                   h166bdaf_2    conda-forge
keyutils                  1.6.1                h166bdaf_0    conda-forge
kiwisolver                1.4.4                    pypi_0    pypi
krb5                      1.20.1               h81ceb04_0    conda-forge
lcms2                     2.14                 hfd0df8a_1    conda-forge
ld_impl_linux-64          2.39                 hcc3a1bd_1    conda-forge
lerc                      4.0.0                h27087fc_0    conda-forge
libblas                   3.9.0           16_linux64_openblas    conda-forge
libbrotlicommon           1.0.9                h166bdaf_8    conda-forge
libbrotlidec              1.0.9                h166bdaf_8    conda-forge
libbrotlienc              1.0.9                h166bdaf_8    conda-forge
libcblas                  3.9.0           16_linux64_openblas    conda-forge
libcups                   2.3.3                h36d4200_3    conda-forge
libdeflate                1.17                 h0b41bf4_0    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc-ng                 12.2.0              h65d4601_19    conda-forge
libgfortran-ng            12.2.0              h69a702a_19    conda-forge
libgfortran5              12.2.0              h337968e_19    conda-forge
libgirepository           1.74.0               ha2a38d2_1    conda-forge
libglib                   2.74.1               h606061b_1    conda-forge
libgomp                   12.2.0              h65d4601_19    conda-forge
libiconv                  1.17                 h166bdaf_0    conda-forge
libjpeg-turbo             2.1.4                h166bdaf_0    conda-forge
liblapack                 3.9.0           16_linux64_openblas    conda-forge
libnsl                    2.0.0                h7f98852_0    conda-forge
libopenblas               0.3.21          pthreads_h78a6416_3    conda-forge
libpng                    1.6.39               h753d276_0    conda-forge
librsvg                   2.54.4               h7abd40a_0    conda-forge
libsqlite                 3.40.0               h753d276_0    conda-forge
libstdcxx-ng              12.2.0              h46fd767_19    conda-forge
libtiff                   4.5.0                h6adf6a1_2    conda-forge
libuuid                   2.32.1            h7f98852_1000    conda-forge
libwebp-base              1.2.4                h166bdaf_0    conda-forge
libxcb                    1.13              h7f98852_1004    conda-forge
libxml2                   2.10.3               h7463322_0    conda-forge
libzlib                   1.2.13               h166bdaf_4    conda-forge
matplotlib                3.6.3                    pypi_0    pypi
matplotlib-base           3.6.3           py311h8597a09_0    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
ncurses                   6.3                  h27087fc_1    conda-forge
networkx                  3.0                pyhd8ed1ab_0    conda-forge
numpy                     1.24.1                   pypi_0    pypi
openjpeg                  2.5.0                hfec8fc6_2    conda-forge
openssl                   3.0.8                h0b41bf4_0    conda-forge
packaging                 23.0               pyhd8ed1ab_0    conda-forge
pandas                    1.5.3                    pypi_0    pypi
pango                     1.50.14              hd33c08f_0    conda-forge
pcre2                     10.40                hc3806b6_0    conda-forge
pillow                    9.4.0                    pypi_0    pypi
pip                       22.3.1             pyhd8ed1ab_0    conda-forge
pixman                    0.40.0               h36c2ea0_0    conda-forge
platformdirs              3.1.0              pyhd8ed1ab_0    conda-forge
pooch                     1.7.0              pyhd8ed1ab_0    conda-forge
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
pycairo                   1.23.0                   pypi_0    pypi
pycparser                 2.21               pyhd8ed1ab_0    conda-forge
pygobject                 3.42.2                   pypi_0    pypi
pyopenssl                 23.0.0             pyhd8ed1ab_0    conda-forge
pyparsing                 3.0.9              pyhd8ed1ab_0    conda-forge
pysocks                   1.7.1              pyha2e5f31_6    conda-forge
python                    3.11.0          he550d4f_1_cpython    conda-forge
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python_abi                3.11                    3_cp311    conda-forge
pytz                      2022.7.1           pyhd8ed1ab_0    conda-forge
rdkit                     2022.09.4       py311h277d903_1    conda-forge
readline                  8.1.2                h0f457ee_0    conda-forge
reportlab                 3.6.12                   pypi_0    pypi
requests                  2.28.2             pyhd8ed1ab_0    conda-forge
scipy                     1.10.1                   pypi_0    pypi
setuptools                66.1.1             pyhd8ed1ab_0    conda-forge
sigcpp-2.0                2.10.8               h27087fc_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
sparsehash                2.0.4                h9c3ff4c_0    conda-forge
spyrmsd                   0.5.2              pyhd8ed1ab_0    conda-forge
sqlalchemy                1.4.46                   pypi_0    pypi
tk                        8.6.12               h27826a3_0    conda-forge
tqdm                      4.64.1             pyhd8ed1ab_0    conda-forge
typing-extensions         4.4.0                hd8ed1ab_0    conda-forge
typing_extensions         4.4.0              pyha770c72_0    conda-forge
tzdata                    2022g                h191b570_0    conda-forge
urllib3                   1.26.14            pyhd8ed1ab_0    conda-forge
wheel                     0.38.4             pyhd8ed1ab_0    conda-forge
xorg-fixesproto           5.0               h7f98852_1002    conda-forge
xorg-inputproto           2.3.2             h7f98852_1002    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.0.10               h7f98852_0    conda-forge
xorg-libsm                1.2.3             hd9c2040_1000    conda-forge
xorg-libx11               1.7.2                h7f98852_0    conda-forge
xorg-libxau               1.0.9                h7f98852_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h7f98852_1    conda-forge
xorg-libxfixes            5.0.3             h7f98852_1004    conda-forge
xorg-libxi                1.7.10               h7f98852_0    conda-forge
xorg-libxrender           0.9.10            h7f98852_1003    conda-forge
xorg-libxtst              1.2.3             h7f98852_1002    conda-forge
xorg-recordproto          1.14.2            h7f98852_1002    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h7f98852_1002    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
zlib                      1.2.13               h166bdaf_4    conda-forge
zstandard                 0.19.0                   pypi_0    pypi
zstd                      1.5.2                h3eb15da_6    conda-forge