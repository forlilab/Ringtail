# shell script to install ringtail gui with unnamed viewer
#!/bin/sh
# Define micromamba root (change this if your install location is different)
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

# create and prepare environment
micromamba create -y -n ringtail_test_gui "python>=3.9,<3.12"
micromamba activate ringtail_test_gui

# ringtail dependencies
micromamba install -c conda-forge meeko rdkit multiprocess scipy pandas matplotlib chemicalite pyside6

# saltviewer dependencies
micromamba install -c conda-forge prody pyopengl pylinalg wgpu-py pygltflib pygfx line_profiler

micromamba activate ringtail_test_gui
mkdir ringtail_test_gui
cd ringtail_test_gui

# download and install ringtail
git clone -b performance git@github.com:forlilab/PrivateRingtail.git
cd PrivateRingtail
pip install -e .
micromamba activate ringtail_test_gui

# download molviewer
git clone -b pyside git@github.com:forlilab/salticidae.git




