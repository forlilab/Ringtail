These examples are shows using micromamba as the environment manager,
other environment managers may be used. 

# create and prepare environment
micromamba create -y -n ringtail_test_gui "python>=3.9,<3.12"
micromamba activate ringtail_test_gui

# ringtail dependencies
micromamba install -c conda-forge meeko rdkit scipy pandas matplotlib pyside6

# saltviewer dependencies
micromamba install -c conda-forge pyopengl pylinalg wgpu-py pygltflib pygfx line_profiler tomli
pip install prody==2.6.1

mkdir ringtail_test_gui
cd ringtail_test_gui

# download and install ringtail
git clone -b performance git@github.com:forlilab/PrivateRingtail.git
cd PrivateRingtail
pip install .

# download molviewer
git clone -b pyside git@github.com:forlilab/salticidae.git

# activate the environment after all installs
micromamba activate ringtail_test_gui

# open the ringtail gui
python GUI/ringtail_gui.py

