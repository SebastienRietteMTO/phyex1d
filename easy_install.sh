#/usr/bin/bash

set -e

# Command line arguments
if (( $# == 1 )); then
  directory=$1
  mode='user'
elif (( $# == 2 )); then
  directory=$1
  mode=$2
  if [ "$mode" != 'user' -a "$mode" != 'dev' -a "$mode" != 'maxidev' ]; then
    echo "MODE must be user, dev or maxidev"
    exit 2
  fi
else
  echo "Usage $0 DIRECTORY <MODE>"
  echo "DIRECTORY must not exist and will be created, MODE can be 'user', 'dev' or 'maxidev'"
  exit 1
fi

# Directory
if [ -e "$directory" ]; then
  echo "DIRECTORY must not exist"
  exit 3
fi
mkdir "$directory"
cd "$directory"

# Python env
python3 -m venv phyex1d.env
. phyex1d.env/bin/activate

git clone https://github.com/GdR-DEPHY/DEPHY-SCM.git

git_pip() {
  git clone $1
  cd $(basename $1 .git)
  pip3 install -e .
  cd ..
}

if [ $mode == 'maxidev' ]; then
  git_pip https://github.com/SebastienRietteMTO/pppy.git
  git_pip https://github.com/UMR-CNRM/pyfortool.git
  git_pip https://github.com/SebastienRietteMTO/pyfxtran.git
  git_pip https://github.com/UMR-CNRM/ctypesForFortran.git
fi

if [ $mode == 'user' ]; then
  pip3 install phyex1d
else
  git_pip https://github.com/SebastienRietteMTO/phyex1d.git
fi

git clone https://github.com/UMR-CNRM/PHYEX.git
pip install -r PHYEX/requirements.txt
. PHYEX/tools/env.sh
cd PHYEX/build/with_ecbuild
./make_ecbuild.sh
cd ../../../

cat - <<EOF > execute.sh
cd $PWD
. PHYEX/tools/env.sh
. phyex1d.env/bin/activate
phyex1d DEPHY-SCM/ARMCU/REF/ARMCU_REF_SCM_driver.nc \\
        --exp dt=60 class=PhysicsAromeThetaR grid=L90mesonh \\
        --plot rv.png rv y_var=P \\
        --plot qc.png qc y_var=P vmin=0. vmax=5.E-5 \\
        --plot T.png T y_var=P \\
        --plot Theta.png Theta y_var=P \\
        --force
eog *.png
EOF
chmod +x execute.sh
./execute.sh
