#!/bin/bash

# Required modules to compile on Harvard's Odyssey cluster
if [ $(hostname -f | grep -o "rc.fas.harvard.edu") ]; then
    module load fftw/3.3.4-fasrc09
    module load gcc/5.3.0-fasrc01
    module load openmpi/1.10.4-fasrc01
fi

if [ ! -d s2hat ]; then
    mkdir -p s2hat
fi

if [ ! -e s2hat/s2hat.h ]; then
    if [ ! -e s2hat_v2.55beta_30april2012.tar.gz ]; then
        wget http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat/software/s2hat_v2.55beta_30april2012.tar.gz
    fi
    tar -C s2hat -xvf s2hat_v2.55beta_30april2012.tar.gz
fi

if [ ! -e s2hat/s2hat_pure.h ]; then
    if [ ! -e pureS2HAT_v1.1beta_24feb2010.tar.gz ]; then
        wget http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/software/pures2hat/download/pureS2HAT_v1.1beta_24feb2010.tar.gz
    fi
    tar -C s2hat -xvf pureS2HAT_v1.1beta_24feb2010.tar.gz
fi

cd s2hat
if [ ! -e Makefile ]; then
    ln -s ../s2hat.makefile Makefile
fi
make

