#!/bin/bash -e

# get install location
if [ $# -eq 0 ]
	then
		echo 'No install location defined, using' $HOME'/.local/'
		prefix=$HOME/.local/
	else
		prefix=$1
		echo 'Will install in' $prefix
fi

# make a destination directory for runtime files
export TEMPO2=$prefix/share/tempo2
mkdir -p $TEMPO2

curl -L https://api.github.com/repos/mattpitkin/tempo2/tarball/master -o tempo2.tgz
tar zxvf tempo2.tgz

cd mattpitkin-tempo2-*

./bootstrap
./configure --prefix=$prefix
make && make install
cp -r T2runtime/* $TEMPO2
cd ..

rm -rf mattpitkin-tempo2-*
rm -rf tempo2.tgz
echo "Set TEMPO2 environment variable to ${TEMPO2} to make things run more smoothly."

