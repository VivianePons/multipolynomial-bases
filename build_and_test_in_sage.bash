#!/bin/bash

set -e

cp -r /multipolynomial-basis .
cd multipolynomial-basis

echo
echo
echo "==================================================================="
echo "SageMath version"
sage --version


echo
echo
echo "==================================================================="
echo "Installing Multipolynomial basis"
make install

echo
echo
echo "==================================================================="
echo "Testing SnapPy"
make test

