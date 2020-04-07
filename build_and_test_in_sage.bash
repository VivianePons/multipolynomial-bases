#!/bin/bash

set -e

cp -r /multipolynomial-bases .
cd multipolynomial-bases

echo
echo
echo "==================================================================="
echo "SageMath version"
sage --version


echo
echo
echo "==================================================================="
echo "Installing Multipolynomial bases"
sage -pip install --upgrade --no-index -v .

echo
echo
echo "==================================================================="
echo "Testing Myltupolynomial bases"
sage setup.py test

