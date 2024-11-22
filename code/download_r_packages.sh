#!/bin/bash

echo "Installing R packages..."
cd variants
mkdir .Rlibs
./ces_installation.R
echo ""
echo "...done installing R packages."
echo ""
