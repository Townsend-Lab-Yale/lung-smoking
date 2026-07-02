#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

echo "Checking for python version..."
python_version=`python --version`
if [ "$python_version" != "Python 3.9.5" ]
then
    echo -e "Warning:\r\n We are working with python version 3.9.5, "`
    `"it has to be installed before running this setup.\r\n\r\n"`
    `"A good way to have multiple python versions on your system is "`
    `"by installing pyenv (https://github.com/pyenv/pyenv).\r\n\r\n"`
    `"If pyenv is installed, run:\r\n\r\n"`
    `"env PYTHON_CONFIGURE_OPTS=\"--enable-shared\" pyenv install 3.9.5\r\n"
    read -p "If you still want to continue with a different python version, "`
    `"press enter to continue or Cntrl-C to abort setup."
fi
echo $python_version
echo ""


echo "Creating a python virtual environment..."
python -m venv .venv
echo "done."
echo ""


echo "Installing all python packages..."
echo ""
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
echo "done installing python packages."
echo ""
echo ""


echo "Configuring the macOS PyTensor compiler..."
PYTENSOR_HOOK='# lung-smoking macOS PyTensor configuration'
if ! grep -Fq "$PYTENSOR_HOOK" .venv/bin/activate
then
    {
        echo ""
        echo "$PYTENSOR_HOOK"
        echo '_LUNG_SMOKING_ACTIVATION_HOOK=1'
        echo 'if ! . "$VIRTUAL_ENV/../configure_macos_pytensor.sh"; then'
        echo '    unset _LUNG_SMOKING_ACTIVATION_HOOK'
        echo '    deactivate'
        echo '    return 1 2>/dev/null || exit 1'
        echo 'fi'
        echo 'unset _LUNG_SMOKING_ACTIVATION_HOOK'
    } >> .venv/bin/activate
fi

if [ "$(uname -s)" = "Darwin" ]
then
    (
        . .venv/bin/activate
        lung_smoking_verify_macos_toolchain

        echo "Running a two-gene PyTensor compilation smoke test..."
        PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}" \
            "$VIRTUAL_ENV/bin/python" -c \
            'import numpy as np
from cancer_epistasis import estimate_lambdas

samples = np.array([100, 30, 20, 5], dtype=int)
result = estimate_lambdas(
    samples,
    upper_bound_prior=np.ones(4, dtype=float),
    draws=1,
    kwargs={"progressbar": False},
)
lambdas = np.asarray(result["lambdas"], dtype=float)
assert lambdas.shape == (4,), lambdas.shape
assert np.all(np.isfinite(lambdas)), lambdas
assert np.all(lambdas >= 0), lambdas
print("Verified two-gene flux estimation:", lambdas)'
    )
fi
echo "done configuring the PyTensor compiler."
echo ""
