#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"$SCRIPT_DIR/download_python_packages.sh"
"$SCRIPT_DIR/download_data.sh"
"$SCRIPT_DIR/merge_data.sh"
"$SCRIPT_DIR/download_r_packages.sh"
