#!/bin/bash

echo "Merging MAF and clinical files..."
echo ""
source .venv/bin/activate
python prepare_maf_data.py
python merge_MAF_clinical.py
python store_panel_info.py
deactivate
echo "...done."
echo ""
echo ""
