import os
import requests, sys
import pandas as pd

from locations import gene_list_file
from locations import gene_coordinates_file

coordinates = pd.read_csv(gene_list_file,
                          names=["gene", "chromosome", "start", "end"])

coordinates = coordinates.set_index("gene")

for name in coordinates.index:
    server = "http://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens/" + name + "?content-type=application/json"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    coordinates.loc[name] = [decoded.get(key) for key in ('seq_region_name','start','end')]


no_coordinate = list(coordinates[coordinates['chromosome'].isna()].index)
no_coordinate += list(coordinates[coordinates['start'].isna()].index)
no_coordinate += list(coordinates[coordinates['end'].isna()].index)
no_coordinate = set(no_coordinate)
if no_coordinate != set([]):
    print("No coordinate found for the following genes:")
    for gene in no_coordinate:
        print(gene)

coordinates = coordinates.dropna()
coordinates = coordinates.astype({"chromosome":str, "start":int, "end":int})

coordinates.to_csv(gene_coordinates_file)
