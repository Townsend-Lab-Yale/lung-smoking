import pandas as pd
from lift_coords import lift_over
## The package lift_coords requires to have the LiftOver command-line
## program. Obtain it from https://genome-store.ucsc.edu/

from locations import full_maf_file_names

destination_build = 'grch38'


builds = {'TCGA':'grch38',
          'FM-AD':'grch38',
          'OncoSG':'grch37',
          'Broad':'grch37',
          'MSK2015':'grch37',
          'TSP':'grch37',
          'MSK2017':'grch37',
          'MSK2018':'grch37',
          'TracerX':'grch37',
          'Genie':'grch37'
          }


print("Importing original data frames...")
original_dfs = {db:pd.read_csv(full_maf_file_names[db],
                               sep="\t",
                               comment="#")
                for db in builds.keys()}
print("...done.")
print("")


## Fix typo in TracerX
original_dfs['TracerX'].End_Position = original_dfs['TracerX'].End_Position.map(
    lambda v: 178035542 if v == 'DEC1178035542' else v)


lifted_dfs = {}
unlifted_dfs = {db:[] for db in builds.keys()}

for db, build in builds.items():
    if build != destination_build:
        print(f"Lifting coordinates of {db} from {build} to {destination_build}...")
        lifted_dfs[db], unlifted_dfs[db] = lift_over(original_dfs[db],
                                                     build,
                                                     destination_build)
        print("...done.")
        print("")
    else:
        lifted_dfs[db] = original_dfs[db]

    lifted_dfs[db].to_csv(full_maf_file_names[db].replace(".txt", "_lifted.txt"))
