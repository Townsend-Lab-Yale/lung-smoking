import pandas as pd
import numpy as np
from count_combinations import compute_samples
from cancer_epistasis import numbers_positive_lambdas
from cancer_epistasis import estimate_lambdas
from locations import merged_maf_file_name

## We fix M=3 so
numbers_positive_lambdas = numbers_positive_lambdas[3]

db = pd.read_csv(merged_maf_file_name)
db = db[db['Variant_Classification'] != 'Silent']


genes = ['TP53', 'KRAS', 'LRP1B']
samples = compute_samples(db, mutations=genes)


bounds = 1
print(f"Bounds: {bounds}")
MAP = estimate_lambdas(samples, draws=1,
                       upper_bound_prior=bounds,
                       kwargs={'return_raw':True})
print(f"MAP: {MAP[0]['lambdas']}")
hit_bound = (MAP[0]['lambdas'] == bounds)

while np.sum(hit_bound) > 0:
    prop_bounds = np.array([1.1 if x else 1 for x in hit_bound])*bounds
    print(f"Proposed bounds: {prop_bounds}")
    prop_MAP = estimate_lambdas(samples, draws=1,
                                upper_bound_prior=prop_bounds,
                                kwargs={'return_raw':True,
                                        'start':MAP[0]})
    print(f"Proposed MAP: {prop_MAP[0]['lambdas']}")
    if prop_MAP[1].fun > MAP[1].fun:
        MAP = prop_MAP
        bounds = prop_bounds
        hit_bound = (MAP[0]['lambdas'] == bounds)
    else:
        print("Proposed bounds have lower likelihood")
