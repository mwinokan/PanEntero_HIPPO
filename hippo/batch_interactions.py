
import hippo
import os
from tqdm import tqdm

os.system('cp -v aug15_rgen.sqlite aug23_interaction_bench.sqlite')

animal = hippo.HIPPO('aug23_interaction_bench', 'aug23_interaction_bench.sqlite')

### Fragments

for i,pose in tqdm(enumerate(animal.poses(tag='hits'))):
    pose.calculate_interactions()

recipe = hippo.Recipe.from_json(animal.db, 'recipe_jun3_80bases_from35instock.json', allow_db_mismatch=True)

# ### Recipe products

# gen = hippo.RandomRecipeGenerator(db=animal.db, start_with=recipe, suppliers=['Stock', 'Enamine'])

# comps = animal.compounds[gen.route_pool.product_ids]
# poses = comps.poses

# ref = animal.poses['Ax0310a']
# for i,pose in tqdm(enumerate(poses)):
#     pose.reference = ref
#     pose.calculate_interactions()

animal.db.close()
