
import hippo
import os
from tqdm import tqdm

os.system('cp -v aug7_rgen.sqlite aug14_batch_fp.sqlite')

animal = hippo.HIPPO('aug14_batch_fp', 'aug14_batch_fp.sqlite')

recipe = hippo.Recipe.from_json(animal.db, 'recipe_jun3_80bases_from35instock.json', allow_db_mismatch=True)

gen = hippo.RandomRecipeGenerator(db=animal.db, start_with=recipe, suppliers=['Stock', 'Enamine'])

comps = animal.compounds[gen.route_pool.product_ids]
poses = comps.poses

ref = animal.poses['Ax0310a']
for pose in tqdm(poses):
    pose.reference = ref
    pose.calculate_fingerprint()

animal.db.close()
