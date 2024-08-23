
import hippo
import os

import logging
logger = logging.getLogger('HIPPO')

os.system('cp -v aug15_rgen.sqlite aug15_batch.sqlite')

animal = hippo.HIPPO('aug15_batch', 'aug15_batch.sqlite')

recipe = hippo.Recipe.from_json(animal.db, 'recipe_jun3_80bases_from35instock.json', allow_db_mismatch=True)

gen = hippo.RandomRecipeGenerator(db=animal.db, start_with=recipe, suppliers=['Stock', 'Enamine'])

for i in range(2000):
    logger.header(f'{i=}')
    r = gen.generate(shuffle=True, currency='USD', debug=False)

animal.db.close()
