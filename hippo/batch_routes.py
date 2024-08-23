
import os
os.system('cp -v june7_knitwork.sqlite batch_routes_jun13.sqlite ')

import hippo
from sqlite3 import DatabaseError

from mlog import setup_logger
logger = setup_logger('batch_routes_jun13')

animal = hippo.HIPPO('batch_routes_jun13', 'batch_routes_jun13.sqlite')

animal.db.create_table_route()
animal.db.create_table_component()

comps = animal.add_enamine_quote(path='../enamine_quotes/Q1870545_USD_20mg.xlsx', orig_name_col='Customer Code')
comps = animal.add_enamine_quote(path='../enamine_quotes/Q1870545_USD_50mg.xlsx', orig_name_col='Customer Code')
comps = animal.add_enamine_quote(path='../enamine_quotes/Q1870545_USD_100mg.xlsx', orig_name_col='Customer Code')
comps = animal.add_enamine_quote(path='../enamine_quotes/Q1870545_USD_200mg.xlsx', orig_name_col='Customer Code')

recipe = hippo.Recipe.from_json(animal.db, 'recipe_jun3_80bases_from35instock.json', allow_db_mismatch=True)

gen = hippo.RandomRecipeGenerator(db=animal.db, start_with=recipe, suppliers=['Stock', 'Enamine'])

from tqdm import tqdm

for i,c in tqdm(enumerate(gen.product_pool), total=len(gen.product_pool)):
	
	try:
		reactions = c.reactions
	except DatabaseError:
		logger.error(f"Error getting {c}'s reactions")
		continue

	for reaction in reactions:

		try:
			recipes = reaction.get_recipes()
		except DatabaseError:
			logger.error(f"Error getting {reaction}'s ({c}) recipes")
			continue

		for recipe in recipes:

			route = animal.register_route(recipe=recipe)

			# logger.info(f'registered {route=}')

	if i%100 == 0:
		logger.success('Committing...')
		animal.db.commit()

animal.db.close()
