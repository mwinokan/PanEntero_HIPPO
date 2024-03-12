#!/usr/bin/env python3

import hippo3 as hippo
from pathlib import Path
import pandas as pd
import numpy as np
import molparse as mp
import os
import json
from tqdm import tqdm

import logging
logger = logging.getLogger('load_syndirella_outputs')

os.system('cp -v 2a_hits.sqlite batch0.sqlite')
# os.system('cp -v 2a_hits.sqlite full_batch0.sqlite')

logger.debug(f'{hippo.__file__=}')

animal = hippo.HIPPO('PanEntero_HIPPO','batch0.sqlite')

inspiration_lookup = {
	"A71EV2A-x0310_0A": "A71x0310a",
	"A71EV2A-x0416_0A": "A71x0416a",
}

def main():
	
	print(animal.compounds)
	print(animal.poses)

	animal.alias_dict = {}
	
	batch0_root = Path('/data/xchem-fragalysis/kfieseler/A71EV2A_run2/batch0')
	parse_syndirella_products_batch(batch0_root)

	print(animal.compounds)
	print(animal.poses)

def parse_syndirella_products_csv(df: pd.DataFrame,
								  product_csv: bool,
								  base_output_path: str):
	"""
	Parse a products csv within a base compound. 
	:param df: dataframe of the products csv (could be internal step or final step). 
	:param product_csv: bool if it is the final products csv or not.
	:param base_output_path: path to output folder for the base compound. 
	:return to_hippo_df: pd.DataFrame with all of the reactants and products within that batch folder.
	"""

	for i, row in tqdm(df.iterrows()):
		comp_metadata: dict = {}
		pose_metadata: dict = {}

		tags = []
		
		if product_csv:

			if np.isnan(row['∆∆G']):
				continue # there was no successful placement
			else:
				pose_metadata['placement_ddG'] = row['∆∆G']
				pose_metadata['placement_mRMSD'] = row.comRMSD

				if pose_metadata['placement_ddG'] > 0:
					continue
				
				if pose_metadata['placement_mRMSD'] > 2:
					continue
				
			# inspiration
			inspiration1: str = inspiration_lookup[eval(row.regarded)[0]]
			inspiration2: str = inspiration_lookup[eval(row.regarded)[1]]
			
			# pose_path
			name = row['name']
			pose_path = os.path.join(str(base_output_path), name, f'{name}.minimised.mol')

		else:

			logger.debug('not product_csv')
			
			inspiration1 = None
			inspiration2 = None
			pose_path = None

		if i%5000 == 0:
			logger.debug(f'row {i=}')
			
		# register the reactants
		reactant1 = animal.register_compound(smiles=row.r1_smiles, return_compound=False, commit=False)
		
		# r2 reactant could not exist yet
		if 'r2_smiles' in row:
			reactant2 = animal.register_compound(smiles=row.r2_smiles, return_compound=False, commit=False)
			reactants = [reactant1, reactant2]
		else:
			reactants = [reactant1]
			
		# register the product
		product_smiles: str = row.smiles
		
		pose_name = row['name']

		split_pose = pose_name.split('-')
		comp_name = '-'.join(split_pose[:-3]+[split_pose[-2]])
		comp_metadata['syndirella_name'] = comp_name


		product = animal.register_compound(smiles=product_smiles, return_compound=True, metadata=comp_metadata, commit=False)

		animal.alias_dict[comp_name] = product.id

		### CHANGE TO DIRECT DB MANIPULATION?

		if 'base' in pose_name:
			product.tags.add('base', commit=False)

		elif 'base_compound' in row:
			base = animal.alias_dict[row.base_compound]
			product.set_base(base, commit=False)
			product.tags.add('elab', commit=False)
		
		### END CHANGE
		
		# register the reaction
		if 'flag' in row:
			if not isinstance(row['flag'],float):
				comp_metadata['flag'] = row.flag

		# register the reaction
		reaction = animal.register_reaction(type=row.reaction, product=product, reactants=reactants, commit=False)

		inspirations = []
		
		if inspiration1:
			inspirations.append(animal.poses[inspiration1])

		if inspiration2:
			inspirations.append(animal.poses[inspiration2])

		if not pose_path:
			print(row)
			raise Exception

		# register the pose
		pose = animal.register_pose(compound=product, name='A', target='A71EV2A', path=pose_path, metadata=pose_metadata, inspirations=inspirations, commit=False)

		# if i > 1000:
		# 	break

	logger.title('committing...')
	animal.db.commit()

	return

def parse_syndirella_products_batch(path: str):
	"""
	Parse all base compounds within a batch. 
	:param path: path to batch folder. 
	:return batch_df with all of the reactants and products within that batch folder. 
	"""

	logger.var('path',path)

	files = list(path.iterdir())

	logger.var('#subdirectories to check', len(files))

	# loop over subdirectories
	for i,subdirectory in enumerate(files):
		if not subdirectory.is_dir():
			continue

		# if i < 3:
		# 	continue
		
		logger.title(f'{subdirectory.name} #{i+1}/{len(files)}')

		product_csvs = list(subdirectory.glob('**/*products*.csv'))

		product_csvs = [c for c in product_csvs if '.ipynb_checkpoints' not in str(c)]
		
		product_csvs = [c for c in product_csvs if str(c).endswith('_placements.csv') or not is_last_step(c)]

		logger.var('#product_csvs', len(product_csvs))

		for j,product_csv in enumerate(product_csvs):

			logger.header(f'{product_csv.name} #{j+1}/{ len(product_csvs)}')

			df = pd.read_csv(product_csv, low_memory=False)

			has_regarded = 'regarded' in df.columns

			base_output_dir= subdirectory / 'output'
			parse_syndirella_products_csv(df, 
										   has_regarded, 
										   base_output_dir)

		if i > 1:
			break

	return

def is_last_step(csv):
	a,b = csv.name.split('_')[-1].removesuffix('.csv').split('of')
	return int(a) == int(b)

if __name__ == "__main__":
	main()
