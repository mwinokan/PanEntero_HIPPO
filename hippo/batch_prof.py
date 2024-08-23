
import hippo
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import os

os.system('cp -v 2a_hits.sqlite apr29_prof.sqlite')

animal = hippo.HIPPO('apr29_prof', 'apr29_prof.sqlite')

map = {pose.metadata['observation_longname']:pose.id for pose in animal.poses.get_by_tag('hits') if ('observation_longname' in pose.metadata)}

# find all elab pickles:
for i,file in enumerate(Path('/data/xchem-fragalysis/kfieseler/A71EV2A_run4').glob('??????????????-??????????-?/*_to_hippo.pkl.gz')):
    # print(i, file)
    animal.add_syndirella_elabs(file, inspiration_map=map, stop_after=1000)
    break

animal.db.close()
