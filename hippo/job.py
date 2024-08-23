
import hippo
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import os

os.system('cp -v 2a_hits.sqlite may2_job.sqlite')

animal = hippo.HIPPO('may2_job', 'may2_job.sqlite')

map = {pose.metadata['observation_longname']:pose.id for pose in animal.poses.get_by_tag('hits') if ('observation_longname' in pose.metadata)}

for i,file in enumerate(Path('/data/xchem-fragalysis/kfieseler/A71EV2A_run4').glob('??????????????-??????????-?/*_to_hippo.pkl.gz')):
    print(i, file)
    animal.add_syndirella_elabs(file, inspiration_map=map)

animal.db.close()
