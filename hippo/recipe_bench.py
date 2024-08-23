
import hippo

animal = hippo.HIPPO('may2_job_quote', 'may2_job_quote.sqlite')

bases = animal.compounds(tag='Syndirella base')[:10]

import cProfile, pstats, io
from pstats import SortKey

with cProfile.Profile() as pr:
	recipe = bases.get_recipe()

	# p = pstats.Stats(OUTPUT_FILE)
	pr.dump_stats('recipe_bench.prof')

	s = io.StringIO()
	sortby = SortKey.CUMULATIVE
	ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
	ps.print_stats()
	print(s.getvalue())

	# pr.print_stats()
	# pr.sort_stats('time').print_stats(10)

