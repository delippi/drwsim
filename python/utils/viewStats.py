import pstats
out='../out/'
stats = pstats.Stats(out+'main.profile')
stats.strip_dirs().sort_stats('time').print_stats()
