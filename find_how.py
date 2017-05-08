import json, os, sys

sp = 'C2H4'
reactor = 'autoignition'

for fat in os.listdir(reactor):
	if 'phi' in fat:
		fd = os.path.join(reactor, fat)
		for pnt in os.listdir(fd):
			if 'point' in pnt:
				sp_kept = json.load(open(os.path.join(fd, pnt)))['species']
				if sp in sp_kept.keys():
					print fd
					print sp_kept[sp]['by_K']
					sys.exit()