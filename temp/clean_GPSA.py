import os, shutil

prj_name = os.path.join('prj_mix','prj_mix_dec+tol')
#prj_name = 'prj_H2'
name2del = 'GPSA'

fld_root = os.path.join('prj',prj_name,'detailed','raw')

for fo in os.listdir(fld_root):
	fld_fo = os.path.join(fld_root,fo)
	if not os.path.isdir(fld_fo): continue
	for reactor in os.listdir(fld_fo):
		fld_reactor = os.path.join(fld_fo,reactor)
		if not os.path.isdir(fld_reactor): continue
		for fat in os.listdir(fld_reactor):
			fld_fat = os.path.join(fld_reactor, fat)
			if not os.path.isdir(fld_fat): continue
			
			f_json = os.path.join(fld_fat, name2del+'.json')
			if os.path.exists(f_json):
				os.remove(f_json)
				print 'deleted '+str(f_json)
			
			fld = os.path.join(fld_fat, name2del)
			if os.path.exists(fld):
				shutil.rmtree(fld)
				print 'deleted '+str(fld)