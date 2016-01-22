#also see temp.py for work on summarizing the output
import csv, glob, os
import numpy as np
from concurrency_figures import DegreeTab

def simav2csv(dirpath, fpattern, outname, force=False):
	"""Compute  average of the replicates
	(where the replicates are in a given `dirpath` and share `fpattern`).
	The results files are assumed to have a single header line.
	Past computation of simav used if possible,
	unless `force` is True.
	"""
	assert outname.endswith('csv')
	#full path to results files
	patternpath = os.path.join(dirpath, fpattern)
	outpath = os.path.join(dirpath, outname) 
	#get all completed replicates of this simulation
	fnames = glob.glob(patternpath)
	n_sims = len(fnames)
	if n_sims == 0:
		print 'WARN: no data for ', patternpath
		return None
	elif n_sims < 100:
		print 'WARN: incomplete data for ', patternpath
		print 'Only {0:d} simulations run. (Computing average anyway!)'.format(n_sims)
	#check if simav needs updating
	ftimes = [os.path.getmtime(fname) for fname in fnames]
	use_existing = False
	try:
		csvtime = os.path.getmtime(outpath)
		if any((csvtime < ftime) for ftime in ftimes):
			print 'Updating ' + outpath
		else:
			use_existing = not force
	except os.error:
		print 'Creating ' + outpath

	#compute average across simulations
	if use_existing:
		print 'reuse existing simav'
		simav = np.loadtxt(outpath, delimiter=',', skiprows=1)
		point_prevalence = DegreeTab(simav[:,8:20] + simav[:,20:32]).point_prevalence()
		kappa3 = None  #must compute within loop (see below)
		with open(outpath, 'r') as temp:
			csv_header = next(temp).strip()
			headers = tuple(hdr for hdr in csv_header.split(','))
	else:  # compute the average of the replicates
		asum = 0
		kappa3 = 0  #must be computed in the loop!
		point_prevalence = 0 #compute in the loop *if* population changes
		#get the header from the first file (assumed to be identical across results files)
		with open(fnames[0], 'r') as temp:
			csv_header = next(temp).strip()
			headers = tuple(hdr for hdr in csv_header.split(','))
		for fname in fnames: #get the individual replicates one at a time
			anew = np.loadtxt(fname, delimiter=',', skiprows=1)
			asum += anew
			male_partner_dist = anew[:,8:20]
			female_partner_dist = anew[:,20:32]
			all_tab = DegreeTab(male_partner_dist + female_partner_dist)
			assert (all_tab._data==(anew[:,8:20]+anew[:,20:32])).all()
			#get kappa3 for every period: array shape (t,)
			kappa3 += all_tab.kappa3
			#get point prevalence for every period: array shape (t,)
			point_prevalence += all_tab.point_prevalence()
		assert asum.shape in [(250, 2+6+4*12), (250, 2+6+4*12+2)]
		pptest = DegreeTab(asum[:,8:20]+asum[:,20:32]).point_prevalence()
		simav = asum / float(n_sims)             #deflate sum to average
		"""
		Write `a` to dirpath, where `a` is simulation average.
		(Average of 100 replications of simulation for this `eps`.)
		"""
		#this data set is not huge, so make and write a string representing all of it
		output = [csv_header]
		output.extend((','.join(str(d) for d in row) for row in simav))
		print 'writing ' + outpath
		with open(outpath, 'w') as fout:
			fout.write('\n'.join(output))
		#create other data used for graphs
		kappa3 /= float(n_sims)                   #deflate sum to average
		point_prevalence /= float(n_sims)         #deflate sum to average
		assert np.allclose(point_prevalence, pptest)
	return dict(headers=headers, simav=simav, point_prevalence=point_prevalence, kappa3=kappa3)


import xlwt #http://www.python-excel.org/



def summarize(outdir):
	"""Return None.
	Summarize a concurrency experiment.
	"""
	#first, turn all our csv files into a workbook
	xls_header = [
	'n_infected_males', 'n_infected_females',
	'malePrimaryTransToday', 'femalePrimaryTransToday',
	'maleAsymptomaticTransToday', 'femaleAsymptomaticTransToday',
	'maleSymptomaticTransToday', 'femaleSymptomaticTransToday',
	'maleDistPartnerships', '', '', '', '', '', '', '', '', '', '', '',
	'femaleDistPartnerships', '', '', '', '', '', '', '', '', '', '', '',
	'maleDistHIVp', '', '', '', '', '', '', '', '', '', '', '',
	'femaleDistHIVp', '', '', '', '', '', '', '', '', '', '', '',
	]
	#a list of csv files and their sheet names
	sheets = ['eps{0:02d}simav'.format(i) for i in range(11)]
	fnames = [os.path.join(outdir,s+'.csv') for s in sheets]
	outname = os.path.basename(outdir) + '.xls'
	xlspth = os.path.join(outdir, outname)
	csv2wb(xlspth, fnames, sheets, header=xls_header, skiplines=1)


def csv2wb(xlsout, fpaths, sheetnames=None, header=None, skiplines=0):
	"""Return None. Merge CSV files in fpaths to a single .xls file,
	whose path is given by xlsout. (Silently overwritten!!)
	Use `skiplines` to skip header lines in the CSV files."""
	wbk = xlwt.Workbook()
	if sheetnames is None: #use root filename
		sheetnames = [os.path.basename(f) for f in fpaths]
		sheetnames = [os.path.splitext(f)[0] for f in sheetnames]
	else:
		assert len(fpaths)==len(sheetnames)
	for pth,shtnm in zip(fpaths, sheetnames):
		sheet = wbk.add_sheet(shtnm)
		hrows = 0 #number of header rows
		if header is not None:
			for (c,hdr) in enumerate(header):
				sheet.write(0,c,hdr)
			hrows += 1
		with open(pth,'r') as fin:
			indata = csv.reader(fin)
			for _ in range(skiplines):
				next(indata)
			for rowct, line in enumerate(indata):
				for colnum, datum in enumerate(line):
					sheet.write(hrows+rowct, colnum, float(datum))
	print 'writing summary to ', xlsout
	wbk.save(xlsout)

####### code for table 2
def table4simav(dirpath, title):
	"""Return None. Provides a summary table for one simulation.
	"""
	table = [title]
	table_data = dict()  #map eps -> y2pp (a dict)
	years = range(249)
	for eps in range(11):
		sheetname = 'eps{0:02d}simav'.format(eps)
		fname = os.path.join(dirpath , sheetname + '.csv')
		fname = os.path.normpath(fname)
		y2pp = dict()  #map year -> data for this eps
		table_data[eps] = y2pp
		with open(fname, 'r') as fin:
			next(fin)  #discard header line
			for rownum, line in enumerate(fin):
				year = rownum
				if year in years:
					line = line.split(',')

					male_partner_dist = np.array([map(float, line[8:20])])
					female_partner_dist = np.array([map(float, line[20:32])])
					all_tab = DegreeTab(male_partner_dist + female_partner_dist)
					#get point prevalence (ok to compute from averaged data)
					concurrency_point_prevalence = all_tab.point_prevalence()

					total_inf = float(line[0]) + float(line[1])
					pp = total_inf / 200  #200=20000/100  -> point prevalance as pct
					y2pp[year] = pp, round(100*concurrency_point_prevalence[0])

	header = 'year  '
	for eps in range(11):
		header += '{0:10}'.format('eps{0:02d}'.format(eps) )
	table.append(header)
	for year in years:
		data = '{0:4d}'.format(year)
		for eps in range(11):
			data += '{0:6.2f}({1:2.0f})'.format(*table_data[eps][year])
		table.append(data)
	return '\n'.join(table)


if __name__ == '__main__':
	for duration in 1,2,3,5, 10, 15: #partnership duration in years
		outdir = "C:/temp/out/dilute25dur/duration{:02d}".format(duration)
		summarize(outdir)
		title = '\n\nEHG (sexfreq=100), but with average partnership duration of {0:02d} years'.format(duration)
		print table4simav(outdir, title)

