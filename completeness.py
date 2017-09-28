'''
This is program written by
	Anis Mohammed Vengasseri
	anis.mhd@gmail.com
	reacamv.com

'''
from numpy import linspace, ceil, zeros, cumsum, loadtxt, sqrt
import matplotlib.pyplot as plt
import sys
def num_items_array(in_array, num_yr, max_yr, dT):
	freq_bin = zeros(num_yr, dtype=int)
	for i in range(num_yr):
		min_val = max_yr - (i+1)*dT
		max_val = max_yr - i*dT
		freq_bin[i] = len( in_array[ (in_array>min_val)*(in_array<=max_val) ] )
	return freq_bin

def completeness(year, magnitude, Mmin, Mmax, dT, dM):
	if not (len(year) == len(magnitude)):
		print "number of year and magnitude are not matching... exiting..."
		return None
	max_yr = int(max(year))
	num_yr = int(ceil((max(year)-min(year))/dT))
	min_yr = int(max_yr - num_yr*dT)
	dT = int(dT)
	num_bin = int(ceil((Mmax-Mmin)/dM))
	mag_bin = {}
	mag_bin[0] = {}
	mag_bin[0]['min magnitude'] = 0.0
	mag_bin[0]['max magnitude'] = Mmin
	cum_T = cumsum(zeros(num_yr) + dT)
	mag_bin[0]['year bin'] = num_items_array( year[magnitude<=Mmin] , num_yr, max_yr, dT)
	mag_bin[0]['year bin cumsum'] = cumsum(mag_bin[0]['year bin'])
	mag_bin[0]['standard deviation'] = sqrt(cumsum(mag_bin[0]['year bin']))/cum_T
	for i in range(1,num_bin+1):
		mag_bin[i] = {}
		mag_bin[i]['min magnitude'] = Mmin + (i-1)*dM
		mag_bin[i]['max magnitude'] = min(Mmin + i*dM,Mmax)
		mag_bin[i]['year bin'] = num_items_array( year[(magnitude>mag_bin[i]['min magnitude'])*(magnitude<=mag_bin[i]['max magnitude'])] , num_yr, max_yr, dT)
		mag_bin[i]['year bin cumsum'] = cumsum(mag_bin[i]['year bin'])
		mag_bin[i]['standard deviation'] = sqrt(cumsum(mag_bin[i]['year bin']))/cum_T
	header = '{0:12s},<={1:10.2f}'.format('Year',Mmin)
	for i in range(1,num_bin+1):
		header = header + ',M{0:5.2f}-{1:5.2f}'.format(mag_bin[i]['min magnitude'],mag_bin[i]['max magnitude'])
	print header
	for i in range(num_yr):
		line_str = 'Y{0:5d}-{1:5d}'.format(max_yr-(i+1)*dT, max_yr-i*dT)
		for j in range(0,num_bin+1):
			line_str = line_str + ',{0:12d}'.format(mag_bin[j]['year bin'][i])
		print line_str
	print "Cumulative Number of Earthquakes"
	print header
	for i in range(num_yr):
		line_str = 'Y{0:5d}-{1:5d}'.format(max_yr-(i+1)*dT, max_yr-i*dT)
		for j in range(0,num_bin+1):
			line_str = line_str + ',{0:12d}'.format(mag_bin[j]['year bin cumsum'][i])
		print line_str
	print "Data to be used for completeness check"
	header = '{0:12s},{2:12s},{3:12s},<={1:10.2f}'.format('Year',Mmin,'cumulT','1/SQRT(T)')
	for i in range(1,num_bin+1):
		header = header + ',M{0:5.2f}-{1:5.2f}'.format(mag_bin[i]['min magnitude'],mag_bin[i]['max magnitude'])
	print header
	for i in range(num_yr):
		line_str = 'Y{0:5d}-{1:5d},{2:12d},{3:12.10f}'.format(max_yr-(i+1)*dT, max_yr, int(cum_T[i]), 1/sqrt(cum_T[i]))
		for j in range(0,num_bin+1):
			line_str = line_str + ',{0:12.7f}'.format(sqrt(mag_bin[j]['year bin cumsum'][i])/cum_T[i])
		print line_str
	'''
	for i in range(0,num_bin+1):
		plt.close('all')
		plt.figure()
		plt.loglog(cum_T, 1/sqrt(cum_T))
		plt.scatter(cum_T, mag_bin[i]['standard deviation'], 'r')
		plt.grid(b=True, which='both', color='0.65',linestyle='--')
		plt.show()
		fname = 'completeness{0:d}.pdf'.format(i)
	'''
if __name__ == '__main__':
	A = loadtxt(sys.argv[1])
	completeness(A[:,0], A[:,3], 4.0, 8.0, 10, 1.0)