import lc_read

data_fname = 'tests/data/BUGGY.fix.lch'
data_start, data_end = lc_read.get_data_timelimits(data_fname)
print('{} has data from {} to {}'.format(data_fname,data_start, data_end))
stream = lc_read.read(data_fname,chan_map='SPOBS2')
print(stream)
#print(stream[0].data[:10])
stream.plot(size=(800, 600), equal_scale='False', method='full')