from obspy import read, read_inventory
import spectral as sp
import sys

## Set variables
# DOCTAR data
datafile_base='../0_Data/DOCTAR/D01.DO.???..D.2011.236.000000';
metafile='../0_Data/DOCTAR/DOCTAR_inventory.xml';
length=12*3600; # 12 hours, avoids EQ at 2011-08-24T18
drivechan='3'
respchan='4'
noisechan='equal'
# Marianas data
datafile_base='../Marianas/XF.B12.2012_10_01-05_2Hz.mseed';
metafile='../Marianas/XF.B12.station.xml';
length=None  # read all
drivechan='1'
respchan='Z'
noisechan='response'

def main():
    ## Read in data and responses
    print('Reading channels ',end='')
    st=read(datafile_base)
    for tr in st:
        print(tr.stats.channel,end=', ')
    print('and prepping meta/data')
    inv=read_inventory(metafile)
    st.attach_response(inv)
    # Cut data to have start at same time on all channels
    dt=st[0].stats.starttime+60
    if length:
        st.trim(dt,dt+length)

    ## Calculate spectra and transfer functions
    print('Calculating PSDs...')
    s=sp.calc_PSDs(st,window_length=1000);
    print('Calculating coherences...')
    c=sp.calc_cohers(st,window_length=1000);
    print('Calculating transfer functions...')
    x=sp.calc_XF(s,c,drivechan,respchan,noisechan)

    ## Plot all
    print('Plotting...')
    header='PLOTS/{}.{}'.format(st[0].stats.network,st[0].stats.station)
    st.plot(size=(800,800),outfile='{}_time_obspy.png'.format(header))
    sp.plot_spects(s,outfile='{}_spect_spectral.png'.format(header))
    sp.plot_cohers(c,outfile='{}_coher_spectral.png'.format(header))
    sp.plot_XF(x,outfile='{}_XF_spectral.png'.format(header))

    ## Plot PPSD
#     from obspy.signal import PPSD
#     ppsd = PPSD(st[1].stats, inv)
#     print(ppsd.id)
#     ppsd.add(st[1]) 
#     ppsd.plot() 

    
if __name__ == "__main__": main()