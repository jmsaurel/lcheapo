XX.TEST.2019-11-07.mseed was created using:
    s = lcread(20191107T14_SPOBS09_F02.raw.lch, station='TEST', network='XX',
               obs_type='SPOBS2')
    s.write('test.mseed', 'MSEED', encoding='STEIM1', byteorder='<')
