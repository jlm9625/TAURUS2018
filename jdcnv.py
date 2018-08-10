from astropy.time import Time

def jdcnv(year,month,day,hour,min,sec):
    
    mm = str(month)
    if len(mm) != 2: mm = '0'+mm
    dd = str(day)
    if len(dd) != 2: dd = '0'+dd
    hh = str(hour)
    if len(hh) != 2: hh = '0'+hh

    mn = str(min)
    if len(mn) != 2: mn = '0'+mn
    
    if sec < 10.0:
        ss = '0'+str(sec)
    else: ss = str(sec)
        
    timestamp = str(year) + '-'+mm+'-'+dd+'T'+hh+':'+mn+':'+ss
    #pdb.set_trace()

            
    times = [timestamp]
    t = Time(times,format='isot',scale='utc')
    thisjd = t.jd.copy()
    
    return thisjd[0]
    

#def k2t0(cnum,t0):
    
#    if cnum == 4: timestamp_start =     
    