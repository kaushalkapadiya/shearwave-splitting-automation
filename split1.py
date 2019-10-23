from obspy import geodetics
from obspy import read
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.io.sac.util import get_sac_reftime
import splitwavepy as sw
from obspy.taup import TauPyModel
import os

#window size
windowsize=10

#filters f1=min f2=max
f1=[0.01,0.04,0.02,0.01,0.01,0.02,0.01,0.02,0.01,0.05]
f2=[0.1,0.125,0.3,0.3,0.4,1.0,0.15,0.25,1.0,0.2]

#parameters to trim around skstime
minsks=15
maxsks=30



def bestvalue(t):

    d=1000
    td=1000
    index=1000
    ti=-1
    for i in range(len(t)):
        if(t[i][1]>=10):
            if(d>t[i][1]):
                d=t[i][1]
                td=t[i][3]
                index=t[i][4]
                ti=i
            elif(d==t[i][1] and t[i][2]<td):
                td=t[i][3]
                index=t[i][4]
                ti=i
            else:
                continue
    if(index==1000):
        raise Exception("not possible")
    return index,ti



#read data

def read_data(path):
    try:
        dfile=os.listdir('C:/Users/kkapa/Desktop/RAJ-SKS')
        loop=len(dfile)
        inc=0
    except:
        raise Exception("file not found")

    while(inc!=loop-1 and inc!=loop-2 and inc!=loop-3):

        s1,s2,s3='','',''
        refile=''

        while(True):

            try:
                ext=dfile[inc][-3::1]
                if(ext=='sac' or ext=='SAC'):
                    s1+=dfile[inc]
                    refile+=dfile[inc]
                    inc+=1
                    # print(s1)
                    break
                else:
                    inc+=1
            except:
                inc+=1

        while(True):
            try:
                ext=dfile[inc][-3::1]
                if(ext=='sac' or ext=='SAC'):
                    s2+=dfile[inc]
                    inc+=1
                    # print(s2)
                    break
                else:
                    inc+=1
            except:
                inc+=1

        while(True):
            try:
                ext=dfile[inc][-3::1]
                if(ext=='sac' or ext=='SAC'):
                    s3+=dfile[inc]
                    inc+=1
                    # print(s3)
                    break
                else:
                    inc+=1
            except:
                inc+=1

        st=read(path+'/'+s1, debug_headers=True)+read(path+'/'+s2, debug_headers=True)+read(path+'/'+s3, debug_headers=True)


        # st = read('2011.052.10.57.52.4000.XX.KTL.00.BHE.M.sac', debug_headers=True)
        # st += read('2011.052.10.57.52.4000.XX.KTL.00.BHN.M.sac', debug_headers=True)
        # st += read('2011.052.10.57.52.4000.XX.KTL.00.BHZ.M.sac', debug_headers=True)
        #2011.052.10.57.52.4000.XX.KTL.00.BHZ.M   2011.052.10.57.52.4000.XX.KTL.00.BHN.M    2011.052.10.57.52.4000.XX.KTL.00.BHE.M
        # for
        #extarcting the required information form the data
        tr=st[0]
        # print(tr.stats)
        # st.plot()
        evtime=tr.stats['starttime']
        endtime=tr.stats['endtime']
        tr=tr.stats['sac']
        tr=dict(tr)
        b=tr['b']
        evla=tr['evla']
        evlo=tr['evlo']
        stla=tr['stla']
        stlo=tr['stlo']
        evdp=tr['evdp']
        model = TauPyModel('iasp91')
        arrivals = model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,phase_list=['SKS'])
        skstime = evtime+arrivals[0].time-b
        # print(skstime)




        #applying filters

        dist, az, baz = geodetics.base.gps2dist_azimuth(evla,evlo,stla,stlo)

        figurefile=path+'/'+refile[0:30]
        resultfile=refile[0:30]+'_results.txt'
        resultfile=path+'/'+resultfile
        f=open(resultfile,'w')
        # f=open('2011.052.10.57.52'+'_result1.txt','w')
        f.write('  EventId'+'\t    '+'Baz'+'\t\t'+' filter'+'\t\t'+'SI'+'\t'+'Split/Null'+'\t\t\t'+'EigenM'+'\t\t\t'+' TransM'+'\t\t\t'+' CrossM'+'\t\t\t\t'+'\n')
        f.write('2011.052.10.57.52  '+str(round(baz,2))+'\t'+'f1'+'\t'+'f2'+'\t'+'\t'+'\t\t'+'|'+'phi'+'\t'+'dev'+'\t'+'t'+'\t'+'dt'+'\t'+'|'+'phi'+'\t'+'dev'+'\t'+'t'+'\t'+'dt'+'\t'+'|'+'phi'+'\t'+'dev'+'\t'+'t'+'\t'+'dt'+'\t'+'\n')
        f.close()


        for j in range(len(f1)):

            st.filter("bandpass",freqmin=f1[j],freqmax=f2[j])
            # st.plot()
            # trim around SKS
            st.trim(skstime-minsks,skstime+maxsks)
            # st.plot()

            #creating pair
            north = st[1].data
            east = st[0].data
            sample_interval = st[0].stats.delta
            # print(sample_interval)
            realdata = sw.Pair(north, east, delta=sample_interval)
            si=realdata.splitting_intensity()
            x,y=realdata.cordinatewindow()
            diff=int(y-x)-windowsize
            # realdata.plot()

            try:
                #initial EigenM
                measure = sw.EigenM(realdata)
                temp=measure.measurements()
                print(-1,temp)
                temp=list(temp)
                temp.append(-1)
                m=[]
                m.append(temp)
                #initial TransM
                measure1 = sw.TransM(realdata, pol=baz)
                temp=measure1.measurements()
                print(-1,temp)
                temp=list(temp)
                temp.append(-1)
                m1=[]
                m1.append(temp)

                #initial CrossM
                measure2 = sw.CrossM(realdata)
                temp=measure2.measurements()
                print(-1,temp)
                temp=list(temp)
                temp.append(-1)
                m2=[]
                m2.append(temp)
            except:
                print("please check this data manually canot apply filter ",j+1)
                print("for file names")
                print(s1,s2,s3)
                continue




            #setting windows

            for i in range(diff):

                a=realdata
                # a.plot()

                try:
                    a.set_window(x+i,x+windowsize+i)
                    # a.plot()

                    try:
                        #EigenM
                        measure = sw.EigenM(a)
                        temp=measure.measurements()
                        print(i,measure.measurements())
                        temp=list(temp)
                        temp.append(i)
                        m.append(temp)

                        #TransM

                        measure1 = sw.TransM(realdata, pol=baz)
                        temp=measure1.measurements()
                        print(i,measure1.measurements())
                        temp=list(temp)
                        temp.append(i)
                        m1.append(temp)

                        #CrossM

                        measure2 = sw.CrossM(a)
                        temp=measure2.measurements()
                        print(i,measure2.measurements())
                        temp=list(temp)
                        temp.append(i)
                        m2.append(temp)
                    except:
                        continue
                except:
                    continue
            try:
                index,ti=bestvalue(m)
                print(index)
                phi,dev,t,dt=m[ti][0],m[ti][1],m[ti][2],m[ti][3]
                a=realdata
                if(index==-1):
                    a.set_window(x,y)
                else:
                    a.set_window(x+index,x+windowsize+index)
                # a.plot()
                measure = sw.EigenM(a)
                # measure.plot()
                fname=figurefile+'_EigenM_'+str(j)+'.pdf'
                # to save the plot pass 'save' and file name
                # to save and show the plot pass 'showandsave' and file name
                #to show the plot pass nothing

                measure.plot('save',fname)

                # f=open('2011.052.10.57.52'+'_result.txt','a')
                # f.write('eigenm'+'\t'+str(f1[j])+'\t'+str(f2[j])+'\t'+str(phi)+'\t'+str(dev)+'\t'+str(t)+'\t'+str(dt)+'\n')
                # f.close()
            except:
                print("not saved")
                continue


            try:
                index,ti=bestvalue(m1)
                print(index)
                phi1,dev1,t1,dt1=m1[ti][0],m1[ti][1],m1[ti][2],m1[ti][3]

                a=realdata
                if(index==-1):
                    a.set_window(x,y)
                else:
                    a.set_window(x+index,x+windowsize+index)
                # a.plot()
                measure1 = sw.TransM(a, pol=baz)
                # measure1.plot()
                fname=figurefile+'_TransM_'+str(j)+'.pdf'
                measure1.plot('save',fname)

                # f=open('2011.052.10.57.52'+'_result.txt','a')
                # f.write('transm'+'\t'+str(f1[j])+'\t'+str(f2[j])+'\t'+str(phi)+'\t'+str(dev)+'\t'+str(t)+'\t'+str(dt)+'\n')
                # f.close()
            except:
                continue


            try:
                index,ti=bestvalue(m2)
                print(index)

                phi2,dev2,t2,dt2=m2[ti][0],m2[ti][1],m2[ti][2],m2[ti][3]

                a=realdata
                if(index==-1):
                    a.set_window(x,y)

                else:
                    a.set_window(x+index,x+windowsize+index)
                # a.plot()
                measure2 = sw.CrossM(a)
                # measure2.plot()
                fname=figurefile+'_CrossM_'+str(j)+'.pdf'
                measure2.plot('save',fname)

                # f=open('2011.052.10.57.52'+'_result.txt','a')
                # f.write('CrossM'+'\t'+str(f1[j])+'\t'+str(f2[j])+'\t'+str(phi)+'\t'+str(dev)+'\t'+str(t)+'\t'+str(dt)+'\n')
                # f.write('\n')
                # f.close()
            except:
                continue
            try:
                f=open(resultfile,'a')
                f.write('\t'+'\t\t\t'+str(f1[j])+'\t'+str(f2[j])+'\t'+str(round(si,2))+'\t\t\t'+'|'+str(round(phi,3))+'\t'+str(round(dev,3))+'\t'+str(round(t,3))+'\t'+str(round(dt,3))+'\t'+'|'+str(round(phi1,3))+'\t'+str(round(dev1,3))+'\t'+str(round(t1,3))+'\t'+str(round(dt1,3))+'\t'+'|'+str(round(phi2,3))+'\t'+str(round(dev2,3))+'\t'+str(round(t2,3))+'\t'+str(round(dt2,3))+'\t'+'\n')

                f.close()
            except:
                print("canot write")
                continue






#enter data

# data=str(input("enter the path of the dir"))
read_data('C:/Users/kkapa/Desktop/RAJ-SKS')
