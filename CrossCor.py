import numpy as np
import glob
import sys
import getopt
import os
from astropy.io import fits
from scipy.interpolate import interp1d,splrep,splev,sproot
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,SpanSelector,CheckButtons
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
#from stsci.convolve import boxcar 
from astropy.convolution import convolve, Box1DKernel

#Usage: run CrossCor.py imagename ('objid')
#Usage: run CrossCor.py imagename ('unsure')
#Usage: run CrossCor.py imagename ('inter')
'''
Inputs:
imagename - string -  set to be your big.fits file, including the path.
objid - optional string -  set to the 6-letter object id of a single object to open only that one
usure - optional string -  set to the string 'unsure' to open only those objects which the user flagged as unsure
inter - optional string -  set to the string 'inter' to run interactively; this is automatically on if looking at 1 objid
emiss - optional string -  set to the string 'emiss' to run through a list of emission-line only galaxies

Change to your dataset:
wave1,wave2 - (int,int) - the wavelength range over which to correlate the data, can be changed by the user in the GUI
zrange - (int,int) - the range of redshifts over which to search
emclip - boolean - set True to clip emission lines by default, or False to keep them in. 
                   This setting can always be toggled within the GUI, but this value is what displays first.
specplot - boolean - set to 1 if you want a second plot showing common spectral lines overlaid on the galaxy at the currect redshift. 
                     Can be kind of clunky, so default is off, but could be useful to enable for 'unsure' galaxies. 
                     Can be toggled in the GUI. Requires galaxlylines.dat in the same folder as your image

Other:
outloc - string - where your files will be output, defaults to the same as input
temploc - string - where you store your templates. This defaults to you imagename location
outname - string - name of your output files. Defaults to 'cc_' + (your big.fits name). e.g. 'a' creates a.txt and verb_a.txt.
imloc - string - this will be set to the path to the big.fits file automatically
imname - string - this will be set to just the big.fits file automatically
tcorr - boolean - will be automatically set to 0 if the correction code has not been run.


NOTE: The outfile must correspond to the image that it was generated from. 
      The file is generated once with all objects in the mask if it does not exist, and then is modified from there.
'''


wave1,wave2 = (3900,8200)
zrange = (0.01,0.4)
emclip = True
markbads = False #Read in a separate file and auto-mark bads
haselist = False #Read in list of objects to assume emission
specplot = 0

outloc = 0
outname = 0
temploc = sys.path[0]+'/Templates/'
readunsure = 0
interactive = 0
objid = 0
tcorr = 1 

fullCmdArguments = sys.argv
argumentList = fullCmdArguments[2:]
imagename = fullCmdArguments[1]
unixOptions = "iueo:b:l:t:"  
gnuOptions = ["inter", "unsure","emiss","objid=","bads=","elist=","temploc="]  

try:  
    arguments, values = getopt.getopt(argumentList, unixOptions, gnuOptions)
except getopt.error as err:  
    # output error, and return with an error code
    print (str(err))
    sys.exit(2)

for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--inter"):
        interactive = 1
    elif currentArgument in ("-u", "--unsure"):
        readunsure = 1
    elif currentArgument in ("-o", "--objid"):
        objid = currentValue
        interactive = 1
    elif currentArgument in ("-e", "--emiss"):
        emclip = False
    elif currentArgument in ("-b", "--bads"):
        badfile = currentValue
        markbads = True
    elif currentArgument in ("-l", "--elist"):
        emissfile = currentValue
        haselist = True
    elif currentArgument in ("-t", "--temploc"):
        temploc = currentValue

imname = imagename
while imname.find('/') != -1:
    imname = imname[imname.find('/')+1:]
imloc = imagename.replace(imname,'')
if not temploc: temploc = imloc
if not outloc: outloc = imloc
if not outname: outname = 'cc_' + imname.replace('.fits','')


outfile = outloc + outname + '.txt'
outfilev = outloc + 'verb_' + outname + '.txt'

outfile = outloc + outname + '.txt'
#Trying to fix chanigng to a new computer, set to the mask you are doing
#outfilev = outfile.replace('cc_feb16_abig','verb_a6')
#outfilev = outfile.replace('cc_feb16_bbig','verb_b6')
#outfilev = outfile.replace('cc_feb17_bbig','verb_b7')
#outfilev = outfile.replace('cc_feb16_dbig','verb_d6')
#outfilev = outfile.replace('cc_feb16_ebig','verb_e6')
#outfilev = outfile.replace('cc_feb17_ebig','verb_e7')
#outfilev = outfile.replace('cc_feb16_fbig','verb_f6')
#outfilev = outfile.replace('cc_feb16_gbig','verb_g6')
#outfilev = outfile.replace('cc_feb16_hbig','verb_h6')
#outfilev = outfile.replace('cc_feb16_ibig','verb_i6')
outfilev = outfile.replace('cc_feb17_jbig','verb_j7')

#Check if telluric correction is done
if not glob.glob(imloc+'cor_??????_'+imname):
    tcorr = 0

#Check if the text files is set up. If not, create one
if not os.path.exists(outfile):
    print 'No output text files exists. Creating ' + outfile + ' and ' + outfilev
    if tcorr: imarr = glob.glob(imloc + 'cor_??????_' + imname)
    else: imarr = glob.glob(imloc + '??????_' + imname)
    f = open(outfile,'w+')
    f2 = open(outfilev,'w+')
    f.write('#OBJID  temp    z           dzhi        dzlo        ccmax     chi2    rchi2  eclip')
    f.write('     S/N      Star Bad  Unsure    ImageName            Revisit  Note  Unusable\n')
    f2.write('#OBJID z23 dzhi23 dzlo23 ccmax23 chi223 rchi223 z24 dzhi24 dzlo24 ccmax24 chi224 ')
    f2.write('rchi224 z25 dzhi25 dzlo25 ccmax25 chi225 rchi225 z26 dzhi26 dzlo26 ccmax26 chi226 ')
    f2.write('rchi226 z27 dzhi27 dzlo27 ccmax27 chi227 rchi227 temp eclip S/N Star Bad Unsure ')
    f2.write('ImageName Revisit Note Unusable\n')
    for i in range(len(imarr)):
        objname = imarr[i].replace(imloc,'')
        f.write(('%06d   %d    %.6f    %.6f   ' %(0,0,0,0)+
                 '%.6f    %7.2f   %7.2f   %2.2f     ' %(0,0,0,0)+
                 '%d      %2.2f      %d    %d     ' %(0,0,0,0)+
                 '%d ' % (0)) + objname + 
                ('      %d       %d       %d' % (0,0,0)) + '\n')
        f2.write(('%06d %.6f %.6f %.6f ' %(0,0,0,0)+
                  '%7.3f %7.3f %2.2f %.6f ' %(0,0,0,0)+
                  '%.6f %.6f %7.3f %7.3f ' %(0,0,0,0)+
                  '%2.2f %.6f %.6f %.6f ' %(0,0,0,0)+
                  '%7.3f %7.3f %2.2f %.6f ' %(0,0,0,0)+
                  '%.6f %.6f %7.3f %7.3f ' %(0,0,0,0)+
                  '%2.2f %.6f %.6f %.6f ' %(0,0,0,0)+
                  '%7.3f %7.3f %2.2f %d ' %(0,0,0,0)+
                  '%d %3.2f %d %d ' %(0,0,0,0)+
                  '%d' % (0)) + ' '+ objname + 
                 (' %d %d %d' % (0,0,0)) + '\n')
    f.close()
    f2.close()


skylines = [5577.0,5890.0,6300,6364]
telluric = [(7584,7650)]

def divz(X,Y):
    return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

def svdfit(b,y):
    decomp = np.linalg.svd(b,full_matrices=False)
    sol1 = np.transpose(decomp[2])
    sol2 = divz(1.0,decomp[1])
    sol3 = np.dot(np.transpose(decomp[0]),y)
    if np.sometrue(sol3):
        solr = (sol2*sol3)
        soll = np.dot(sol1,solr)
    else:
        soll = np.zeros(sol3.shape)
    return soll

def legendre(x,nl):
    x = np.asarray(x,dtype='float64')
    l = (nl+1,)+x.shape
    p = np.zeros(l,dtype='float64')
    p[0] = np.ones(x.shape,dtype='float64')
    if nl > 0: p[1] = x
    if nl > 1:
        for j in range(1,nl):
            p[j+1] = ((2*j+1)*x*p[j] - j*p[j-1])/(j+1)
    return p

def outGraph(xdata,ydata,objid,kind,lim=0):
    #kind is a short string describing the graph
    filepath = outloc
    fig, ax1 = plt.subplots(figsize=(10,8))
    ax1.plot(xdata,ydata,color='red')
    ax1.set_title(('%06d' % objid) + ', ' + kind)
    if lim != 0: ax1.set_ylim(lim)
    fig.savefig(filepath + ('%06d' % objid) + '_' + kind)
    plt.close(fig)
    
def FitPoly(s,n,g,w,o,w1,w2):
    u = g*np.greater_equal(w,w1)*np.less_equal(w,w2)
    x = np.arange(len(s))
    xu,su,nu,gu,wu = np.compress(u,[x,s,n,g,w],1)
    xmin,xmax = xu.min(),xu.max()
    xc = 0.5*(xmax+xmin)
    xs = 0.5*(xmax-xmin)
    xl = (xu-xc)/xs
    basis = legendre(xl,o)
    wt = divz(gu,nu)
    co = svdfit(np.transpose(basis)*wt[::,np.newaxis],su*wt)
    xl = (x-xc)/xs
    basis = legendre(xl,o)
    pl = np.dot(co,basis) * np.greater_equal(x,xmin) * np.less_equal(x,xmax) * np.greater_equal(w,w1) * np.less_equal(w,w2)
    return pl

def Smooth(y,good,p=50,h=25):
    m = np.zeros(y.shape,y.dtype)
    for j in range(len(y)):
        a,b = np.clip([j-h,j+h+1],0,len(y)-1)
        u = np.compress(good[a:b],y[a:b])
        if len(u):
            if p == 50:
                m[j] = np.median(u)
            else:
                m[j] = np.sort(u)[len(u)*p/100]
    return m
    
class Spectrum():
    def __init__(self,file,clip=0,verb=0,good=0,changeclip=0):
        self.file = file
        self.bigfile = imloc + imname
        self.verb = verb
        self.clipemission = clip
        self.Read()
        self.bighdu = fits.open(self.bigfile)[0]
        self.image_data = self.bighdu.data[self.head['csecta']:(self.head['csectb']+1)]
        self.lobar = self.head['lowbound']-self.head['csecta']-1
        self.hibar = self.head['upbound']-self.head['csecta']+1
        self.bighead = self.bighdu.header
        if tcorr:
            self.spec = self.data[0]
            self.orig = self.data[1]
            self.noise = self.data[2]
            self.flat = self.data[3]
            self.sky = self.data[4]
            self.tc = self.data[5]
            self.mask = self.data[6]
            self.fit = self.data[7]
        else:
            self.spec = self.data[0]
            self.noise = self.data[1]
            self.flat = self.data[2]
            self.sky = self.data[3]
        self.objid = self.head['OBJID']
        self.findSN()
        self.GetBadPixels(good=good,changeclip=changeclip)
		
    def Read(self):
        self.hdu = fits.open(self.file)[0]
        self.head = self.hdu.header
        self.data = self.hdu.data.astype('float32')
        self.crval1 = self.head["crval1"]
        self.crpix1 = self.head["crpix1"]
        self.cdelt1 = self.head["cdelt1"]
        self.naxis1 = self.head["naxis1"]
        self.dcflag = self.head["dc-flag"]
        self.wavelength = (1.0+np.arange(self.naxis1)-self.crpix1)*self.cdelt1 + self.crval1
        if self.dcflag: self.wavelength = np.power(10.0,self.wavelength)

    def findSN(self):
        sig = np.median(self.spec[450:650])
        noi = np.median(self.noise[450:650])
        if noi:
            self.signoise = sig/noi
        else:
            self.signoise = 0
		    
    def GetBadPixels(self,good=0,changeclip=0):
        if hasattr(self,"noise"):
#            bad = 1.*np.logical_or(np.less_equal(self.spec,1e-4),np.less_equal(self.noise,1e-4))
            bad = 1.*np.logical_or(np.less_equal(self.spec,-3*self.noise),np.less_equal(self.noise,1e-4))
            for skyline in skylines:
                bad = 1.*np.logical_or(bad,np.greater(self.wavelength,skyline-3*self.cdelt1)*np.less(self.wavelength,skyline+3*self.cdelt1))
            for absorb in telluric:
                bad = 1.*np.logical_or(bad,np.greater(self.wavelength,absorb[0])*np.less(self.wavelength,absorb[-1]))
            #bad = 1.*np.greater(boxcar(bad,(8,)),0)       #stsci boxcar
            bad = 1.*np.greater(convolve(bad,Box1DKernel(8)),0) #astropy boxcar
        else:
            bad = np.zeros(self.spec.shape,dtype='float')
        sm = Smooth(self.spec,np.logical_not(bad))
        if self.verb:
            lim = (np.median(self.spec)*0.25,np.median(self.spec)*1.75)
            outGraph(self.wavelength,self.spec,self.objid,'spec',lim=lim)
            outGraph(self.wavelength,sm,self.objid,'sm50',lim=lim)
        sr = self.spec-sm
        ss = 1.49*Smooth(abs(sr),np.logical_not(bad),h=50)
        if self.verb: outGraph(self.wavelength,ss,self.objid,'ss50')
        self.good2 = np.ones(sr.shape,dtype='float32')  # KEEPS EMISSION LINES IN
        if changeclip:
            if self.clipemission:
                self.good = good*np.less(sr,2.5*ss)
            else:
                self.good = good*np.ones(sr.shape,dtype='float32')*np.less(sr,2.5*ss)+np.greater(sr,2.5*ss)
            self.continuum = 1.*np.less(abs(sr),2.5*ss)  # NOTE THAT POLYNOMIAL FIT TO CONTINUUM MASKS OUT EMISSION/STRONG ABSORPTION
        else:
            self.good = np.logical_not(bad)
            if self.clipemission:
                self.good = self.good*np.less(sr,2.5*ss) # CLIPS OUT EMISSION LINES AND BAD POSITIVE SKY LINE RESIDUALS
            else:
                self.good = self.good*np.ones(sr.shape,dtype='float32')  # KEEPS EMISSION LINES IN
            self.continuum = 1.*np.less(abs(sr),2.5*ss)  # NOTE THAT POLYNOMIAL FIT TO CONTINUUM MASKS OUT EMISSION/STRONG ABSORPTION
            if self.verb: outGraph(self.wavelength,self.spec*self.good,self.objid,'sg50',lim=lim)
        self.sm = sm
        self.ss = ss
        
class Template(Spectrum):
    def __init__(self,n,clip=0,verb=0):
        self.file = temploc + "spDR2-%03d.fit" % (n)
        self.verb = verb
        self.clipemission = clip
        self.Read()
        self.spec = self.data[0]
        self.verb = 0
        self.GetBadPixels()
        wu,su = np.compress(self.good,[self.wavelength,self.spec],1)
        self.tck = splrep(wu,su,s=0,k=3)
        self.interp = lambda w: splev(w,self.tck)*np.greater_equal(w,wu.min())*np.less_equal(w,wu.max())
        #self.interp = interp1d(wu,su,fill_value=0,bounds_error=0)
    def Redshift(self,z,w):
        return self.interp(w)/(1.0+z)
        
class Plot():
    def __init__(self,g):
        '''
        Variable initialization
        All of these should be set automatically, but for an explanation:
        image - image location + file, will be read as the spectrum
        objimage - just the image name
        verb - set to 1 to display a few more outputs
        eclip - whether to keep or clip emission lines
        wavechange - used when the user changes the wavelegnth from wave1 or wave2
        zbound - float,float - stores the lower and upper bound for the redshift range
        lowave and hiwave - the starting values of wave1 and wave2, but they get modified as the user changes them
        create - sets to 0 after creating the axes for the first time
        add - toggles between 1 (add) and 0 (subtract) for modifying the mask
        star,baddata,unsure,flag1,flag2,flag3 - 0 or 1, flags that the user can set
        failure - boolean - is set to 1 if the cross correlation fails, then the code outputs only zeros
        mkspec - boolean - toggles whether or not to display the spec plot
        '''
        self.image = g
        self.objimage = self.image
        while self.objimage.find('/') != -1:
            self.objimage = self.objimage[self.objimage.find('/')+1:]
        verb = 0
        eclip = emclip
        self.eclip = eclip
        self.wavechange = 0
        self.zbound = zrange
        self.lowave = wave1
        self.hiwave = wave2
        self.create = 1
        self.range = 100
        self.add = 0
        self.star = 0
        self.baddata = 0
        self.unsure = 0
        self.interactive = interactive
        self.failure = 0
        self.flag1 = 0
        self.flag2 = 0
        self.flag3 = 0
        self.mkspec = specplot

    def doCC(self,g,eclip,verb,newmask=0,wavechange=0,good=0,changeclip=0):
        '''
        Manages the cross correlation, sending the templates through one at a time and checking for failures.
        '''
        print 'Computing Redshift for ' + self.image
        if newmask:
            G = g
        else:
            G = Spectrum(g,clip=eclip,verb=verb,good=good,changeclip=changeclip)
            self.G = G

        self.t23 = CCcalc(23,G,self.range,eclip,verb,self.zbound,wavechange=wavechange)
        self.t24 = CCcalc(24,G,self.range,eclip,verb,self.zbound,wavechange=wavechange)
        self.t25 = CCcalc(25,G,self.range,eclip,verb,self.zbound,wavechange=wavechange)
        self.t26 = CCcalc(26,G,self.range,eclip,verb,self.zbound,wavechange=wavechange)
        self.t27 = CCcalc(27,G,self.range,eclip,verb,self.zbound,wavechange=wavechange)

        if self.t23.failure or self.t24.failure or self.t25.failure or self.t26.failure or self.t27.failure:
            self.failure = 1
            return
            
        maxarr = np.asarray((self.t23.ccmax,self.t24.ccmax,self.t25.ccmax,self.t26.ccmax,self.t27.ccmax))
        chisqarr = np.asarray((self.t23.rchi2,self.t24.rchi2,self.t25.rchi2,self.t26.rchi2,self.t27.rchi2))
        maxindex = np.argmax(maxarr)
        minindex = np.argmin(chisqarr)
        #If best chisq and cc_max at same point, it's easy
        if maxindex == minindex:
            self.tempid = maxindex + 23
        else: #It just couldn't be easy, could it?
            bestzs = np.asarray((self.t23.zsmax,self.t24.zsmax,self.t25.zsmax,self.t26.zsmax,self.t27.zsmax))
            z_o_a = 0.85 * (self.zbound[1] - self.zbound[0]) + self.zbound[0]
            zneighbors = [len(np.where(bestzs - bestzs[ii] <= 0.015)[0]) for ii in range(5)]
            #First: is one best guess at the upper-z bounds? It's probably worse
            if bestzs[maxindex] >= z_o_a and bestzs[minindex] < z_o_a:
                self.tempid = minindex + 23
            elif bestzs[maxindex] < z_o_a and bestzs[minindex] >= z_o_a:
                self.tempid = maxindex + 23
            #Second: is one a 3+ hit mode? Go with that.
            elif zneighbors[maxindex] >= 3:
                self.tempid = maxindex + 23
            elif zneighbors[minindex] >= 3:
                self.tempid = minindex + 23
            #Third: Is one a 2-hit mode and the other isn't? Try that
            elif zneighbors[maxindex] >= 2 and zneighbors[minindex] == 1:
                self.tempid = maxindex + 23
            elif zneighbors[minindex] >= 2 and zneighbors[maxindex] == 1:
                self.tempid = minindex + 23
            #Fourth: Okay, this is probably not too hot. Just go with the rchisq
            else:
                self.tempid = minindex + 23
            
        self.findTemp()

    def findTemp(self):
        if self.tempid == 23: self.temp = self.t23
        elif self.tempid == 24: self.temp = self.t24
        elif self.tempid == 25: self.temp = self.t25
        elif self.tempid == 26: self.temp = self.t26
        elif self.tempid == 27: self.temp = self.t27

    def mkspecplot(self):
        self.lfile = np.genfromtxt(imloc + 'galaxylines.dat',skip_header=11,dtype=None,names=('loc','type','name'))
        self.fig2 = plt.figure(figsize=(10,6))
        self.a2x = self.fig2.add_axes([0, 0, 1, 1])
        G = self.G
        self.galplot = self.a2x.plot(G.wavelength,G.spec,color='black')
        xmin,xmax = self.a2x.get_xlim()
        self.a2x.set_xlim(xmin,xmax)
        self.a2x.set_ylim(np.median(G.spec)*-0.5,np.median(G.spec)*3.5)
        ymin,ymax = self.a2x.get_ylim()
        for i in range(len(self.lfile)):
            loc,ctype,name = self.lfile[i]
            shiftloc = loc*(1+self.temp.zsmax)
            if ctype == 2: color = 'indianred'
            else: color = 'cornflowerblue'
            if name == 'Break': color = 'mediumseagreen'
            elif name == 'Hbeta': color = 'cornflowerblue'
            self.a2x.plot((shiftloc,shiftloc),(-100000,100000),color=color,linestyle='--',alpha=0.75)
        
    
    def createPlot(self):
        '''
        Sets up the GUI for the first time 
        '''
        if self.failure:
            G = self.G
            f3 = open(outfile,"r+")
            f4 = open(outfilev,'r+')
            d = f3.readlines()
            e = f4.readlines()
            f3.seek(0)
            f4.seek(0)
            for i in d:
                if self.objimage in i:
                    f3.write(('%06d   %d    %.6f    %.6f   ' % (G.objid,0,0,0)+
                              '%.6f    %7.2f   %7.2f   %2.2f     ' % (0,0,0,0)+
                              '%d      %2.2f      %d    %d     ' % (0,0,0,1)+
                              '%d ' % (1)) + self.objimage + 
                             ('      %d      %d      %d' % (0,0,1)) + '\n')

                else:
                    f3.write(i)
            for i in e:
                if self.objimage in i:
                    f4.write(('%06d %.6f %.6f %.6f '%(G.objid,0,0,0)+
                              '%7.3f %7.3f %2.2f %.6f '%(0,0,0,0)+
                              '%.6f %.6f %7.3f %7.3f '%(0,0,0,0)+
                              '%2.2f %.6f %.6f %.6f '%(0,0,0,0)+
                              '%7.3f %7.3f %2.2f %.6f '%(0,0,0,0)+
                              '%.6f %.6f %7.3f %7.3f '%(0,0,0,0)+
                              '%2.2f %.6f %.6f %.6f '%(0,0,0,0)+
                              '%7.3f %7.3f %2.2f %d '%(0,0,0,0)+
                              '%d %3.2f %d %d '%(0,0,0,1)+
                              '%d' % (1)) + ' '+ self.objimage + 
                             (' %d %d %d' % (0,0,1)) + '\n')
                else:
                    f4.write(i)
            f3.close()
            f4.close()
            return
                
        if self.create:
            G = self.G
            self.fig = plt.figure(figsize=(14,8))
            self.ax1 = self.fig.add_axes([0.05, 0.7, 0.43, 0.25])
            self.ax2 = self.fig.add_axes([0.53, 0.7, 0.43, 0.20])
            self.ax2b = self.fig.add_axes([0.65, 0.40, 0.3, 0.2])
            self.ax2c = self.fig.add_axes([0.65, 0.1, 0.3, 0.2])
            self.ax3 = self.fig.add_axes([0.19, 0.1, 0.45, 0.50])
            self.ax5 = self.fig.add_axes([0.53, 0.9, 0.43, 0.08])
            self.ax1.set_title('Spectrum of ' + self.objimage + ', S/N = ' + str(G.signoise))
            self.ax2.set_title('Galaxy Template Comparison')
            self.ax2b.set_title('H and K Line Zoom')
            self.ax2c.set_title('O3 and H-Beta Line Zoom')
            self.ax3.set_title('Cross Correlation                        ')
            self.refreshText()
            self.l1, = self.ax1.plot(G.wavelength,G.spec,color='cornflowerblue',label='Spectrum')
            self.l2, = self.ax1.plot(G.wavelength,G.good*G.spec,color='orange',label='Mask')
            self.l3, = self.ax1.plot(G.wavelength,G.sm,color='indianred',label='Smooth',visible=False)
            self.l4, = self.ax1.plot(G.wavelength,self.temp.Gp,color='darkgreen',label='Fit')
            self.l5, = self.ax1.plot((self.lowave,self.lowave),(-100000,100000),color='black')
            self.l6, = self.ax1.plot((self.hiwave,self.hiwave),(-100000,100000),color='black')
            self.ax1.set_ylim(np.median(G.spec)*-0.5,np.median(G.spec)*4)
            self.ax1.set_xlabel('$\lambda$')
            self.ax1.set_ylabel('Counts')
            self.ax1.legend(loc=1)
            self.k1, = self.ax2.plot(G.wavelength,self.temp.Gr,color='cornflowerblue')
            self.k2, = self.ax2.plot(G.wavelength,self.temp.Trs[self.temp.zp][self.temp.zmax],color='black',alpha=0.5)
            self.k1b, = self.ax2b.plot(G.wavelength,self.temp.Gr2,color='cornflowerblue')
            self.k2b, = self.ax2b.plot(G.wavelength,self.temp.Trs[self.temp.zp][self.temp.zmax],color='black',alpha=0.5)
            self.k1c, = self.ax2c.plot(G.wavelength,self.temp.Gr2,color='cornflowerblue')
            self.k2c, = self.ax2c.plot(G.wavelength,self.temp.Trs[self.temp.zp][self.temp.zmax],color='black',alpha=0.5)
            l1bloc = 3933.7*(1+self.temp.zsmax)
            l2bloc = 3968.5*(1+self.temp.zsmax)
            self.l1b, = self.ax2b.plot((l1bloc,l1bloc),(-10,10),color='mediumseagreen',linestyle='--')
            self.l2b, = self.ax2b.plot((l2bloc,l2bloc),(-10,10),color='mediumseagreen',linestyle='--')
            l1cloc = 4861.3*(1+self.temp.zsmax)
            l2cloc = 4959*(1+self.temp.zsmax)
            l3cloc = 5007*(1+self.temp.zsmax)
            self.l1c, = self.ax2c.plot((l1cloc,l1cloc),(-10,10),color='mediumseagreen',linestyle='--')
            self.l2c, = self.ax2c.plot((l2cloc,l2cloc),(-10,10),color='mediumseagreen',linestyle='--')
            self.l3c, = self.ax2c.plot((l3cloc,l3cloc),(-10,10),color='mediumseagreen',linestyle='--')
            self.ax2.set_xlim(self.temp.w1-100,self.temp.w2+100)
            self.ax2.set_ylim(-4,5.5)
            self.ax2b.set_xlim(3800*(1+self.temp.zsmax),4100*(1+self.temp.zsmax))
            self.ax2b.set_ylim(-4,5.5)
            self.ax2c.set_xlim(4750*(1+self.temp.zsmax),5050*(1+self.temp.zsmax))
            self.ax2c.set_ylim(-4,5.5)
            self.ax2b.yaxis.tick_right()
            self.ax2c.yaxis.tick_right()
            self.ax2.set_xlabel('$\lambda$                                 ')
            self.ax2.set_ylabel('(g-f)/$\sigma$')
            self.ax2b.set_xlabel('$\lambda$')
            
            self.ax2c.set_xlabel('$\lambda$')
            self.m1, = self.ax3.plot(self.temp.zs[0],self.temp.cc[0],color='indianred',zorder=1)
            self.m1_23, = self.ax3.plot(self.t23.zs[0],self.t23.cc[0]*(max(self.temp.cc[0])/max(self.t23.cc[0])),
                                        color='#ff9933',alpha=0.35,zorder=0)
            self.m1_24, = self.ax3.plot(self.t24.zs[0],self.t24.cc[0]*(max(self.temp.cc[0])/max(self.t24.cc[0])),
                                        color='Yellow',alpha=0.35,zorder=0)
            self.m1_25, = self.ax3.plot(self.t25.zs[0],self.t25.cc[0]*(max(self.temp.cc[0])/max(self.t25.cc[0])),
                                        color='#33cc33',alpha=0.35,zorder=0)
            self.m1_26, = self.ax3.plot(self.t26.zs[0],self.t26.cc[0]*(max(self.temp.cc[0])/max(self.t26.cc[0])),
                                        color='#0066cc',alpha=0.35,zorder=0)
            self.m1_27, = self.ax3.plot(self.t27.zs[0],self.t27.cc[0]*(max(self.temp.cc[0])/max(self.t27.cc[0])),
                                        color='#9933ff',alpha=0.35,zorder=0)
            if self.temp.zp>0: self.m2, = self.ax3.plot(self.temp.zs[1],self.temp.cc[1],color='firebrick')
            else: self.m2, = self.ax3.plot((0,0),(0,0),color='firebrick',zorder=2)
            if self.temp.zp>1: self.m3, = self.ax3.plot(self.temp.zs[2],self.temp.cc[2],color='darkred')
            else: self.m3, = self.ax3.plot((0,0),(0,0),color='darkred',zorder=2)
            if self.temp.zp>2: self.m4, = self.ax3.plot(self.temp.zs[3],self.temp.cc[3],color='darkred')
            else: self.m4, = self.ax3.plot((0,0),(0,0),color='darkred',zorder=2)
            self.g1, = self.ax3.plot((self.temp.zs[self.temp.zp][self.temp.zmax],
                                      self.temp.zs[self.temp.zp][self.temp.zmax]),
                                     (-10000,10000),color='mediumseagreen',zorder=1)
            self.ax3.set_xlim(self.zbound[0],self.zbound[1])
            self.ax3.set_ylim(min(self.temp.cc[0])*0.95,max(self.temp.cc[0]*1.05))
            self.ax3.set_xlabel('z')
            self.ax3.set_ylabel('Correlation Coeffecient')
            self.ax5.imshow(G.image_data, cmap='gray',clim=(-240.0, 240.0),aspect='auto')
            self.ax5.plot((0,12000),(G.lobar,G.lobar),color='mediumseagreen')
            self.ax5.plot((0,12000),(G.hibar,G.hibar),color='mediumseagreen')
            self.ax5.get_xaxis().set_visible(False)
            self.ax5.get_yaxis().set_visible(False)
            ax2min,ax2max = self.ax2.get_xlim()
            self.ax5.set_xlim((ax2min-4900)/2,(ax2max-4900)/2)
            rax = plt.axes([0.01, 0.34, 0.08, 0.12])
            check = CheckButtons(rax, ('Spectrum', 'Mask','Smooth','Fit'), (True, True, False, True))
            self.create = 0
                
#==================================================================
#Buttons and sliders
#==================================================================
        def checkButton(label):
            if label == 'Spectrum':
                self.l1.set_visible(not self.l1.get_visible())
            elif label == 'Mask':
                self.l2.set_visible(not self.l2.get_visible())
            elif label == 'Smooth':
                self.l3.set_visible(not self.l3.get_visible())
            elif label == 'Fit':
                self.l4.set_visible(not self.l4.get_visible())
            self.ax1.legend(loc=1)
            plt.draw()
        check.on_clicked(checkButton)

        def ccslider(xmin, xmax):
            indmin, indmax = np.searchsorted(self.temp.zs[0], (xmin, xmax))
            indmax = min(len(self.temp.zs[0]) - 1, indmax)
            self.newcc = np.mean((indmin,indmax))
            self.temp.rerunCC(self.newcc,self.range,self.zbound)
            self.updatePlot()
        span3 = SpanSelector(self.ax3, ccslider, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))

        def waveChangeHiF(event):
            self.hiwave += 50
            updateWave()
        wavechangehif = plt.axes([0.298, 0.643, 0.015, 0.02])
        bwavechangehif = Button(wavechangehif, '>')
        bwavechangehif.on_clicked(waveChangeHiF)
        
        def waveChangeHiB(event):
            self.hiwave -= 50
            updateWave()
        wavechangehib = plt.axes([0.282, 0.643, 0.015, 0.02])
        bwavechangehib = Button(wavechangehib, '<')
        bwavechangehib.on_clicked(waveChangeHiB)

        def waveChangeLoF(event):
            self.lowave += 50
            updateWave()
        wavechangelof = plt.axes([0.232, 0.643, 0.015, 0.02])
        bwavechangelof = Button(wavechangelof, '>')
        bwavechangelof.on_clicked(waveChangeLoF)
        
        def waveChangeLoB(event):
            self.lowave -= 50
            updateWave()
        wavechangelob = plt.axes([0.216, 0.643, 0.015, 0.02])
        bwavechangelob = Button(wavechangelob, '<')
        bwavechangelob.on_clicked(waveChangeLoB)

        def waveChangeReset(event):
            self.lowave = wave1
            self.hiwave = wave2
            updateWave()
        wavechangereset = plt.axes([0.176, 0.643, 0.036, 0.02])
        bwavechangereset = Button(wavechangereset, 'Reset')
        bwavechangereset.on_clicked(waveChangeReset)

        def updateWave():
            self.wavechange = (self.lowave,self.hiwave)
            self.l5.set_xdata((self.lowave,self.lowave))
            self.l6.set_xdata((self.hiwave,self.hiwave))
            plt.draw()

        def adjRange(event):
            self.range -= 40
            if self.range == -20:
                self.range = 100
            badjrange.label.set_text('Range: ' + ('%3d' % self.range))
            plt.draw()
        adjrange = plt.axes([0.01, 0.11, 0.12, 0.035])
        badjrange = Button(adjrange, 'Range: ' + ('%3d' % self.range))
        badjrange.on_clicked(adjRange)

        def switchTempF(event):
            self.tempid += 1
            if self.tempid == 28: self.tempid = 23
            bswitchtemp.label.set_text('Template = ' + ('%2d' % self.tempid))
            self.findTemp()
            self.updatePlot()
        def switchTempB(event):
            self.tempid -= 1
            if self.tempid == 22: self.tempid = 27
            bswitchtemp.label.set_text('Template = ' + ('%2d' % self.tempid))
            self.findTemp()
            self.updatePlot()
        switchtemp = plt.axes([0.025, 0.15, 0.09, 0.055])
        bswitchtemp = Button(switchtemp, 'Template = ' + ('%2d' % self.tempid))
        switchtempF = plt.axes([0.118, 0.15, 0.012, 0.055])
        bswitchtempF = Button(switchtempF, '>')
        bswitchtempF.on_clicked(switchTempF)
        switchtempB = plt.axes([0.010, 0.15, 0.012, 0.055])
        bswitchtempB = Button(switchtempB, '<')
        bswitchtempB.on_clicked(switchTempB)

        def toggleClip(event):
            self.eclip = not self.eclip
            btoggleclip.label.set_text('Clip = ' + str(self.eclip))
            self.worktext = self.fig.text(0.5,0.01,'Working...',fontsize=24)
            plt.pause(.1)
            self.doCC(self.image,self.eclip,0,wavechange=self.wavechange,good=self.G.good,changeclip=1)
            #self.doCC(self.image,self.eclip,0)
            plt.gcf().texts.remove(self.worktext)
            bswitchtemp.label.set_text('Template = ' + ('%2d' % self.tempid))
            G = self.G
            self.l2.set_ydata(G.good*G.spec)
            self.l3.set_ydata(G.sm)
            self.l4.set_ydata(self.temp.Gp)
            self.k1.set_ydata(self.temp.Gr)
            self.k1b.set_ydata(self.temp.Gr2)
            self.k1c.set_ydata(self.temp.Gr2)
            self.updatePlot()
        toggleclip = plt.axes([0.010, 0.21, 0.12, 0.035])
        btoggleclip = Button(toggleclip, 'Clip = ' + str(self.eclip))
        btoggleclip.on_clicked(toggleClip)

        def maskMode(event):
            self.add = not self.add
            if self.add: bmode.label.set_text('Mode: Add')
            else: bmode.label.set_text('Mode: Sub')
            plt.draw()
        mode = plt.axes([0.010, 0.31, 0.06, 0.035])
        bmode = Button(mode, 'Mode: Sub')
        bmode.on_clicked(maskMode)

        def maskslider(xmin, xmax):
            indmin, indmax = np.searchsorted(G.wavelength, (xmin, xmax))
            indmax = min(len(G.wavelength) - 1, indmax)
            if indmin == indmax:
                hidiff = abs(self.hiwave - G.wavelength[indmin])
                lodiff = abs(self.lowave - G.wavelength[indmin])
                if hidiff > lodiff:
                    self.lowave = G.wavelength[indmin]
                else:
                    self.hiwave = G.wavelength[indmin]
                updateWave()
            else:
                if self.add == 1:
                    self.G.good = np.concatenate((self.G.good[:indmin],np.ones(indmax-indmin),self.G.good[indmax:]))
                    self.G.continuum = np.concatenate((self.G.continuum[:indmin],np.ones(indmax-indmin),self.G.continuum[indmax:]))
                if self.add == 0:
                    self.G.good = np.concatenate((self.G.good[:indmin],np.zeros(indmax-indmin),self.G.good[indmax:]))
                    self.G.continuum = np.concatenate((self.G.continuum[:indmin],np.zeros(indmax-indmin),self.G.continuum[indmax:]))
                self.l2.set_ydata(self.G.good*self.G.spec)
                plt.draw()
        span = SpanSelector(self.ax1, maskslider, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))

        def buttonGo(event):
            self.worktext = self.fig.text(0.5,0.01,'Working...',fontsize=24)
            plt.pause(.1)
            self.doCC(self.G,self.eclip,0,newmask=1,wavechange=self.wavechange)
            plt.gcf().texts.remove(self.worktext)
            bswitchtemp.label.set_text('Template = ' + ('%2d' % self.tempid))
            self.k1.set_ydata(self.temp.Gr)
            self.k1b.set_ydata(self.temp.Gr2)
            self.k1c.set_ydata(self.temp.Gr2)
            self.l4.set_ydata(self.temp.Gp)
            self.updatePlot()
        go = plt.axes([0.070, 0.31, 0.02, 0.035])
        bgo = Button(go, 'Go')
        bgo.on_clicked(buttonGo)

        def Bsave(event):
            self.savetext = self.fig.text(0.5,0.01,'Saving...',fontsize=24)
            plt.pause(.1)
            plt.gcf().texts.remove(self.savetext)
            G = self.G
            f3 = open(outfile,"r+")
            f4 = open(outfilev,'r+')
            d = f3.readlines()
            e = f4.readlines()
            f3.seek(0)
            f4.seek(0)
            for i in d:
                if self.objimage in i:
                    f3.write(('%06d   %d    %.6f    %.6f   ' % (G.objid,self.tempid,self.temp.zsmax,self.temp.dzhi)+
                              '%.6f    %7.2f   %7.2f   %2.2f     ' % (self.temp.dzlo, self.temp.ccmax, self.temp.chi2,self.temp.rchi2)+
                              '%d      %2.2f      %d    %d     ' % (self.eclip,G.signoise,self.star,self.baddata)+
                              '%d ' % (self.unsure)) + self.objimage + 
                             ('      %d      %d      %d' % (self.flag1, self.flag2, self.flag3)) + '\n')
                else:
                    f3.write(i)
            for i in e:
                if self.objimage in i:
                    f4.write(('%06d %.6f %.6f %.6f ' %(G.objid,self.t23.zsmax,self.t23.dzhi,self.t23.dzlo) +
                              '%7.3f %7.3f %2.2f %.6f ' %( self.t23.ccmax, self.t23.chi2,self.t23.rchi2,self.t24.zsmax) +
                              '%.6f %.6f %7.3f %7.3f ' %(self.t24.dzhi,self.t24.dzlo, self.t24.ccmax, self.t24.chi2) +
                              '%2.2f %.6f %.6f %.6f ' %(self.t24.rchi2,self.t25.zsmax,self.t25.dzhi,self.t25.dzlo) +
                              '%7.3f %7.3f %2.2f %.6f ' %( self.t25.ccmax, self.t25.chi2,self.t25.rchi2,self.t26.zsmax) +
                              '%.6f %.6f %7.3f %7.3f ' %(self.t26.dzhi,self.t26.dzlo, self.t26.ccmax, self.t26.chi2) +
                              '%2.2f %.6f %.6f %.6f ' %(self.t26.rchi2,self.t27.zsmax,self.t27.dzhi,self.t27.dzlo) +
                              '%7.3f %7.3f %2.2f %d ' %( self.t27.ccmax, self.t27.chi2,self.t27.rchi2,self.tempid) +
                              '%d %3.2f %d %d ' %(self.eclip,G.signoise,self.star,self.baddata) +
                              '%d' % (self.unsure)) + ' '+ self.objimage + 
                             (' %d %d %d' % (self.flag1, self.flag2, self.flag3)) + '\n')

                else:
                    f4.write(i)
            f3.close()
            f4.close()
            filelocation2 = outloc + 'cc_' + self.objimage.replace('.fits','.png')
            filelocation3 = outloc + 'cc_' + self.objimage.replace('.fits','2.png')
            self.fig.savefig(filelocation2)
            if self.mkspec: self.fig2.savefig(filelocation3)
            plt.clf()
            plt.close(self.fig)
            if self.mkspec: plt.close(self.fig2)

        def checkStarFlag():
            if self.star == 1: 
                bsflag.label.set_text('Unmark Star')
                bsflag.color = 'IndianRed'
                bsflag.hovercolor='Red'
            else: 
                bsflag.label.set_text('Mark Star')
                bsflag.color = '0.85'
                bsflag.hovercolor='0.95'
            plt.draw()

        def starFlag(event):
            if self.star == 1: self.star = 0
            else: self.star = 1
            checkStarFlag()
        sflag = plt.axes([0.15, 0.0, 0.08, 0.05])
        if G.spec[1000:1500].mean() > 10000:
            self.star = 1
            bsflag = Button(sflag, 'Mark Star',color = 'IndianRed',hovercolor = 'Red')
        else:
            self.star = 0
            bsflag = Button(sflag, 'Mark Star',color = '0.85',hovercolor = '0.95')
        bsflag.on_clicked(starFlag)
        checkStarFlag()

        def dataFlag(event):
            if self.baddata == 1:
                self.baddata = 0
                bdflag.label.set_text('Mark Bad')
                bdflag.color = '0.85'
                bdflag.hovercolor='0.95'
            else:
                self.baddata = 1
                bdflag.label.set_text('Unmark Bad')
                bdflag.color = 'IndianRed'
                bdflag.hovercolor='Red'
            plt.draw()
        dflag = plt.axes([0.25, 0.0, 0.08, 0.05])
        if self.baddata == 0: bdflag = Button(dflag, 'Mark Bad',color='0.85',hovercolor='0.95')
        else: bdflag = Button(dflag, 'Unmark Bad',color='IndianRed',hovercolor='Red')
        bdflag.on_clicked(dataFlag)
        

        def unsureFlag(event):
            if self.unsure == 1:
                self.unsure = 0
                buflag.label.set_text('Mark Unsure')
            else:
                self.unsure = 1
                buflag.label.set_text('Unmark Unsure')
            plt.draw()
        uflag = plt.axes([0.35, 0.0, 0.08, 0.05])
        buflag = Button(uflag, 'Mark Unsure')
        buflag.on_clicked(unsureFlag)

        def eFlag1(event):
            if self.flag1 == 1:
                self.flag1 = 0
                bflag1.label.set_text('Mark Revisit')
            else:
                self.flag1 = 1
                bflag1.label.set_text('Unmark Revisit')
            plt.draw()
        flag1 = plt.axes([0.65, 0.0, 0.08, 0.05])
        bflag1 = Button(flag1, 'Mark Revisit')
        bflag1.on_clicked(eFlag1)
        
        def eFlag2(event):
            if self.flag2 == 1:
                self.flag2 = 0
                bflag2.label.set_text('Mark Note')
            else:
                self.flag2 = 1
                bflag2.label.set_text('Unmark Range')
            plt.draw()
        flag2 = plt.axes([0.75, 0.0, 0.08, 0.05])
        bflag2 = Button(flag2, 'Mark Range')
        bflag2.on_clicked(eFlag2)
        
        def eFlag3(event):
            if self.flag3 == 1:
                self.flag3 = 0
                bflag3.label.set_text('Mark Unusable')
            else:
                self.flag3 = 1
                bflag3.label.set_text('Unmark Unusable')
            plt.draw()
        flag3 = plt.axes([0.85, 0.0, 0.08, 0.05])
        bflag3 = Button(flag3, 'Mark Unusable')
        bflag3.on_clicked(eFlag3)

        def changeZ(event):
            self.inputtext = self.fig.text(0.45,0.01,'Input Lower Bound (float)',fontsize=14)
            plt.pause(.1)
            loz = input('Input lower bound (float): ')
            plt.gcf().texts.remove(self.inputtext)
            self.inputtext = self.fig.text(0.45,0.01,'Input Upper Bound (float)',fontsize=14)
            plt.pause(.1)
            hiz = input('Input upper bound (float): ')
            plt.gcf().texts.remove(self.inputtext)
            self.zbound = loz,hiz
            self.worktext = self.fig.text(0.5,0.01,'Working...',fontsize=24)
            plt.pause(.1)
            self.doCC(self.G,self.eclip,0,newmask=1,wavechange=self.wavechange)
            plt.gcf().texts.remove(self.worktext)
            bswitchtemp.label.set_text('Template = ' + ('%2d' % self.tempid))
            self.k1.set_ydata(self.temp.Gr)
            self.k1b.set_ydata(self.temp.Gr2)
            self.k1c.set_ydata(self.temp.Gr2)
            self.l4.set_ydata(self.temp.Gp)
            self.updatePlot()
        changez = plt.axes([0.0, 0.06, 0.1, 0.03])
        bchangez = Button(changez, 'Change z Range')
        bchangez.on_clicked(changeZ)

        def checkSpecFlag():
            if self.star == 1: 
                bsflag.label.set_text('Unmark Star')
                bsflag.color = 'Red'
            else: 
                bsflag.label.set_text('Mark Star')
                bsflag.color = 0.85
            plt.draw()
        
        def mkSpec(event):
            if self.mkspec == 1:
                self.mkspec = 0
                plt.close(self.fig2)
                bmkspec.label.set_text('Show Spec Plot')
            else:
                self.mkspec = 1
                self.mkspecplot()
                plt.show(block=False)
                plt.figure(1)
                bmkspec.label.set_text('Hide Spec Plot')
            plt.draw()
        mkspec = plt.axes([0.01, 0.63, 0.08, 0.03])
        if specplot: bmkspec = Button(mkspec, 'Hide Spec Plot')
        else: bmkspec = Button(mkspec, 'Show Spec Plot')
        bmkspec.on_clicked(mkSpec)

        save = plt.axes([0.0, 0.0, 0.1, 0.05])
        if self.interactive:
            bsave = Button(save, 'Save')
            bsave.on_clicked(Bsave)
            if self.mkspec: self.mkspecplot()
            plt.show()
        else:
            Bsave(self)
        

#==============================================================
            
    def updatePlot(self):
        '''
        Updates the plot as changes are made (e.g. toggling the emission clip or changing the wavelength range)
        '''
        self.k2.set_ydata(self.temp.Trs[self.temp.zp][self.temp.zmax])
        self.k2b.set_ydata(self.temp.Trs[self.temp.zp][self.temp.zmax])
        self.k2c.set_ydata(self.temp.Trs[self.temp.zp][self.temp.zmax])
        self.m1.set_ydata(self.temp.cc[0])
        self.m1.set_xdata(self.temp.zs[0])
        self.m1_23.set_ydata(self.t23.cc[0]*(max(self.temp.cc[0])/max(self.t23.cc[0])))
        self.m1_23.set_xdata(self.t23.zs[0])
        self.m1_24.set_ydata(self.t24.cc[0]*(max(self.temp.cc[0])/max(self.t24.cc[0])))
        self.m1_24.set_xdata(self.t24.zs[0])
        self.m1_25.set_ydata(self.t25.cc[0]*(max(self.temp.cc[0])/max(self.t25.cc[0])))
        self.m1_25.set_xdata(self.t25.zs[0])
        self.m1_26.set_ydata(self.t26.cc[0]*(max(self.temp.cc[0])/max(self.t26.cc[0])))
        self.m1_26.set_xdata(self.t26.zs[0])
        self.m1_27.set_ydata(self.t27.cc[0]*(max(self.temp.cc[0])/max(self.t27.cc[0])))
        self.m1_27.set_xdata(self.t27.zs[0])
        if self.temp.zp>0:
            self.m2.set_ydata(self.temp.cc[1])
            self.m2.set_xdata(self.temp.zs[1])
        else:
            self.m2.set_ydata((0,0))
            self.m2.set_xdata((0,0))
        if self.temp.zp>1:
            self.m3.set_ydata(self.temp.cc[2])
            self.m3.set_xdata(self.temp.zs[2])
        else:
            self.m3.set_ydata((0,0))
            self.m3.set_xdata((0,0))
        if self.temp.zp>2:
            self.m4.set_ydata(self.temp.cc[3])
            self.m4.set_xdata(self.temp.zs[3])
        else:
            self.m4.set_ydata((0,0))
            self.m4.set_xdata((0,0))
        self.g1.set_xdata((self.temp.zs[self.temp.zp][self.temp.zmax],self.temp.zs[self.temp.zp][self.temp.zmax]))
        self.ax3.set_xlim(self.zbound[0],self.zbound[1])
        self.ax3.set_ylim(min(self.temp.cc[0])*0.95,max(self.temp.cc[0]*1.05))
        self.ax2.set_xlim(self.temp.w1-100,self.temp.w2+100)
        l1bloc = 3933.7*(1+self.temp.zsmax)
        l2bloc = 3968.5*(1+self.temp.zsmax)
        self.l1b.set_xdata((l1bloc,l1bloc))
        self.l2b.set_xdata((l2bloc,l2bloc))
        l1cloc = 4861.3*(1+self.temp.zsmax)
        l2cloc = 4959*(1+self.temp.zsmax)
        l3cloc = 5007*(1+self.temp.zsmax)
        self.l1c.set_xdata((l1cloc,l1cloc))
        self.l2c.set_xdata((l2cloc,l2cloc))
        self.l3c.set_xdata((l3cloc,l3cloc))
        self.ax2b.set_xlim(3800*(1+self.temp.zsmax),4100*(1+self.temp.zsmax))
        self.ax2c.set_xlim(4750*(1+self.temp.zsmax),5050*(1+self.temp.zsmax))
        ax2min,ax2max = self.ax2.get_xlim()
        self.ax5.set_xlim((ax2min-4900)/2,(ax2max-4900)/2)
        self.refreshText(remove = 1)
        if self.mkspec:
            plt.close(self.fig2)
            self.mkspecplot()
            plt.show(block=False)
            plt.figure(1)
        plt.draw()

    def refreshText(self,remove = 0):
        if remove:
            plt.gcf().texts.remove(self.ztext)
            plt.gcf().texts.remove(self.zhitext)
            plt.gcf().texts.remove(self.zlotext)
            plt.gcf().texts.remove(self.temptext)
            plt.gcf().texts.remove(self.linetext)
            plt.gcf().texts.remove(self.t23text)
            plt.gcf().texts.remove(self.t24text)
            plt.gcf().texts.remove(self.t25text)
            plt.gcf().texts.remove(self.t26text)
            plt.gcf().texts.remove(self.t27text)
        self.ztext = self.fig.text(0.4285,0.610,"z = %.6f" % (self.temp.zsmax),fontsize = 12)
        self.zhitext = self.fig.text(0.5095,0.624,"+%.6f" % (self.temp.dzhi),fontsize = 8)
        self.zlotext = self.fig.text(0.513,0.604,"-%.6f" % (self.temp.dzlo),fontsize = 8)
        self.temptext = self.fig.text(0.005,0.58,'$T$     $ccmax$         $z$                 $\chi^2$',fontsize = 9)
        self.linetext = self.fig.text(0.005,0.58,'__   _______    _________     ______',fontsize=9)
        self.t23text = self.fig.text(0.005,0.56,"23: %7.2f   %.6f   %7.2f" % (self.t23.ccmax,self.t23.zsmax,self.t23.chi2),fontsize = 9)
        self.t24text = self.fig.text(0.005,0.545,"24: %7.2f   %.6f   %7.2f" % (self.t24.ccmax,self.t24.zsmax,self.t24.chi2),fontsize = 9)
        self.t25text = self.fig.text(0.005,0.53,"25: %7.2f   %.6f   %7.2f" % (self.t25.ccmax,self.t25.zsmax,self.t25.chi2),fontsize = 9)
        self.t26text = self.fig.text(0.005,0.515,"26: %7.2f   %.6f   %7.2f" % (self.t26.ccmax,self.t26.zsmax,self.t26.chi2),fontsize = 9)
        self.t27text = self.fig.text(0.005,0.50,"27: %7.2f   %.6f   %7.2f" % (self.t27.ccmax,self.t27.zsmax,self.t27.chi2),fontsize = 9)
       
        
class CCcalc:
    '''
    Performs the cross correlation
    G - galaxy spectrum after going the class Spectrum()
    range - int - number of points per iteration. 
            Default 100, user can set to 60 or 20 in the GUI if they want to zoom in on a small peak that is beside a large one
    t - int - template number, either 23,24,25,26,27
    failure - boolean - set to 1 if the correlation fails
    z1, z2 - float,float - bounds over which to fit redshift
    dz0 - float - starting stepsize for the first iteration
    nz - int - number of steps
    zp - int - iteration number (first iteration is 0)
    keepgoing - boolean - true if another iteration should be done, 
                automatically set to false if enough iterations have been done to accurately compute errors
    w1,w2 - int,int - wavelength range over which to correlate. Data outside of this range is ignored.
    wavechange - int,int - if the user has modified the wavelength in the GUI, wavechange contains the new lower and upper bound
    o - int - order of the polynomial to fit to the galaxy and templates. Scales with the length of data over which we are fitting
    cc - array of correlation coefficients
    zs - array of redshifts corresponding to those coefficients
    Trs - array of normalized templates shifted to the redshifts in zs
    '''
    def __init__(self,t,G,inrange,eclip,verb,zbound,wavechange=0):
        self.G = G
        self.range = inrange
        self.t = t
        self.failure = 0
    
        T = Template(t,clip=eclip,verb=verb)
        self.T = T
        
        z1,z2 = zbound
        dz0 = 4*G.cdelt1/6500.0
        nz = int((z2-z1)/dz0+1)
        zp = 0
        keepgoing = True
        w1,w2 = wave1,wave2
        if wavechange:
            self.w1 = wavechange[0]
            self.w2 = wavechange[1]
        else:
            self.w1 = w1
            self.w2 = w2
        o = int((w2-w1)/250)
        self.o = o

        self.cc = []
        self.zs = []
        self.Trs = []
        
        dz = (z2-z1)/(nz-1)

        self.computeCC(z1,z2,dz,zp,keepgoing,self.w1,self.w2,o,G,T,self.t,verb,intro=1)

    def rerunCC(self,newcc,inrange,zbound):
        '''
        This is done if the user clicks on a different peak - 
        instead of rerunning the whole template, we can skip the first run through and start at the second, near where the user clicked.
        See above docstring for explanation of variables
        '''
        self.range = inrange
        self.cc = [self.cc[0]]
        self.zs = [self.zs[0]]
        self.Trs = [self.Trs[0]]
        zmax = self.zs[0][int(newcc)]
        z1,z2 = zbound
        nz = 100
        dz = (z2-z1)/(nz-1)
        if self.range == 100:
            z1,z2 = zmax-5*dz,zmax+4.9*dz
        elif self.range == 60:
            z1,z2 = zmax-3*dz,zmax+2.9*dz
        elif self.range == 20:
            z1,z2 = zmax-1*dz,zmax+0.9*dz
        dz = dz/10.0
        zp = 1
        keepgoing = True
        self.zp = zp
        self.computeCC(z1,z2,dz,zp,keepgoing,self.w1,self.w2,self.o,self.G,self.T,self.t,0)
            

    def computeCC(self,z1,z2,dz,zp,keepgoing,w1,w2,o,G,T,t,verb,intro=0):
        '''
        Main iterative calculations are done here
        '''

        dGw = G.wavelength[1:]-G.wavelength[:-1]
        dGw = np.compress(np.greater(G.wavelength[1:],w1)*np.greater(G.wavelength[:-1],w2),dGw)
        Gw = np.compress(np.greater(G.wavelength[1:],w1)*np.greater(G.wavelength[:-1],w2),(G.wavelength[1:]+G.wavelength[:-1])/2.)
        dzfloor = dGw.mean()/Gw.mean()

        while keepgoing:
            if intro:
                #print z1,z2,dzfloor,dz
                zs = np.arange(z1,z2+dz/2,dz,dtype='float64')
            else:
                #print zsmax,z1,z2,dzfloor,dz
                zs = np.sort(np.concatenate([np.arange(zsmax,z2+dz/2,dz,dtype='float64'), np.arange(zsmax-dz,z1-dz/2,-dz,dtype='float64')]))
            self.zs.append(zs)
                
            if verb: print "Using order=",o

            try:
                Gp = FitPoly(G.spec,G.noise,G.continuum,G.wavelength,o,w1,w2)
            except(ValueError):
                print 'Could not fit polynomial'
                keepgoing = 0
                self.failure = 1
            if not self.failure:
                self.Gp = Gp
                Gm = 1.*G.good*np.greater_equal(G.wavelength,w1)*np.less_equal(G.wavelength,w2)
                Gw = divz(Gm,G.noise)
                Gn = divz(np.sum(G.spec*Gw),np.sum(Gw))
                if verb:
                    lim = (np.median(G.spec)*0.25,np.median(G.spec)*1.75)
                    outGraph(G.wavelength,Gp,G.objid,'pl',lim=lim)

                Gp2 = FitPoly(G.spec,G.noise,G.continuum,G.wavelength,int((1e4-min([w1,4900]))/250.),min([w1,4900]),10000)
                Gm2 = 1.*G.good2
                Gw2 = divz(Gm2,G.noise)
                Gr2 = (G.spec-Gp2)*Gw2
                self.Gr2 = Gr2

                Tzs = np.asarray([T.Redshift(z,G.wavelength/(1+z)) for z in zs])
                Tzm = np.add.reduce(Tzs*Gw[np.newaxis,::],1)/np.sum(Gw)
                Tzs = (Tzs/Tzm[::,np.newaxis]) * Gn

                Tps = np.asarray([FitPoly(Tz,G.noise,G.continuum,G.wavelength,o,w1,w2) for Tz in Tzs])

                Gr = (G.spec-Gp)*Gw
                Trs = (Tzs-Tps)*Gw
                self.Gr = Gr
                self.Trs.append(Trs)
                
                cc = np.add.reduce(Trs*Gr[np.newaxis,::],1)
                self.cc.append(cc)
                if verb: outGraph(zs,cc,G.objid,"c_temp%d_%d" % (t,zp))
                
                zmax = np.argmax(cc)
                self.zmax = zmax
                #print "(%d,%d,%.5f)" % (t,zp,zs[zmax])
                try:
                    dz1 = cc[zmax]-cc[zmax-1]
                    dz2 = cc[zmax]-cc[zmax+1]
                except(IndexError):
                    dz1 = 1
                    dz2 = 1
                self.dz1 = dz1
                self.dz2 = dz2
                if dz1 < 0.5 and dz2 < 0.5: keepgoing = False
    
                if intro:
#                    dz = 0.3/99.0
                    intro = 0
                    
                if keepgoing:
                    cctck = splrep(zs,cc,s=0,k=3)
                    zsmax = fsolve(lambda z: splev(z,cctck,1),(zs[zmax],))[0]

                    dz = dz/10.0
                    dz = max([dz,dzfloor])
                    z1,z2 = zs[zmax]-10*dz,zs[zmax]+9.9*dz
#                    if self.range == 100:
#                        z1,z2 = zs[zmax]-5*dz,zs[zmax]+4.9*dz
#                    elif self.range == 60:
#                        z1,z2 = zs[zmax]-3*dz,zs[zmax]+2.9*dz
#                    elif self.range == 20:
#                        z1,z2 = zs[zmax]-1*dz,zs[zmax]+0.9*dz
#                    dz = dz/2.5
#                    dz = max([dz,dzfloor])
                    if zp < 10: zp += 1
                    else: keepgoing = False

#                    zp += 1
                        
        if not self.failure:
                '''
                Our results of zsmax and ccmax are calculated here, after keepgoing has been set to False
                '''
                cctck = splrep(zs,cc,s=0,k=3)
                zsmax = fsolve(lambda z: splev(z,cctck,1),(zs[zmax],))[0]
                ccmax = splev(zsmax,cctck)

                #An alternative way to compute the uncertainly, the the one below is preferred
                #zlo = fsolve(lambda z: splev(z,cctck)-(ccmax-0.5),(zs[zmax-3],))[0]
                #zhi = fsolve(lambda z: splev(z,cctck)-(ccmax-0.5),(zs[zmax+3],))[0]
                #dzlo = zs[zmax]-zlo
                #dzhi = zhi-zs[zmax]
                #print "         +%.6f" % (dzhi)
                #print "z = %.6f" % (zs[zmax])
                #print "         -%.6f" % (dzlo)
        
                cc2tck = splrep(zs,cc-(ccmax-0.5),s=0,k=3)
                try: 
                    ccroots = sproot(cc2tck,mest=2*len(cc))
                    if len(ccroots)<2:
                        shrnk_zs = zs[zmax - 10: zmax + 10]
                        shrnk_cc = cc[zmax - 10: zmax + 10] - ccmax
                        shrnk_cc2tck = splrep(shrnk_zs,shrnk_cc, s=0, k=3)
                        ccroots = sproot(shrnk_cc2tck,mest=2*len(shrnk_cc))
                except(TypeError):
                    ccroots = [0]
                    
                self.zp = zp
                #Here we estimate our uncertainly in z
                if len(ccroots)>=2:
                    zlo,zhi = np.sort(np.take(ccroots,np.argsort(abs(ccroots-zsmax)))[:2])
                    dzlo,dzhi = zsmax-zlo,np.abs(zhi-zsmax)
                    print "         +%.6f" % (dzhi)
                    print "z = %.6f" % (zsmax)
                    print "         -%.6f" % (dzlo)
                    self.zsmax = zsmax
                    self.ccmax = ccmax
                    self.dzhi = dzhi
                    self.dzlo = dzlo
                else:
                    zlo = fsolve(lambda z: splev(z,cctck)-(ccmax-0.5),(zs[zmax-3],))[0]
                    zhi = fsolve(lambda z: splev(z,cctck)-(ccmax-0.5),(zs[zmax+3],))[0]
                    dzlo = zs[zmax]-zlo
                    dzhi = zhi-zs[zmax]
                    print "Fewer than two roots: solved independently"
                    print "         +%.6f" % (dzhi)
                    print "z = %.6f" % (zsmax)
                    print "         -%.6f" % (dzlo)
                    self.zsmax = zsmax
                    self.ccmax = ccmax
                    self.dzhi = dzhi
                    self.dzlo = dzlo

                    #self.zsmax = 1
                    #self.ccmax = 1
                    #self.dzhi = 1
                    #self.dzlo = 1
                    #print "Error, fewer than 2 roots"

                #Here we compute chi2 and reduce chi2
                Tzs = T.Redshift(zsmax,G.wavelength/(1+z))
                Tzm = np.add.reduce(Tzs*Gw)/np.sum(Gw)
                Tzs = Tzs/Tzm * Gn
                Tps = FitPoly(Tzs,G.noise,G.continuum,G.wavelength,o,w1,w2)
                Gm = 1.*G.good*np.greater_equal(G.wavelength,w1)*np.less_equal(G.wavelength,w2)
                Gr = (G.spec-Gp)*Gm
                Trs = (Tzs-Tps)*Gm
                self.chi2 = np.sum(divz(Trs-Gr,G.noise)**2*Gm)
                self.rchi2 = self.chi2/(np.sum(np.greater(Gm,0))-1)
                print "Template %d,  Chi^2=%.2f   RChi^2=%.2f  CCmax=%.2f  zsmax=%.6f" % (t,self.chi2,self.rchi2,ccmax,zsmax)
                
    

'''
Shell of the code - this begins the cross correlation object by object
Plot() - sets up the object as a class that can be correlated
doCC() - performs the cross correlation
createPlot() - opens the GUI, which allows the user to modify and save the results of the correlation
'''
if tcorr: tcstr = 'cor_'
else: tcstr = ''
if markbads: bads = np.loadtxt(badfile,unpack=True)
if haselist: emitters = np.loadtxt(emissfile,unpack=True)

if readunsure == 1:
    data = np.genfromtxt(outfile,dtype=None,names=True)
    dataunsure = data[data['Unsure'] == 1]
    dataunsure = dataunsure[dataunsure['Star'] == 0]
    dataunsure = dataunsure[dataunsure['Revisit'] == 0]
    dataunsure = dataunsure[dataunsure['Note'] == 0]
    dataunsure = dataunsure[dataunsure['Unusable'] == 0]
    #dataunsure = dataunsure[dataunsure['SN'] > 5]
    tarr = [i for i in reversed(dataunsure[np.argsort(dataunsure['SN'])])]
    imarr = []
    for i in range(len(tarr)):
        imarr.append(imloc+tarr[i]['ImageName'])
    rerun = True
elif objid:
    imarr = glob.glob(imloc + tcstr + objid + '_' + imname)
    if not imarr:
        sys.exit('Cannot find ' + imloc + tcstr + objid + '_' + imname)
    rerun = True
else:
    rerun = False
    outdata = np.genfromtxt(outfile,dtype=None,names=True)
    imarr = glob.glob(imloc+tcstr+'??????_' + imname)
    if not imarr:
        sys.exit('No image files found of form ' + imloc +tcstr+'??????_' + imname)

max_k = len(imarr)
for k in range(max_k):
    j = imarr[k]
    if not rerun:
        if outdata[k]['OBJID'] != 0:
            continue #Already been worked
    im = Plot(j)
    if markbads:
        if float(im.image[:6]) in bads:
            print 'Marked Bad!'
            im.baddata = 1
    if haselist:
        if float(im.image[:6]) in emitters:
            print 'Assuming Emission'
            im.eclip = False
    im.doCC(im.image,im.eclip,0)
    print 'Image {0} / {1}'.format(k+1,max_k)
    print imarr[k]
    im.createPlot()
    
 

