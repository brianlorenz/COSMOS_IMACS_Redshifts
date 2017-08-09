import numpy as np
import glob
from sys import *
from astropy.io import fits
from scipy.interpolate import interp1d,splrep,splev,sproot
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,SpanSelector,CheckButtons
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from stsci.convolve import boxcar

skylines = [5577.0,5890.0,6300,6364]
telluric = [(7584,7650)]

def divz(X,Y):
    return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

def svdfit(b,y):
    decomp = np.linalg.svd(b,full_matrices=False)
    #print decomp[2].shape
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
    filepath = '/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOUT/verb/'
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
    #a = np.compress(u,[x,s,n,g,w],1)
    #print a
    #print a.shape
    #jepwfjweo =fewqpfie
    #xu = np.compress(u,x)
    #su = np.compress(u,s)
    #nu = np.compress(u,n)
    #gu = np.compress(u,g)
    #wu = np.compress(u,w)
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
    def __init__(self,file,clip=0,verb=0):
        self.file = file
        self.bigfile = self.file[11:]
        self.verb = verb
        self.clipemission = clip
        self.Read()
        self.bighdu = fits.open(self.bigfile)[0]
        self.image_data = self.bighdu.data[self.head['csecta']:(self.head['csectb']+1)]
        self.lobar = self.head['lowbound']-self.head['csecta']-1
        self.hibar = self.head['upbound']-self.head['csecta']+1
        self.bighead = self.bighdu.header
       # self.spec = self.data[0]
        self.orig = self.data[1]
        self.noise = self.data[2]
        self.flat = self.data[3]
        self.sky = self.data[4]
        #self.tc = self.data[5]
       # self.mask = self.data[6]
       # self.fit = self.data[7]
        self.objid = self.head['OBJID']
        self.findSN()
        self.GetBadPixels()
		
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
		    
    def GetBadPixels(self):
        if hasattr(self,"noise"):
            bad = 1.*np.logical_or(np.less_equal(self.spec,1e-4),np.less_equal(self.noise,1e-4))
            for skyline in skylines:
                bad = 1.*np.logical_or(bad,np.greater(self.wavelength,skyline-3*self.cdelt1)*np.less(self.wavelength,skyline+3*self.cdelt1))
            for absorb in telluric:
                bad = 1.*np.logical_or(bad,np.greater(self.wavelength,absorb[0])*np.less(self.wavelength,absorb[-1]))
            bad = 1.*np.greater(boxcar(bad,(8,)),0)
        else:
            bad = np.zeros(self.spec.shape,dtype='float')
        sm = Smooth(self.spec,np.logical_not(bad))
        if self.verb: lim = (np.median(self.spec)*0.25,np.median(self.spec)*1.75)
        if self.verb: outGraph(self.wavelength,self.spec,self.objid,'spec',lim=lim)
        if self.verb: outGraph(self.wavelength,sm,self.objid,'sm50',lim=lim)
        #if self.verb: qdump("spec.fits",self.spec,t=self.file)
        #if self.verb: qdump("sm50.fits",sm,t=self.file)
        sr = self.spec-sm
        ss = 1.49*Smooth(abs(sr),np.logical_not(bad),h=50)
        if self.verb: outGraph(self.wavelength,ss,self.objid,'ss50')
        #if self.verb: qdump("ss50.fits",ss,t=self.file)
        #### ADDED CODE HERE - multipled self.good by _not(bad)
        self.good = np.logical_not(bad)
        self.good2 = np.ones(sr.shape,dtype='float32')  # KEEPS EMISSION LINES IN
        if self.clipemission:
            self.good = self.good*np.less(sr,2.5*ss) # CLIPS OUT EMISSION LINES AND BAD POSITIVE SKY LINE RESIDUALS
        else:
            self.good = self.good*np.ones(sr.shape,dtype='float32')  # KEEPS EMISSION LINES IN
        self.continuum = 1.*np.less(abs(sr),2.5*ss)  # NOTE THAT POLYNOMIAL FIT TO CONTINUUM MASKS OUT EMISSION/STRONG ABSORPTION
        if self.verb: outGraph(self.wavelength,self.spec*self.good,self.objid,'sg50',lim=lim)
        #if self.verb: qdump("sg50.fits",self.spec*self.good,t=self.file)
        self.sm = sm
        self.ss = ss
        
class Template(Spectrum):
    def __init__(self,n,clip=0,verb=0):
        self.file = "spDR2-%03d.fit" % (n)
        self.verb = verb
        self.clipemission = clip
        self.Read()
        self.spec = self.data[0]
        self.verb = 0
        self.GetBadPixels()
        wu,su = np.compress(self.good,[self.wavelength,self.spec],1)
        #wu = np.compress(self.good,self.wavelength)
        #su = np.compress(self.good,self.spec)
        self.interp = interp1d(wu,su,fill_value=0,bounds_error=0)
    def Redshift(self,z,w):
        return self.interp(w)/(1.0+z)
        
class Plot():
    def __init__(self,g,m=0):
        #Set parameters here:
        self.image = g
        verb = 0
        self.m = m
        if m:
            #eclip = True
            eclip = False
        else:
            eclip = False
        self.eclip = eclip
        self.wavechange = 0
        self.lowave = 5250
        self.hiwave = 8000
        self.create = 1
        self.range = 100
        self.add = 0
        self.star = 0
        self.baddata = 0
        self.unsure = 0
        self.failure = 0
        self.flag1 = 0
        self.flag2 = 0
        self.flag3 = 0
        #self.doCC(self.image,eclip,verb)
        #self.createPlot(hide=1)

    def doCC(self,g,eclip,verb,newmask=0,wavechange=0):
        print 'Computing Redshift for ' + self.image
        if newmask:
            G = g
        else:
            G = Spectrum(g,clip=eclip,verb=verb)
            self.G = G

        self.t23 = CCcalc(23,G,self.range,eclip,verb,wavechange=wavechange)
        self.t24 = CCcalc(24,G,self.range,eclip,verb,wavechange=wavechange)
        self.t25 = CCcalc(25,G,self.range,eclip,verb,wavechange=wavechange)
        self.t26 = CCcalc(26,G,self.range,eclip,verb,wavechange=wavechange)
        self.t27 = CCcalc(27,G,self.range,eclip,verb,wavechange=wavechange)

        if self.t23.failure or self.t24.failure or self.t25.failure or self.t26.failure or self.t27.failure:
            self.failure = 1
            return
            
        maxarr = np.asarray((self.t23.ccmax,self.t24.ccmax,self.t25.ccmax,self.t26.ccmax,self.t27.ccmax))
        maxindex = np.argmax(maxarr)
        self.tempid = maxindex + 23

        self.findTemp()

    def findTemp(self):
        if self.tempid == 23: self.temp = self.t23
        elif self.tempid == 24: self.temp = self.t24
        elif self.tempid == 25: self.temp = self.t25
        elif self.tempid == 26: self.temp = self.t26
        elif self.tempid == 27: self.temp = self.t27

    def Save(self):
        G = self.G
        f.write(('%06d   %d    %.6f    %.6f   %.6f    %7.2f   %7.2f   %2.2f     %d      %2.2f      %d    %d     %d ' % (G.objid,self.tempid,self.temp.zsmax,self.temp.dzhi,self.temp.dzlo, self.temp.ccmax, self.temp.chi2,self.temp.rchi2,self.eclip,G.signoise,self.star,self.baddata,self.unsure)) + self.image + ('      %d      %d      %d' % (self.flag1, self.flag2, self.flag3)) + '\n')
        f2.write(('%06d %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %d %d %3.2f %d %d %d' % (G.objid,self.t23.zsmax,self.t23.dzhi,self.t23.dzlo, self.t23.ccmax, self.t23.chi2,self.t23.rchi2,self.t24.zsmax,self.t24.dzhi,self.t24.dzlo, self.t24.ccmax, self.t24.chi2,self.t24.rchi2,self.t25.zsmax,self.t25.dzhi,self.t25.dzlo, self.t25.ccmax, self.t25.chi2,self.t25.rchi2,self.t26.zsmax,self.t26.dzhi,self.t26.dzlo, self.t26.ccmax, self.t26.chi2,self.t26.rchi2,self.t27.zsmax,self.t27.dzhi,self.t27.dzlo, self.t27.ccmax, self.t27.chi2,self.t27.rchi2,self.tempid,self.eclip,G.signoise,self.star,self.baddata,self.unsure)) + ' '+ self.image + (' %d %d %d' % (self.flag1, self.flag2, self.flag3)) + '\n')
        filelocation2 = '/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '/' + 'cor_' + self.image.replace('.fits','.png')
        self.fig.savefig(filelocation2)
        plt.clf()
        plt.close()

    def mkspecplot(self):
        self.lfile = np.genfromtxt('galaxylines.dat',dtype=None,names=True)
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
            elif ctype == 3: color = 'mediumseagreen'
            else: color = 'cornflowerblue'
            self.a2x.plot((shiftloc,shiftloc),(-100000,100000),color=color,linestyle='--',alpha=0.75)
        
    
    def createPlot(self,m=0):
        #if self.m: self.mkspecplot()
        if self.failure:
            if not m:
                G = self.G
                f.write(('%06d   %d    %.6f    %.6f   %.6f    %7.2f   %7.2f   %2.2f     %d      %2.2f      %d    %d     %d ' % (G.objid,0,0,0,0,0,0,0,0,0,0,1,1)) + self.image + ('      %d       %d       %d' % (0,0,1)) + '\n')
                f2.write(('%06d %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %d %d %3.2f %d %d %d' % (G.objid,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1)) + ' '+ self.image + (' %d %d %d' % (0,0,1)) + '\n')
                return
            else:
                G = self.G
                f3 = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '.txt',"r+")
                f4 = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + 'verb_' + letnum + '.txt','r+')
                d = f3.readlines()
                e = f4.readlines()
                f3.seek(0)
                f4.seek(0)
                for i in d:
                    if self.image in i:
                        f3.write(('%06d   %d    %.6f    %.6f   %.6f    %7.2f   %7.2f   %2.2f     %d      %2.2f      %d    %d     %d ' % (G.objid,0,0,0,0,0,0,0,0,0,0,1,1)) + self.image + ('      %d      %d      %d' % (0,0,1)) + '\n')
                    else:
                        f3.write(i)
                for i in e:
                    if self.image in i:
                        f4.write(('%06d %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %d %d %3.2f %d %d %d' % (G.objid,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1)) + ' '+ self.image + (' %d %d %d' % (0,0,1)) + '\n')
                    else:
                        f4.write(i)
                f3.close()
                f4.close()
                return
                
        if self.create:
            G = self.G
            self.fig = plt.figure(figsize=(14,8))
            #self.ax1 = plt.subplot2grid((5,8),(0,0),colspan=4,rowspan=2)
            self.ax1 = self.fig.add_axes([0.05, 0.7, 0.43, 0.25])
            #self.ax2 = plt.subplot2grid((5,8),(0,4),colspan=4,rowspan=2)
            self.ax2 = self.fig.add_axes([0.53, 0.7, 0.43, 0.20])
            #self.ax2b = self.fig.add_axes([0.75, 0.43, 0.1, 0.2])
            #self.ax2c = self.fig.add_axes([0.87, 0.43, 0.1, 0.2])
            self.ax2b = self.fig.add_axes([0.65, 0.40, 0.3, 0.2])
            self.ax2c = self.fig.add_axes([0.65, 0.1, 0.3, 0.2])
            #self.ax3 = plt.subplot2grid((5,8),(2,1),colspan=5,rowspan=3)
            self.ax3 = self.fig.add_axes([0.19, 0.1, 0.45, 0.50])
            #self.ax4 = self.fig.add_axes([0.75, 0.05, 0.2, 0.10])
            #self.ax4 = plt.subplot2grid((5,8),(3,6),rowspan=2,colspan=2)
            self.ax5 = self.fig.add_axes([0.53, 0.9, 0.43, 0.08])
            self.ax1.set_title('Spectrum of ' + self.image + ', S/N = ' + str(G.signoise))
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
            #self.ax4.yaxis.tick_right()
            self.ax2.set_xlabel('$\lambda$                                 ')
            self.ax2.set_ylabel('(g-f)/$\sigma$')
            self.ax2b.set_xlabel('$\lambda$')
            
            self.ax2c.set_xlabel('$\lambda$')
            self.m1, = self.ax3.plot(self.temp.zs[0],self.temp.cc[0],color='indianred')
            if self.temp.zp>0: self.m2, = self.ax3.plot(self.temp.zs[1],self.temp.cc[1],color='firebrick')
            else: self.m2, = self.ax3.plot((0,0),(0,0),color='firebrick')
            if self.temp.zp>1: self.m3, = self.ax3.plot(self.temp.zs[2],self.temp.cc[2],color='darkred')
            else: self.m3, = self.ax3.plot((0,0),(0,0),color='darkred')
            if self.temp.zp>2: self.m4, = self.ax3.plot(self.temp.zs[3],self.temp.cc[3],color='darkred')
            else: self.m4, = self.ax3.plot((0,0),(0,0),color='darkred')
            self.g1, = self.ax3.plot((self.temp.zs[self.temp.zp][self.temp.zmax],self.temp.zs[self.temp.zp][self.temp.zmax]),(-10000,10000),color='mediumseagreen')
            self.ax3.set_xlim(0.19,0.51)
            self.ax3.set_ylim(min(self.temp.cc[0])*0.95,max(self.temp.cc[0]*1.05))
            self.ax3.set_xlabel('z')
            self.ax3.set_ylabel('Correlation Coeffecient')
            #self.h1, = self.ax4.plot(self.temp.zs[self.temp.zp],self.temp.cc[self.temp.zp],color='firebrick')
            #self.g2, = self.ax4.plot((self.temp.zs[self.temp.zp][self.temp.zmax],self.temp.zs[self.temp.zp][self.temp.zmax]),(-10000,10000),color='mediumseagreen')
            #self.ax4.set_xlim(min(self.temp.zs[self.temp.zp]),max(self.temp.zs[self.temp.zp]))
            #self.ax4.set_ylim(max(0,self.temp.ccmax-10),self.temp.ccmax)
            #self.ax4.set_title('CC Peak Zoom')
            #self.ax4.set_xlabel('z')
            #self.ax4.set_ylabel('Correlation Coeffecient')
            #self.ax4.tick_params(axis='both',labelsize='small')
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
            self.temp.rerunCC(self.newcc,self.range)
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
            self.lowave = 5250
            self.hiwave = 8000
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
            self.doCC(self.image,self.eclip,0)
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
            f3 = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '.txt',"r+")
            f4 = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + 'verb_' + letnum + '.txt','r+')
            d = f3.readlines()
            e = f4.readlines()
            f3.seek(0)
            f4.seek(0)
            for i in d:
                if self.image in i:
                    f3.write(('%06d   %d    %.6f    %.6f   %.6f    %7.2f   %7.2f   %2.2f     %d      %2.2f      %d    %d     %d ' % (G.objid,self.tempid,self.temp.zsmax,self.temp.dzhi,self.temp.dzlo, self.temp.ccmax, self.temp.chi2,self.temp.rchi2,self.eclip,G.signoise,self.star,self.baddata,self.unsure)) + self.image + ('      %d      %d      %d' % (self.flag1, self.flag2, self.flag3)) + '\n')
                else:
                    f3.write(i)
            for i in e:
                if self.image in i:
                    f4.write(('%06d %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %.6f %.6f %.6f %7.3f %7.3f %2.2f %d %d %3.2f %d %d %d' % (G.objid,self.t23.zsmax,self.t23.dzhi,self.t23.dzlo, self.t23.ccmax, self.t23.chi2,self.t23.rchi2,self.t24.zsmax,self.t24.dzhi,self.t24.dzlo, self.t24.ccmax, self.t24.chi2,self.t24.rchi2,self.t25.zsmax,self.t25.dzhi,self.t25.dzlo, self.t25.ccmax, self.t25.chi2,self.t25.rchi2,self.t26.zsmax,self.t26.dzhi,self.t26.dzlo, self.t26.ccmax, self.t26.chi2,self.t26.rchi2,self.t27.zsmax,self.t27.dzhi,self.t27.dzlo, self.t27.ccmax, self.t27.chi2,self.t27.rchi2,self.tempid,self.eclip,G.signoise,self.star,self.baddata,self.unsure)) + ' '+ self.image + (' %d %d %d' % (self.flag1, self.flag2, self.flag3)) + '\n')
                else:
                    f4.write(i)
            f3.close()
            f4.close()
            filelocation2 = '/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '/' + 'cor_' + self.image.replace('.fits','.png')
            filelocation3 = '/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '/' + 'cor_' + self.image.replace('.fits','2.png')
            self.fig.savefig(filelocation2)
            #self.fig2.savefig(filelocation3)
            plt.clf()
            plt.close(self.fig)
            #plt.close(self.fig2)

        def checkStarFlag():
            if self.star == 1: bsflag.label.set_text('Unmark Star')
            else: bsflag.label.set_text('Mark Star')
            plt.draw()

        def starFlag(event):
            if self.star == 1: self.star = 0
            else: self.star = 1
            checkStarFlag()
        sflag = plt.axes([0.15, 0.0, 0.08, 0.05])
        bsflag = Button(sflag, 'Mark Star')
        bsflag.on_clicked(starFlag)
        if G.spec[1000:1500].mean() > 10000:
            self.star = 1
        checkStarFlag()

        def dataFlag(event):
            if self.baddata == 1:
                self.baddata = 0
                bdflag.label.set_text('Mark Bad')
            else:
                self.baddata = 1
                bdflag.label.set_text('Unmark Bad')
            plt.draw()
        dflag = plt.axes([0.25, 0.0, 0.08, 0.05])
        bdflag = Button(dflag, 'Mark Bad')
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
                bflag2.label.set_text('Mark Range')
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

        if m:
            save = plt.axes([0.0, 0.0, 0.1, 0.05])
            bsave = Button(save, 'Save')
            bsave.on_clicked(Bsave)
        
        #plt.tight_layout(pad=3.0)
        if m:
            #self.mkspecplot()
            plt.show()
        else:
            plt.show(block=False)
            move = raw_input('Enter when done')
            self.savetext = self.fig.text(0.5,0.01,'Saving...',fontsize=24)
            plt.pause(.1)
            plt.gcf().texts.remove(self.savetext)
            self.Save()
    
        
    
        
            
    def updatePlot(self):
        self.k2.set_ydata(self.temp.Trs[self.temp.zp][self.temp.zmax])
        self.k2b.set_ydata(self.temp.Trs[self.temp.zp][self.temp.zmax])
        self.k2c.set_ydata(self.temp.Trs[self.temp.zp][self.temp.zmax])
        self.m1.set_ydata(self.temp.cc[0])
        self.m1.set_xdata(self.temp.zs[0])
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
        #self.g2.set_xdata((self.temp.zs[self.temp.zp][self.temp.zmax],self.temp.zs[self.temp.zp][self.temp.zmax]))
        #self.h1.set_xdata(self.temp.zs[self.temp.zp])
        #self.h1.set_ydata(self.temp.cc[self.temp.zp])
        self.ax3.set_ylim(min(self.temp.cc[0])*0.95,max(self.temp.cc[0]*1.05))
        #self.ax4.set_xlim(min(self.temp.zs[self.temp.zp]),max(self.temp.zs[self.temp.zp]))
        #self.ax4.set_ylim(max(0,self.temp.ccmax-10),self.temp.ccmax)
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
        #self.r2.set_ydata((self.cc[self.zp][self.zmax]-0.5,self.cc[self.zp][self.zmax]-0.5))
        self.refreshText(remove = 1)
        if self.m:
            #plt.close(self.fig2)
            #self.mkspecplot()
            #plt.show(block=False)
            #plt.figure(2)
            pass
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
           # plt.gcf().texts.remove(self.SNtext)
        #self.ztext = self.fig.text(0.476,0.552,"z = %.6f" % (self.temp.zsmax),fontsize = 12)
        #self.zhitext = self.fig.text(0.557,0.566,"+%.6f" % (self.temp.dzhi),fontsize = 8)
        #self.zlotext = self.fig.text(0.5605,0.546,"-%.6f" % (self.temp.dzlo),fontsize = 8)
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
        #self.SNtext = self.fig.text(0.47,0.93,'S/N = %4.2f'%self.G.signoise,fontsize = 12)
       
        
class CCcalc:
    def __init__(self,t,G,inrange,eclip,verb,wavechange=0):
        self.G = G
        self.range = inrange
        self.t = t
        self.failure = 0
    
        T = Template(t,clip=eclip,verb=verb)
        self.T = T
        
        z1,z2 = 0.2,0.5
        dz0 = 4*G.cdelt1/6500.0
        nz = int((z2-z1)/dz0+1)
        zp = 0
        keepgoing = True
        w1,w2 = (3500.0*(1+z2)),8000.0
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

    def rerunCC(self,newcc,inrange):
        self.range = inrange
        self.cc = [self.cc[0]]
        self.zs = [self.zs[0]]
        self.Trs = [self.Trs[0]]
        zmax = self.zs[0][int(newcc)]
        z1,z2,nz = 0.2,0.5,100
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
            while keepgoing:
                zs = np.arange(z1,z2+dz/2,dz,dtype='float64')
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
                    #if verb: qdump("pl.fits",Gp,t=g)

                    Gp2 = FitPoly(G.spec,G.noise,G.continuum,G.wavelength,20,4900,10000)
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
                    #qdump("cc%d%d.fits"%(t,zp),cc,extras=[("crval1",z1),("cdelt1",dz),("crpix1",1),("dc-flag",0)],verb=0)

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
                        dz = 0.3/99.0
                        intro = 0
                    
                    if keepgoing:
                        if self.range == 100:
                            z1,z2 = zs[zmax]-5*dz,zs[zmax]+4.9*dz
                        elif self.range == 60:
                            z1,z2 = zs[zmax]-3*dz,zs[zmax]+2.9*dz
                        elif self.range == 20:
                            z1,z2 = zs[zmax]-1*dz,zs[zmax]+0.9*dz
                        dz = dz/10.0
                        zp += 1
                        
            if not self.failure:
                cctck = splrep(zs,cc,s=0,k=3)
                zsmax = fsolve(lambda z: splev(z,cctck,1),(zs[zmax],))[0]
                ccmax = splev(zsmax,cctck)
                #print "Template = %d" % (t)
        
                #print zsmax,ccmax
                
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
                except(TypeError):
                    ccroots = [0]
                self.zp = zp
                if len(ccroots)>=2:
                    zlo,zhi = np.sort(np.take(ccroots,np.argsort(abs(ccroots-zsmax)))[:2])
                    dzlo,dzhi = zsmax-zlo,zhi-zsmax
                    print "         +%.6f" % (dzhi)
                    print "z = %.6f" % (zsmax)
                    print "         -%.6f" % (dzlo)
                    self.zsmax = zsmax
                    self.ccmax = ccmax
                    self.dzhi = dzhi
                    self.dzlo = dzlo
                else:
                    self.zsmax = 1
                    self.ccmax = 1
                    self.dzhi = 1
                    self.dzlo = 1
                    print "Error, less than 2 roots"
                
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
                
    
objid = 0
#if len(argv[1]) == 2:
#    letnum = argv[1]
elif len(argv[1]) == 9:
    objid = argv[1][:6]
    letnum = argv[1][7:]
elif  argv[1] == 'unsure':
    a6 = np.genfromtxt('crossCorOut/a6.txt',dtype=None,names=True)
    b6 = np.genfromtxt('crossCorOut/b6.txt',dtype=None,names=True)
    d6 = np.genfromtxt('crossCorOut/d6.txt',dtype=None,names=True)
    e6 = np.genfromtxt('crossCorOut/e6.txt',dtype=None,names=True)
    f6 = np.genfromtxt('crossCorOut/f6.txt',dtype=None,names=True)
    g6 = np.genfromtxt('crossCorOut/g6.txt',dtype=None,names=True)
    h6 = np.genfromtxt('crossCorOut/h6.txt',dtype=None,names=True)
    i6 = np.genfromtxt('crossCorOut/i6.txt',dtype=None,names=True)
    b7 = np.genfromtxt('crossCorOut/e7.txt',dtype=None,names=True)
    e7 = np.genfromtxt('crossCorOut/b7.txt',dtype=None,names=True)

    data = np.concatenate((a6,b6,d6,e6,f6,g6,h6,i6,b7,e7))
    dataunsure = data[data['Unsure'] == 1]
    dataunsure = dataunsure[dataunsure['Star'] == 0]
    dataunsure = dataunsure[dataunsure['Flag1'] == 0]
    dataunsure = dataunsure[dataunsure['Flag2'] == 0]
    dataunsure = dataunsure[dataunsure['Flag3'] == 0]
    #dataunsure = dataunsure[dataunsure['SN'] > 5]
    #unsureims = [(j['ImageName'],j['SN']) for j in dataunsure]
    #print unsureims
    tarr = [i for i in reversed(dataunsure[np.argsort(dataunsure['SN'])])]
    imarr = []
    for i in range(len(tarr)):
        imarr.append(tarr[i]['ImageName'])
    #imarr = [j for j in dataunsure['ImageName'] if j[17]+j[15] == 'g6']
    im = [Plot(i,m=1) for i in imarr]
    #for j in range(len(im)):
    #    im[j].doCC(im[j].image,im[j].eclip,0)
    for k in range(len(im)):
        im[k].doCC(im[k].image,im[k].eclip,0)
        print 'Image ', k+1, ' / ', len(im)
        print imarr[k]
        letnum = imarr[k][17]+imarr[k][15]
        im[k].createPlot(m=1)
    exit('Finished unsure objects')
else:
     exit('Enter a mask, eg "a6", to run full mask. \nEnter objid and mask, eg "123456 a6", to open just that object')
exten = 'feb1' + letnum[1] + '_' + letnum[0] + 'big.fits'
if objid:
    im = Plot(glob.glob('cor_' + objid + '_' + exten)[0],m=1)
    im.doCC(im.image,im.eclip,0)
    #f = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '.txt','a+')
    #outpkl = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + 'cc_' + objid + '.pkl','wb')
    im.createPlot(m=1)
    #dilltest = dill.dumps(im)
    #readpkl = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + 'cc_' + objid + '.pkl','rb')
    #testpkl = dill.loads(dilltest)
    #testpkl.createPlot(m=1)
    
else:
    yesno = raw_input('Are you sure you want to overwrite template ' + letnum + '? (y/n) ')
    if yesno == 'y':
        imarr = glob.glob('cor_??????_' + exten)
        f = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + letnum + '.txt','w+')
        f2 = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/crossCorOut/' + 'verb_' + letnum + '.txt','w+')
        f.write('#OBJID  temp    z           dzhi        dzlo        ccmax     chi2    rchi2  eclip     S/N      Star Bad  Unsure    ImageName            Flag1  Flag2  Flag3\n')
        f2.write('#OBJID z23 dzhi23 dzlo23 ccmax23 chi223 rchi223 z24 dzhi24 dzlo24 ccmax24 chi224 rchi224 z25 dzhi25 dzlo25 ccmax25 chi225 rchi225 z26 dzhi26 dzlo26 ccmax26 chi226 rchi226 z27 dzhi27 dzlo27 ccmax27 chi227 rchi227 temp eclip S/N Star Bad Unsure ImageName Flag1 Flag2 Flag3\n')
        im = [Plot(i) for i in imarr]
        for j in range(len(im)):
            im[j].doCC(im[j].image,im[j].eclip,0)
        #for k in range(len(im)):
            print 'Image ', j+1, ' / ', len(im)
            print imarr[j]
            im[j].createPlot()
        f.close()
        f2.close()
    else:
        pass
    
 

