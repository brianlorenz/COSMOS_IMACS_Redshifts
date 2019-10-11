from sys import argv
import glob 
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,SpanSelector,CheckButtons
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy import interpolate


class Plot():
    def __init__(self,letnum,imnum=0):
        self.order = 15
        self.new = 1
        self.letnum = letnum
        self.create = 1
        if not imnum:
            self.imnum = 44
            #self.imnum=65
            #self.imnum=81
            self.openFiles()
        else:
            self.imnum=imnum
            self.openFiles(m=1)

    def divz(self,X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

    def openFiles(self,m=0):
        exten = 'feb1' + self.letnum[1] + '_' + self.letnum[0] + 'big.fits'
        self.imagearr = glob.glob('??????_' + exten)
        self.tcfile = glob.glob('tc_' + '*' + exten)[0]
        hdutc = fits.open(self.tcfile)[0]
        self.headertc = hdutc.header
        self.tc,spec,mask,fit = hdutc.data.astype(float)
        if not m:
            for i in range(len(self.imagearr)):
                self.imnum = i
                self.openImage()
        else:
            self.openImage(m=1)

    def openImage(self,m=0):
        self.image = self.imagearr[self.imnum]
        imagefile = fits.open(self.image)
        hdu = imagefile[0]
        self.header = hdu.header
        self.obj,self.noise,self.flat,self.sky = hdu.data.astype(float)
        imagefile.close()
        
        
        crval1 = self.header["crval1"]
        self.cdelt1 = self.header["cdelt1"]
        crpix1 = self.header["crpix1"]
        dcflag = self.header["dc-flag"]
        self.x = np.arange(0,self.obj.shape[0],1,float)
        self.wavelength = crval1 + (1.+self.x-crpix1)*self.cdelt1
        if dcflag: self.wavelength = np.power(10.0,w)
        self.runCalc()
        if not m:
            self.createPlot()
            self.saveSpec()
        else:
            self.createPlot(m=1)
       

    def Smooth(self,y,p=50):
        m = np.zeros(y.shape,y.dtype)
        for j in range(len(y)):
            a,b = np.clip([j-25,j+26],0,len(y)-1)
            if p == 50:
                m[j] = np.median(y[a:b])
            else:
                m[j] = np.sort(y[a:b])[len(y[a:b])*p/100]
        return m

    def svdfit(self,b,y):
        decomp = np.linalg.svd(b,full_matrices=False)
        #print decomp[2].shape
        sol1 = np.transpose(decomp[2])
        sol2 = self.divz(1.0,decomp[1])
        sol3 = np.dot(np.transpose(decomp[0]),y)
        if np.sometrue(sol3):
            solr = (sol2*sol3)
            soll = np.dot(sol1,solr)
        else:
            soll = np.zeros(sol3.shape)
        return soll

    def runCalc(self,redo=0):
        self.m = self.Smooth(self.obj,p=95)
        self.m = self.m - np.median(self.m-self.obj)

        self.diff = self.obj-self.m
        q = 1.49*self.Smooth(abs(self.diff),p=50)
        mq = np.median(q)

        if not redo:
            self.mask = np.less(q,mq)*np.less(abs(self.diff),5*q)*np.less(self.wavelength,8500.0)

        self.norm = (self.x-(self.x[0]+self.x[-1])/2)/ ((self.x[-1]-self.x[0])/2)
        basis = np.asarray([np.power(self.norm,i) for i in range(self.order)])
        sol = self.svdfit(np.transpose(basis)*self.mask[::,np.newaxis],self.obj*self.mask)

        self.fit = np.dot(sol,basis)
        self.norm = self.divz(self.obj,self.fit)

        self.normgal = self.norm-1
        self.tc1 = self.tc-1

        mask = np.greater(self.obj,1e-8).astype(float)
        mask = mask*np.less(self.norm,np.sort(self.norm)[int(round(0.95*len(self.norm)))])
        self.mask2 = mask

        #fig,axarr=plt.subplots(3,1,figsize=(10,8))
        #ax = axarr[0]
        #ax2 = axarr[1]
        #ax3 = axarr[2]
        
        self.divisions = 4.0
        self.cr = np.zeros((40*int(self.divisions)+1,),float)
        tss = []
        self.plotarr = []
        for j in range(len(self.cr)):
            k = j/self.divisions-(len(self.cr)-1)/(self.divisions*2)
            telspeci = interpolate.interp1d(self.wavelength,self.tc1,fill_value = 0,bounds_error=False)#,kind='cubic')
            telspec = telspeci(self.wavelength+self.cdelt1*k)
            telmaski = interpolate.interp1d(self.wavelength,mask,fill_value = 0,bounds_error=False)#,kind='cubic')
            telmask = telmaski(self.wavelength+self.cdelt1*k)
            #for n in range(len(self.wavelength)):
                #telspec = np.interp(self.wavelength[n]+k,self.wavelength,self.tc1)
                #telmask = np.interp(self.wavelength[n]+k,self.wavelength,mask)
            #telspec = np.interp(self.wavelength+self.cdelt1*k,self.wavelength,self.tc1)
            #telmask = np.interp(self.wavelength+self.cdelt1*k,self.wavelength,mask)
            #if k == 0: telspec = self.tc1
            #elif k < 0: telspec = np.concatenate([abs(k)*[0.0],self.tc1[:k]])
            #elif k > 0: telspec = np.concatenate([self.tc1[k:],k*[0.0]])
            #if k == 0: telmask = mask
            #elif k < 0: telmask = np.concatenate([abs(k)*[0.0],mask[:k]])
            #elif k > 0: telmask = np.concatenate([mask[k:],k*[0.0]])
            self.cr[j] = np.sum(telspec*self.normgal*mask*telmask)
            product = telspec*self.normgal*mask*telmask
            if k == 30:
                ax.plot(self.wavelength,telspec,color='red')
                ax.plot(self.wavelength,telmask,color='black')
                ax.plot(self.wavelength,self.normgal,color='blue')
                ax.plot(self.wavelength,mask,color='green')
                ax.plot(self.wavelength,product,color='orange')
                ax.set_title(str(self.cr[j]))
            if k == 30.5:
                ax2.plot(self.wavelength,telspec,color='red')
                ax2.plot(self.wavelength,telmask,color='black')
                ax2.plot(self.wavelength,self.normgal,color='blue')
                ax2.plot(self.wavelength,mask,color='green')
                ax2.plot(self.wavelength,product,color='orange')
                ax2.set_title(str(self.cr[j]))
            if k == 31:
                ax3.plot(self.wavelength,telspec,color='red')
                ax3.plot(self.wavelength,telmask,color='black')
                ax3.plot(self.wavelength,self.normgal,color='blue')
                ax3.plot(self.wavelength,mask,color='green')
                ax3.plot(self.wavelength,product,color='orange')
                ax3.set_title(str(self.cr[j]))
                ax.set_ylim(-1.1,1.1)
                ax.set_xlim(7550,7800)
                ax2.set_ylim(-1.1,1.1)
                ax2.set_xlim(7550,7800)
                ax3.set_ylim(-1.1,1.1)
                ax3.set_xlim(7550,7800)
                plt.tight_layout()
                plt.show()
            tss.append(telspec)
            self.plotarr.append(k)

        self.peak = (np.argmax(self.cr)/self.divisions)-(len(self.cr)-1)/(2*self.divisions)
        print 'Peak =', self.peak

        telspec = np.interp(self.wavelength+self.cdelt1*self.peak,self.wavelength,self.tc)
        #if self.peak == 0: telspec = self.tc
        #elif self.peak < 0: telspec = np.concatenate([np.abs(self.peak)*[1.0],self.tc[:self.peak]])
        #elif self.peak > 0: telspec = np.concatenate([self.tc[self.peak:],self.peak*[1.0]])

        self.corrgal = self.divz(self.obj,telspec)
        self.corrnorm = self.divz(self.norm,telspec)
        self.telspec = telspec
        

    def createPlot(self,m=0):
        if self.create or not m:
            self.fig = plt.figure(figsize=(10,8))
            self.ax1 = plt.subplot2grid((3,3),(0,0),colspan=3)
            self.ax2 = plt.subplot2grid((3,3),(1,0),colspan=1)
            self.ax3 = plt.subplot2grid((3,3),(1,1),colspan=1)
            self.ax4 = plt.subplot2grid((3,3),(1,2),colspan=1)
            self.ax5 = plt.subplot2grid((3,3),(2,0),colspan=3)
            self.create = 0
        self.updatePlot()
        if m:
            def Right(event):
                self.peak -= 1
                self.updatePlot(redo=1)
                plt.draw()

            def Right5(event):
                self.peak -= 5
                self.updatePlot(redo=1)
                plt.draw()
                
            def Left(event):
                self.peak += 1
                self.updatePlot(redo=1)
                plt.draw()

            def Left5(event):
                self.peak += 5
                self.updatePlot(redo=1)
                plt.draw()

            def Save(event):
                self.saveSpec(m=1)
                        
            right = plt.axes([0.56, 0.337, 0.03, 0.03])
            bright = Button(right, '>')
            bright.on_clicked(Right)
            right5 = plt.axes([0.60, 0.337, 0.03, 0.03])
            bright5 = Button(right5, '>>')
            bright5.on_clicked(Right5)
            left = plt.axes([0.47, 0.337, 0.03, 0.03])
            bleft = Button(left, '<')
            bleft.on_clicked(Left)
            left5 = plt.axes([0.43, 0.337, 0.03, 0.03])
            bleft5 = Button(left5, '<<')
            bleft5.on_clicked(Left5)
            save = plt.axes([0, 0, 0.05, 0.05])
            bsave = Button(save, 'Save')
            bsave.on_clicked(Save)
            plt.show()

    def updatePlot(self,redo=0):
        if redo:
            telspec = np.interp(self.wavelength+self.cdelt1*self.peak,self.wavelength,self.tc)
            self.corrgal = self.divz(self.obj,telspec)
            self.corrnorm = self.divz(self.norm,telspec)
            self.telspec = telspec
            self.m1.set_ydata(self.corrnorm)
            self.m2.set_ydata(self.telspec)
            self.j2.set_xdata((self.peak,self.peak))
            self.i1.set_ydata(self.corrgal)
            self.ax3.set_title('Corrected, shift = ' + str(int(self.peak)))
        else:
            self.ax1.set_title('Spectrum of ' + self.image + ', imnum = ' + str(self.imnum))
            self.ax2.set_title('Uncorrected')
            self.ax3.set_title('Corrected, shift = ' + str(int(self.peak)))
            self.ax4.set_title('Cross Correlation')
            self.ax5.set_title('Full Corrected Spectrum')
            l1, = self.ax1.plot(self.wavelength,self.obj,color='cornflowerblue')
            l2, = self.ax1.plot(self.wavelength,self.obj*self.mask,color='orange')
            l3, = self.ax1.plot(self.wavelength,self.fit,color='darkgreen')
            self.ax1.set_ylim(np.median(self.obj)*-0.5,np.median(self.obj)*4)
            self.ax1.set_xlabel('$\lambda$')
            self.ax1.set_ylabel('Counts')
            k1, = self.ax2.plot(self.wavelength,self.normgal,color='cornflowerblue')
            k2, = self.ax2.plot(self.wavelength,self.tc1,color='black')
            k2, = self.ax2.plot(self.wavelength,self.mask2,color='orange')
            self.ax2.set_xlim(7100,8500)
            self.ax2.set_ylim(-0.1,1.1)
            self.ax2.set_xlabel('$\lambda$')
            self.ax2.set_ylabel('Normalized Counts')
            self.m1, = self.ax3.plot(self.wavelength,self.corrnorm,color='cornflowerblue')
            self.m2, = self.ax3.plot(self.wavelength,self.telspec,color='black')
            self.ax3.set_xlim(7100,8500)
            self.ax3.set_ylim(0,2)
            self.ax3.set_xlabel('$\lambda$')
            self.ax3.set_ylabel('Normalized Counts')
            j1, = self.ax4.plot(self.plotarr,self.cr,color='indianred')#,ls='none',marker='.')
            self.j2, = self.ax4.plot((self.peak,self.peak),(-100,100),color='mediumseagreen')
            self.ax4.set_xlabel('Pixels to shift')
            self.ax4.set_ylabel('Correlation Coefficient')
            self.ax4.set_ylim(min(self.cr)*.95,max(self.cr)*1.05)
            self.i1, = self.ax5.plot(self.wavelength,self.corrgal,color='cornflowerblue')
            self.ax5.set_ylim(np.median(self.obj)*-0.5,np.median(self.obj)*4)
            self.ax5.set_xlabel('$\lambda$')
            self.ax5.set_ylabel('Counts')
            plt.tight_layout()      
            
    
            

    def saveSpec(self,m=0):
        headerout = self.header
        [headerout.pop(key) for key in headerout.keys() if key[:3] == "ROW"]
        headerout['ROW1'] = 'Corrected Spectrum'
        headerout['ROW2'] = 'Original Spectrum'
        headerout['ROW3'] = 'Noise'
        headerout['ROW4'] = 'Flat'
        headerout['ROW5'] = 'Sky'
        headerout['ROW6'] = 'Shifted Telluric Correction'
        headerout['ROW7'] = 'Pixel Mask'
        headerout['ROW8'] = 'Polynomial Fit'
        headerout['NTERMS'] = self.order
        headerout['TC_FILE'] = self.tcfile
        headerout['TC_SHIFT'] = str(int(self.peak)) + ' pixels'
        dataout = np.concatenate((self.corrgal,self.obj,self.noise,self.flat,self.sky,self.telspec,self.mask,self.fit))
        dataout = dataout.reshape(8,len(self.wavelength))      
        hdu = fits.PrimaryHDU(header = headerout, data = dataout)
        filelocation = '/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/' + 'cor_' + self.image
        hdu.writeto(filelocation,overwrite=True)
        filelocation2 = '/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/corImageOut/' + self.letnum + '/' + 'cor_' + self.image.replace('.fits','.png')
        if not m:
            self.fig.savefig(filelocation2)
        else:
            self.fig.savefig(filelocation2.replace('cor_','cor_h_'))
        print('Correction saved to ' + filelocation)
        if self.new and not m:
            f = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/' + 'cor_out_' + self.letnum + '.txt','w+')
            f.write('#Objid Objnum tc_shift imnum\n')
            self.new = 0;
        else:
            f = open('/Users/blorenz/Desktop/COSMOSData/corFitsFileOut/' + 'cor_out_' + self.letnum + '.txt','a+')
        if not m:
            f.write('%06d %3d    %4d     %3d\n' % (self.header['OBJID'],self.header['OBJNUM'],self.peak,self.imnum))
        else:
            f.write('%06d %3d    %4d     %3d Fixed by hand\n' % (self.header['OBJID'],self.header['OBJNUM'],self.peak,self.imnum))
        f.close()
        plt.close(self.fig)

        

