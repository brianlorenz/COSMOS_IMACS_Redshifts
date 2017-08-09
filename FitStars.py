from sys import argv
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,SpanSelector,CheckButtons
from matplotlib import gridspec

class Plot():
    def __init__(self):
        self.image = argv[1]
        self.order = int(argv[2])
        self.add = 1

        hdu = fits.open(self.image)[0]
        self.header = hdu.header
        data = hdu.data
        self.obj,noise,flat,sky = data.astype(float)
    
        crval1 = self.header["crval1"]
        cdelt1 = self.header["cdelt1"]
        crpix1 = self.header["crpix1"]
        dcflag = self.header["dc-flag"]
        self.x = np.arange(0,self.obj.shape[0],1,float)
        self.wavelength = crval1 + (1.+self.x-crpix1)*cdelt1
        if dcflag: self.wavelength = np.power(10.0,w)
        self.runCalc()
        self.createPlot()

        
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
  
    def divz(self,X,Y):
        return X/np.where(Y,Y,Y+1)*np.not_equal(Y,0)

    def updatePlot(self):
        self.runCalc(redo=1)
        self.l3.set_ydata(self.m)
        self.l4.set_ydata(self.diff)
        self.l5.set_ydata(self.fit)
        self.k1.set_ydata(self.corr)
        self.ax1.legend(loc=1)
        plt.draw()

    def runCalc(self,redo=0,fixtel=0):
        self.m = self.Smooth(self.obj,p=95)
        self.m = self.m - np.median(self.m-self.obj)

        self.diff = self.obj-self.m
        q = 1.49*self.Smooth(abs(self.diff),p=50)
        mq = np.median(q)

        if not redo:
            self.mask = np.less(q,mq)*np.less(abs(self.diff),5*q)*np.less(self.wavelength,8500.0)

        norm = (self.x-(self.x[0]+self.x[-1])/2)/ ((self.x[-1]-self.x[0])/2)
        basis = np.asarray([np.power(norm,i) for i in range(self.order)])
        sol = self.svdfit(np.transpose(basis)*self.mask[::,np.newaxis],self.obj*self.mask)

        self.fit = np.dot(sol,basis)
        n = self.divz(self.obj,self.fit)
        
        Bband = 6865,7060
        H20band = 7148,7340
        Aband = 7580,7730
        H20band2 = 8120,8410
        between = lambda X,Y,Z: np.greater(X,Y)*np.less(X,Z)
        atm = between(self.wavelength,Bband[0],Bband[1]) + between(self.wavelength,Aband[0],Aband[1]) + between(self.wavelength,H20band[0],H20band[1]) + between(self.wavelength,H20band2[0],H20band2[1])
        if not fixtel:
            self.good = np.logical_not(atm)
        self.corr = np.where(self.good,1.0,n)
        

    def createPlot(self):
        #Plot for poster/presentation/writeup
        self.fig2 , self.a2x = plt.subplots(figsize=(10,4))
        self.a2x.plot(self.wavelength,self.obj,label='OBJID ' + str(self.header['OBJID']),color='cornflowerblue')
        self.a2x.plot(self.wavelength,self.obj*self.mask,label='Mask',color='orange')
        self.a2x.set_xlabel('$\lambda$')
        self.a2x.set_ylabel('Counts')
        self.a2x.set_title('Spectrum of ' + self.image)
        self.a2x.set_ylim(-1000,max(self.obj)*1.1)
        #plt.axvspan(6865,7060, color='grey', alpha=0.35,label='Atmospheric Absorption Band')
        #plt.axvspan(7148,7340, color='grey', alpha=0.35)
        #plt.axvspan(7580,7730, color='grey', alpha=0.35)
        #plt.axvspan(8120,8410, color='grey', alpha=0.35)
        self.a2x.legend()
        plt.show()
        #Plot for poster/presentation/writeup
        #self.fig2 , self.a2x = plt.subplots(figsize=(10,4))
        #self.a2x.plot(self.wavelength,self.corr,color='black')
        #self.a2x.set_xlabel('$\lambda$')
        #self.a2x.set_ylabel('Telluric Transmission')
        #self.a2x.set_title('Tellluric Correction for ' + self.image)
        #self.a2x.set_ylim(-0.1,1.1)
        #plt.show()
        fweoj= fqwejiw
        self.fig = plt.figure(figsize=(10,8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        self.ax1 = plt.subplot(gs[0])
        self.ax2 = plt.subplot(gs[1])
        self.l1, = self.ax1.plot(self.wavelength, self.obj, label='Data', color='cornflowerblue')
        self.l2, = self.ax1.plot(self.wavelength,self.mask*self.obj, label='Mask', color='orange')
        self.l3, = self.ax1.plot(self.wavelength,self.m, label='Smooth', color='indianred',visible=False)
        self.l4, = self.ax1.plot(self.wavelength,self.diff, label='Difference', color='violet',visible=False)
        self.l5, = self.ax1.plot(self.wavelength,self.fit, label=(str(self.order) + '-Term Fit'), color='darkgreen',visible=True)
        self.ax1.set_xlabel('$\lambda$')
        self.ax1.set_ylabel('Counts')
        self.ax1.set_title('Spectrum of ' + self.image)
        self.ax1.set_ylim(min(self.diff)*1.1,max(self.obj)*1.1)
    
        self.k1, = self.ax2.plot(self.wavelength,self.corr, label='Correction', color='black')
        self.k3, = self.ax2.plot(self.wavelength,self.good,color='orange')
        self.ax2.set_xlabel('$\lambda$')
        self.ax2.set_ylabel('Transmission Coefficient')
        self.ax2.set_title('Telluric Correction')
        self.ax2.set_ylim(0,1.2)
    
        rax = plt.axes([0.105, 0.376, 0.15, 0.15])
        check = CheckButtons(rax, ('Obj', 'Masked Obj','Smooth Obj','Difference', 'Fitted Function'), (True, True, False, False, True))

        self.text = self.fig.text(0.922,0.433,'Add')
        self.ax1.legend(loc=1)

        def maskMode(event):
            self.add = not self.add
            plt.gcf().texts.remove(self.text)
            if self.add:
                self.text = self.fig.text(0.922,0.433,'Add')
            else:
                self.text = self.fig.text(0.922,0.433,'Sub')
            plt.draw()

        def saveSpec(event):
            headerout = self.header
            #[headerout.pop(key) for key in headerout.keys() if key[:3] == "ROW"]
            headerout['ROW1'] = 'Telluric Correction'
            headerout['ROW2'] = 'Spectrum'
            headerout['ROW3'] = 'Pixel Mask'
            headerout['ROW4'] = 'Polynomial Fit'
            headerout['NTERMS'] = self.order
            dataout = np.concatenate((self.corr,self.obj,self.mask,self.fit))
            dataout = dataout.reshape(4,len(self.wavelength))      
            hdu = fits.PrimaryHDU(header = headerout, data = dataout)
            filelocation = '/Users/blorenz/Desktop/COSMOSData/FitsFileOut/' + 'tc_' + self.image
            hdu.writeto(filelocation,overwrite=True)
            filelocation2 = '/Users/blorenz/Desktop/COSMOSData/FitsFileOut/ImageOut/' + 'tc_' + self.image.replace('.fits','.png')
            self.fig.savefig(filelocation2)
            print('Correction saved to ' + filelocation)

        def UP(event):
            self.order += 1
            self.l5.set_label((str(self.order) + '-Term Fit'))
            self.updatePlot()

        def UP5(event):
            self.order += 5
            self.l5.set_label((str(self.order) + '-Term Fit'))
            self.updatePlot()

        def DOWN(event):
            self.order -= 1
            if self.order == 0:
                print 'Order cannot be less than 1, setting to 1'
                self.order = 1
            self.l5.set_label((str(self.order) + '-Term Fit'))
            self.updatePlot()

        def DOWN5(event):
            self.order -= 5
            if self.order <= 0:
                print 'Order cannot be less than 1, setting to 1'
                self.order = 1
            self.l5.set_label((str(self.order) + '-Term Fit'))
            self.updatePlot()
    
        def checkButton(label):
            if label == 'Obj':
                self.l1.set_visible(not self.l1.get_visible())
            elif label == 'Masked Obj':
                self.l2.set_visible(not self.l2.get_visible())
            elif label == 'Smooth Obj':
                self.l3.set_visible(not self.l3.get_visible())
            elif label == 'Difference':
                self.l4.set_visible(not self.l4.get_visible())
            elif label == ('Fitted Function'):
                self.l5.set_visible(not self.l5.get_visible())
            self.ax1.legend(loc=1)
            plt.draw()
        check.on_clicked(checkButton)
    
        def maskslider(xmin, xmax):
            indmin, indmax = np.searchsorted(self.wavelength, (xmin, xmax))
            indmax = min(len(self.wavelength) - 1, indmax)
            if self.add == 1:
                self.mask = np.concatenate((self.mask[:indmin],np.ones(indmax-indmin),self.mask[indmax:]))
            if self.add == 0:
                self.mask = np.concatenate((self.mask[:indmin],np.zeros(indmax-indmin),self.mask[indmax:]))
            self.l2.set_ydata(self.mask*self.obj)
            self.updatePlot()
        span = SpanSelector(self.ax1, maskslider, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))

        def fixtellslider(xmin,xmax):
            indmin, indmax = np.searchsorted(self.wavelength, (xmin, xmax))
            indmax = min(len(self.wavelength) - 1, indmax)
            indmincomp = indmin+self.wavelength[0]
            indmaxcomp = indmax+self.wavelength[0]
            self.good = np.concatenate((self.good[:indmin],(indmax-indmin)*[1.0],self.good[indmax:]))
            self.k3.set_ydata(self.good)
            self.runCalc(redo=1,fixtel=1)
            self.k1.set_ydata(self.corr)
            plt.draw()
        span2 = SpanSelector(self.ax2, fixtellslider, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
        

        #update = plt.axes([0.804, 0.376, 0.085, 0.05])
        #bupdate = Button(update, 'Update Fit')
        #bupdate.on_clicked(updatePlot)
        mode = plt.axes([0.893, 0.376, 0.085, 0.05])
        bmode = Button(mode, 'Mask Mode')
        bmode.on_clicked(maskMode)
        savespec = plt.axes([0, 0, 0.05, 0.05])
        bsavespec = Button(savespec, 'Save')
        bsavespec.on_clicked(saveSpec)
        up = plt.axes([0.26, 0.401, 0.02, 0.02])
        bup = Button(up, '+')
        bup.on_clicked(UP)
        down = plt.axes([0.26, 0.376, 0.02, 0.02])
        bdown = Button(down, '-')
        bdown.on_clicked(DOWN)
        up5 = plt.axes([0.285, 0.401, 0.02, 0.02])
        bup5 = Button(up5, '+5')
        bup5.on_clicked(UP5)
        down5 = plt.axes([0.285, 0.376, 0.02, 0.02])
        bdown5 = Button(down5, '-5')
        bdown5.on_clicked(DOWN5)
        plt.tight_layout()
        plt.show()

    

Plot()
