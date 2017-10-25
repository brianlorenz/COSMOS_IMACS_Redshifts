import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.widgets import Button,SpanSelector,CheckButtons
from matplotlib import gridspec

#Usage: run FitStars.py 'imagename' (order) 
'''
imagename - string - name of your image, including file path
order - int - optional - order of the polynomial to be fit. Default is 15, and can be changed easily in the GUI

outloc - string - where your files will be output, defaults to the same as input

imname - string - automatically set to just the image name
imloc - string - automatically set to just the image location
'''
outloc = 0
order = 15

imagename = sys.argv[1]
imname = imagename
while imname.find('/') != -1:
    imname = imname[imname.find('/')+1:]
imloc = imagename.replace(imname,'')
if not outloc: outloc = imloc

class Plot():
    def __init__(self):
        '''
        Main routine that sets up the calculations and the GUI

        order - int - order of the polynomial fit
        add - boolean - 1 to add to mask or correction, 0 to subtract.
        
        
        '''
        if len(sys.argv) == 3: self.order = int(sys.argv[2])
        else: self.order = order
        self.add = 0

        hdu = fits.open(imagename)[0]
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
        '''
        Smooths the spectrum by comparing a given pixel to the 50 pixels beside it
        y - the spectrum to smooth
        p - float, 0 to 1 - the percentage to smooth to. i.e. 0.95 takes the 95% highest pixel at each step. 0.5 is just a median.
        '''
        m = np.zeros(y.shape,y.dtype)
        for j in range(len(y)):
            a,b = np.clip([j-25,j+26],0,len(y)-1)
            if p == 50:
                m[j] = np.median(y[a:b])
            else:
                m[j] = np.sort(y[a:b])[len(y[a:b])*p/100]
        return m

    def svdfit(self,b,y):
        '''
        Fits the polynomial with single value decomposition
        b - basis function
        y - spectrum
        '''
        decomp = np.linalg.svd(b,full_matrices=False)
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
        self.k3.set_ydata(self.good)
        self.ax1.legend(loc=1)
        plt.draw()

    def runCalc(self,redo=0,fixtel=0):
        '''
        Main calculations are here
        redo - boolean - set to 1 when the mask has been changed by the user
        fixtel - boolean - set to 1 when the telluric spectrum mask has been changed by the user
        obj - spectrum
        m - smoothed spectrum
        diff - difference between object and smoothed spectrum
        mask - the mask applied to to polynomial, cutting bad pixels
        fit - the polynomial fit
        corr - the telluric correction spectrum (primary output)
        '''
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

        #Atmospheric absorption lines
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
        '''
        Sets up the plot
        '''
        #Plot for poster/presentation/writeup
        #self.fig2 , self.a2x = plt.subplots(figsize=(10,4))
        #self.a2x.plot(self.wavelength,self.obj,label='OBJID ' + str(self.header['OBJID']),color='cornflowerblue')
        #self.a2x.plot(self.wavelength,self.obj*self.mask,label='Mask',color='orange')
        #self.a2x.set_xlabel('$\lambda$')
        #self.a2x.set_ylabel('Counts')
        #self.a2x.set_title('Spectrum of ' + imname)
        #self.a2x.set_ylim(-1000,max(self.obj)*1.1)
        #plt.axvspan(6865,7060, color='grey', alpha=0.35,label='Atmospheric Absorption Band')
        #plt.axvspan(7148,7340, color='grey', alpha=0.35)
        #plt.axvspan(7580,7730, color='grey', alpha=0.35)
        #plt.axvspan(8120,8410, color='grey', alpha=0.35)
        #self.a2x.legend()
        #plt.show()
        #Second Plot
        #self.fig2 , self.a2x = plt.subplots(figsize=(10,4))
        #self.a2x.plot(self.wavelength,self.corr,color='black')
        #self.a2x.set_xlabel('$\lambda$')
        #self.a2x.set_ylabel('Telluric Transmission')
        #self.a2x.set_title('Tellluric Correction for ' + imname)
        #self.a2x.set_ylim(-0.1,1.1)
        #plt.show()
        self.fig = plt.figure(figsize=(10,8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        self.ax1 = plt.subplot(gs[0])
        #Shows atmospheric absorption regions
        plt.axvspan(6865,7060, color='grey', alpha=0.08,label='Atmospheric Absorption Band')
        plt.axvspan(7148,7340, color='grey', alpha=0.08)
        plt.axvspan(7580,7730, color='grey', alpha=0.08)
        plt.axvspan(8120,8410, color='grey', alpha=0.08)
        self.ax2 = plt.subplot(gs[1])
        #Shows atmospheric absorption regions
        plt.axvspan(6865,7060, color='grey', alpha=0.08,label='Atmospheric Absorption Band')
        plt.axvspan(7148,7340, color='grey', alpha=0.08)
        plt.axvspan(7580,7730, color='grey', alpha=0.08)
        plt.axvspan(8120,8410, color='grey', alpha=0.08)
        self.l1, = self.ax1.plot(self.wavelength, self.obj, label='Data', color='cornflowerblue')
        self.l2, = self.ax1.plot(self.wavelength,self.mask*self.obj, label='Mask', color='orange')
        self.l3, = self.ax1.plot(self.wavelength,self.m, label='Smooth', color='indianred',visible=False)
        self.l4, = self.ax1.plot(self.wavelength,self.diff, label='Difference', color='violet',visible=False)
        self.l5, = self.ax1.plot(self.wavelength,self.fit, label=(str(self.order) + '-Term Fit'), color='darkgreen',visible=True)
        self.ax1.set_xlabel('$\lambda$')
        self.ax1.set_ylabel('Counts')
        self.ax1.set_title('Spectrum of ' + imname)
        self.ax1.set_ylim(min(self.diff)*1.1,max(self.obj)*1.1)
    
        self.k1, = self.ax2.plot(self.wavelength,self.corr, label='Correction', color='black')
        self.k3, = self.ax2.plot(self.wavelength,self.good,color='orange')
        self.ax2.set_xlabel('$\lambda$')
        self.ax2.set_ylabel('Transmission Coefficient')
        self.ax2.set_title('Telluric Correction')
        self.ax2.set_ylim(0,1.2)
    
        rax = plt.axes([0.105, 0.376, 0.15, 0.15])
        check = CheckButtons(rax, ('Obj', 'Masked Obj','Smooth Obj','Difference', 'Fitted Function'), (True, True, False, False, True))

        self.ax1.legend(loc=1)

        print 'Drag over the top plot to add or subtract regions of the mask, toggling your mode with the button on the right. Drag over the bottom plot to add or remove regions from the telluric correction file once you are happy witht he top plot. You can change the order of the polynomial with the +,-,+5,-5 buttons. Click Save when you are finished or close the window to end without saving.'

#=================================================
#Buttons and Sliders
#=================================================

#Toggle whether or not to add or subtract
        def maskMode(event):
            self.add = not self.add
            if self.add:
                bmode.label.set_text('Mode: Add')
            else:
                bmode.label.set_text('Mode: Sub')
            plt.draw()

#Save the spectrum when finished
        def saveSpec(event):
            headerout = self.header
            headerout['ROW1'] = 'Telluric Correction'
            headerout['ROW2'] = 'Spectrum'
            headerout['ROW3'] = 'Pixel Mask'
            headerout['ROW4'] = 'Polynomial Fit'
            headerout['NTERMS'] = self.order
            dataout = np.concatenate((self.corr,self.obj,self.mask,self.fit))
            dataout = dataout.reshape(4,len(self.wavelength))      
            hdu = fits.PrimaryHDU(header = headerout, data = dataout)
            filelocation = outloc + 'tc_' + imname
            try: hdu.writeto(filelocation,overwrite=True)
            except TypeError: hdu.writeto(filelocation,clobber=True)
            filelocation2 = outloc + 'tc_' + imname.replace('.fits','.png')
            self.fig.savefig(filelocation2)
            print('Correction saved to ' + filelocation)
            plt.close(self.fig)

#Reset the plots if you messed something up
        def bReset(event):
            self.runCalc()
            self.l2.set_ydata(self.mask*self.obj)
            self.updatePlot()
            
#Buttons for changing the order
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

#Checks to display different parts of the star spectrum
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

#Slider that adds or subtracts from the mask
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

#Slider that adds or subtractions from the telluric correction
        def fixtellslider(xmin,xmax):
            indmin, indmax = np.searchsorted(self.wavelength, (xmin, xmax))
            indmax = min(len(self.wavelength) - 1, indmax)
            indmincomp = indmin+self.wavelength[0]
            indmaxcomp = indmax+self.wavelength[0]
            if not self.add: self.good = np.concatenate((self.good[:indmin],(indmax-indmin)*[1.0],self.good[indmax:]))
            else: self.good = np.concatenate((self.good[:indmin],(indmax-indmin)*[0.0],self.good[indmax:]))
            self.runCalc(redo=1,fixtel=1)
            self.k1.set_ydata(self.corr)
            self.k3.set_ydata(self.good)
            plt.draw()
        span2 = SpanSelector(self.ax2, fixtellslider, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
        

#Creation of the buttons
        mode = plt.axes([0.878, 0.376, 0.10, 0.05])
        bmode = Button(mode, 'Mode: Sub')
        bmode.on_clicked(maskMode)
        savespec = plt.axes([0, 0, 0.1, 0.05])
        bsavespec = Button(savespec, 'Save')
        bsavespec.on_clicked(saveSpec)
        reset = plt.axes([0.5, 0.01, 0.1, 0.03])
        breset = Button(reset, 'Reset')
        breset.on_clicked(bReset)
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

#=================================================

#Run the main routine
Plot()
