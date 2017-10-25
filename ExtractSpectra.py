import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from astropy.io import fits
from scipy.optimize import curve_fit
from matplotlib.widgets import Slider, Button, SpanSelector


#Usage: run ExtractSpectra.py imloc (objnum)
#Example: run ExtractSpectra.py "feb16_abig.fits"
    #This extracts all objects from the mask
#Example: run ExtractSpectra.py "feb16_abig.fits" 45
    #This pulls up only object 45, then you can move through the rest from there. 
'''
Extract a .FITS file for each object with spectrum, noise, flat, and sky

Takes big.fits, bigsig.fits, bigext.fits, bigsky.fits files and outputs a 
.fits file for each object, containing one row each of spectrum, noise, flat, and sky. 
objnum - the objects' number in the mask (see header)
imloc - location of the big.fits file. Input to the code
outloc - location to output all of the extracted .fits files
showbad - set to 1 if you want to immediately fix bad fits by hand, 0 if you do not
'''
class FunPlot:
    def __init__(self,image):
        self.image = image
        while image.find('/') != -1:
            image = image[image.find('/')+1:]
        self.imname = image
        self.readFits()
    #Compute the gaussians to be used later to determine which rows to plot
        self.gaussians = np.add.reduce(self.hdu.data,1)
        self.gaussians = np.nan_to_num(self.gaussians)
    #Variable for plotting the gaussian 
        self.redo = 0
        self.toggle = 0
        self.outfile = outloc + 'out_' + self.imname.replace('.fits','.txt')
    #Uncomment these to display the gaussians
        #fig, ax = plt.subplots()
        #ax.plot(np.arange(len(self.gaussians)),self.gaussians)
        #plt.show()

    #Plot for presentation, not necessary for code
    def gaussPlot(self,x1,x2):
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_axes([0.1, 0.1, .85, 0.35])
        ax.plot(np.arange(x1,x2),self.gaussians[x1:x2],color='cornflowerblue',label='Data')
        amin,amax = ax.get_ylim()
        ax.plot((self.csecta[21],self.csecta[21]),(-1000000,1000000),color='grey',label='Slit Division')
        ax.plot((self.csecta[22],self.csecta[22]),(-1000000,1000000),color='grey')
        ax.plot((self.csecta[23],self.csecta[23]),(-1000000,1000000),color='grey')
        ax.plot((self.csecta[24],self.csecta[24]),(-1000000,1000000),color='grey')
        ax.plot((self.csecta[25],self.csecta[25]),(-1000000,1000000),color='grey')
        #ax.plot((self.csecta[26],self.csecta[26]),(-1000000,1000000),color='grey')
        ax.set_title('Data Compressed Along Row Axis')
        ax.set_xlabel('Row Number')
        ax.set_ylabel('Counts')
        ax.set_ylim(amin,amax)
        ax.set_xlim(self.csecta[21]-0.25,self.csecta[25]+0.25)
        ax.set_xticks([self.csecta[21],self.csecta[22],self.csecta[23],self.csecta[24],self.csecta[25]])
        ax.legend(bbox_to_anchor=(0.5, 1.05))
        fig2 = plt.figure(figsize=(10,8))
        ax2 = fig2.add_axes([0.1, 0.1, .85, 0.35])
        image_data = self.hdu.data[self.csecta[21]:self.csecta[25]+1]
        ax2.imshow(image_data, cmap='gray',clim=(-240.0, 240.0),aspect='auto')
        sub = self.csecta[21]
        ax2.set_yticks([self.csecta[21]-sub,self.csecta[22]-sub,self.csecta[23]-sub,self.csecta[24]-sub,self.csecta[25]-sub])
        ax2.set_yticklabels((self.csecta[21],self.csecta[22],self.csecta[23],self.csecta[24],self.csecta[25]))
        ax2.set_xticklabels((0,4900,4900+500*2,4900+1000*2,4900+1500*2,4900+2000*2,4900+2500*2,4900+3000*2))
        ax2.set_xlabel('$\lambda$ ($\AA$)')
        ax2.set_ylabel('Row Number')
        plt.show()
    

        
#Function to read the .fits file and compute header information, done during initialization.
    def readFits(self):
        fitsfile = fits.open(self.image)
        self.hdu = fitsfile[0]
        fitsfile2 = fits.open(self.image.replace(".fits","sig.fits"))
        self.hdu2 = fitsfile2[0]
        fitsfile3 = fits.open(self.image.replace(".fits","ext.fits"))
        self.hdu3 = fitsfile3[0]
        fitsfile4 = fits.open(self.image.replace(".fits","sky.fits"))
        self.hdu4 = fitsfile4[0]
        #Compute values for use later
        self.nslits = self.hdu.header['nslits']
        self.csecta = np.array([self.hdu.header["csect%da"%(j+1)] for j in range(self.nslits)])
        self.csectb = np.array([self.hdu.header["csect%db"%(j+1)] for j in range(self.nslits)])
        print('Read ' + self.image)

        
#Set up Gaussian Function
    def gauss2(self, x, mu, sigma, y0):
        g = self.amp2(x,mu,sigma,y0)*np.exp(-(x-mu)**2/(2.*sigma**2))+y0
        return g

    def amp2(self, x, mu, sigma,y0):
        g = np.exp(-(x-mu)**2/(2.*sigma**2))
        t = np.add.reduce((self.ydata-y0)*g)
        b = np.add.reduce(g*g)
        if b: A = t/b
        else: A = 0
        return A

    #Funciton to compute the gaussians
    def getRows(self):
        datapeak = np.nanmax(self.gaussians[self.csecta[self.objnum-1]+8:self.csectb[self.objnum-1]-8])
        peakindex = np.argwhere(self.gaussians == datapeak)
        threshold = datapeak*0.2
        try:
            start = next(i[0] for i in reversed(list(enumerate(list(self.gaussians)))) if i[1] < threshold and i[0] < peakindex)
            end = next(i[0] for i in enumerate(list(self.gaussians)) if i[1] < threshold and i[0] > peakindex)
            if self.redo:
                pass
            else:
                self.ydata = self.gaussians[start-5:end+5]
                self.xstart = start-5
                self.xdata = np.arange(len(self.ydata))
            p0 = (3,len(self.xdata)/2,10000)
            coeff2, var_matrix2 = curve_fit(self.gauss2, self.xdata, self.ydata,p0=p0)
            self.gausscurve2 = self.gauss2(self.xdata,coeff2[0],coeff2[1],coeff2[2])
            print coeff2
            amp2 = self.amp2(self.xdata,coeff2[0],coeff2[1],coeff2[2])
            self.mu2 = coeff2[0] 
            self.stddev2 = np.abs(coeff2[1])
            self.flag = 0
        #As for as I know, all of these errors are due to the fitting of Gaussians - when it cannot fit, an error returns and it will just leave the stddev at 0.
        except (StopIteration,ValueError,TypeError,RuntimeError):
            self.stddev2 = 1
            if len(peakindex) != 1:
                self.mu2 = 0
            else:
                self.mu2 = 0
            self.gausscurve2 = np.zeros(len(self.xdata))
            print("Could not compute Gaussian Fit")
            self.flag = 1


    
#Plotting function - should be run after findRows(). Give this an object number, and it will plot an appropriate amount of rows 
#determined by the standard deviations of the gaussian fits.
    def plotObj(self,objnum,close=0):
        self.redo = 0
        self.close = close
        self.instruct = 1
        self.objnum=objnum

        #The function called whenever the figure needs to be updated. This will set the xxdata, ydata, scales, and titles for each of the plots
        def updatePlot(addrow=0):
            '''
            leftbound and rightbound - this is the region plotted that will be output as a .fits file
            flag - 0 if fit was good, 1 if fit failed
            xstart - beginning of the region that the Gaussian fit over
            mu2 - mean of the Gaussian fit measured in pixels from xstart.
            stddev2 - the standard deviation
            numrows - the determined number of rows to sum over, 1.4*stddev2
            gausscurve2 - the y_data of the Gaussian Fit
            '''
            self.image_data = self.hdu.data[self.csecta[self.objnum-1]:self.csectb[self.objnum-1]+1]
            self.fitsplot.set_data(self.image_data)
            if addrow:
                pass
            else:
                if self.flag == 0:
                    self.leftbound = self.xstart+self.mu2-self.numrows
                    self.rightbound = self.xstart+self.mu2+self.numrows+1
                else:
                    self.leftbound = 0
                    self.rightbound = 0
            self.l1.set_ydata(np.add.reduce(self.hdu.data[self.leftbound:self.rightbound]))
            self.l2.set_ydata(np.sqrt(np.add.reduce(self.hdu2.data[self.leftbound:self.rightbound]**2)))
            self.ax3.set_title('Gaussian Fit, Rows ' + str(self.leftbound) + ' to ' + str(self.rightbound-1) + ', ' + self.imname)
            #self.zoomplot.set_xdata(self.wavelength)
            #self.zoomplot2.set_xdata(self.wavelength)
            #self.zoomplot.set_ydata(np.add.reduce(self.hdu.data[self.leftbound:self.rightbound]))
            #self.zoomplot2.set_ydata(np.sqrt(np.add.reduce(self.hdu2.data[self.leftbound:self.rightbound]**2)))
            self.k3.set_xdata((self.leftbound,self.leftbound))
            self.k4.set_xdata((self.rightbound-1,self.rightbound-1))
            if self.redo:
                self.xstart = self.csecta[self.objnum-1]
            else:
                self.k2.set_xdata(self.xdata+self.xstart)
                self.k2.set_ydata(self.gausscurve2)
                self.x3data = self.xdata
                self.x3start = self.xstart
                self.y3data = self.ydata
            self.k1.set_xdata(self.xdata+self.xstart)
            self.k1.set_ydata(self.ydata)
            self.ax.set_xlim(self.wavelength[0],self.wavelength[-1])
            self.ax.set_ylim(min(self.dataplot),max(self.dataplot))
            self.h1.set_ydata((self.leftbound-self.csecta[self.objnum-1],self.leftbound-self.csecta[self.objnum-1]))
            self.h2.set_ydata((self.rightbound-1-self.csecta[self.objnum-1],self.rightbound-1-self.csecta[self.objnum-1]))
            self.ax4.set_yticks((0,25))
            self.ax4.set_yticklabels([self.csecta[self.objnum-1],self.csecta[self.objnum-1]+25])
            if addrow:
                pass
            else:
                #self.ax2.set_xlim(self.wavelength[0],self.wavelength[-1])
                #self.ax2.set_ylim(min(self.dataplot),max(self.dataplot))
                try:
                    self.ax3.set_xlim(self.xdata[0]+self.xstart,self.xdata[-1]+self.xstart)
                    self.ax3.set_ylim(min(min(self.ydata),min(self.gausscurve2)),max(max(self.ydata),max(self.gausscurve2))*1.1)
                except IndexError:
                    print self.xstart
                    self.ax3.set_xlim(self.xstart,self.xstart + 30)
                    self.ax3.set_ylim(0,10000)
                except ValueError:
                    print self.xstart
                    self.ax3.set_xlim(self.xstart,self.xstart + 30)
                    self.ax3.set_ylim(0,10000)
            self.ax.legend([self.l1,self.l2],[('OBJID ' + str(self.objids[self.objnum-1])),"Noise"],bbox_to_anchor=(1.22, 0.22))
            self.ax3.legend(bbox_to_anchor=(1.22, 1.0))
            plt.draw()

        #The calculates the gaussian fit for the given object number, and creates the plots if this is the first time. 
        #This does not need ot be called everytime a plot needs to be updated (e.g. zoomzlider) but is called whenever a 
        #gaussian needs to be refit (e.g. zoomslider3 or moving to the next objnum). 
        def calcPlot(objnum,createplot=1):
            '''
            ax1 - contains the plot of the data and noise (noise is added in quadriture)
            ax3 - contains the Gaussian plot and the zoomslider to refit the Gaussian
            ax4 - contains the display of the relevent section of the big.fits file.
            '''
            self.getRows()
            self.redo = 0
            self.wavelength = (1.+np.arange(self.hdu.header["naxis1"])-self.hdu.header["crpix1"])*self.hdu.header["cdelt1"] + self.hdu.header["crval1"]
            self.objids = [self.hdu.header["objid%d"%(j+1)]
                           .replace('ob','').replace('ref','') 
                           for j in range(self.nslits)]
            self.numrows = int(round(1.4*self.stddev2))
            self.mu2 = int(round(self.mu2))
            self.dataplot = np.add.reduce(self.hdu.data[self.xstart+self.mu2-self.numrows:self.xstart+self.mu2+self.numrows+1])
            self.dataplot2 = np.add.reduce(self.hdu2.data[self.xstart+self.mu2-self.numrows:self.xstart+self.mu2+self.numrows+1]**2)
            self.dataplot2 = np.sqrt(self.dataplot2)
            if createplot == 1:
                self.fig = plt.figure(figsize=(10,8))
                #self.ax = plt.subplot2grid((2,2),(0,0),colspan=2)
                #self.ax2 = plt.subplot2grid((2,2),(1,0))
                #self.ax3 = plt.subplot2grid((2,2),(1,0),colspan=2)
                self.ax = self.fig.add_axes([0.12, 0.1, 0.7, 0.3])
                self.ax3 = self.fig.add_axes([0.12, 0.6, 0.7, 0.3])
                self.ax4 = self.fig.add_axes([0.12, 0.4, 0.7, 0.1])
                #self.ax4.xaxis.tick_top()
                self.ax4.get_xaxis().set_visible(False)
                self.ax4.set_ylabel('Row Number')
                self.image_data = self.hdu.data[self.csecta[self.objnum-1]:self.csectb[self.objnum-1]+1]
                self.fitsplot = self.ax4.imshow(self.image_data, cmap='gray',clim=(-240.0, 240.0),aspect='auto')
                self.h1, = self.ax4.plot((0,len(self.wavelength)),(0,0),color='mediumseagreen')
                self.h2, = self.ax4.plot((0,len(self.wavelength)),(0,0),color='mediumseagreen')
                self.ax4.set_xlim(0,len(self.wavelength))
                self.l1, = self.ax.plot(self.wavelength,self.dataplot,label=('OBJID ' + str(self.objids[self.objnum-1])), color='cornflowerblue')
                self.l2, = self.ax.plot(self.wavelength,self.dataplot2, label='Noise', color='orange')
                #self.zoomplot, = self.ax2.plot(self.wavelength,self.dataplot,label=('OBJID ' + str(self.objids[self.objnum-1])), color='cornflowerblue')
                #self.zoomplot2, = self.ax2.plot(self.wavelength,self.dataplot2,label='Noise', color='orange')
                self.k1, = self.ax3.plot(self.xdata,self.ydata, color = 'cornflowerblue',label='Data')
                self.k2, = self.ax3.plot(self.xdata,self.gausscurve2, color = 'red',label='Gaussian Fit')
                self.k3, = self.ax3.plot((self.mu2-self.numrows,self.mu2-self.numrows),(-1000000,1000000),color='mediumseagreen',label='Plotted Rows')
                self.k4, = self.ax3.plot((self.mu2+self.numrows+1,self.mu2+self.numrows+1),(-1000000,1000000),color='mediumseagreen')
                self.ax.set_xlabel('$\lambda$')
                self.ax.set_ylabel('Counts')
                #self.ax2.set_xlabel('$\lambda$')
                #self.ax2.set_ylabel('Counts')
                self.ax3.set_xlabel('Pixel Index')
                self.ax3.set_ylabel('Counts')
                #self.ax2.set_title('Zoom')
                #self.ax3.set_title('Gaussian Fit')
                if self.instruct:
                    print "Drag over the Gaussian plot to attempt to refit a Guassian over the given region."
                    print "If the fitting isn't working, click on the Gaussian plot to plot a single row, then use the arrow buttons to adjust the range."
                    print "When finished, press savespec before moving on."
                    self.instruct = 0
            updatePlot()
            redoGauss()


        def redoGauss():
            self.xdata = np.arange(self.csectb[self.objnum-1]-self.csecta[self.objnum-1])
            self.ydata = self.gaussians[self.csecta[self.objnum-1]:self.csectb[self.objnum-1]]
            #self.k1.set_xdata(self.xdata)
            #self.k1.set_ydata(self.ydata)
            self.ax3.set_xlim(self.xdata[0],self.xdata[-1])
            self.ax3.set_ylim(min(self.ydata),max(self.ydata))
            self.redo = 1
            updatePlot()

        #Funcitons for the Buttons
        def next(event):
            self.objnum += 1
            self.redo = 0
            calcPlot(self.objnum,createplot=0)
        
        def prev(event):
            self.objnum -= 1
            self.redo = 0
            calcPlot(self.objnum,createplot=0)

        def addRowRight(event):
            self.rightbound += 1
            updatePlot(addrow=1)

        def rmRowRight(event):
            self.rightbound -= 1
            updatePlot(addrow=1)

        def addRowLeft(event):
            self.leftbound += 1
            updatePlot(addrow=1)

        def rmRowLeft(event):
            self.leftbound -= 1
            updatePlot(addrow=1)

        def saveSpectrum(event):
            hduarray = np.add.reduce(self.hdu.data[self.leftbound:self.rightbound])
            hdu2array = np.add.reduce(self.hdu2.data[self.leftbound:self.rightbound])
            hdu3array = np.add.reduce(self.hdu3.data[self.leftbound:self.rightbound])
            hdu4array = np.add.reduce(self.hdu4.data[self.leftbound:self.rightbound]) 
            dataout = np.concatenate((hduarray,hdu2array,hdu3array,hdu4array))
            dataout = dataout.reshape(4,len(self.wavelength))
            headerout = self.hdu.header[:156]
            headerout['ROW1'] = 'big'
            headerout['ROW2'] = 'bigsig'
            headerout['ROW3'] = 'bigext'
            headerout['ROW4'] = 'bigsky'
            headerout['OBJNUM'] = self.objnum
            headerout['OBJID'] = int(self.objids[self.objnum-1])
            headerout['CSECTA'] = self.csecta[self.objnum-1]
            headerout['CSECTB'] = self.csectb[self.objnum-1]
            headerout['SUMROWS'] = self.rightbound-self.leftbound
            headerout['LOWBOUND'] = self.leftbound
            headerout['UPBOUND'] = self.rightbound-1
            hdu = fits.PrimaryHDU(header = headerout, data = dataout)
            filelocation2 = outloc + ('%06d' % int(self.objids[self.objnum-1])) + '_' + self.imname
            self.fig.savefig(outloc + ('%06d' % int(self.objids[self.objnum-1])) + '_' + self.imname.replace('.fits','.png'))
            hdu.writeto(filelocation2,clobber=True,output_verify='ignore')#overwrite=True)
            print('Spectrum saved to ' + filelocation2)
            f3 = open(self.outfile,"r+")
            d = f3.readlines()
            f3.seek(0)
            f3.write(d[0])
            for i in d[1:]:
                if self.objnum == int(i[7:10]):
                    f3.write('%06d %3d    %4d %4d   %4d     %4d    %d\n' % (int(self.objids[self.objnum-1]),self.objnum,
                                                                            self.mu2+self.xstart,self.stddev2,
                                                                            self.leftbound,self.rightbound-1,self.flag))
                else:
                    f3.write(i)
            f3.close()
            savetext = self.fig.text(0.85,0.3,'Saved!',fontsize=24)
            plt.pause(0.1)
            plt.gcf().texts.remove(savetext)
            if self.close:
                plt.close(self.fig)

        def Refresh(event):
            plt.close(self.fig)
            if self.close:
                self.plotObj(self.objnum,close=1)
            else:
                self.plotObj(self.objnum)

        def zoomGauss(event):
            self.toggle = not self.toggle
            if self.toggle:
                self.ax3.set_xlim(self.x3data[0]+self.x3start,self.x3data[-1]+self.x3start)
                self.ax3.set_ylim(min(min(self.y3data),min(self.gausscurve2)),max(max(self.y3data),max(self.gausscurve2))*1.1)
            else:
                self.ax3.set_xlim(self.xdata[0]+self.xstart,self.xdata[-1]+self.xstart)
                self.ax3.set_ylim(min(min(self.ydata),min(self.gausscurve2)),max(max(self.ydata),max(self.gausscurve2))*1.1)
            plt.draw()
            


  #Functions for zooming in on the first plot and refitting on the gauss                  
        #def zoomslider(xmin, xmax):
            #indmin, indmax = np.searchsorted(self.wavelength, (xmin, xmax))
            #indmax = min(len(self.wavelength) - 1, indmax)
            #thisx = self.wavelength[indmin:indmax]
            #thisy = self.dataplot[indmin:indmax]
            #thisy2 = self.dataplot2[indmin:indmax]
            #if indmin==indmax:
                #thisx = self.wavelength
                #thisy = self.dataplot
                #thisy2 = self.dataplot2
            #self.zoomplot.set_data(thisx, thisy)
            #self.zoomplot2.set_data(thisx, thisy2)
            
            
        def zoomslider3(xmin, xmax):
            indmin, indmax = np.searchsorted(self.xdata+self.xstart, (xmin, xmax))
            indmax = min(len(self.xdata) - 1, indmax)
            if indmin == indmax:
                self.leftbound = indmin+self.csecta[self.objnum-1]
                self.rightbound = indmin+self.csecta[self.objnum-1]+1
                updatePlot(addrow=1)
            else:
                self.xstart = self.xdata[indmin]+self.xstart
                self.xdata = self.xdata[indmin:indmax]-self.xdata[indmin]
                self.ydata = self.ydata[indmin:indmax]
                self.redo = 1
                calcPlot(self.objnum,createplot=0)
            #fig.canvas.draw()

        #Immediately run calcplot after calling plotObj()
        calcPlot(self.objnum)
        
        #Define the button locations
        if not self.close:
            axnext = plt.axes([0.85, 0.575, 0.10, 0.03])
            axprev = plt.axes([0.85, 0.535, 0.10, 0.03])
            bnext = Button(axnext, 'Next OBJID')
            bprev = Button(axprev, 'Prev OBJID')
            bnext.on_clicked(next)
            bprev.on_clicked(prev)
        zoomgauss = plt.axes([0.72, 0.535, 0.1, 0.03])
        savespec = plt.axes([0.85, 0.40, 0.1, 0.1])
        refresh = plt.axes([0.03, 0.535, 0.075, 0.03])
        addrowright = plt.axes([0.58, 0.535, 0.03, 0.03])
        rmrowright = plt.axes([0.545, 0.535, 0.03, 0.03])
        addrowleft = plt.axes([0.355, 0.535, 0.03, 0.03])
        rmrowleft = plt.axes([0.32, 0.535, 0.03, 0.03])
        #Turn them into buttons
        bzoomgauss = Button(zoomgauss, 'Zoom Gauss')
        bsavespec = Button(savespec, 'Save Spec')
        brefresh = Button(refresh, 'Refresh')
        baddrowright = Button(addrowright, '>')
        brmrowright = Button(rmrowright, '<')
        baddrowleft = Button(addrowleft, '>')
        brmrowleft = Button(rmrowleft, '<')
        #Point to the relevent function when clicked
        bzoomgauss.on_clicked(zoomGauss)
        bsavespec.on_clicked(saveSpectrum)
        brefresh.on_clicked(Refresh)
        baddrowright.on_clicked(addRowRight)
        brmrowright.on_clicked(rmRowRight)
        baddrowleft.on_clicked(addRowLeft)
        brmrowleft.on_clicked(rmRowLeft)
        #Define the sliders
        #span = SpanSelector(self.ax, zoomslider, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
        span3 = SpanSelector(self.ax3, zoomslider3, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
        #Display the plots       
        #plt.tight_layout(pad=5.0)
        plt.show()


    def outputData(self):
        new = 1
        for i in range (0,self.nslits):
        #for i in range(110,115):
            self.objnum = i+1
            #Compute the Gaussian Fit
            self.getRows()
            #Setup the variables
            numrows = int(round(1.4*self.stddev2))
            mu2 = int(round(self.mu2))
            leftbound = self.xstart+mu2-numrows
            rightbound = self.xstart+mu2+numrows+1
            wavelength = (1.+np.arange(self.hdu.header["naxis1"])-self.hdu.header["crpix1"])*self.hdu.header["cdelt1"] + self.hdu.header["crval1"]
            objids = [self.hdu.header["objid%d"%(j+1)]
                      .replace('ob','').replace('ref','') 
                      for j in range(self.nslits)]
            #Plot then save the plot
            #fig, axarr = plt.subplots(3,1,figsize=(10,8))
            #ax1 = axarr[0]
            #ax2 = axarr[1]
            #ax4 = axarr[2]
            fig = plt.figure(figsize=(10,8))
            ax1 = fig.add_axes([0.12, 0.1, 0.7, 0.3])
            ax2 = fig.add_axes([0.12, 0.6, 0.7, 0.3])
            ax4 = fig.add_axes([0.12, 0.4, 0.7, 0.1])
            dataplot = np.add.reduce(self.hdu.data[leftbound:rightbound])
            l1, = ax1.plot(wavelength,dataplot,label=('OBJID ' + str(objids[self.objnum-1])), color='cornflowerblue')
            k2, = ax2.plot(self.xdata+self.xstart,self.gausscurve2,color = 'red',label='Gaussian Fit')
            xdata2 = np.arange(self.csectb[self.objnum-1]-self.csecta[self.objnum-1])+self.csecta[self.objnum-1]
            ydata2 = self.gaussians[self.csecta[self.objnum-1]:self.csectb[self.objnum-1]]
            k1, = ax2.plot(xdata2,ydata2,color = 'cornflowerblue',label='Data')
            k3, = ax2.plot((leftbound,leftbound),(-1000000,1000000),color='mediumseagreen',label='Plotted Rows')
            k4, = ax2.plot((rightbound-1,rightbound-1),(-1000000,1000000),color='mediumseagreen')
            ax2.set_xlim(xdata2[0],xdata2[-1])
            try:
                ax2.set_ylim(max(min(min(ydata2),min(self.gausscurve2)),-100000),max(max(ydata2),max(self.gausscurve2))*1.1)
            except ValueError: ax2.set_ylim(-10000,10000)
            ax1.set_xlabel('$\lambda$')
            ax1.set_ylabel('Counts')
            ax2.set_xlabel('Pixel Index')
            ax2.set_ylabel('Counts')
            ax1.set_title('Spectrum of objnum ' + str(self.objnum) + ' in ' + self.imname)
            ax2.set_title('Gaussian Fit, Plotted Rows ' + str(leftbound) + ' to ' + str(rightbound-1))
            ax1.legend()
            ax2.legend()
            ax4.get_xaxis().set_visible(False)
            ax4.set_ylabel('Row Number')
            image_data = self.hdu.data[self.csecta[self.objnum-1]:self.csectb[self.objnum-1]+1]
            fitsplot = ax4.imshow(image_data, cmap='gray',clim=(-240.0, 240.0),aspect='auto')
            h1, = ax4.plot((0,len(wavelength)),
                           (leftbound-self.csecta[self.objnum-1],
                            leftbound-self.csecta[self.objnum-1]),
                           color='mediumseagreen')
            h2, = ax4.plot((0,len(wavelength)),
                           (rightbound-1-self.csecta[self.objnum-1],
                            rightbound-1-self.csecta[self.objnum-1]),
                           color='mediumseagreen')
            ax4.set_xlim(0,len(wavelength))
            ax4.set_yticks((0,25))
            ax4.set_yticklabels([self.csecta[self.objnum-1],self.csecta[self.objnum-1]+25])
            #plt.tight_layout()
            print objids[self.objnum-1]
            fig.savefig(outloc + ('%06d' % int(objids[self.objnum-1])) + '_' + self.imname.replace('.fits','.png'))
            #Write the .fits file
            hduarray = np.add.reduce(self.hdu.data[leftbound:rightbound])
            hdu2array = np.sqrt(np.add.reduce(self.hdu2.data[leftbound:rightbound]**2))
            hdu3array = np.add.reduce(self.hdu3.data[leftbound:rightbound])
            hdu4array = np.add.reduce(self.hdu4.data[leftbound:rightbound]) 
            dataout = np.concatenate((hduarray,hdu2array,hdu3array,hdu4array))
            dataout = dataout.reshape(4,len(wavelength))
            headerout = self.hdu.header[:156]
            headerout['ROW1'] = 'big'
            headerout['ROW2'] = 'bigsig'
            headerout['ROW3'] = 'bigext'
            headerout['ROW4'] = 'bigsky'
            headerout['OBJNUM'] = self.objnum
            headerout['OBJID'] = int(objids[self.objnum-1])
            headerout['CSECTA'] = self.csecta[self.objnum-1]
            headerout['CSECTB'] = self.csectb[self.objnum-1]
            headerout['SUMROWS'] = rightbound-leftbound
            headerout['LOWBOUND'] = leftbound
            headerout['UPBOUND'] = rightbound-1
            hdu = fits.PrimaryHDU(header = headerout, data = dataout)
            filelocation2 = outloc + ('%06d' % int(objids[self.objnum-1])) + '_' + self.imname
            hdu.writeto(filelocation2,clobber=True,output_verify='ignore')#overwrite=True)
            print 'Spectrum saved to ' + filelocation2 
            #Add the info to the text file output
            if new:
                f = open(self.outfile,'w+')
                f.write('#Objid Objnum Mean Stddev Lowbound Upbound Flag\n')
                new = 0;
            else:
                f = open(self.outfile,'a+')
            f.write('%06d %3d    %4d %4d   %4d     %4d    %d\n' % (int(objids[self.objnum-1]),self.objnum,
                                                                   mu2+self.xstart,self.stddev2,leftbound,
                                                                   rightbound-1,self.flag))
            f.close()
            plt.close(fig)

def objid_to_objnum(datafile,objid):
    '''Take an objid and return the objnum stored in datafile'''
    try:
        with open(datafile,'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                oid,onum = line.split()[:2]
                if int(oid) == int(objid):
                    return int(onum)
        return KeyError("No such objid in database")
    except IOError:
        print "You cannot send a six-digit objid if you don't"
        print "have a datafile (e.g., out_maskbig.txt) to query"
        print "for the objnum."
        return IOError

if __name__ == '__main__':

    try: imloc = sys.argv[1]
    except: sys.exit('Usage: run ExtractSpectra.py imagelocation (objnum)')
    outloc = './'
    showbad = 1 #set to 1 if you want to immediately fix bad fits by hand, 0 if you do not
    
    a = FunPlot(imloc)
    
    if len(sys.argv) == 3:
        if len(sys.argv[2]) == 6:
            #User typed in 6-digit ID, not slit index
            objnum =  objid_to_objnum(a.outfile,sys.argv[2])                   
            a.plotObj(objnum)
        else:
            a.plotObj(int(sys.argv[2]))
    elif len(sys.argv) == 2:
        yesno = raw_input('Are you sure you want to overwrite output from ' + sys.argv[1] + '? (y/n) ').lower()
        if yesno == 'y':
            a.outputData()
            if showbad:
                outdata = np.genfromtxt(a.outfile,dtype=None,names=True)
                outdata = outdata[outdata['Flag'] == 1]
                for i in range(len(outdata)):
                    a.plotObj(outdata[i]['Objnum'],close=1)            
            

