Period_D&P is a program which calculates parameters of flux changes for different types of periodical variable stars. 
It uses photometric observations from text file, plots light curve, makes fitting by Fourier series, then calculates 
value of the period, plots phase curve and cuts all points which have too large deviation between fit and observations.

For using this code you need either install Spyder with Anaconda or just include next libraries:
scipy.optimize     #for the method of LS
numpy              #for math stuff
matplotlib.pyplot  #for plotting
time               #to know time of calculations
tkinter            #graphic interface
os      	         #to create directories
