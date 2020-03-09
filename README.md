Please use periodDP_V1.0.0.py as last and stable version

Period_D&P is a program which calculates parameters of flux changes for different types of periodical variable stars. 
It uses photometric observations from text file, plots light curve, makes fitting by Fourier series, then calculates 
value of the period, plots phase curve and cuts all points which have too large deviation between fit and observations.

For using this code you need either install Spyder with Anaconda or just download next libraries:scipy.optimize, numpy   , matplotlib.pyplot, time, tkinter, os   

In linux, all you need - to run the .sh file (these commands are included  there)

Authors:

  Pavlo Kashko (kashko.pavlo@gmail.com) - main author;
  
  Dmytro Tvardovskyi (dmytro.tvardovskyi@gmail.com) - graphic interface
  
  Viktor Khalack - scientific advisor
  
  ©Kashko2019 (see LICENSE)
  
  DOI: 10.5281/zenodo.3257785
 
If use, please cite:
  @misc{kashko2020revealing,
      title={Revealing the nature of HD63401},
      author={Pavlo Kashko and Viktor Khalack and Oleksandr Kobzar and Dmytro Tvardovskyi and Mathieu Perron-Cormier},
      year={2020},
      eprint={2003.02925},
      archivePrefix={arXiv},
      primaryClass={astro-ph.SR}
  }
