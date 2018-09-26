Read Me for Lab 0 and my many files

Directories beginning with 3.1 through 3.3 contain data from the data taking
  portion of the lab.

DataAnalysis4p1-m10.py   
  Contains the data analysis script required for section 4.1 when the CCD was set to
  -10 degrees C
DataAnalysis4p1-p10.py   
  Contains the data analysis script required for section 4.1 when the CCD was set to 
  +10 degrees C
DataAnalysis4p2-m10.py   
  Contains the data analysis script required for section 4.2 when the CCD was set to
  -10 degrees C
DataAnalysis4p2-p10.py   
  Contains the data analysis script required for section 4.2 when the CCD was set to
  +10 degrees C
DataAnalysis4p3.py       
  Contains the data analysis script required for section 4.3 (The CCD was set to -5
  degrees C)
  Note: It takes a few minutes to run.
DataAnalysis4p3.py       
  Contains the data analysis script required for section 4.3 - Flat fields 
DataAnalysis4p4.py       
  Contains the data analysis script required for section 4.4 - Bad Pixel Maps
DataAnalysis4p5-crop.py       
  Contains the data analysis script required for section 4.5 - cropping 50 micrometer
  line from the master spectrum flat.
DataAnalysis4p5-lamp.py       
  Contains the data analysis script required for section 4.5 - plotting average
  pixel value for every column in the lamp spectrum.
DataAnalysis4p5.py       
  Contains the data analysis script required for section 4.5 - analyzing and 
  normalizing data from 50 micrometer master spectrum flat then applying it
  to the lamp spectrum

sensitivity_data_0deg.txt
  Tracks the x, y, and brightness values of flat_master.fits from the center to the
  upper right corner of the fits file
sensitivity_data_90deg.txt
  Tracks the x, y, and brightness values of 
  3.2_CCDFlats/3.2_flat_real_rotated.00000010.FIT
  from the center to the upper right corner of the fits file

flat_master.fits 
  The master flat file from section 4.3
flat_master_r90.fits 
  The flat file from section 4.3 after the CCD was rotated 90 degrees.
neg10BIAS_master.fits
  The master bias file when the CCD temperature was -10 degrees C. Section 4.1
pos10BIAS_master.fits
  The master bias file when the CCD temperature was +10 degrees C. Section 4.1
neg10dark_adjusted.fits
  The Adjusted (master dark minus the master bias) master dark file when the 
  CCD's temperature was -10 degrees C. Section 4.2
pos10dark_adjusted.fits
  The Adjusted (master dark minus the master bias) master dark file when the 
  CCD's temperature was +10 degrees C. Section 4.2
bad_pixel_map.fits
  The bad pixel map created from the bias and dark frames. Section 4.4
masterflat_50E-6m.fits
  The fits file of the master spectrum flat cropped to only the 50 micrometer slit

Important Terminal Outputs 
  A file of information I was asked to calculate that I keep mostly as my 
  own log book

lab0_intro.pdf
  The Lab 0 booklet from the course GitHub page
quick-start_guide_to_CCDSoft.pdf
  The quick start guide to CCD Soft from the course GitHub page
spectrograph_instructions.pdf
  The spectrograph instructions from the course GitHub page

neg10BIAS_clipped.pdf
  The count distributions for the BIAS frames, excluding outliers, when the 
  CCD was set to -10 degrees C. Section 4.1
neg10BIAS_raw.pdf
  The count distributions for the BIAS frames, including outliers, when the 
  CCD was set to -10 degrees C. Section 4.1
neg10DARK_adjusted.pdf
  The count distributions for the DARK frames, excluding outliers, 
  minus the master bias (for -10 degrees C) when the CCD was set to 
  -10 degrees C. Section 4.2
neg10DARK_cut.pdf
  The count distributions for the DARK frames, excluding outliers, 
  when the CCD was set to -10 degrees C. Section 4.2
neg10DARK_darkcurrent.pdf
  The count distributions for the DARK current, excluding outliers,
  minus the master bias, when the CCD was set to -10 degrees C. Section 4.2
neg10DARK_raw.pdf
  The count distributions for the DARK frames, including outliers, 
  when the CCD was set to -10 degrees C. Section 4.2
neg10DARK_typical-exposure.pdf
  The graph of number of typical counts (above the bias) as a function of 
  the exposure time when the CCD was set to -10 degrees C

pos10BIAS_clipped.pdf
  The count distributions for the BIAS frames, excluding outliers, when the 
  CCD was set to +10 degrees C. Section 4.1
pos10BIAS_raw.pdf
  The count distributions for the BIAS frames, including outliers, when the 
  CCD was set to +10 degrees C. Section 4.1
pos10DARK_adjusted.pdf
  The count distributions for the DARK frames, excluding outliers, 
  minus the master bias (for +10 degrees C) when the CCD was set to 
  -10 degrees C. Section 4.2
pos10DARK_cut.pdf
  The count distributions for the DARK frames, excluding outliers, 
  when the CCD was set to +10 degrees C. Section 4.2
pos10DARK_darkcurrent.pdf
  The count distributions for the DARK current, excluding outliers,
  minus the master bias, when the CCD was set to +10 degrees C. Section 4.2
pos10DARK_raw.pdf
  The count distributions for the DARK frames, including outliers, 
  when the CCD was set to +10 degrees C. Section 4.2
pos10DARK_typical-exposure.pdf
  The graph of number of typical counts (above the bias) as a function of 
  the exposure time when the CCD was set to +10 degrees C

flat_master_raw.pdf
  The histogram showing the distribution of counts in the pixels of the flat frames.
flat_master_raw_r90.pdf
  The histogram showing the distribution of counts in the pixels of the flat frames
  when the CCD is rotated by 90 degrees.

brightness-distance.pdf
  The relative brightness plotted as a function of relative distance for the flats 
  at both 0 degrees and 90 degrees

spectrograph_crop_fit.pdf
  Average count value per column of pixels for the dome flat with quadratic fit. 
  Section 4.5
spectrograph_crop_values-Lamp.pdf
  Average count value per column of pixels in the arc lamp spectrum.
  Section 4.5
spectrograph_cropnorm_fit.pdf
  Normalized average count value per column of pixels for the dome flat.
  Section 4.5
spectrograph_lampadj_fit.pdf
  Average count value per column of pixels in the arc lamp spectrum adjusted by the
  the quadratic fit from spectrograph_crop_fit.pdf
  Section 4.5
spectrograph_lampcrop_fit.pdf
  Same as spectrograph_crop_values-Lamp.pdf
  Section 4.5
