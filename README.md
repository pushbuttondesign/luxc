# luxC

#### A C toolbox for lighting and colour science
#### Many of these functions are validated against the equilavent python functions in: https://github.com/pushbuttondesign/luxpy

### Usage:
- include .h in your file and add the .c at compile time
- fill the SPD_data structure with a list of spectral data values and corrosponding wavelength values
- pass a pointer to the structure to the various operations you want to perform on the data
- most operations pass back a pointer to same data structure, with additinal fields filled in

### TODO
- Validate all the functions, implement ceedling unit testing
- Optimise calculation speed, current CCT trial and error calculation in colourtemp() is very slow and duplicated, swap to alternative equation, several exist
- Stop usage of global structures and pass each variable individually to each function, this will avoid duplicate calculations or needing to add 'bit set' flags & checks

*********************************************

## Colour space termonology used in this program
#### Tristumulus values  -  Associated chromaticity coordinates
- CIE 1931 (XYZ)  -  CIE 1931 (xyY)
- CIE 1964 (X10Y10Z10)  -  CIE 1931 (x10y10Y10)
- CIE 1976 (L*u*v*)  -  CIE 1976 (u’v’L*) CIE 1976 UCS Diagram
- CIE 1976 (L*10u*10v*10) -  CIE 1976 (u’10v’10L*)
- CIE 1976 (L*a*b*) -  none
- CIE 1976 (L*10a*10b*10)  -  none
- CIE 1964 (U*V*W*)  -  CIE 1960 (uvW*), CIE 1960 UCS Diagram (u = u', v = 2/3v', W* = L*)
- CIE 2002 (CIECAM02)  -  CIE 2002 (J'a'b') CIECAM02 UCS Diagram

**************************************

## The following functions are included:

- blackbody - Computes radiance for every wavelength of blackbody of given temprature between given wavelength range

- colourtemp - Caculate colour temprature in Kelvin of normalied spectral irradiance curve

- illuminantD - Calculates CIE Illuminant D as normalised spectral radiance at 5nm intervals between 300nm - 830nm. CIE guidelines state that other intervals should be caculated by liner interpolation

- TM_30_15 - Calculates TM-30-15 Rf and Rg values

- normaliseSPDtoONE - normalise Spectral Power Distribution between 0 and 1

- normaliseSPDtoYonehundred - normalise SPD so that Y = 100 (this is for CIE 1931 (XYZ) tristimulus values caculation)

- stderror - Caculates standard error of the estimate (liner regression) of two curves

- SPDtoCoord10 - Calculates CIE 1976 (u’v’L*) CIE 1976 UCS Diagram, CIE 1931 (xyY) colour co-ordinates and CIE 1931 (XYZ) tristimulus values from either the spectral irradiance or spectral radiance data. 10 degree values used. No luminance normalisation takes place

- SPDtoCoord2 - Calculates CIE 1976 (u’v’L*) CIE 1976 UCS Diagram, CIE 1931 (xyY) colour co-ordinates and CIE 1931 (XYZ) tristimulus values from either the spectral irradiance or spectral radiance data. 2 degree values used. No luminance normalisation takes place

- CIECAM02 - Caculate CIECAM02 Cordinates

- CIECAM02forTM30 - Caculate CIECAM02 Cordinates as required by TM-30-15 specification from 10 deg XYZ values

- delta - Calculates chromaticity euclidean difference, a delta between two points in the same CIE colour space

- ajust_wavlength_range - Ajust the wavelength range of an SPD data vector such that data between the selected range is returned.

- ajust_wavelength_interval - Ajusts the wavelength interval of data supplied, starting from the same wavelength

- gaussian - Returns a gaussian curve around a center

- melanopic_luminance - Calculates the equivalent melanopic luminance as defined by Lucas et al., "Measuring and using light in the melanopsin age." Trends in Neuroscience, Jan 2014. Referenced by Well Building Standard v1 2016 which uses it to calculate equilivant melanopic lux (EML) as defined by that standard, Appendix C Table L1 & L2

- circadian_luminance - Calculates circadian luminance as defined by http://lucasgroup.lab.manchester.ac.uk/research/measuringmelanopicilluminance/

***********************************************************

## These caculations were pulled from a number of refferences:

- The CIE standards - [http://cie.co.at/publications/international-standards]

- The IES TM-30-15 standard - [https://www.ies.org/product/ies-method-for-evaluating-light-source-color-rendition/]

- A toolbox for caculating melanopic_luminance - [http://lucasgroup.lab.manchester.ac.uk/measuringmelanopicilluminance/]

#### Where I needed additional detail I turned to the following detailed books:

- Colour calculations with the CIE system - [https://www.amazon.co.uk/Color-Science-Concepts-Quantitative-Formulae/dp/0471399183]

- The SLL lighting handbook - [https://www.cibse.org/knowledge/knowledge-items/detail?id=a0q20000008I6xeAAC]

- The IES lighting handbook - [https://www.ies.org/product/lighting-handbook-10th-edition/]

- Zumtobel lighting handbook - [https://www.zumtobel.com/PDB/teaser/EN/lichthandbuch.pdf]

#### Other refferences providing useful background to many of these caculations and their application to circadian or human centric lighting:

- Effects of light on human circadian rhythms, our understanding & its limits summerised in a roadmap produced by the CIE in 2016 - [http://www.cie.co.at/publications/research-roadmap-healthful-interior-lighting-applications]

- Another summary of light and its effects on people is 'Measuring and using light in the Melanopsin age' - [https://personalpages.manchester.ac.uk/staff/robert.lucas/Lucas%20et%20al%2014.pdf]

- Experiential and psyco-physics based effects of lighting touching on some of the limitations of the CIE photometric system to be aware of and proposes solutions and some new ways of thinking. - [https://www.amazon.co.uk/Value-Metrics-Better-Lighting-Monograph/dp/0819493228]

***************************

## Units used in this program:

- wavelengths                           //range [nm]
- spectral_radiance                     //range [W·sr−1·m−2·nm-1]
- spectral_rad_normalised               //range [W·sr−1·m−2·nm-1]
- spectral_irradiance                   //range [W·m2·nm-1]
- spectral_irrad_normalised             //range [W·m2·nm-1]
- radiance                              //single number [W·sr−1·m−2]
- luminance                             //single number [lm]
- irradiance                            //single number [W·m2]
- illuminance                           //single number [lux]
- colour_temp                           //single number [K]
- radiant_intensity                     //single number [W·sr-1]
- radiant_flux                          //single number [W]
- array_len                             //single number [count]
- spectral_power_distribution_interval  //single number, wavelength interals
- uprime_10deg;
- vprime_10deg;
- twothirds_vprime_10deg;
- uprime_2deg;
- vprime_2deg;
- twothirds_vprime_2deg;
- x_2deg;
- y_2deg;
- z_2deg;
- X_2deg;
- Y_2deg;
- Z_2deg;
- astar_2deg;
- bstar_2deg;
- Lstar_2deg;
- ustar_2deg;
- vstar_2deg;
- Jprime_2deg;
- aprime_2deg;
- bprime_2deg;
- x_10deg;
- y_10deg;
- z_10deg;
- X_10deg;
- Y_10deg;
- Z_10deg;
- astar_10deg;
- bstar_10deg;
- Lstar_10deg;
- ustar_10deg;
- vstar_10deg;
- Jprime_10deg;
- aprime_10deg;
- bprime_10deg;
- hueangle;
- TM_Rf;
- TM_Rg;
