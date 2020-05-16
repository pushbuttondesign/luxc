//
//  SPDlib.h
//  SPD library
//
//  Libary contains functions for lighting and colour science
//  They operate on spectral power distribution data
//

#ifndef SPDlib_h
#define SPDlib_h

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//macros
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#define M_E (2.71828)
#endif
#define MINmacro(a, b) (((a) + (b) - fabsl((a) - (b))) * 0.5)
#define MAXmacro(a, b) (((a) + (b) + fabsl((a) - (b))) * 0.5)
#define CIE_W_MIN 360    //360nm is bottom of visible range as defined by CIE
#define CIE_W_MAX 830    //830nm is top of visible range as defined by CIE
#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))
#define tand(x) (tan(fmod((x),360) * M_PI / 180))
#define asind(x) (asin(fmod((x),360) * M_PI / 180))
#define acosd(x) (acos(fmod((x),360) * M_PI / 180))
#define atand(x) (atan(fmod((x),360) * M_PI / 180))
#define DEGREES_TO_RADIANS(degrees) ((M_PI * degrees)/180)
#define RADIANS_TO_DEGREES(radians) ((radians) * (180.0 / M_PI))

//global variables
struct SPD_data{
    long double *wavelength;             //range [nm]
    long double *spec_radiance;          //range [W·sr−1·m−2·nm-1]
    long double *spec_rad_normalised;    //range [W·sr−1·m−2·nm-1]
    long double *spec_irradiance;        //range [W·m2·nm-1]
    long double *spec_irrad_normalised;  //range [W·m2·nm-1]
    long double radiance;               //single number [W·sr−1·m−2]
    long double luminance;              //single number [lm]
    long double irradiance;             //single number [W·m2]
    long double illuminance;            //single number [lux]
    long double colour_temp;            //single number [K]
    long double radiant_intensity;      //single number [W·sr-1]
    long double radiant_flux;           //single number [W]
    int arr_len;                        //single number [count]
    int wav_interval;                   //single number, wavelength interals
};

struct CIE_coordinates{
    long double uprime_10deg;
    long double vprime_10deg;
    long double twothirds_vprime_10deg;
    long double uprime_2deg;
    long double vprime_2deg;
    long double twothirds_vprime_2deg;
    long double x_2deg;
    long double y_2deg;
    long double z_2deg;
    long double X_2deg;
    long double Y_2deg;
    long double Z_2deg;
    long double astar_2deg;
    long double bstar_2deg;
    long double Lstar_2deg;
    long double ustar_2deg;
    long double vstar_2deg;
    long double Jprime_2deg;
    long double aprime_2deg;
    long double bprime_2deg;
    long double x_10deg;
    long double y_10deg;
    long double z_10deg;
    long double X_10deg;
    long double Y_10deg;
    long double Z_10deg;
    long double astar_10deg;
    long double bstar_10deg;
    long double Lstar_10deg;
    long double ustar_10deg;
    long double vstar_10deg;
    long double Jprime_10deg;
    long double aprime_10deg;
    long double bprime_10deg;
    long double hueangle;
    long double TM_Rf;
    long double TM_Rg;
};

struct SPDandcoord {
    struct SPD_data SPD_data;
    struct CIE_coordinates CIEcoord;
};

extern int length_of_illuminantD;

//Colour space termonology used in this program
//Tristumulus values – associated chromaticity coordinates
//CIE 1931 (XYZ) - CIE 1931 (xyY)
//CIE 1964 (X10Y10Z10) - CIE 1931 (x10y10Y10)
//CIE 1976 (L*u*v*) - CIE 1976 (u’v’L*) CIE 1976 UCS Diagram
//CIE 1976 (L*10u*10v*10) - CIE 1976 (u’10v’10L*)
//CIE 1976 (L*a*b*) - none
//CIE 1976 (L*10a*10b*10) - none
//CIE 1964 (U*V*W*) – CIE 1960 (uvW*) CIE 1960 UCS Diagram (u = u', v = 2/3v', W* = L*)
//CIE 2002 (CIECAM02) - CIE 2002 (J'a'b') CIECAM02 UCS Diagram

//function declarations
struct SPD_data* blackbody(int min, int max, int temprature, int interval);
long double colourtemp(int w_min, int w_max, int w_int, int cct_min, int cct_max, int cct_inter, struct SPD_data* test, int vector_flag);
struct SPD_data* illuminantD(int w_min, int w_max, int interval, int temprature, struct SPD_data* SPD_out);
struct CIE_coordinates* TM_30_15(struct SPD_data* SPD_in, struct CIE_coordinates* results, int vector_flag);
long double* normaliseSPDtoONE(long double* SPD_in, long double* SPD_out, int len);
struct SPD_data* normaliseSPDtoYonehundred(struct SPD_data* SPD_in, int vector_flag);
long double stderror(long double* target, long double* comparison, int len);
struct CIE_coordinates* SPDtoCoord10(struct SPD_data* SPD_in, struct CIE_coordinates* results, long double k, int vector_flag);
struct CIE_coordinates* SPDtoCoord2(struct SPD_data* SPD_in, struct CIE_coordinates* results, long double k, int vector_flag);
struct CIE_coordinates* CIECAM02(struct CIE_coordinates*, int flag);
struct CIE_coordinates* CIECAM02forTM30(struct CIE_coordinates* CEScoordREF, struct CIE_coordinates* adapting_illuminant_coord, int flag);
long double delta(struct CIE_coordinates *coord1, struct CIE_coordinates *coord2, int mode);
struct SPD_data* ajust_wavlength_range(struct SPD_data* SPD_in, struct SPD_data* SPD_out, int minwave, int maxwave, int vector_flag);
struct SPD_data* ajust_wavelength_interval(struct SPD_data* SPD_in, struct SPD_data* SPD_out, int interval, int vector_flag);
int gaussian(long double *results, long double *wavelength, const int wstart, const int wstop, const int wint, const double peak_height, const double width, const double peak_center);

#endif /* SPDlib_h */
