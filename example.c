//tests SPDlib

#include <stdio.h>
#include "SPDlib.h"
#include <time.h>

int main() {
    
    printf("\ntesting SPDlib\n\n");
    
    //setup wavelength
    int wmin = 360;
    int wmax = 800;
    
    //build a test LED file
    struct SPD_data *tmtest = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    tmtest->arr_len = 7 + 1;
    tmtest->wavelength = (long double*)malloc(sizeof(long double) * tmtest->arr_len);
    tmtest->spec_radiance = (long double*)malloc(sizeof(long double) * tmtest->arr_len);
    tmtest->spec_rad_normalised = (long double*)malloc(sizeof(long double) * tmtest->arr_len);
    tmtest->wavelength[0] = 200;
    tmtest->wavelength[1] = 300;
    tmtest->wavelength[2] = 400;
    tmtest->wavelength[3] = 500;
    tmtest->wavelength[4] = 600;
    tmtest->wavelength[5] = 700;
    tmtest->wavelength[6] = 800;
    tmtest->wavelength[7] = 900;
    tmtest->spec_radiance[0] = 0;
    tmtest->spec_radiance[1] = 1;
    tmtest->spec_radiance[2] = 0.6;
    tmtest->spec_radiance[3] = 0.1;
    tmtest->spec_radiance[4] = 0.5;
    tmtest->spec_radiance[5] = 0.3;
    tmtest->spec_radiance[6] = 0.3;
    tmtest->spec_radiance[7] = 0;
    struct SPD_data *tmtest2 = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    tmtest2->arr_len = 900 - 200 + 1;
    tmtest2->wavelength = (long double*)malloc(sizeof(long double) * tmtest2->arr_len);
    tmtest2->spec_radiance = (long double*)malloc(sizeof(long double) * tmtest2->arr_len);
    tmtest2->spec_rad_normalised = (long double*)malloc(sizeof(long double) * tmtest2->arr_len);
    ajust_wavelength_interval(tmtest, tmtest2, 1, 3);
    
    printf("\nsimulated LED spectrum loaded into TM-30 reads:\n");
    for (int i = 0; i < tmtest2->arr_len; i++) {
        printf("%Lf %Lf\n", tmtest2->wavelength[i], tmtest2->spec_radiance[i]);
    }
    
    //build test blackbody input
    struct SPD_data *illA;
    illA = blackbody(wmin, wmax, 4800, 1);
    printf("\n\nsimulated LED spectrum loaded into TM-30 reads:\n");
    for (int i = 0; i < illA->arr_len; i++) {
        printf("%Lf %Lf\n", illA->wavelength[i], illA->spec_radiance[i]);
    }
    
    
    //test TM-30-15()
    struct CIE_coordinates *tmval = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    tmval = TM_30_15(illA, tmval, 3);
    //tmval = TM_30_15(tmtest2, tmval, 3);
    printf("\nTM-30-15 Rf: %Lf\n", tmval->TM_Rf);
    printf("TM-30-15 Rg: %Lf\n", tmval->TM_Rg);
    
    free(tmval);
    free(tmtest->wavelength);
    free(tmtest->spec_radiance);
    free(tmtest->spec_rad_normalised);
    free(tmtest);
    free(tmtest2->wavelength);
    free(tmtest2->spec_radiance);
    free(tmtest2->spec_rad_normalised);
    free(tmtest2);
    
    
    //test blackbody()
    struct SPD_data *TESTBB1, *TESTBB2, *TESTBB3;
    TESTBB1 = blackbody(wmin, wmax, 5000, 1);
    TESTBB2 = blackbody(wmin, wmax, 2856, 1);
    
    
    //test stderror()
    long double COMPAREBB1, COMPAREBB2;
    COMPAREBB1 = stderror(TESTBB1->spec_rad_normalised, TESTBB1->spec_rad_normalised, (wmax - wmin));
    COMPAREBB2 = stderror(TESTBB1->spec_rad_normalised, TESTBB2->spec_rad_normalised, (wmax - wmin));
    printf("standard error of difference between 5000K and 5000K blackbodies is: %Lf\n", COMPAREBB1);
    printf("standard error of difference between 5000K and 2000K blackbodies is: %Lf\n\n", COMPAREBB2);
    
    
    //test ajust_wavlength_range()
    printf("\nBlackbody values are:\n");
    for (int i = 0; i < TESTBB2->arr_len; i += 1) {
        printf("%.0Lf %Lf %Lf\n", TESTBB2->wavelength[i], TESTBB2->spec_radiance[i], TESTBB2->spec_rad_normalised[i]);
    }
    struct SPD_data *devided_TESTBB2 = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    devided_TESTBB2->arr_len = 500 - 200 + 1;
    devided_TESTBB2->wavelength = (long double*)malloc(sizeof(long double) * devided_TESTBB2->arr_len);
    devided_TESTBB2->spec_radiance = (long double*)malloc(sizeof(long double) * devided_TESTBB2->arr_len);
    
    ajust_wavlength_range(TESTBB2, devided_TESTBB2, 200, 500, 3);
    
    printf("\nTruncated Blackbody values are:\n");
    for (int i = 0; i < devided_TESTBB2->arr_len; i += 1) {
        printf("%.0Lf %Lf\n", devided_TESTBB2->wavelength[i], devided_TESTBB2->spec_radiance[i]);
    }
    
    free(TESTBB2->spec_rad_normalised);
    free(TESTBB2->spec_radiance);
    free(TESTBB2->wavelength);
    free(TESTBB2);
    
    free(devided_TESTBB2->spec_radiance);
    free(devided_TESTBB2->wavelength);
    free(devided_TESTBB2);
    
    
    //test ajust_wavelength_interval() upsample
    TESTBB3 = blackbody(wmin, wmax, 2856, 10);
    printf("\nBlackbody values are:\n");
    for (int i = 0; i < TESTBB3->arr_len; i++) {
        printf("%.0Lf %.0Lf\n", TESTBB3->wavelength[i], TESTBB3->spec_radiance[i]);
    }
    
    struct SPD_data intup;
    intup.wav_interval = 3;
    intup.arr_len = floor((wmax - wmin) / (double)intup.wav_interval) + 1;
    intup.wavelength = (long double*)calloc(intup.arr_len, sizeof(long double));
    intup.spec_radiance = (long double*)calloc(intup.arr_len, sizeof(long double));
    
    ajust_wavelength_interval(TESTBB3, &intup, 3, 3);

    printf("\najusted interval is:\n");
    for (int i = 0; i < intup.arr_len; i++) {
        printf("%.0Lf %.0Lf\n", intup.wavelength[i], intup.spec_radiance[i]);
    }
    
    free(intup.wavelength);
	free(intup.spec_radiance);
    
    //test ajust_wavelength_interval() downsample
    struct SPD_data intdown;
    intdown.wav_interval = 13;
    intdown.arr_len = floor((wmax - wmin) / (double)intdown.wav_interval) + 1;
    intdown.wavelength = (long double*)calloc(intdown.arr_len, sizeof(long double));
    intdown.spec_radiance = (long double*)calloc(intdown.arr_len, sizeof(long double));
    
    ajust_wavelength_interval(TESTBB3, &intdown, 13, 3);
    
    printf("\najusted interval is:\n");
    for (int i = 0; i < intdown.arr_len; i++) {
        printf("%.0Lf %.0Lf\n", intdown.wavelength[i], intdown.spec_radiance[i]);
    }
	
	free(intdown.wavelength);
	free(intdown.spec_radiance);
	free(TESTBB3->spec_rad_normalised);
    free(TESTBB3->spec_radiance);
    free(TESTBB3->wavelength);
    free(TESTBB3);
	
    
    //test illumiantD()
    struct SPD_data *illD = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    illD->arr_len = (800 - 300) * 1 + 1;
    illD->wavelength = (long double*)malloc(sizeof(long double) * illD->arr_len);
    illD->spec_radiance = (long double*)malloc(sizeof(long double) * illD->arr_len);
    illD->spec_rad_normalised = (long double*)malloc(sizeof(long double) * illD->arr_len);
    illD = illuminantD(300, 800, 1, 6500, illD);
    printf("illuminant D values are:\n");
    for (int i = 0; i < illD->arr_len; i++) {
        printf("%Lf %Lf\n", illD->wavelength[i], illD->spec_rad_normalised[i]);
    }
    
    
    //setup timer
    clock_t start, end;
    long double duration;
    
    
    //test SPDtoCoord()
    struct CIE_coordinates *CIE10, *CIE2;
    CIE10 = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    CIE2 = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    start = clock();
    SPDtoCoord10(TESTBB1, CIE10, 100, 3);
    end = clock();
    printf("CIE u'10 v'10 cordinates of 5000K BB are: %Lf %Lf\n", CIE10->uprime_10deg, CIE10->vprime_10deg);
    printf("CIE x10 y10 cordinates of 5000K BB are: %Lf %Lf\n", CIE10->x_10deg, CIE10->y_10deg);
    printf("CIE X10 Y10 Z10 cordinates of 5000K BB are: %Lf %Lf %Lf\n", CIE10->X_10deg, CIE10->Y_10deg, CIE10->Z_10deg);
    SPDtoCoord2(TESTBB1, CIE2, 100, 3);
    printf("CIE u'2 v'2  cordinates of 5000K BB are: %Lf %Lf\n", CIE2->uprime_2deg, CIE2->vprime_2deg);
    printf("CIE x'2 y'2 cordinates of 5000K BB are: %Lf %Lf\n", CIE2->x_2deg, CIE2->x_2deg);
    printf("CIE X2 Y2 Z2 cordinates of 5000K BB are: %Lf %Lf %Lf\n\n", CIE2->X_2deg, CIE2->Y_2deg, CIE2->Z_2deg);
    SPDtoCoord10(illD, CIE10, 100, 3);
    printf("CIE u'10 v'10 cordinates of 6500K Illuminant D are: %Lf %Lf\n", CIE10->uprime_10deg, CIE10->vprime_10deg);
    printf("CIE x'10 y'10 cordinates of 6500K Illuminant D are: %Lf %Lf\n", CIE10->x_10deg, CIE10->y_10deg);
    printf("CIE X10 Y10 Z10 cordinates of 6500K Illuminant D are: %Lf %Lf %Lf\n", CIE10->X_10deg, CIE10->Y_10deg, CIE10->Z_10deg);
    SPDtoCoord2(illD, CIE2, 100, 3);
    printf("CIE u'2 v'2 cordinates of 6500K Illuminant D are: %Lf %Lf\n", CIE2->uprime_2deg, CIE2->vprime_2deg);
    printf("CIE x'2 y'2 cordinates of 6500K Illuminant D are: %Lf %Lf\n", CIE2->x_2deg, CIE2->y_2deg);
    printf("CIE X2 Y2 Z2 cordinates of 6500K Illuminant D are: %Lf %Lf %Lf\n\n", CIE2->X_2deg, CIE2->Y_2deg, CIE2->Z_2deg);
    free(CIE2);
    free(CIE10);
    
    
    //print time results for SPDtoCoord()
    duration = end - start;
    duration /= CLOCKS_PER_SEC;
    printf("\nDuration of one SPDtoCoord caculation is: %Lf in seconds\n\n", duration);
    
    
    //test colourtemp()
    start = clock();
    long double CCT2 = colourtemp(wmin, wmax, 1, 1000, 10000, 10, TESTBB1, 3);
    end = clock();
    printf("colour temprature of a 5000K BB is: %.0Lf\n", CCT2);
    long double CCT3 = colourtemp(wmin, wmax, 1, 1000, 10000, 10, illD, 3);
    printf("colour temprature of Illuminant D equasion ran at 6500K is: %.0Lf\n", CCT3);
    
    free(TESTBB1->spec_rad_normalised);
    free(TESTBB1->spec_radiance);
    free(TESTBB1->wavelength);
    free(TESTBB1);
    
    free(illD->spec_radiance);
    free(illD->spec_rad_normalised);
    free(illD->wavelength);
    free(illD);
    
    
    //print time results for colourtemp()
    duration = end - start;
    duration /= CLOCKS_PER_SEC;
    printf("\nDuration of one colourtemp caculation is: %Lf in seconds\n\n", duration);
    
    
    //test gaussian()
    struct SPD_data g;
    g.arr_len = 780 - 380 + 1;
    g.wavelength = (long double*)malloc(sizeof(long double) * g.arr_len);
    g.spec_radiance = (long double*)malloc(sizeof(long double) * g.arr_len);
    
    gaussian(g.spec_radiance, g.wavelength, 380, 780, 1, 1532, 300, 555);
    
    printf("generated gaussian:\n")
    for (int i = 0; i < g.arr_len; i++) {
        printf("%Lf %Lf\n", g.wavelength[i], g.spec_radiance[i]);
    }
    
    
    //test ajust_wavinterval()
	struct SPD_data gint;
	gint.arr_len = (780 - 380) / 5 + 1;
    gint.wavelength = (long double*)calloc(gint.arr_len, sizeof(long double));
    gint.spec_radiance = (long double*)malloc(sizeof(long double) * gint.arr_len);
    gint.spec_irradiance = NULL;
    
	ajust_wavelength_interval(&g, &gint, 5, 3);
    
    
    //test melanopic_luminance()
	printf("input to eml:\n");
    for (int i = 0; i < gint.arr_len; i++) {
        printf("%Lf %Lf\n", gint.wavelength[i], gint.spec_radiance[i]);
    }
    
    long double ML = 5;
    ML = melanopic_luminance(&gint, 3);
    
    printf("\nMelanopic Luminance is: %Lf\n\n", ML);
    
	free(g.wavelength);
    free(g.spec_radiance);
    free(gint.wavelength);
    free(gint.spec_radiance);
    
    return 0;
}
