//
//  SPDlib.h
//  SPD library
//
//  Libary contains functions for lighting and colour science
//  They operate on spectral power distribution data
//
//  TODO:
//      - Validate all functions, there is a small error in CCT calculation
//      - Optimise calculation speed, CCT trial and error calculation is very slow, swap to alternative equation
//      - Stop usage of 3x global structures and pass each variable individually to each function
//        This will avoid duplicate calculations or needing to add 'bit set' flags & checks
//

//include & define
#include "SPDlib.h"

//global variables

//
//  Computes radiance for every wavelength of blackbody of given temprature between given wavelength range
//	lowest temp with values in range 360nm - 830nm is 25K
//	lowest temp with values in range 380nm - 780nm is 520K
//
//  EXTERNAL REQUIREMENTS: must free(): SPD, SPD->wavelength, SPD->spec_radiance, SPD->spec_rad_normalised
//  INPUTS: - min wavelength to begin caculation from (nm), int
//          - max wavelength to end caculation at (nm), int
//          - temprature (kelvin), int
//          - wavelength interval to caculate blackbody at (nm), int
//  OUTPUT: - pointer to malloced SPD_data structure, SPD_data*, containing:
//              - spectral radiance (watt per steradian per square metre (W·sr−1·m−2)) per wavelength at 1nm intervals (W·sr−1·m−2·nm-1))
//              - normalised spectral radiance radiance
//              - total radiance over wavelength range
//
struct SPD_data* blackbody(int min, int max, int temprature, int interval) {
    //Check arguments
    if (min > max || min == max){
        fprintf(stderr, "ERROR: minimum wavelength must be smaller than maximum wavelength for blackbody compute\n");
        return NULL;
    } else if (temprature < 1) {
        fprintf(stderr, "ERROR: temprature must be positive interger for blackbody compute\n");
        return NULL;
    } else if (interval < 1) {
        fprintf(stderr, "ERROR: interval cannot be less than 1nm\n");
        return NULL;
    } else if (temprature < 550) {
        fprintf(stderr, "WARNING: temperature cannot be less than 550K for use with CIE functions\n");
    }
    
    const double H = 6.626070040e-34;   //Planck's constant (Joule-seconds) 6.626070040e-34
    const double C = 299792458;         //speed of light in vacume (meters per second)
    const double K = 1.3806488e-23;     //Boltzmann's constant (Joules per Kelvin) 1.3806488e-23
    const double nm_to_m = 1e-9;        //conversion between nm and m
    
    struct SPD_data *SPD;               //allocate memory based on wavelength bounds
    SPD = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    SPD->arr_len = ((max - min) / interval + 1);
    SPD->wavelength = (long double*)malloc(sizeof(long double) * SPD->arr_len);
    SPD->spec_radiance = (long double*)malloc(sizeof(long double) * SPD->arr_len);
    SPD->spec_rad_normalised = (long double*)malloc(sizeof(long double) * SPD->arr_len);
    SPD->spec_irradiance = NULL;        //specifically NULL to avoid colour temp looking at bad values (other functions look at irrad by default)
    SPD->spec_irrad_normalised = NULL;  //debug - remove this
    
    SPD->radiance = 0;
    for (int i = 0; i < SPD->arr_len; i++) {
        //Computes radiance for every wavelength of blackbody of given temprature
        SPD->wavelength[i] = min + (interval * i);
        SPD->spec_radiance[i] = ((2 * H * pow(C, 2)) / (pow((SPD->wavelength[i] * nm_to_m), 5))) * (1 / (expm1((H * C) / ((SPD->wavelength[i] * nm_to_m) * K * temprature))));
        //Caculate total radiance
        SPD->radiance += SPD->spec_radiance[i] * interval;
    }
    
    //normalise SPD
    normaliseSPDtoONE(SPD->spec_radiance, SPD->spec_rad_normalised, SPD->arr_len);
    
    return SPD;
}

//
//  Caculates standard error of the estimate (liner regression) of two curves
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - Two long double* pointers to two vectors describing values of a curve. Vectors intervals, units and length must match, SPD_data*
//          - Length of each vector, int
//  OUTPUT: - Standard error of estimate, long double
//
long double stderror(long double *target, long double *comparison, int len) {
    //check arguments
    if (len < 1) {
        fprintf(stderr, "ERROR: length of vectors for standard error computation must be a positive interger\n");
        return 1;
    }
    
    long double sum = 0, div, error;
    struct {
        long double *difference;
        long double *square;
    } standard_error;
    
    standard_error.difference = (long double*)calloc(len, sizeof(long double));
    standard_error.square = (long double*)calloc(len, sizeof(long double));
    
    for (int i = 0; i < len; i++) {
        //compute difference and square of difference
        standard_error.difference[i] = target[i] - comparison[i];
        standard_error.square[i] = standard_error.difference[i] * standard_error.difference[i];
        //sum squares
        sum += standard_error.square[i];
    }
    
    //compute error
    div = sum / (len - 2);
    error = sqrt(div);
    
    free(standard_error.difference);
    free(standard_error.square);
    
    return error;
}

//
//  Caculate colour temprature in Kelvin of normalied spectral irradiance curve
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - minimum wavelength, int
//          - maximum wavelength, int
//          - wavelength interval, int
//          - minimum CCT to test, int
//          - maximum CCT to test, int
//          - inverval between CCTs to test at, between stated min and max, int
//          - SPD_data containing normalised irradiance wavelength data to test, SPD_data*
//          - vector flag that tells the program which SPD vector to run caculation on, int
//          - 1 = spectral_irradiance, 2 = spec_irrad_normalised, 3 = spectral_radiance, 4 = spec_rad_normalised
//  OUTPUT: - the temprature in K of the nearest blackbody to the test spectrum as a long double, 0 on error
//
long double colourtemp(int w_min, int w_max, int w_int, int cct_min, int cct_max, int cct_inter, struct SPD_data *test, int vector_flag){
    //Check arguments
    if (w_min > w_max || w_min == w_max){
        fprintf(stderr, "ERROR: minimum wavelength must be smaller than maximum wavelength to caculate colour temprature\n");
        return 0;
    } else if (cct_min > cct_max || cct_min == cct_max) {
        fprintf(stderr, "ERROR: minimum cct must be smaller than maximum cct to caculate colour temprature\n");
        return 0;
    } else if (test == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to compute colour temprature cannot be NULL\n");
        return 0;
    }
    
    //caculate CIE 1976 (u’v’L*) values of test SPD
    struct CIE_coordinates *test_coord;
    test_coord = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    SPDtoCoord2(test, test_coord, 100, vector_flag);
    
    //find blackbody with closest CIE 1960 (uv) UCS diagram match - equilivent to CIE 1976 u’ 2/3v’
    struct SPD_data *bb;
    long double difference1, difference2, temp = 0.0;
    struct CIE_coordinates *bb_coord;
    bb_coord = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    bb = blackbody(w_min, w_max, cct_min, w_int);
    SPDtoCoord2(bb, bb_coord, 100, 3);
    difference1 = sqrtl((powl((test_coord->uprime_2deg - bb_coord->uprime_2deg), 2) + ((4 / 9.0) * powl((test_coord->vprime_2deg - bb_coord->vprime_2deg), 2))));	//euclidian distance
    
    free(bb->spec_radiance);
    free(bb->spec_rad_normalised);
    free(bb->wavelength);
    free(bb);
    for (int i = (cct_min + cct_inter); i <= cct_max; i += cct_inter) {
        bb = blackbody(w_min, w_max, i, w_int);
        SPDtoCoord2(bb, bb_coord, 100, 3);
        difference2 = sqrtl((powl((test_coord->uprime_2deg - bb_coord->uprime_2deg), 2) + ((4 / 9.0) * powl((test_coord->vprime_2deg - bb_coord->vprime_2deg), 2))));
        if (difference2 < difference1) {
            difference1 = difference2;
            temp = i;
        }
        free(bb->spec_radiance);
        free(bb->spec_rad_normalised);
        free(bb->wavelength);
        free(bb);
    }
    free(bb_coord);
    free(test_coord);

    //confirm value is within CIE CCT vallidity tolerance
    if (difference1 <= 0.05) {
        return temp;
    } else {
        fprintf(stderr, "ERROR: no valid CCT value for given SPD\n");
        return 0;
    }
}

//
//  normalise SPD between 0 and 1
//	assumes positive string of numbers
//
//  EXTERNAL REQUIREMENTS: must free() SPD_out
//  INPUTS: - pointer to an SPD data double vector (for instance spec_radiance or spec_irradiance), SPD_data*
//          - long double *SPD_out, pointer to vector to store normalised data
//          - length of data vector, 1 indexed, int
//  OUTPUT: - pointer to SPD data double vector containing normalised values, double*
//
long double* normaliseSPDtoONE(long double *SPD_in, long double *SPD_out, int len){
    //check arguments
    if (len < 1) {
        fprintf(stderr, "ERROR: length of vector to normalise must be a positive interger\n");
        return NULL;
    } else if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL\n");
        return NULL;
    }
    
    long double max_valu = 0, min_valu = 0;
    
    for (int i = 0; i < len; i++) {
        
        //Find largest value
        if (i <= 0) {
            max_valu = SPD_in[0];
            min_valu = SPD_in[0];
        } else if (i > 0){
            if (max_valu < SPD_in[i]) {
                max_valu = SPD_in[i];
            } else if (min_valu > SPD_in[i]) {
                min_valu = SPD_in[i];
            }
        }
    }

	int test = 0;
	//catch devide by 0 error for strings with no
	if (max_valu > 0 && min_valu >= 0)
	{
    	//Normalise SPD
    	for (int i = 0; i < len; i++) {
        	SPD_out[i] = (SPD_in[i] - min_valu) / (max_valu - min_valu);
        	test += SPD_out[i];
    	}
	} else
	{
		memset(SPD_out, '0', len);
	}
	
	if (test == NAN)
	{
		fprintf(stderr, "ERROR: normalisation has non-int value\n");
		return NULL;
	}
	
    return SPD_out;
}

//
//  normalise SPD so that Y = 100 when CIE 1931 (XYZ) tristimulus values are caculated
//
//  EXTERNAL REQUIREMENTS: must free() SPD_out
//  INPUTS: - SPD_data *SPD_in, with array length, wavelength and either spectral irradiance or spectral radiance data filled and corrosponding spec normalised vector malloc'ed. Wavelength interval must be 1, 2 or 5 nm
//          - int vector_flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - pointer to SPD data double vector containing normalised values, double*
//
struct SPD_data* normaliseSPDtoYonehundred(struct SPD_data *SPD_in, int vector_flag){
    //check arguments
    if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL\n");
        return NULL;
    } else if (SPD_in->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0\n");
        return NULL;
    } else if (SPD_in->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL\n");
        return NULL;
    } else if (SPD_in->spec_irradiance == NULL && SPD_in->spec_radiance == NULL) {
        fprintf(stderr, "ERROR: pointers to spectral data vectors cannot all be NULL\n");
        return NULL;
    }
    
    //open CIE XYZ weighting data of same wavelength interval (1, 2 or 5nm) as SPD_in (xBar, yBar, zBar)
    int data_samples;
    SPD_in->wav_interval = SPD_in->wavelength[1] - SPD_in->wavelength[0];
    FILE *fp;
    switch (SPD_in->wav_interval) {
        case 1:
            data_samples = 401;
            fp = fopen("./data/XYZweightings_data_1nm_intervals_2deg_standard_observer.txt", "r");
            break;
        case 2:
            data_samples = 201;
            fp = fopen("./data/XYZweightings_data_2nm_intervals_2deg_standard_observer.txt", "r");
            break;
        case 5:
            data_samples = 81;
            fp = fopen("./data/XYZweightings_data_5nm_intervals_2deg_standard_observer.txt", "r");
            break;
        default:
            fprintf(stderr, "ERROR: Wavelength interval must be 1, 2 or 5 nm\n");
            return NULL;
    }
    if (fp == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/XYZweightings_data...2deg_standard_observer.txt file.\n");
        exit(EXIT_FAILURE);
    }
    
    //ajust range of data to 380nm - 780nm to match CIE recomendations
    //cannot be done on normalised data as normalisation will no longer make sense
    struct SPD_data *correct_range_data;
    correct_range_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    correct_range_data->wavelength = (long double*)malloc(sizeof(long double) * data_samples);
    long double *invector1, *invector2, *outvector1, *outvector2;
    switch (vector_flag) {
        case 1:
            correct_range_data->spec_irradiance = (long double*)malloc(sizeof(long double) * data_samples);
            correct_range_data->spec_irrad_normalised = (long double*)malloc(sizeof(long double) * data_samples);
            outvector2 = correct_range_data->spec_irrad_normalised;
            outvector1 = correct_range_data->spec_irradiance;
            invector1 = SPD_in->spec_irradiance;
            invector2 = SPD_in->spec_irrad_normalised;
            break;
        case 3:
            correct_range_data->spec_radiance = (long double*)malloc(sizeof(long double) * data_samples);
            correct_range_data->spec_rad_normalised = (long double*)malloc(sizeof(long double) * data_samples);
            outvector2 = correct_range_data->spec_rad_normalised;
            outvector1 = correct_range_data->spec_radiance;
            invector1 = SPD_in->spec_radiance;
            invector2 = SPD_in->spec_rad_normalised;
            break;
        default:
            free(correct_range_data->wavelength);
            free(correct_range_data);
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            return NULL;
    }
    ajust_wavlength_range(SPD_in, correct_range_data, 380, 780, vector_flag);
    free(invector1);
    free(invector2);
    free(SPD_in->wavelength);
    free(SPD_in);
    
    //read in CIE XYZ weighting data from selected data file (xBar, yBar, zBar)
    struct colourWeighting{
        long double *xBar;
        long double *yBar;
        long double *zBar;
    } *weightings;
    weightings = (struct colourWeighting*)calloc(1, sizeof(struct colourWeighting));
    weightings->xBar = (long double*)calloc(data_samples, sizeof(long double));
    weightings->yBar = (long double*)calloc(data_samples, sizeof(long double));
    weightings->zBar = (long double*)calloc(data_samples, sizeof(long double));
    for (int i = 0; i < data_samples; i++) {
        fscanf(fp, "%Lf %Lf %Lf", &weightings->xBar[i], &weightings->yBar[i], &weightings->zBar[i]);
    }
    fclose(fp);
    
    //caculate normalisation factor k such that Y = 100
    long double k = 0;
    for (int i = 0; i < data_samples; i++) {
        //k += weightings->yBar[i] * outvector1[i];
        k += weightings->yBar[i] * outvector1[i] * correct_range_data->wav_interval;
    }
    
    //use normalisation factor k to normalise SPD
    for (int i = 0; i < correct_range_data->arr_len; i++) {
        outvector2[i] = outvector1[i] / k * 100;
    }
    
    //free
    free(weightings->xBar);
    free(weightings->yBar);
    free(weightings->zBar);
    free(weightings);
    
    return correct_range_data;
}

//
//  Calculates chromaticity euclidean difference, a delta between two points in the same CIE colour space.
//
//  EXTERNAL REQUIREMENTS: must free() the results->deta variable of selected mode
//  INPUTS: - two CIE_cordinates pointers with u'v' values filled
//          - Mode, int
//          - 0 is u'10 v'10
//  OUTPUT: - returns a long double containing the delta between the two coordinates tested in the space of the mode selected
//

long double delta(struct CIE_coordinates *coord1, struct CIE_coordinates *coord2, int mode){
    //setup memory for results
    long double result_delta;
    
    switch (mode) {
        case 0:
            //check arguments
            if (coord1->uprime_10deg == 0 || coord2->vprime_10deg == 0) {
                fprintf(stderr, "ERROR: uprime_10deg or vprime_10deg values not provided\n");
                return -1;
            }
            result_delta = sqrtl(pow((coord2->uprime_10deg - coord1->uprime_10deg), 2) + pow((coord2->vprime_10deg - coord1->vprime_10deg), 2));
            break;
        default:
            fprintf(stderr, "ERROR: mode must be 0\n");
            return -1;
            break;
    }
    
    return result_delta;
}

//
//  Calculates CIE 1976 (u’v’L*) CIE 1976 UCS Diagram, CIE 1931 (xyY) colour co-ordinates and CIE 1931 (XYZ) tristimulus values from either the spectral irradiance or spectral radiance data
//  No luminance normalisation takes place. 10 degree.
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - SPD_data* with array length, wavelength and either spectral irradiance or spectral radiance data filled. Wavelength interval must be 1, 2 or 5 nm
//          - CIE_coordinates* with memory allocated for results
//          - int, k normalising constant for caculting CIE 1931 (XYZ) tristimulus values, 1 to perform absolute caculation, 100 to caculate k such that Y = 100 (standard for relitive cases), otherwise enter value such that Y is equal to any convenent number. Note k must be constant when compairing tirstimulus values or chromaticty cordinates.
//          - int flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - CIE_coordinates* with CIE 1976 (u’v’L*) CIE 1976 UCS Diagram, CIE 1931 (xyY) colour co-ordinates and CIE 1931 (XYZ) tristimulus values filled in
//

struct CIE_coordinates* SPDtoCoord10(struct SPD_data *SPD_in, struct CIE_coordinates *results, long double k, int vector_flag) {
    //check arguments
    if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL\n");
        return NULL;
    } else if (SPD_in->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0\n");
        return NULL;
    } else if (SPD_in->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL\n");
        return NULL;
    } else if (SPD_in->spec_irradiance == NULL && SPD_in->spec_radiance == NULL) {
        fprintf(stderr, "ERROR: pointers to spectral data vectors cannot all be NULL\n");
        return NULL;
    } else if (vector_flag != 1 && vector_flag != 3) {
    	fprintf(stderr, "ERROR: vector flag must be 1 or 3\n");
        return NULL;
    } else if (k <= 0) {
    	fprintf(stderr, "ERROR: normalising constant cannot be less than 1\n");
        return NULL;
    }
    
    long double ztest = 0;
    switch (vector_flag)
    {
    	case 1:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_irradiance[i];
    		}
    		break;
    	case 3:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_radiance[i];
    		}
    		break;
    }
    
    if (ztest <= 0) {
    	fprintf(stderr, "ERROR: SPD data must contain values >0\n");
        return NULL;
    }
    
    //open CIE XYZ weighting data of same wavelength interval (1, 2 or 5nm) as SPD_in (xBar, yBar, zBar)
    int data_samples;
    SPD_in->wav_interval = SPD_in->wavelength[1] - SPD_in->wavelength[0];
    FILE *fp;
    switch (SPD_in->wav_interval) {
        case 1:
            data_samples = 401;
            fp = fopen("./data/XYZweightings_data_1nm_intervals_10deg_standard_observer.txt", "r");
            break;
        case 2:
            data_samples = 201;
            fp = fopen("./data/XYZweightings_data_2nm_intervals_10deg_standard_observer.txt", "r");
            break;
        case 5:
            data_samples = 81;
            fp = fopen("./data/XYZweightings_data_5nm_intervals_10deg_standard_observer.txt", "r");
            break;
        default:
            fprintf(stderr, "ERROR: Wavelength interval must be 1, 2 or 5 nm\n");
            return NULL;
    }
    if (fp == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/XYZweightings_data...10deg_standard_observer.txt file.\n");
        exit(EXIT_FAILURE);
    }
    
    //ajust range of data to 380nm - 780nm to match CIE recomendations
    //cannot be done on normalised data as normalisation will no longer make sense
    struct SPD_data *correct_range_data;
    correct_range_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    correct_range_data->wavelength = (long double*)calloc(data_samples, sizeof(long double));
    long double *vector;
    switch (vector_flag) {
        case 1:
            correct_range_data->spec_irradiance = (long double*)calloc(data_samples, sizeof(long double));
            vector = correct_range_data->spec_irradiance;
            break;
        case 3:
            correct_range_data->spec_radiance = (long double*)calloc(data_samples, sizeof(long double));
            vector = correct_range_data->spec_radiance;
            break;
        default:
            free(correct_range_data->wavelength);
            free(correct_range_data);
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            return NULL;
    }
    ajust_wavlength_range(SPD_in, correct_range_data, 380, 780, vector_flag);
    
    //test new SPD data for values (ensure all data was not cropped off)
    ztest = 0;
    switch (vector_flag)
    {
    	case 1:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_irradiance[i];
    		}
    		break;
    	case 3:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_radiance[i];
    		}
    		break;
    }
    
    if (ztest <= 0) {
    	fprintf(stderr, "ERROR: SPD data must contain values >0\n");
        return NULL;
    }
    
    //read in CIE XYZ weighting data from selected data file (xBar, yBar, zBar)
    struct colourWeighting{
        long double *xBar;
        long double *yBar;
        long double *zBar;
    } *weightings;
    weightings = (struct colourWeighting*)calloc(1, sizeof(struct colourWeighting));
    weightings->xBar = (long double*)calloc(data_samples, sizeof(long double));
    weightings->yBar = (long double*)calloc(data_samples, sizeof(long double));
    weightings->zBar = (long double*)calloc(data_samples, sizeof(long double));
    for (int i = 0; i < data_samples; i++) {
        fscanf(fp, "%Lf %Lf %Lf", &weightings->xBar[i], &weightings->yBar[i], &weightings->zBar[i]);
    }
    fclose(fp);
    
    //caculate XYZ values using CIE XYZ weightings (xBar, yBar, zBar)
    results->X_10deg = 0;
    results->Y_10deg = 0;
    results->Z_10deg = 0;
    for (int i = 0; i < data_samples; i++) {
        results->X_10deg += weightings->xBar[i] * vector[i] * SPD_in->wav_interval;
        results->Y_10deg += weightings->yBar[i] * vector[i] * SPD_in->wav_interval;
        results->Z_10deg += weightings->zBar[i] * vector[i] * SPD_in->wav_interval;
    }
    
    //fill in luminance
    switch (vector_flag)
    {
    	case 1:
    		SPD_in->illuminance = results->Y_10deg * 683.002;
    		break;
    	case 3:
    		SPD_in->luminance = results->Y_10deg * 683.002;
    		break;
    }
    
    if (k == 100) {
        //caculate normalisation factor k such that Y = 100
        k = 100 / results->Y_10deg;
    }

    //use normalisation factor k
    results->X_10deg *= k;
    results->Y_10deg *= k;
    results->Z_10deg *= k;
    
    //caculate CIE 1931 (xyY) values
    results->x_10deg = results->X_10deg / (results->X_10deg + results->Y_10deg + results->Z_10deg);
    results->y_10deg = results->Y_10deg / (results->X_10deg + results->Y_10deg + results->Z_10deg);
    results->z_10deg = results->Z_10deg / (results->X_10deg + results->Y_10deg + results->Z_10deg);
    
    //caculate CIE 1976 (u’v’L*) values
    results->uprime_10deg = (4 * results->X_10deg) / (results->X_10deg + (15 * results->Y_10deg) + (3 * results->Z_10deg));
    results->vprime_10deg = (9 * results->Y_10deg) / (results->X_10deg + (15 * results->Y_10deg) + (3 * results->Z_10deg));
    results->twothirds_vprime_10deg = (2 / 3.0) * results->vprime_10deg;
    if ((results->Y_10deg / 100) > 0.01) {
        results->Lstar_10deg = 116 * pow((results->Y_10deg / 100), ((long double)1/3)) - 16;
    } else {
        results->Lstar_10deg = 903.3 * (results->Y_10deg / 100);
    }

    //free memory
    free(vector);
    free(correct_range_data->wavelength);
    free(correct_range_data);
    free(weightings->xBar);
    free(weightings->yBar);
    free(weightings->zBar);
    free(weightings);
    
    return results;
}

//
//  Calculates CIE 1976 (u’v’L*) CIE 1976 UCS Diagram, CIE 1931 (xyY) colour co-ordinates and CIE 1931 (XYZ) tristimulus values from either the spectral irradiance or spectral radiance data
//  No luminance normalisation takes place. 2 degree.
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - SPD_data* with array length, wavelength and either spectral irradiance or spectral radiance data filled. Wavelength interval must be 1, 2 or 5 nm
//          - CIE_cordinates* for results with memory allocated
//          - int, k normalising constant for caculting CIE 1931 (XYZ) tristimulus values, 1 to perform absolute caculation, 100 to caculate k such that Y = 100 (standard for relitive cases), otherwise enter value such that Y is equal to any convenent number. Note k must be constant when compairing tirstimulus values or chromaticty cordinates.
//          - int flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - CIE_coordinates* with CIE 1976 (u’v’L*) CIE 1976 UCS Diagram, CIE 1931 (xyY) colour co-ordinates and CIE 1931 (XYZ) tristimulus values filled in
//

struct CIE_coordinates* SPDtoCoord2(struct SPD_data *SPD_in, struct CIE_coordinates *results, long double k, int vector_flag) {
    //check arguments
    if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL\n");
        return NULL;
    } else if (SPD_in->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0\n");
        return NULL;
    } else if (SPD_in->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL\n");
        return NULL;
    } else if (SPD_in->spec_irradiance == NULL && SPD_in->spec_radiance == NULL) {
        fprintf(stderr, "ERROR: pointers to spectral data vectors cannot all be NULL\n");
        return NULL;
    } else if (vector_flag != 1 && vector_flag != 3) {
    	fprintf(stderr, "ERROR: vector flag must be 1 or 3\n");
        return NULL;
    } else if (k <= 0) {
    	fprintf(stderr, "ERROR: normalising constant cannot be less than 1\n");
        return NULL;
    }
    
    long double ztest = 0;
    switch (vector_flag)
    {
    	case 1:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_irradiance[i];
    		}
    		break;
    	case 3:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_radiance[i];
    		}
    		break;
    }
    
    if (ztest <= 0) {
    	fprintf(stderr, "ERROR: SPD data must contain values >0\n");
        return NULL;
    }
    
    //open CIE XYZ weighting data of same wavelength interval (1, 2 or 5nm) as SPD_in (xBar, yBar, zBar)
    int data_samples;
    SPD_in->wav_interval = SPD_in->wavelength[1] - SPD_in->wavelength[0];
    FILE *fp;
    switch (SPD_in->wav_interval) {
        case 1:
            data_samples = 401;
            fp = fopen("./data/XYZweightings_data_1nm_intervals_2deg_standard_observer.txt", "r");
            break;
        case 2:
            data_samples = 201;
            fp = fopen("./data/XYZweightings_data_2nm_intervals_2deg_standard_observer.txt", "r");
            break;
        case 5:
            data_samples = 81;
            fp = fopen("./data/XYZweightings_data_5nm_intervals_2deg_standard_observer.txt", "r");
            break;
        default:
            fprintf(stderr, "ERROR: Wavelength interval must be 1, 2 or 5 nm\n");
            return NULL;
    }
    if (fp == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/XYZweightings_data...2deg_standard_observer.txt file.\n");
        exit(EXIT_FAILURE);
    }
    
    //ajust range of data to 380nm - 780nm to match CIE XYZ weightings data
    //cannot be done on normalised data as normalisation will no longer make sense
    struct SPD_data *correct_range_data;
    correct_range_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    correct_range_data->wavelength = (long double*)calloc(data_samples, sizeof(long double));
    long double *vector;
    switch (vector_flag) {
        case 1:
            correct_range_data->spec_irradiance = (long double*)calloc(data_samples, sizeof(long double));
            vector = correct_range_data->spec_irradiance;
            break;
        case 3:
            correct_range_data->spec_radiance = (long double*)calloc(data_samples, sizeof(long double));
            vector = correct_range_data->spec_radiance;
            break;
        default:
            free(correct_range_data->wavelength);
            free(correct_range_data);
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            return NULL;
    }
    ajust_wavlength_range(SPD_in, correct_range_data, 380, 780, vector_flag);
    
    //test new SPD data for values (ensure all data was not cropped off)
    ztest = 0;
    switch (vector_flag)
    {
    	case 1:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_irradiance[i];
    		}
    		break;
    	case 3:
    		for (int i = 0; i < SPD_in->arr_len; i++)
    		{
    			ztest += SPD_in->spec_radiance[i];
    		}
    		break;
    }
    
    if (ztest <= 0) {
    	fprintf(stderr, "ERROR: SPD data must contain values >0\n");
        return NULL;
    }
    
    //read in CIE XYZ weighting data from selected data file (xBar, yBar, zBar)
    struct colourWeighting{
        long double *xBar;
        long double *yBar;
        long double *zBar;
    } *weightings;
    weightings = (struct colourWeighting*)calloc(1, sizeof(struct colourWeighting));
    weightings->xBar = (long double*)calloc(data_samples, sizeof(long double));
    weightings->yBar = (long double*)calloc(data_samples, sizeof(long double));
    weightings->zBar = (long double*)calloc(data_samples, sizeof(long double));
    for (int i = 0; i < data_samples; i++) {
        fscanf(fp, "%Lf %Lf %Lf", &weightings->xBar[i], &weightings->yBar[i], &weightings->zBar[i]);
    }
    fclose(fp);
    
    //caculate XYZ values using CIE XYZ weightings (xBar, yBar, zBar)
    results->X_2deg = 0;
    results->Y_2deg = 0;
    results->Z_2deg = 0;
    for (int i = 0; i < data_samples; i++ ) {
        results->X_2deg += weightings->xBar[i] * vector[i] * SPD_in->wav_interval;
        results->Y_2deg += weightings->yBar[i] * vector[i] * SPD_in->wav_interval;
        results->Z_2deg += weightings->zBar[i] * vector[i] * SPD_in->wav_interval;
    }
    
    //fill in luminance
    switch (vector_flag)
    {
    	case 1:
    		SPD_in->illuminance = results->Y_2deg * 683.002;
    		break;
    	case 3:
    		SPD_in->luminance = results->Y_2deg * 683.002;
    		break;
    }
    
    
    
    if (k == 100) {
        //caculate normalisation factor k such that Y = 100
        k = 100 / results->Y_2deg;
    }
    
    //caculate using normalisation factor k
    results->X_2deg *= k;
    results->Y_2deg *= k;
    results->Z_2deg *= k;
    
    //caculate CIE 1931 (xyY) values
    results->x_2deg = results->X_2deg / (results->X_2deg + results->Y_2deg + results->Z_2deg);
    results->y_2deg = results->Y_2deg / (results->X_2deg + results->Y_2deg + results->Z_2deg);
    results->z_2deg = results->Z_2deg / (results->X_2deg + results->Y_2deg + results->Z_2deg);
    
    //caculate CIE 1976 (u’v’L*) values
    results->uprime_2deg = (4 * results->X_2deg) / (results->X_2deg + (15 * results->Y_2deg) + (3 * results->Z_2deg));
    results->vprime_2deg = (9 * results->Y_2deg) / (results->X_2deg + (15 * results->Y_2deg) + (3 * results->Z_2deg));
    results->twothirds_vprime_2deg = (2 / 3.0) * results->vprime_2deg;
    if ((results->Y_2deg / 100) > 0.01) {
        results->Lstar_2deg = 116 * pow((results->Y_2deg / 100), ((long double)1/3)) - 16;
    } else {
        results->Lstar_2deg = 903.3 * (results->Y_2deg / 100);
    }
    
    //free memory
    free(vector);
    free(correct_range_data->wavelength);
    free(correct_range_data);
    free(weightings->xBar);
    free(weightings->yBar);
    free(weightings->zBar);
    free(weightings);
    
    return results;
}

//
//  Calculates TM-30-15 Rf and Rg values
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - SPD_data* with array length, wavelength and either spectral irradiance or spectral radiance data. Wavelength interval must be 1, 2 or 5 nm
//          - CIE_coordinates* for results with memory allocated
//          - int flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - CIE_coordinates* with TM_Rf and TM_Rg values fileld
//

struct CIE_coordinates* TM_30_15(struct SPD_data *SPD_test, struct CIE_coordinates *results, int vector_flag) {
    //check arguments
    if (SPD_test == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data cannot be NULL for TM-30-15 caculation\n");
        return results;
    } else if (SPD_test->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0 for TM-30-15 caculation\n");
        results->TM_Rf = 0;
        results->TM_Rg = 0;
        return results;
    } else if (SPD_test->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL for TM-30-15 caculation\n");
        results->TM_Rf = 0;
        results->TM_Rg = 0;
        return results;
    } else if (SPD_test->spec_irradiance == NULL && SPD_test->spec_radiance == NULL) {
        fprintf(stderr, "ERROR: pointers to spectral irradiance and spectral radiance vectors cannot both be NULL for TM-30-15 caculation\n");
        results->TM_Rf = 0;
        results->TM_Rg = 0;
        return results;
    } else if (SPD_test->wavelength[0] > 380 || SPD_test->wavelength[SPD_test->arr_len - 1] < 780) {
        printf("%Lf\n", SPD_test->wavelength[700]);
        fprintf(stderr, "ERROR: TM-30-15 caculation not possible on wavelength data without minimum 380nm - 780nm range\n");
        results->TM_Rf = 0;
        results->TM_Rg = 0;
        return results;
    }
    
    //check wavelength interval is 1, 2 or 5nm to match CIE 1964 10deg standard observer
    //open CES file
    SPD_test->wav_interval = (SPD_test->wavelength[1] - SPD_test->wavelength[0]);
    int data_samples = 0;
    FILE *fpCES;
    switch (SPD_test->wav_interval) {
        case 1:
            data_samples = 401;
            fpCES = fopen("./data/CES1nm.csv", "r");
            break;
        case 2:
            data_samples = 201;
            fpCES = fopen("./data/CES2nm.csv", "r");
            break;
        case 5:
            data_samples = 81;
            fpCES = fopen("./data/CES5nm.csv", "r");
            break;
        default:
            fprintf(stderr, "ERROR: Wavelength interval must be 1, 2 or 5 nm\n");
            results->TM_Rf = 0;
            results->TM_Rg = 0;
            return results;
    }
    if (fpCES == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/CESXnm.csv file.\n");
        exit(EXIT_FAILURE);
    }
    
    //count down 2 lines in CES file
    char c = '0';
    int lines = 0;
    while (lines < 2) {
        while (c != '\r') {
            c = fgetc(fpCES);
        }
        c = '0';    //reset char
        lines++;    //count a line
    }
    
    //ajust wavelength range of test source to 380nm - 780nm in order to match CIE 1964 10deg standard observer
    struct SPD_data *correct_range_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    correct_range_data->wavelength = (long double*)calloc(data_samples, sizeof(long double));
    long double *vector;
    switch (vector_flag) {
        case 1:
            correct_range_data->spec_irradiance = (long double*)malloc(sizeof(long double) * data_samples);
            vector = correct_range_data->spec_irradiance;
            break;
        case 3:
            correct_range_data->spec_radiance = (long double*)malloc(sizeof(long double) * data_samples);
            vector = correct_range_data->spec_radiance;
            break;
        default:
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            results->TM_Rf = 0;
            results->TM_Rg = 0;
            return results;
    }
    ajust_wavlength_range(SPD_test, correct_range_data, 380, 780, vector_flag);
    correct_range_data->arr_len = data_samples;
    correct_range_data->wav_interval = correct_range_data->wavelength[1] - correct_range_data->wavelength[0];
    
    //normalise test source to y = 100
    correct_range_data = normaliseSPDtoYonehundred(correct_range_data, vector_flag);
    
    //caculate CCT of test source
    correct_range_data->colour_temp = colourtemp(correct_range_data->wavelength[0], correct_range_data->wavelength[correct_range_data->arr_len - 1], correct_range_data->wav_interval, 1000, 25000, 1, correct_range_data, vector_flag);
    if (correct_range_data->colour_temp == 0) {
        fprintf(stderr, "ERROR: TM-30-15 caculation not possible on SPD without valid CCT <= 25000k\n");
        return NULL;
    }
    
    //caculate illuminant D or blackbody or mixed D & BB refference source based on CCT of test source
    struct SPD_data *SPD_reff;
    if (correct_range_data->colour_temp <= 4500) {
        //refference is blackbody
        SPD_reff = blackbody(380, 780, correct_range_data->colour_temp, correct_range_data->wav_interval);
        //normalise refference SPD so that Y = 100
        SPD_reff = normaliseSPDtoYonehundred(SPD_reff, 3);
    } else if (correct_range_data->colour_temp >= 5500){
        //refference SPD is illuminant D
        SPD_reff = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
        SPD_reff->arr_len = (780 - 380) / correct_range_data->wav_interval + 1;
        SPD_reff->wavelength = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        SPD_reff->spec_radiance = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        SPD_reff->spec_rad_normalised = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        SPD_reff = illuminantD(380, 780, correct_range_data->wav_interval, correct_range_data->colour_temp, SPD_reff);
        //normalise refference SPD so that Y = 100
        SPD_reff = normaliseSPDtoYonehundred(SPD_reff, 3);
    } else {
        //refference SPD is a mixture of blackbody and illuminant D
        //allocate memory for combined result
        SPD_reff = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
        SPD_reff->arr_len = (780 - 380) / correct_range_data->wav_interval + 1;
        SPD_reff->wavelength = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        SPD_reff->spec_radiance = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        SPD_reff->spec_rad_normalised = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        //illuminat D
        struct SPD_data *SPD_reffD = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
        SPD_reffD->arr_len = (780 - 380) / correct_range_data->wav_interval + 1;
        SPD_reffD->wavelength = (long double*)malloc(sizeof(long double) * SPD_reffD->arr_len);
        SPD_reffD->spec_radiance = (long double*)malloc(sizeof(long double) * SPD_reffD->arr_len);
        SPD_reffD->spec_rad_normalised = (long double*)malloc(sizeof(long double) * SPD_reff->arr_len);
        SPD_reffD = illuminantD(380, 780, correct_range_data->wav_interval, correct_range_data->colour_temp, SPD_reffD);
        //blackbody
        struct SPD_data *SPD_reffP;
        SPD_reffP = blackbody(380, 780, correct_range_data->colour_temp, correct_range_data->wav_interval);
        
        //copy wavelength data from blackbody to combined total
        memcpy(SPD_reff->wavelength, SPD_reffP->wavelength, (sizeof(long double) * data_samples - 1));
        
        //combine according to TM-30 recomendations
        for (int i = 0; i < data_samples; i++) {
            SPD_reff->spec_radiance[i] = (((5500 - correct_range_data->colour_temp) / 1000) * SPD_reffP->spec_radiance[i]) + ((1 - ((5500 - correct_range_data->colour_temp) / 1000)) * SPD_reffD->spec_radiance[i]);
        }
        
        //normalise refference SPD so that Y = 100
        SPD_reff = normaliseSPDtoYonehundred(SPD_reff, 3);
        
        free(SPD_reffP->spec_radiance);
        free(SPD_reffP->spec_rad_normalised);
        free(SPD_reffP->wavelength);
        free(SPD_reffP);
        free(SPD_reffD->spec_radiance);
        free(SPD_reffD->spec_rad_normalised);
        free(SPD_reffD->wavelength);
        free(SPD_reffD);
    }
    
    
    /*
     //debug
     printf("\n\nSPD_reffP\n");
     for (int i = 0; i < SPD_reffP->arr_len; i++) {
     printf("%LF %LF\n", SPD_reffP->wavelength[i], SPD_reffP->spec_rad_normalised[i]);
     }
     
     printf("\n\nSPD_reffD\n");
     for (int i = 0; i < SPD_reffD->arr_len; i++) {
     printf("%LF %LF\n", SPD_reffD->wavelength[i], SPD_reffD->spec_rad_normalised[i]);
     }
    //debug print out refference source info
    printf("\ncolour temprature of input to TM-30: %.0Lf\n", correct_range_data->colour_temp);
    struct CIE_coordinates *CIE2 = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    SPDtoCoord2(correct_range_data, CIE2, 1, 3);
    printf("\nCIE u'2 v'2  cordinates of input to TM-30: %Lf %Lf\n", CIE2->uprime_2deg, CIE2->vprime_2deg);
    printf("CIE x'2 y'2 cordinates of input to TM-30: %Lf %Lf\n", CIE2->x_2deg, CIE2->y_2deg);
    printf("CIE X2 Y2 Z2 cordinates of input to TM-30: %Lf %Lf %Lf\n\n", CIE2->X_2deg, CIE2->Y_2deg, CIE2->Z_2deg);
    printf("colour temprature of refference source is: %.0Lf\n", SPD_reff->colour_temp);
    SPDtoCoord2(SPD_reff, CIE2, 1, 3);
    printf("\nCIE u'2 v'2  cordinates of refference: %Lf %Lf\n", CIE2->uprime_2deg, CIE2->vprime_2deg);
    printf("CIE x'2 y'2 cordinates of refference: %Lf %Lf\n", CIE2->x_2deg, CIE2->y_2deg);
    printf("CIE X2 Y2 Z2 cordinates of refference: %Lf %Lf %Lf\n\n", CIE2->X_2deg, CIE2->Y_2deg, CIE2->Z_2deg);
    free(CIE2);
    
    //debug print
    //printf("\nTM-30 refference is:\n");
    //for (int i = 0; i < SPD_reff->arr_len; i++) {
    //    printf("%LF %LF\n", SPD_reff->wavelength[i], SPD_reff->spec_rad_normalised[i]);
    //}
    
    */
    
    //load the 99 CES values
    struct SPD_data **CESSPD;
    CESSPD = (struct SPD_data**)malloc(sizeof(struct SPD_data*) * 99);
    for (int i = 0; i < 99; i++) {
        CESSPD[i] = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
        CESSPD[i]->spec_irradiance = (long double*)malloc(sizeof(long double) * data_samples);
        if (i < 1) {
            CESSPD[i]->wavelength = (long double*)malloc(sizeof(long double) * data_samples);
        } else {
            CESSPD[i]->wavelength = CESSPD[0]->wavelength;
        }
    }
    
    for (int i = 0; i < data_samples; i++) {
        fscanf(fpCES, "%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf,%Lf\r", &CESSPD[0]->wavelength[i], &CESSPD[0]->spec_irradiance[i], &CESSPD[1]->spec_irradiance[i], &CESSPD[2]->spec_irradiance[i], &CESSPD[3]->spec_irradiance[i], &CESSPD[4]->spec_irradiance[i], &CESSPD[5]->spec_irradiance[i], &CESSPD[6]->spec_irradiance[i], &CESSPD[7]->spec_irradiance[i], &CESSPD[8]->spec_irradiance[i], &CESSPD[9]->spec_irradiance[i], &CESSPD[10]->spec_irradiance[i], &CESSPD[11]->spec_irradiance[i], &CESSPD[12]->spec_irradiance[i], &CESSPD[13]->spec_irradiance[i], &CESSPD[14]->spec_irradiance[i], &CESSPD[15]->spec_irradiance[i], &CESSPD[16]->spec_irradiance[i], &CESSPD[17]->spec_irradiance[i], &CESSPD[18]->spec_irradiance[i], &CESSPD[19]->spec_irradiance[i], &CESSPD[20]->spec_irradiance[i], &CESSPD[21]->spec_irradiance[i], &CESSPD[22]->spec_irradiance[i], &CESSPD[23]->spec_irradiance[i], &CESSPD[24]->spec_irradiance[i], &CESSPD[25]->spec_irradiance[i], &CESSPD[26]->spec_irradiance[i], &CESSPD[27]->spec_irradiance[i], &CESSPD[28]->spec_irradiance[i], &CESSPD[29]->spec_irradiance[i], &CESSPD[30]->spec_irradiance[i], &CESSPD[31]->spec_irradiance[i], &CESSPD[32]->spec_irradiance[i], &CESSPD[33]->spec_irradiance[i], &CESSPD[34]->spec_irradiance[i], &CESSPD[35]->spec_irradiance[i], &CESSPD[36]->spec_irradiance[i], &CESSPD[37]->spec_irradiance[i], &CESSPD[38]->spec_irradiance[i], &CESSPD[39]->spec_irradiance[i], &CESSPD[40]->spec_irradiance[i], &CESSPD[41]->spec_irradiance[i], &CESSPD[42]->spec_irradiance[i], &CESSPD[43]->spec_irradiance[i], &CESSPD[44]->spec_irradiance[i], &CESSPD[45]->spec_irradiance[i], &CESSPD[46]->spec_irradiance[i], &CESSPD[47]->spec_irradiance[i], &CESSPD[48]->spec_irradiance[i], &CESSPD[49]->spec_irradiance[i], &CESSPD[50]->spec_irradiance[i], &CESSPD[51]->spec_irradiance[i], &CESSPD[52]->spec_irradiance[i], &CESSPD[53]->spec_irradiance[i], &CESSPD[54]->spec_irradiance[i], &CESSPD[55]->spec_irradiance[i], &CESSPD[56]->spec_irradiance[i], &CESSPD[57]->spec_irradiance[i], &CESSPD[58]->spec_irradiance[i], &CESSPD[59]->spec_irradiance[i], &CESSPD[60]->spec_irradiance[i], &CESSPD[61]->spec_irradiance[i], &CESSPD[62]->spec_irradiance[i], &CESSPD[63]->spec_irradiance[i], &CESSPD[64]->spec_irradiance[i], &CESSPD[65]->spec_irradiance[i], &CESSPD[66]->spec_irradiance[i], &CESSPD[67]->spec_irradiance[i], &CESSPD[68]->spec_irradiance[i], &CESSPD[69]->spec_irradiance[i], &CESSPD[70]->spec_irradiance[i], &CESSPD[71]->spec_irradiance[i], &CESSPD[72]->spec_irradiance[i], &CESSPD[73]->spec_irradiance[i], &CESSPD[74]->spec_irradiance[i], &CESSPD[75]->spec_irradiance[i], &CESSPD[76]->spec_irradiance[i], &CESSPD[77]->spec_irradiance[i], &CESSPD[78]->spec_irradiance[i], &CESSPD[79]->spec_irradiance[i], &CESSPD[80]->spec_irradiance[i], &CESSPD[81]->spec_irradiance[i], &CESSPD[82]->spec_irradiance[i], &CESSPD[83]->spec_irradiance[i], &CESSPD[84]->spec_irradiance[i], &CESSPD[85]->spec_irradiance[i], &CESSPD[86]->spec_irradiance[i], &CESSPD[87]->spec_irradiance[i], &CESSPD[88]->spec_irradiance[i], &CESSPD[89]->spec_irradiance[i], &CESSPD[90]->spec_irradiance[i], &CESSPD[91]->spec_irradiance[i], &CESSPD[92]->spec_irradiance[i], &CESSPD[93]->spec_irradiance[i], &CESSPD[94]->spec_irradiance[i], &CESSPD[95]->spec_irradiance[i], &CESSPD[96]->spec_irradiance[i], &CESSPD[97]->spec_irradiance[i], &CESSPD[98]->spec_irradiance[i]);
    }
    
    //open CIE XYZ weighting data of same wavelength interval (1, 2 or 5nm) as SPD_in (xBar, yBar, zBar)
    FILE *fp;
    switch (data_samples) {
        case 401:
            fp = fopen("./data/XYZweightings_data_1nm_intervals_10deg_standard_observer.txt", "r");
            break;
        case 201:
            fp = fopen("./data/XYZweightings_data_2nm_intervals_10deg_standard_observer.txt", "r");
            break;
        case 81:
            fp = fopen("./data/XYZweightings_data_5nm_intervals_10deg_standard_observer.txt", "r");
            break;
        default:
            fprintf(stderr, "ERROR: Wavelength interval must be 1, 2 or 5 nm\n");
            return NULL;
    }
    if (fp == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/XYZweightings_data...10deg_standard_observer.txt file.\n");
        exit(EXIT_FAILURE);
    }
    
    //read in CIE XYZ weighting data from selected data file (xBar, yBar, zBar)
    struct colourWeighting {
        long double *xBar;
        long double *yBar;
        long double *zBar;
    } *weightings;
    weightings = (struct colourWeighting*)calloc(1, sizeof(struct colourWeighting));
    weightings->xBar = (long double*)malloc(sizeof(long double) * data_samples);
    weightings->yBar = (long double*)malloc(sizeof(long double) * data_samples);
    weightings->zBar = (long double*)malloc(sizeof(long double) * data_samples);
    for (int i = 0; i < data_samples; i++) {
        fscanf(fp, "%Lf %Lf %Lf", &weightings->xBar[i], &weightings->yBar[i], &weightings->zBar[i]);
    }
    fclose(fp);
    
    //setup memory
    struct CIE_coordinates *CEScoordREF, *CEScoordTEST, *testILLUMINANTcoord, *refferenceILLUMINANTcoord;
    CEScoordREF = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    CEScoordTEST = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    refferenceILLUMINANTcoord = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    testILLUMINANTcoord = (struct CIE_coordinates*)calloc(1, sizeof(struct CIE_coordinates));
    long double colourdiff[99];
    struct SPD_data *CESSPDunderREF = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    CESSPDunderREF->wavelength = CESSPD[0]->wavelength;
    CESSPDunderREF->spec_irradiance = (long double*)malloc(sizeof(long double) * data_samples);
    CESSPDunderREF->arr_len = data_samples;
    struct SPD_data *CESSPDunderTEST = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    CESSPDunderTEST->wavelength = CESSPD[0]->wavelength;
    CESSPDunderTEST->spec_irradiance = (long double*)malloc(sizeof(long double) * data_samples);
    CESSPDunderTEST->arr_len = data_samples;
    long double **binab = (long double**)malloc(sizeof(long double*) * 16);
    for (int i = 0; i < 16; i++) {
        binab[i] = (long double*)calloc(99 * 4, sizeof(long double));  //4 * 99 vales to save a' reff, b' reff, a' test, b' test, max of all CES in one bin
    }
    long double *binabAVG = (long double*)calloc(16 * 4, sizeof(long double));   //average a' reff, b' reff, a' test, b' test
    int bincount[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; //count of number of non-zero a'b' values saved into each bin to allow averaging
    long double kt = 0, Xt = 0, Yt = 0, Zt = 0, kr = 0, Xr = 0, Yr = 0, Zr = 0, theta;
    
    //setup data for adapting conditions to caculate CIECAM02
    SPDtoCoord10(SPD_reff, refferenceILLUMINANTcoord, 100, 3);
    SPDtoCoord10(correct_range_data, testILLUMINANTcoord, 100, vector_flag);
    
    //caculate CIECAM02 chromaticity cordinates of every CES with both test and refference sources
    for (int i = 0; i < 99; i++) {
        //caculate 10deg XYZ with refference source
        kr = 0;
        Xr = 0;
        Yr = 0;
        Zr = 0;
        for (int ii = 0; ii < data_samples; ii++) {
            //CESSPDunderREF->spec_irradiance[ii] = SPD_reff->spec_rad_normalised[ii] * CESSPD[i]->spec_irradiance[ii]; //caculate SPD of each CESSPD illuminated by the refference light
            Xr += SPD_reff->spec_rad_normalised[ii] * weightings->xBar[ii] * CESSPD[i]->spec_irradiance[ii];
            //printf("[%d] %Lf += %Lf * %Lf * %Lf\n", ii, Xr, SPD_reff->spec_rad_normalised[ii], weightings->xBar[ii], CESSPD[i]->spec_irradiance[ii]);
            Yr += SPD_reff->spec_rad_normalised[ii] * weightings->yBar[ii] * CESSPD[i]->spec_irradiance[ii];
            //printf("[%d] %Lf += %Lf * %Lf * %Lf\n", ii, Yr, SPD_reff->spec_rad_normalised[ii], weightings->yBar[ii], CESSPD[i]->spec_irradiance[ii]);
            Zr += SPD_reff->spec_rad_normalised[ii] * weightings->zBar[ii] * CESSPD[i]->spec_irradiance[ii];
            //printf("[%d] %Lf += %Lf * %Lf * %Lf\n", ii, Zr, SPD_reff->spec_rad_normalised[ii], weightings->zBar[ii], CESSPD[i]->spec_irradiance[ii]);
            kr += SPD_reff->spec_rad_normalised[ii] * weightings->yBar[ii];
            //printf("[%d] %Lf += %Lf * %Lf\n", ii, kr, SPD_reff->spec_rad_normalised[ii], weightings->yBar[ii]);
        }
        CEScoordREF->X_10deg = (100 * Xr) / kr;
        CEScoordREF->Y_10deg = (100 * Yr) / kr;
        CEScoordREF->Z_10deg = (100 * Zr) / kr;
        //SPDtoCoord10(CESSPDunderREF, CEScoordREF, 1, 1);  //caculate XYZ values for each CESSPD illuminated by the refference light
        
        //debug print current CESSPD
        //printf("CES SPD [%d] reads\n", i);
        //for (int iii = 0; iii < data_samples; iii++) {
        //    printf("%Lf %Lf\n", CESSPD[i]->wavelength[iii], CESSPD[i]->spec_irradiance[iii]);
        //}
        
        //caculate CIE 2002 (J'a'b') CIECAM02 UCS Diagram coordinates of each CES with refference source
        CIECAM02forTM30(CEScoordREF, refferenceILLUMINANTcoord, 10);
        //debug
        //printf("CESSPD[%d] under refference light:\n%Lf X %Lf Y %Lf Z\n%Lf J' %Lf a' %Lf b'\n%Lf hueangle\n\n", i, CEScoordREF->X_10deg, CEScoordREF->Y_10deg, CEScoordREF->Z_10deg, CEScoordREF->Jprime_10deg, CEScoordREF->aprime_10deg, CEScoordREF->bprime_10deg, CEScoordREF->hueangle);
        //printf("CESSPD[%d] hue angle under refference light: %Lf\n", i, CEScoordREF->hueangle);
        
        //caculate 10deg XYZ with test source
        //long double k = 0;
        kt = 0;
        Xt = 0;
        Yt = 0;
        Zt = 0;
        for (int ii = 0; ii < data_samples; ii++) {
            CESSPDunderTEST->spec_irradiance[ii] = correct_range_data->spec_radiance[ii] * CESSPD[i]->spec_irradiance[ii];
            //debug swap back to this line and change normiliseSPDtoYonehundred() to avoid remalloc
            //CESSPDunderTEST->spec_irradiance[ii] = vector[ii] * CESSPD[i]->spec_irradiance[ii];
            //k += correct_range_data->spec_radiance[ii] * weightings->yBar[ii];
            Xt += correct_range_data->spec_radiance[ii] * weightings->xBar[ii] * CESSPD[i]->spec_irradiance[ii];
            //printf("[%d] %Lf += %Lf * %Lf * %Lf\n", ii, Xt, correct_range_data->spec_radiance[ii], weightings->xBar[ii], CESSPD[i]->spec_irradiance[ii]);
            Yt += correct_range_data->spec_radiance[ii] * weightings->yBar[ii] * CESSPD[i]->spec_irradiance[ii];
            //printf("[%d] %Lf += %Lf * %Lf * %Lf\n", ii, Yt, correct_range_data->spec_radiance[ii], weightings->yBar[ii], CESSPD[i]->spec_irradiance[ii]);
            Zt += correct_range_data->spec_radiance[ii] * weightings->zBar[ii] * CESSPD[i]->spec_irradiance[ii];
            //printf("[%d] %Lf += %Lf * %Lf * %Lf\n", ii, Zt, correct_range_data->spec_radiance[ii], weightings->zBar[ii], CESSPD[i]->spec_irradiance[ii]);
            kt += correct_range_data->spec_radiance[ii] * weightings->yBar[ii];
            //printf("[%d] %Lf += %Lf * %Lf\n", ii, kt, correct_range_data->spec_radiance[ii], weightings->yBar[ii]);
        }
        //k /= 100;
        //SPDtoCoord10(CESSPDunderTEST, CEScoordTEST, k, 1);
        CEScoordTEST->X_10deg = 100 * Xt / kt;
        CEScoordTEST->Y_10deg = 100 * Yt / kt;
        CEScoordTEST->Z_10deg = 100 * Zt / kt;
        
        //caculate CIE 2002 (J'a'b') CIECAM02 UCS Diagram cordinates of each CES with test source
        CIECAM02forTM30(CEScoordTEST, testILLUMINANTcoord, 10);
        //printf("CESSPD[%d] under test light:\n%Lf X %Lf Y %Lf Z\n%Lf J' %Lf a' %Lf b'\n\n", i, CEScoordTEST->X_10deg, CEScoordTEST->Y_10deg, CEScoordTEST->Z_10deg, CEScoordTEST->Jprime_10deg, CEScoordTEST->aprime_10deg, CEScoordTEST->bprime_10deg);
        
        //caculate colour difference for each CES (euclidian distance)
        colourdiff[i] = sqrtl(pow(CEScoordTEST->Jprime_10deg - CEScoordREF->Jprime_10deg, 2) + pow(CEScoordTEST->aprime_10deg - CEScoordREF->aprime_10deg, 2) + pow(CEScoordTEST->bprime_10deg - CEScoordREF->bprime_10deg, 2));
        //printf("colourdiff[%d] = %Lf\n\n", i, colourdiff[i]);
        fflush(stdin);
        
        //assign to hue angle bin based on chromaticity cordinates under refference light
        if (CEScoordREF->hueangle < 22.5) {
            binab[0][bincount[0]] = CEScoordREF->aprime_10deg;
            binab[0][bincount[0]] = CEScoordREF->bprime_10deg;
            binab[0][bincount[0]] = CEScoordTEST->aprime_10deg;
            binab[0][bincount[0]] = CEScoordTEST->bprime_10deg;
            binabAVG[0] += CEScoordREF->aprime_10deg;
            binabAVG[1] += CEScoordREF->bprime_10deg;
            binabAVG[2] += CEScoordTEST->aprime_10deg;
            binabAVG[3] += CEScoordTEST->bprime_10deg;
            bincount[0]++;
        } else if (CEScoordREF->hueangle < 45.0) {
            binab[1][bincount[1]] = CEScoordREF->aprime_10deg;
            binab[1][bincount[1]] = CEScoordREF->bprime_10deg;
            binab[1][bincount[1]] = CEScoordTEST->aprime_10deg;
            binab[1][bincount[1]] = CEScoordTEST->bprime_10deg;
            binabAVG[4] += CEScoordREF->aprime_10deg;
            binabAVG[5] += CEScoordREF->bprime_10deg;
            binabAVG[6] += CEScoordTEST->aprime_10deg;
            binabAVG[7] += CEScoordTEST->bprime_10deg;
            bincount[1]++;
        } else if (CEScoordREF->hueangle < 67.5) {
            binab[2][bincount[2]] = CEScoordREF->aprime_10deg;
            binab[2][bincount[2]] = CEScoordREF->bprime_10deg;
            binab[2][bincount[2]] = CEScoordTEST->aprime_10deg;
            binab[2][bincount[2]] = CEScoordTEST->bprime_10deg;
            binabAVG[8] += CEScoordREF->aprime_10deg;
            binabAVG[9] += CEScoordREF->bprime_10deg;
            binabAVG[10] += CEScoordTEST->aprime_10deg;
            binabAVG[11] += CEScoordTEST->bprime_10deg;
            bincount[2]++;
        } else if (CEScoordREF->hueangle < 90.0) {
            binab[3][bincount[3]] = CEScoordREF->aprime_10deg;
            binab[3][bincount[3]] = CEScoordREF->bprime_10deg;
            binab[3][bincount[3]] = CEScoordTEST->aprime_10deg;
            binab[3][bincount[3]] = CEScoordTEST->bprime_10deg;
            binabAVG[12] += CEScoordREF->aprime_10deg;
            binabAVG[13] += CEScoordREF->bprime_10deg;
            binabAVG[14] += CEScoordTEST->aprime_10deg;
            binabAVG[15] += CEScoordTEST->bprime_10deg;
            bincount[3]++;
        } else if (CEScoordREF->hueangle < 112.5) {
            binab[4][bincount[4]] = CEScoordREF->aprime_10deg;
            binab[4][bincount[4]] = CEScoordREF->bprime_10deg;
            binab[4][bincount[4]] = CEScoordTEST->aprime_10deg;
            binab[4][bincount[4]] = CEScoordTEST->bprime_10deg;
            binabAVG[16] += CEScoordREF->aprime_10deg;
            binabAVG[17] += CEScoordREF->bprime_10deg;
            binabAVG[18] += CEScoordTEST->aprime_10deg;
            binabAVG[19] += CEScoordTEST->bprime_10deg;
            bincount[4]++;
        } else if (CEScoordREF->hueangle < 135.0) {
            binab[5][bincount[5]] = CEScoordREF->aprime_10deg;
            binab[5][bincount[5]] = CEScoordREF->bprime_10deg;
            binab[5][bincount[5]] = CEScoordTEST->aprime_10deg;
            binab[5][bincount[5]] = CEScoordTEST->bprime_10deg;
            binabAVG[20] += CEScoordREF->aprime_10deg;
            binabAVG[21] += CEScoordREF->bprime_10deg;
            binabAVG[22] += CEScoordTEST->aprime_10deg;
            binabAVG[23] += CEScoordTEST->bprime_10deg;
            bincount[5]++;
        } else if (CEScoordREF->hueangle < 157.5) {
            binab[6][bincount[6]] = CEScoordREF->aprime_10deg;
            binab[6][bincount[6]] = CEScoordREF->bprime_10deg;
            binab[6][bincount[6]] = CEScoordTEST->aprime_10deg;
            binab[6][bincount[6]] = CEScoordTEST->bprime_10deg;
            binabAVG[24] += CEScoordREF->aprime_10deg;
            binabAVG[25] += CEScoordREF->bprime_10deg;
            binabAVG[26] += CEScoordTEST->aprime_10deg;
            binabAVG[27] += CEScoordTEST->bprime_10deg;
            bincount[6]++;
        } else if (CEScoordREF->hueangle < 180.0) {
            binab[7][bincount[7]] = CEScoordREF->aprime_10deg;
            binab[7][bincount[7]] = CEScoordREF->bprime_10deg;
            binab[7][bincount[7]] = CEScoordTEST->aprime_10deg;
            binab[7][bincount[7]] = CEScoordTEST->bprime_10deg;
            binabAVG[28] += CEScoordREF->aprime_10deg;
            binabAVG[29] += CEScoordREF->bprime_10deg;
            binabAVG[30] += CEScoordTEST->aprime_10deg;
            binabAVG[31] += CEScoordTEST->bprime_10deg;
            bincount[7]++;
        } else if (CEScoordREF->hueangle < 202.5) {
            binab[8][bincount[8]] = CEScoordREF->aprime_10deg;
            binab[8][bincount[8]] = CEScoordREF->bprime_10deg;
            binab[8][bincount[8]] = CEScoordTEST->aprime_10deg;
            binab[8][bincount[8]] = CEScoordTEST->bprime_10deg;
            binabAVG[32] += CEScoordREF->aprime_10deg;
            binabAVG[33] += CEScoordREF->bprime_10deg;
            binabAVG[34] += CEScoordTEST->aprime_10deg;
            binabAVG[35] += CEScoordTEST->bprime_10deg;
            bincount[8]++;
        } else if (CEScoordREF->hueangle < 225.0) {
            binab[9][bincount[9]] = CEScoordREF->aprime_10deg;
            binab[9][bincount[9]] = CEScoordREF->bprime_10deg;
            binab[9][bincount[9]] = CEScoordTEST->aprime_10deg;
            binab[9][bincount[9]] = CEScoordTEST->bprime_10deg;
            binabAVG[36] += CEScoordREF->aprime_10deg;
            binabAVG[37] += CEScoordREF->bprime_10deg;
            binabAVG[38] += CEScoordTEST->aprime_10deg;
            binabAVG[39] += CEScoordTEST->bprime_10deg;
            bincount[9]++;
        } else if (CEScoordREF->hueangle < 247.5) {
            binab[10][bincount[10]] = CEScoordREF->aprime_10deg;
            binab[10][bincount[10]] = CEScoordREF->bprime_10deg;
            binab[10][bincount[10]] = CEScoordTEST->aprime_10deg;
            binab[10][bincount[10]] = CEScoordTEST->bprime_10deg;
            binabAVG[40] += CEScoordREF->aprime_10deg;
            binabAVG[41] += CEScoordREF->bprime_10deg;
            binabAVG[42] += CEScoordTEST->aprime_10deg;
            binabAVG[43] += CEScoordTEST->bprime_10deg;
            bincount[10]++;
        } else if (CEScoordREF->hueangle < 270.0) {
            binab[11][bincount[11]] = CEScoordREF->aprime_10deg;
            binab[11][bincount[11]] = CEScoordREF->bprime_10deg;
            binab[11][bincount[11]] = CEScoordTEST->aprime_10deg;
            binab[11][bincount[11]] = CEScoordTEST->bprime_10deg;
            binabAVG[44] += CEScoordREF->aprime_10deg;
            binabAVG[45] += CEScoordREF->bprime_10deg;
            binabAVG[46] += CEScoordTEST->aprime_10deg;
            binabAVG[47] += CEScoordTEST->bprime_10deg;
            bincount[11]++;
        } else if (CEScoordREF->hueangle < 292.5) {
            binab[12][bincount[12]] = CEScoordREF->aprime_10deg;
            binab[12][bincount[12]] = CEScoordREF->bprime_10deg;
            binab[12][bincount[12]] = CEScoordTEST->aprime_10deg;
            binab[12][bincount[12]] = CEScoordTEST->bprime_10deg;
            binabAVG[48] += CEScoordREF->aprime_10deg;
            binabAVG[49] += CEScoordREF->bprime_10deg;
            binabAVG[50] += CEScoordTEST->aprime_10deg;
            binabAVG[51] += CEScoordTEST->bprime_10deg;
            bincount[12]++;
        } else if (CEScoordREF->hueangle < 315.0) {
            binab[13][bincount[13]] = CEScoordREF->aprime_10deg;
            binab[13][bincount[13]] = CEScoordREF->bprime_10deg;
            binab[13][bincount[13]] = CEScoordTEST->aprime_10deg;
            binab[13][bincount[13]] = CEScoordTEST->bprime_10deg;
            binabAVG[52] += CEScoordREF->aprime_10deg;
            binabAVG[53] += CEScoordREF->bprime_10deg;
            binabAVG[54] += CEScoordTEST->aprime_10deg;
            binabAVG[55] += CEScoordTEST->bprime_10deg;
            bincount[13]++;
        } else if (CEScoordREF->hueangle < 337.5) {
            binab[14][bincount[14]] = CEScoordREF->aprime_10deg;
            binab[14][bincount[14]] = CEScoordREF->bprime_10deg;
            binab[14][bincount[14]] = CEScoordTEST->aprime_10deg;
            binab[14][bincount[14]] = CEScoordTEST->bprime_10deg;
            binabAVG[56] += CEScoordREF->aprime_10deg;
            binabAVG[57] += CEScoordREF->bprime_10deg;
            binabAVG[58] += CEScoordTEST->aprime_10deg;
            binabAVG[59] += CEScoordTEST->bprime_10deg;
            bincount[14]++;
        } else if (CEScoordREF->hueangle <= 360.0) {
            binab[15][bincount[15]] = CEScoordREF->aprime_10deg;
            binab[15][bincount[15]] = CEScoordREF->bprime_10deg;
            binab[15][bincount[15]] = CEScoordTEST->aprime_10deg;
            binab[15][bincount[15]] = CEScoordTEST->bprime_10deg;
            binabAVG[60] += CEScoordREF->aprime_10deg;
            binabAVG[61] += CEScoordREF->bprime_10deg;
            binabAVG[62] += CEScoordTEST->aprime_10deg;
            binabAVG[63] += CEScoordTEST->bprime_10deg;
            bincount[15]++;
        } else {
            fprintf(stderr, "ERROR: hue angle value not valid\n");
            return NULL;
        }
    }
    free(CESSPDunderREF->spec_irradiance);
    free(CESSPDunderREF);
    free(CESSPDunderTEST->spec_irradiance);
    free(CESSPDunderTEST);
    free(CESSPD[0]->wavelength);
    for (int i = 0; i < 99; i++) {
        free(CESSPD[i]->spec_irradiance);
        free(CESSPD[i]);
    }
    free(CESSPD);
    free(CEScoordREF);
    free(CEScoordTEST);
    free(refferenceILLUMINANTcoord);
    free(testILLUMINANTcoord);
    fclose(fpCES);
    //free(vector);
    free(correct_range_data->spec_rad_normalised);
    free(correct_range_data->spec_radiance);
    free(correct_range_data->wavelength);
    free(correct_range_data);
    free(weightings->xBar);
    free(weightings->yBar);
    free(weightings->zBar);
    free(weightings);
    
    //caculate fidelity index Rf as arithmetic mean of colour differnces for each CES scaled between 0 and 100
    long double Rf = 0;
    for (int i = 0; i < 99; i++) {
        Rf += (1.0 / 99) * colourdiff[i];
    }
    long double Rfprime = 100 - 7.54 * Rf;
    results->TM_Rf = 10 * log(pow(M_E, Rfprime / 10) + 1);
    
    //average a' and b' of test and refference in each bin
    for (int i = 0; i < (16 * 4); i += 4) {
        if (i != 63) {
            for (int ii = 0; ii < 4; ii++) {
                binabAVG[i + ii] /= bincount[(int)floor(i / 4)];
            }
        }
    }
    
    //caculate the area enclosed by the 16 averaged chromaticity cordinates for refference & test source
    long double DIFFaprimeREF[16], MEANbprimeREF[16], DIFFaprimeTEST[16], MEANbprimeTEST[16];
    long double areaREF = 0, areaTEST = 0;
    
    for (int i = 0; i < 16; i++) {
        if (i == 0) {   //first loop
            DIFFaprimeREF[i] = binabAVG[i + 4] - binabAVG[i];
            //printf("DIFFaprimeREF[%d] %Lf = binabAVG[%d] %Lf - binabAVG[%d] %Lf\n", i, DIFFaprimeREF[i], i+4, binabAVG[i+4], i, binabAVG[i]);
            MEANbprimeREF[i] = (binabAVG[i + 1 + 4] + binabAVG[i + 1]) / 2;
            DIFFaprimeTEST[i] = binabAVG[i + 2 + 4] - binabAVG[i + 2];
            MEANbprimeTEST[i] = (binabAVG[i + 3 + 4] + binabAVG[i + 3]) / 2;
        } else if (i != 15) {   //mid loops
            DIFFaprimeREF[i] = binabAVG[(i + 1) * 4] - binabAVG[i * 4];
            //printf("DIFFaprimeREF[%d] %Lf = binabAVG[%d] %Lf - binabAVG[%d] %Lf\n", i, DIFFaprimeREF[i], (i + 1) * 4, binabAVG[(i + 1) * 4], i * 4, binabAVG[i * 4]);
            MEANbprimeREF[i] = (binabAVG[(i + 1) * 4 + 1] + binabAVG[(i * 4) + 1]) / 2;
            DIFFaprimeTEST[i] = binabAVG[(i + 1) * 4 + 2] - binabAVG[(i * 4) + 2];
            MEANbprimeTEST[i] = (binabAVG[(i + 1) * 4 + 3 ] + binabAVG[(i * 4) + 3]) / 2;
        } else {    //final loop
            DIFFaprimeREF[i] = binabAVG[0] - binabAVG[i * 4];
            //printf("DIFFaprimeREF[%d] %Lf = binabAVG[%d] %Lf - binabAVG[%d] %Lf\n", i, DIFFaprimeREF[i], 0, binabAVG[0], i * 4, binabAVG[i * 4]);
            MEANbprimeREF[i] = (binabAVG[1] + binabAVG[(i * 4) + 1]) / 2;
            DIFFaprimeTEST[i] = binabAVG[2] - binabAVG[(i * 4) + 2];
            MEANbprimeTEST[i] = (binabAVG[3] + binabAVG[(i * 4) + 3]) / 2;
        }
        areaREF += DIFFaprimeREF[i] * MEANbprimeREF[i];
        areaTEST += DIFFaprimeTEST[i] * MEANbprimeTEST[i];
    }
    
    //caculate gaumit index Rg
    results->TM_Rg = 100 * (areaTEST / areaREF);
    
    //free
    free(SPD_reff->spec_radiance);
    free(SPD_reff->spec_rad_normalised);
    free(SPD_reff->wavelength);
    free(SPD_reff);
    free(binabAVG);
    for (int i = 0; i < 16; i++) {
        free(binab[i]);
    }
    free(binab);
    return results;
}

//
//  Calculates CIE Illuminant D as normalised spectral radiance at 5nm intervals between 300nm - 830nm
//  CIE guidelines state that other intervals should be caculated by liner interpolation
//
//  EXTERNAL REQUIREMENTS:
//  INPUTS: - int, minimum wavelength to report
//          - int, maximum wavelength to report
//          - int, interval of data entrys
//          - int, temprature of SPD to be returned min 4000 - max 25,000 (kelvin)
//          - SPD_data* SPD_out, with wavelength, radiance and normalised spectral radiance malloc'ed to (w_max - w_min) / interval + 1
//  OUTPUT: - SPD_data* SPD_out with normlaied spectral radiance per wavelength at 5nm intervals filled (W·sr−1·m−2·nm-1))
//

struct SPD_data* illuminantD(int w_min, int w_max, int interval, int temprature, struct SPD_data* SPD_out){
    //check input
    if (temprature < 4000 || temprature > 25000) {
        fprintf(stderr, "ERROR: Illuminant D can only be caculated for 4000K - 25000K\n");
        SPD_out->colour_temp = 0;
        return SPD_out;
    }
    
    //allocate memory for CIE caculation
    struct SPD_data *illD = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    illD->arr_len = (830 - 300) / 5 + 1;
    illD->wavelength = (long double*)malloc(sizeof(long double) * illD->arr_len);
    illD->spec_radiance = (long double*)malloc(sizeof(long double) * illD->arr_len);
    
    //caculate CIE 1931 (xyY) chromaticity cordinates for illuminant D at given CCT
    struct CIE_coordinates CIEcoord;
    if (temprature <= 7000) {
        CIEcoord.x_2deg = (-4.607E9 / pow(temprature, 3)) + (2.9678E6 / pow(temprature, 2)) + (0.09911E3 / temprature) + 0.244063;
        CIEcoord.y_2deg = (-3.0 * pow(CIEcoord.x_2deg, 2)) + (2.87 * CIEcoord.x_2deg) - 0.275;
    } else if (temprature > 7000){
        CIEcoord.x_2deg = (-2.0064E9 / pow(temprature, 3)) + (1.9018E6 / pow(temprature, 2)) + (0.24748E3 / temprature) + 0.23704;
        CIEcoord.y_2deg = (-3.0 * pow(CIEcoord.x_2deg, 2)) + (2.87 * CIEcoord.x_2deg) - 0.275;
    }
    
    //open S components file
    FILE *sval = fopen("./data/illuminantD_Scomponents.csv", "r");
    if (sval == NULL) {
        fprintf(stderr, "ERROR: Could not open ./data/illuminantD_Scomponents2.csv");
        return NULL;
    }
    //skip header
    char c = '0';
    while (c != '\r') {
        c = fgetc(sval);
    }
    
    //caculate SPD
    long double M1 = (-1.3515 - (1.7703 * CIEcoord.x_2deg) + (5.9114 * CIEcoord.y_2deg)) / (0.0241 + (0.2562 * CIEcoord.x_2deg) - (0.7341 * CIEcoord.y_2deg));
    long double M2 = (0.03 - (31.4424 * CIEcoord.x_2deg) + (30.0717 * CIEcoord.y_2deg)) / (0.0241 + (0.2562 * CIEcoord.x_2deg) - (0.7341 * CIEcoord.y_2deg));
    long double *sval_wavelength = (long double*)calloc(illD->arr_len, sizeof(long double));
    long double *s0 = (long double*)malloc(sizeof(long double) * illD->arr_len);
    long double *s1 = (long double*)malloc(sizeof(long double) * illD->arr_len);
    long double *s2 = (long double*)malloc(sizeof(long double) * illD->arr_len);
    for (int i = 0; i < illD->arr_len; i++) {
        fscanf(sval, "%Lf,%Lf,%Lf,%Lf", &sval_wavelength[i], &s0[i], &s1[i], &s2[i]);
        illD->wavelength[i] = sval_wavelength[i];
        illD->spec_radiance[i] = s0[i] + (M1 * s1[i]) + (M2 * s2[i]);
    }
    fclose(sval);
    free(sval_wavelength);
    free(s0);
    free(s1);
    free(s2);
    
    ajust_wavlength_range(illD, SPD_out, w_min, w_max, 3);
    free(illD->spec_radiance);
    free(illD->wavelength);
    free(illD);
    
    //ajust wavelength interval
    if (interval != 5) {
        struct SPD_data *newSPD_out = calloc(1, sizeof(struct SPD_data));
        newSPD_out->arr_len = (w_max - w_min) / interval + 1;
        newSPD_out->wavelength = malloc(sizeof(long double) * newSPD_out->arr_len);
        newSPD_out->spec_radiance = malloc(sizeof(long double) * newSPD_out->arr_len);
        newSPD_out->spec_rad_normalised = malloc(sizeof(long double) * newSPD_out->arr_len);
        ajust_wavelength_interval(SPD_out, newSPD_out, interval, 3);
        free(SPD_out->wavelength);
        free(SPD_out->spec_radiance);
        free(SPD_out->spec_rad_normalised);
        free(SPD_out);
        
        //normalise
        normaliseSPDtoONE(newSPD_out->spec_radiance, newSPD_out->spec_rad_normalised, newSPD_out->arr_len);
        newSPD_out->wav_interval = newSPD_out->wavelength[1] - newSPD_out->wavelength[0];
        
        return newSPD_out;
    }
    
    //normalise
    normaliseSPDtoONE(SPD_out->spec_radiance, SPD_out->spec_rad_normalised, SPD_out->arr_len);
    SPD_out->wav_interval = SPD_out->wavelength[1] - SPD_out->wavelength[0];
    
    return SPD_out;
}

//
//  Ajust the wavelength range of an SPD data vector such that data between the selected range is returned.
//  If range is larger than data supplied, 0s are padded. Only spec_irrad or spec_rad values are returned. Spec_irrad is default if avalible.
//  The wavelength intervals will match.
//
//  EXTERNAL REQUIREMENTS:
//  INPUTS: - struct SPD_data* with wavelength, arr_len and a spectral value filled in
//          - struct SPD_data* with memory allocated for the same values
//          - min wavelength for new range in nm
//          - max wavelength for new range in nm
//          - int flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - struct SPD_data* with wavelength, arr_len and spec_irrad or spec_rad values changed to suit new range
//
struct SPD_data* ajust_wavlength_range(struct SPD_data* SPD_in, struct SPD_data* SPD_out, int minwave, int maxwave, int vector_flag) {
    SPD_out->wav_interval = SPD_in->wavelength[1] - SPD_in->wavelength[0];
    SPD_out->arr_len = ((maxwave - minwave) / SPD_out->wav_interval) + 1;
    //check arguments
    if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL when changing wavelength range\n");
        exit(EXIT_FAILURE);
    } else if (SPD_in->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0 when changing wavelength range\n");
        exit(EXIT_FAILURE);
    } else if (SPD_in->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL when changing wavelength range\n");
        exit(EXIT_FAILURE);
    } else if (SPD_in->spec_irradiance == NULL && SPD_in->spec_radiance == NULL) {
        fprintf(stderr, "ERROR: pointers to spectral irradiance and spectral radiance vectors cannot both be NULL when changing wavelength range\n");
        exit(EXIT_FAILURE);
    } else if (maxwave - minwave <= 0) {
        fprintf(stderr, "ERROR: cannot change wavelength range to range where maximum wavelength is less than minimum wavelength\n");
        exit(EXIT_FAILURE);
    }
    
    long double *invector, *outvector;
    switch (vector_flag) {
        case 1:
            invector = SPD_in->spec_irradiance;
            outvector = SPD_out->spec_irradiance;
            break;
        case 3:
            invector = SPD_in->spec_radiance;
            outvector = SPD_out->spec_radiance;
            break;
        default:
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            exit(EXIT_FAILURE);
    }
    
    //find locations of minimum and maximum wavelengths, if data is shorter, default value is 0
    int read_index = 0, write_index = 0, last_read_index = 0;
    for (int i = 0; i < SPD_in->arr_len; i++) {
        if (SPD_in->wavelength[i] >= minwave) {
            read_index = i;
            write_index = (SPD_in->wavelength[i] - minwave) / SPD_out->wav_interval;
            if (write_index > SPD_out->arr_len) {
                fprintf(stderr, "ERROR: SPD range cannot be changed as wavelength ranges do not overlap\n");
                exit(EXIT_FAILURE);
            }
            goto nestedbreak1;
        } else if (i == (maxwave - minwave)) {
            fprintf(stderr, "ERROR: SPD range cannot be changed as wavelength ranges do not overlap\n");
            exit(EXIT_FAILURE);
        }
    }
    nestedbreak1:;
    for (int i = 0; i < SPD_in->arr_len; i++) {
        if (SPD_in->wavelength[i] > maxwave) {
            last_read_index = i - 1;
            goto nestedbreak2;
        } else if (i == SPD_in->arr_len - 1) {
            last_read_index = SPD_in->arr_len - 1;
        }
    }
    nestedbreak2:;
    
    //pull out the spec_irradiance or spec_radiance values between >=minwave in nm and <=maxwave in nm
    int ii = 0;
    for (int i = 0; i < SPD_out->arr_len; i++) {
        if (i < write_index) {
            outvector[i] = 0;
        } else if (i <= write_index + (last_read_index - read_index)) {
            outvector[i] = invector[read_index + ii];
            ii++;
        } else if (i > write_index + (last_read_index - read_index)) {
            outvector[i] = 0;
        }
        SPD_out->wavelength[i] = minwave + (SPD_out->wav_interval * i);
    }
    SPD_out->spec_radiance = outvector;
    
    return SPD_out;
}

//
//  Ajusts the wavelength interval of data supplied, starting from the same wavelength
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - struct SPD_data* SPD_in, with wavelength, arr_len and a spectral value filled in, wavelength intervals are assumed to be regular
//          - struct SPD_data* SPD_out, with wavelength, arr_len and the same spectral vector malloc'ed for new size
//          - int interval, new wavelength interval
//          - int vector_flag, indicating which vector to ajust. 1 = spectral_irradiance, 2 = spec_irrad_normalised, 3 = spectral_radiance, 4 = spec_irrad_normalised
//  OUTPUT: - struct SPD_DATA*, with wavelength, arr_len and spec_irrad or spec_rad values changed to suit new interval
//
struct SPD_data* ajust_wavelength_interval(struct SPD_data *SPD_in, struct SPD_data *SPD_out, int interval, int vector_flag){
    if (interval < 1) {
        fprintf(stderr, "ERROR: Not possible to have interval of less than 1\n");
        return SPD_out;
    }
    
    //use irradiance or radiance depending on mode
    long double *invector, *outvector;
    switch (vector_flag) {
        case 1:
            invector = SPD_in->spec_irradiance;
            outvector = SPD_out->spec_irradiance;
            break;
        case 2:
            invector = SPD_in->spec_irrad_normalised;
            outvector = SPD_out->spec_irrad_normalised;
            break;
        case 3:
            invector = SPD_in->spec_radiance;
            outvector = SPD_out->spec_radiance;
            break;
        case 4:
            invector = SPD_in->spec_rad_normalised;
            outvector = SPD_out->spec_rad_normalised;
            break;
        default:
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            exit(EXIT_FAILURE);
    }
    
    SPD_in->wav_interval = SPD_in->wavelength[1] - SPD_in->wavelength[0];
    
    //fill in wavelength for SPD_out
    SPD_out->wav_interval = interval;
    SPD_out->wavelength[0] = SPD_in->wavelength[0];		//initial wavelength is the same
    for(int i = 1; i < SPD_out->arr_len; i++) {
    	SPD_out->wavelength[i] = SPD_out->wavelength[i - 1] + interval;
	}
    
    
    //if intervals match, return same values
    if (interval == SPD_in->wav_interval) {
        for (int i = 0; i < SPD_in->arr_len; i++) {
            outvector[i] = invector[i];
        }
        SPD_out->wav_interval = SPD_in->wav_interval;
        return SPD_out;
    //if desired interval is greater, use liner interpolation to downsample
    } else if (interval > SPD_in->wav_interval)
    {	
        int count = 0;		//SPD_out array element index
        for (int i = 0; i < SPD_in->arr_len; i++)		//step through output
        {
        	if (i != SPD_in->arr_len && count != SPD_out->arr_len)				//if not last wavelength
        	{
        		if(SPD_in->wavelength[i] == SPD_out->wavelength[count])		//if wavelengths match
        		{
        			outvector[count] = invector[i];		//copy data from in to out
        			count++;						//increment count
        		} else if(SPD_out->wavelength[count] - SPD_in->wavelength[i] <= SPD_in->wav_interval)		//if SPD_in->wavelength is within one SPD_in interval away from desired SPD_out wavelength
        		{
        			//liner interpolation
                    long double xa = SPD_in->wavelength[i];
                    long double xb = SPD_in->wavelength[i + 1];
                    long double ya = invector[i];
                    long double yb = invector[i + 1];
                    outvector[count] = (ya * (xb - SPD_out->wavelength[count]) + yb * (SPD_out->wavelength[count] - xa)) / (xb - xa);
                    count++;
        		}
        	}
        }
    //if desired interval is smaller, use liner interpolation to upsample
    } else if (interval < SPD_in->wav_interval) {
        int count = 0;		//SPD_out array element index
        for (int i = 0; i < SPD_in->arr_len; i++)		//step through input
        {
        	if (i != SPD_in->arr_len)				//if not last wavelength
        	{
        		if(SPD_in->wavelength[i] == SPD_out->wavelength[count])		//if wavelengths match at start
        		{
        			outvector[count] = invector[i];		//copy data from in to out
        			count++;						//increment count
        		} else
        		{
        			for (int ii = 0; ii < floor(SPD_in->wav_interval / (double)interval); ii++)	//interpolate a number of times
        			{
                		//liner interpolation
                    	long double xa = SPD_in->wavelength[i - 1];
                    	long double xb = SPD_in->wavelength[i];
                    	long double ya = invector[i - 1];
                    	long double yb = invector[i];
                    	outvector[count] = (ya * (xb - SPD_out->wavelength[count]) + yb * (SPD_out->wavelength[count] - xa)) / (xb - xa);
                    	count++;
                	}
                	if (count < SPD_out->arr_len)		//if not last wavelength
                	{
                		//check for matching wavelength at end
                		if(SPD_in->wavelength[i] == SPD_out->wavelength[count])		//if wavelengths match
        				{
        					outvector[count] = invector[i];		//copy data from in to out
        					count++;						//increment count
        				}
        			}
        		}
        	}
        }
    }
    
    return SPD_out;
}

//
//  Caculate CIECAM02 Cordinates as required by TM-30-15 specification from 10 deg XYZ values
//
//  EXTERNAL REQUIREMENTS: NONE
//  INPUTS: - struct CEScoordREF* with X_10deg, Y_10deg, Z_10deg values filled in
//          - int flag, currently unused
//  OUTPUT: - CIE_coordinates* CEScoordREF, with CIE 2002 (J'a'b') 10 deg values filled
//
struct CIE_coordinates* CIECAM02forTM30(struct CIE_coordinates *CEScoordREF, struct CIE_coordinates *adapting_illuminant_coord, int flag) {
    //static variables provided by TM-30
    static float c = 0.69, Fl = 0.7937, n = 0.2, Nbb = 1.0003, z = 1.9272, Ncb = 1.0003, Nc = 1, F = 1, La = 100, D = 1;
    
    //adopted white is test or refference illuminant respectivly
    long double Rw = (0.7328 * adapting_illuminant_coord->X_10deg) + (0.4296 * adapting_illuminant_coord->Y_10deg) + (-0.1624 * adapting_illuminant_coord->Z_10deg);
    long double Gw = (-0.7036 * adapting_illuminant_coord->X_10deg) + (1.6975 * adapting_illuminant_coord->Y_10deg) + (0.0061 * adapting_illuminant_coord->Z_10deg);
    long double Bw = (0.0030 * adapting_illuminant_coord->X_10deg) + (0.0136 * adapting_illuminant_coord->Y_10deg) + (0.9834 * adapting_illuminant_coord->Z_10deg);
    
    //CES sample under test or refference illuminant
    long double R = (0.7328 * CEScoordREF->X_10deg) + (0.4296 * CEScoordREF->Y_10deg) + (-0.1624 * CEScoordREF->Z_10deg);
    long double G = (-0.7036 * CEScoordREF->X_10deg) + (1.6975 * CEScoordREF->Y_10deg) + (0.0061 * CEScoordREF->Z_10deg);
    long double B = (0.0030 * CEScoordREF->X_10deg) + (0.0136 * CEScoordREF->Y_10deg) + (0.9834 * CEScoordREF->Z_10deg);
    
    //applying a chromatic adaption transformation, the corresponding colour under the illuminant is
    long double Rc = (100 / Rw) * R;
    long double Gc = (100 / Gw) * G;
    long double Bc = (100 / Bw) * B;

    long double Rwc = (100 / Rw) * Rw;
    long double Gwc = (100 / Gw) * Gw;
    long double Bwc = (100 / Bw) * Bw;
    
    //the cone responce are converted to XYZ colour space and back
    long double Xc = (1.096124 * Rc) + (-0.278869 * Gc) + (0.182745 * Bc);
    long double Yc = (0.454369 * Rc) + (0.473533 * Gc) + (0.072098 * Bc);
    long double Zc = (-0.009628 * Rc) + (-0.005698 * Gc) + (1.015326 * Bc);
    
    long double Xwc = (1.096124 * Rwc) + (-0.278869 * Gwc) + (0.182745 * Bwc);
    long double Ywc = (0.454369 * Rwc) + (0.473533 * Gwc) + (0.072098 * Bwc);
    long double Zwc = (-0.009628 * Rwc) + (-0.005698 * Gwc) + (1.015326 * Bwc);
    
    long double rho = (0.38971 * Xc) + (0.68898 * Yc) + (-0.07868 * Zc);
    long double gamma = (-0.22981 * Xc) + (1.18340 * Yc) + (0.04641 * Zc);
    long double beta = (0 * Xc) + (0 * Yc) + (1 * Zc);
    
    long double rhoW = (0.38971 * Xwc) + (0.68898 * Ywc) + (-0.07868 * Zwc);
    long double gammaW = (-0.22981 * Xwc) + (1.18340 * Ywc) + (0.04641 * Zwc);
    long double betaW = (0 * Xwc) + (0 * Ywc) + (1 * Zwc);
    
    //luminance level adaption factor is then applied
    long double rhoA;
    if (rho < 0) {
        rho *= -1;
        rhoA = (((400 * pow((Fl * rho / 100), 0.42)) / (27.13 + pow((Fl * rho / 100), 0.42))) * -1) + 0.1;
    } else {
        rhoA = (((400 * pow((Fl * rho / 100), 0.42)) / (27.13 + pow((Fl * rho / 100), 0.42))) + 0.1);
    }
    long double gammaA;
    if (gamma < 0) {
        gamma *= -1;
        gammaA = (((400 * pow((Fl * gamma / 100), 0.42)) / (27.13 + pow((Fl * gamma / 100), 0.42))) * -1) + 0.1;
    } else {
        gammaA = (((400 * pow((Fl * gamma / 100), 0.42)) / (27.13 + pow((Fl * gamma / 100), 0.42))) + 0.1);
    }
    long double betaA;
    if (beta < 0) {
        beta *= -1;
        betaA = (((400 * pow((Fl * beta / 100), 0.42)) / (27.13 + pow((Fl * beta / 100), 0.42))) * -1) + 0.1;
    } else {
        betaA = (((400 * pow((Fl * beta / 100), 0.42)) / (27.13 + pow((Fl * beta / 100), 0.42))) + 0.1);
    }
    
    long double rhoAw;
    if (rhoW < 0) {
        rhoW *= -1;
        rhoAw = (((400 * pow((Fl * rhoW / 100), 0.42)) / (27.13 + pow((Fl * rhoW / 100), 0.42))) * -1) + 0.1;
    } else {
        rhoAw = (((400 * pow((Fl * rhoW / 100), 0.42)) / (27.13 + pow((Fl * rhoW / 100), 0.42))) + 0.1);
    }
    long double gammaAw;
    if (gammaW < 0) {
        gammaW *= -1;
        gammaAw = (((400 * pow((Fl * gammaW / 100), 0.42)) / (27.13 + pow((Fl * gammaW / 100), 0.42))) * -1) + 0.1;
    } else {
        gammaAw = (((400 * pow((Fl * gammaW / 100), 0.42)) / (27.13 + pow((Fl * gammaW / 100), 0.42))) + 0.1);
    }
    long double betaAw;
    if (betaW < 0) {
        betaW *= -1;
        betaAw = (((400 * pow((Fl * betaW / 100), 0.42)) / (27.13 + pow((Fl * betaW / 100), 0.42))) * -1) + 0.1;
    } else {
        betaAw = (((400 * pow((Fl * betaW / 100), 0.42)) / (27.13 + pow((Fl * betaW / 100), 0.42))) + 0.1);
    }
    
    //CIECAM02 correlates for red-green and yellow-blue can be determined
    long double a = rhoA - ((12.0 / 11) * gammaA) + ((1.0 / 11) * betaA);
    long double b = (1.0 / 9) * (rhoA + gammaA - (2 * betaA));
    
    //achromatic responce
    long double A = ((2 * rhoA) + gammaA + ((1.0 / 20) * betaA) - 0.305) * Nbb;
    long double Aw = ((2 * rhoAw) + gammaAw + ((1.0 / 20) * betaAw) - 0.305) * Nbb;
    
    //hue angle in degrees
    long double h = 0;
    if (a == 0) {
        if (b == 0 ) {
            h = 0;
        } else {
            if (b >= 0) {
                h = (360 / (2 * M_PI)) * atan2(b, a);
            } else {
                h = 360 + (360 / (2 * M_PI)) * atan2(b, a);
            }
        }
    } else {
        if (b >= 0) {
            h = (360 / (2 * M_PI)) * atan2(b, a);
        } else {
            h = 360 + (360 / (2 * M_PI)) * atan2(b, a);
        }
    }
    
    CEScoordREF->hueangle = h;
    long double et = (1.0 / 4) * (cos(((h * M_PI) / 180.0) + 2.0) + 3.8);   //ecetrencity
    long double t = ((50000.0 / 13) * Nc * Ncb) * et * sqrtl(pow(a, 2) + pow(b, 2)) / (rhoA + gammaA + ((21.0 / 20) * betaA));
    
    //CIECAM02 apperance correlates can be caculated as
    long double J = 100 * pow(A / Aw, c * z);
    long double C = pow(t, 0.9) * pow(J / 100.0, 0.5) * pow((1.64 - pow(0.29, n)), 0.73);
    long double M = C * pow(Fl, 0.25);
    
    //caculate CIE 2002 (J'a'b') CIECAM02 UCS Diagram coordinates of each CES with refference source
    long double Mprime = (1 / 0.0228) * log(1 + 0.0228 * M);
    CEScoordREF->Jprime_10deg = ((1 + 100 * 0.007) * J) / (1 + 0.007 * J);
    CEScoordREF->aprime_10deg = Mprime * cos((h * M_PI) / 180);
    CEScoordREF->bprime_10deg = Mprime * sin((h * M_PI) / 180);
    
    return CEScoordREF;
}

//
//  Caculate CIECAM02 Cordinates as required by TM-30-15 specification from 10 deg XYZ values
//
//  EXTERNAL REQUIREMENTS: NONE
//  INPUTS: - struct CEScoordREF* with X_10deg, Y_10deg, Z_10deg values filled in
//          - int flag, currently unused
//  OUTPUT: - CIE_coordinates* CEScoordREF, with CIE 2002 (J'a'b') 10 deg values filled
//
struct CIE_coordinates* CIECAM02(struct CIE_coordinates *CEScoordREF, int flag) {
    
    //For DEBUG only, varibales reporduce the worked example in "measuring colour 4th edition"
    CEScoordREF->X_10deg = 19.31;
    CEScoordREF->Y_10deg = 23.93;
    CEScoordREF->Z_10deg = 10.14;                                   //input - sample under test conditions
    long double Xw = 98.88, Yw = 90, Zw = 32.03;                    //input - adopted white in test conditions
    long double c = 0.69, Nc = 1, F = 1, Yb = 18, La = 200;         //input - test condition parameters
    long double xwr = 1.0 / 3, ywr = 1.0 / 3, Ywr = 100;            //constant - reference whie in refrence conditions
    long double Xwr = ((xwr * Ywr) / ywr), Zwr = (1 - xwr - ywr) * (Ywr / ywr);
    
    long double n = Yb / Yw;                                        //caculated values
    long double Nbb = 0.725 * pow((1 / n), 0.2);
    long double Ncb = Nbb;
    long double z = 1.48 + pow(n, 0.5);
    
    //compute D factor
    long double D = F * (1 - (1 / 3.6) * pow(M_E, ((-La-42.0) / 92)));
    if (D > 1) {
        D = 1;
    } else if (D < 0){
        D = 0;
    }
    
    //caculate RGB tristimulus values for adapting conditions first
    long double Rw = (0.7328 * Xw) + (0.4296 * Yw) + (-0.1624 * Zw);
    long double Gw = (-0.7036 * Xw) + (1.6975 * Yw) + (0.0061 * Zw);
    long double Bw = (0.0030 * Xw) + (0.0136 * Yw) + (0.9834 * Zw);
    
    //refference white in reference conditions
    long double Rwr = (0.7328 * Xwr) + (0.4296 * Ywr) + (-0.1624 * Zwr);
    long double Gwr = (-0.7036 * Xwr) + (1.6975 * Ywr) + (0.0061 * Zwr);
    long double Bwr = (0.0030 * Xwr) + (0.0136 * Ywr) + (0.9834 * Zwr);
    
    //CES values
    long double R = (0.7328 * CEScoordREF->X_10deg) + (0.4296 * CEScoordREF->Y_10deg) + (-0.1624 * CEScoordREF->Z_10deg);
    long double G = (-0.7036 * CEScoordREF->X_10deg) + (1.6975 * CEScoordREF->Y_10deg) + (0.0061 * CEScoordREF->Z_10deg);
    long double B = (0.0030 * CEScoordREF->X_10deg) + (0.0136 * CEScoordREF->Y_10deg) + (0.9834 * CEScoordREF->Z_10deg);
    
    long double Dr = ((Yw / Ywr) * (Rwr / Rw) * D) + (1 - D);
    long double Dg = ((Yw / Ywr) * (Gwr / Gw) * D) + (1 - D);
    long double Db = ((Yw / Ywr) * (Bwr / Bw) * D) + (1 - D);
    
    long double Rc = Dr * R;
    long double Gc = Dg * G;
    long double Bc = Db * B;
    
    long double Rwc = Dr * Rw;
    long double Gwc = Dg * Gw;
    long double Bwc = Db * Bw;
    
    long double k = 1 / (5 * La + 1);
    long double Fl = pow((0.2 * k), 4) * (5 * La) + 0.1 * pow((1 - pow(k, 4)), 2) * pow((5 * La), 1.0 / 3);
    
    long double Xc = (1.096124 * Rc) + (-0.278869 * Gc) + (0.182745 * Bc);
    long double Yc = (0.454369 * Rc) + (0.473533 * Gc) + (0.072098 * Bc);
    long double Zc = (-0.009628 * Rc) + (-0.005698 * Gc) + (1.015326 * Bc);
    
    long double Xwc = (1.096124 * Rwc) + (-0.278869 * Gwc) + (0.182745 * Bwc);
    long double Ywc = (0.454369 * Rwc) + (0.473533 * Gwc) + (0.072098 * Bwc);
    long double Zwc = (-0.009628 * Rwc) + (-0.005698 * Gwc) + (1.015326 * Bwc);
    
    long double rho = (0.38971 * Xc) + (0.68898 * Yc) + (-0.07868 * Zc);
    long double gamma = (-0.22981 * Xc) + (1.18340 * Yc) + (0.04641 * Zc);
    long double beta = (0 * Xc) + (0 * Yc) + (1 * Zc);
    
    long double rhoW = (0.38971 * Xwc) + (0.68898 * Ywc) + (-0.07868 * Zwc);
    long double gammaW = (-0.22981 * Xwc) + (1.18340 * Ywc) + (0.04641 * Zwc);
    long double betaW = (0 * Xwc) + (0 * Ywc) + (1 * Zwc);
    
    long double rhoA;
    if (rho < 0) {
        rho *= -1;
        rhoA = (((400 * pow((Fl * rho / 100), 0.42)) / (27.13 + pow((Fl * rho / 100), 0.42))) * -1) + 0.1;
    } else {
        rhoA = (((400 * pow((Fl * rho / 100), 0.42)) / (27.13 + pow((Fl * rho / 100), 0.42))) + 0.1);
    }
    long double gammaA;
    if (gamma < 0) {
        gamma *= -1;
        gammaA = (((400 * pow((Fl * gamma / 100), 0.42)) / (27.13 + pow((Fl * gamma / 100), 0.42))) * -1) + 0.1;
    } else {
        gammaA = (((400 * pow((Fl * gamma / 100), 0.42)) / (27.13 + pow((Fl * gamma / 100), 0.42))) + 0.1);
    }
    long double betaA;
    if (beta < 0) {
        beta *= -1;
        betaA = (((400 * pow((Fl * beta / 100), 0.42)) / (27.13 + pow((Fl * beta / 100), 0.42))) * -1) + 0.1;
    } else {
        betaA = (((400 * pow((Fl * beta / 100), 0.42)) / (27.13 + pow((Fl * beta / 100), 0.42))) + 0.1);
    }
    
    long double rhoAw;
    if (rhoW < 0) {
        rhoW *= -1;
        rhoAw = (((400 * pow((Fl * rhoW / 100), 0.42)) / (27.13 + pow((Fl * rhoW / 100), 0.42))) * -1) + 0.1;
    } else {
        rhoAw = (((400 * pow((Fl * rhoW / 100), 0.42)) / (27.13 + pow((Fl * rhoW / 100), 0.42))) + 0.1);
    }
    long double gammaAw;
    if (gammaW < 0) {
        gammaW *= -1;
        gammaAw = (((400 * pow((Fl * gammaW / 100), 0.42)) / (27.13 + pow((Fl * gammaW / 100), 0.42))) * -1) + 0.1;
    } else {
        gammaAw = (((400 * pow((Fl * gammaW / 100), 0.42)) / (27.13 + pow((Fl * gammaW / 100), 0.42))) + 0.1);
    }
    long double betaAw;
    if (betaW < 0) {
        betaW *= -1;
        betaAw = (((400 * pow((Fl * betaW / 100), 0.42)) / (27.13 + pow((Fl * betaW / 100), 0.42))) * -1) + 0.1;
    } else {
        betaAw = (((400 * pow((Fl * betaW / 100), 0.42)) / (27.13 + pow((Fl * betaW / 100), 0.42))) + 0.1);
    }
    
    long double a = rhoA - ((12.0 / 11) * gammaA) + ((1.0 / 11) * betaA);
    long double b = (1.0 / 9) * (rhoA + gammaA - (2 * betaA));
    long double h = atan(b / a);
    long double hr = h * (180 / M_PI);
    if (a < 0 && b >= 0) {
        h = hr + 180;
    } else if (a < 0 && b < 0) {
        h = hr + 180;
    } else if (a >= 0 && b < 0) {
        h = hr + 360;
    } else {
        h = hr;
    }
    
    long double hi[5] = {20.14, 90.00, 164.25, 237.53, 380.14};
    long double ei[5] = {0.8, 0.7, 1.0, 1.2, 0.8};
    long double hprime = h;
    if (h < hi[0]) {
        hprime += 360;
    }
    
    long double Hi = 0;
    int ivalue = 0;
    for (int i = 0; i < 5; i++) {
        if (hi[i] <= hprime) {
            if (hprime < hi[i + 1]) {
                switch (i) {
                    case 0:
                        Hi = 0;
                        ivalue = 0;
                        break;
                    case 1:
                        Hi = 100;
                        ivalue = 1;
                        break;
                    case 2:
                        Hi = 200;
                        ivalue = 2;
                        break;
                    case 3:
                        Hi = 300;
                        ivalue = 3;
                        break;
                    case 4:
                        Hi = 400;
                        ivalue = 4;
                        break;
                }
                goto nestedbreak;
            }
        }
    }
nestedbreak:;
    
    long double H = Hi + ((100 * ((hprime - hi[ivalue]) / ei[ivalue])) / (((hprime - hi[ivalue]) / ei[ivalue]) + ((hi[ivalue + 1] - hprime) / (ei[ivalue + 1]))));
    long double et = (1.0 / 4) * (cos(((h * M_PI) / 180.0) + 2.0) + 3.8);
    
    long double A = ((2 * rhoA) + gammaA + ((1.0 / 20) * betaA) - 0.305) * Nbb;
    long double Aw = ((2 * rhoAw) + gammaAw + ((1.0 / 20) * betaAw) - 0.305) * Nbb;
    
    long double J = 100 * pow(A / Aw, c * z);
    long double Q = (4 / c) * pow(J / 100, 0.5) * (Aw + 4) * pow(Fl, 0.25);
    long double t = ((50000.0 / 13) * Nc * Ncb) * et * sqrtl(pow(a, 2) + pow(b, 2)) / (rhoA + gammaA + ((21.0 / 20) * betaA));
    long double C = pow(t, 0.9) * pow(J / 100.0, 0.5) * pow((1.64 - pow(0.29, n)), 0.73);
    long double M = C * pow(Fl, 0.25);
    long double s = 100 * pow(M / Q, 0.5);
    long double Mprime = (1 / 0.0228) * log(1 + 0.0228 * M);
    
    //caculate CIE 2002 (J'a'b') CIECAM02 UCS Diagram coordinates of each CES with refference source
    CEScoordREF->Jprime_10deg = ((1 + 100 * 0.007) * J) / (1 + 0.007 * J);
    CEScoordREF->aprime_10deg = Mprime * cos((h * M_PI) / 180);
    CEScoordREF->bprime_10deg = Mprime * sin((h * M_PI) / 180);
    
    return CEScoordREF;
}

//
//  Returns a gaussian curve around a center
//
//  EXTERNAL REQUIREMENTS: must free(): SPD, SPD->wavelength, SPD->spec_radiance, SPD->spec_rad_normalised
//  INPUTS: - pointer to array with memory allocated to store results in, long dounle*
//          - pointer to array with memory allocated to store wavelength in, long double*
//          - start wavelength (nm), int
//          - stop wavelenght (nm), int
//          - wavelength interval (nm), int
//          - intensity of gaussian peak (W·m2·nm-1) or (W·sr−1·m−2·nm-1), double
//          - width of gaussian (nm), double
//          - position of gaussian peak (nm), double
//  OUTPUT: - 0 on scucess
//
int gaussian(long double *results, long double *wavelength, const int wstart, const int wstop, const int wint, const double peak_height, const double width, const double peak_center)
{
    int count = 0;
    for (int i = wstart; i <= wstop; i += wint)
    {
        if (i == wstart) {
            wavelength[0] = wstart;
        } else {
            wavelength[count] = wavelength[count - 1] + wint;
        }
        results[count] = pow(1 / (width * sqrt(2 * M_PI)), exp( -(1.0 / 2)) * pow((i - peak_center) / width, 2)) * peak_height;
        count++;
    }
    
    return 0;
}

//
//  Calculates the equivalent melanopic luminance
//	as defined by Lucas et al., "Measuring and using light in the melanopsin age." Trends in Neuroscience, Jan 2014
//	http://lucasgroup.lab.manchester.ac.uk/research/measuringmelanopicilluminance/
//	Referenced by Well Building Standard v1 2016 which uses it to calculate equilivant melanopic lux (EML)
//	as defined by that standard, Appendix C Table L1 & L2
//	requires corneal irradiance as an input
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - struct SPD_data* with wavelength, arr_len and a spectral value filled in
//          - int flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - equivalent melanopic luminance, long double
//			- -1 on error
//
long double melanopic_luminance(struct SPD_data *SPD_in, int vector_flag)
{	
	
	//check arguments
    if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL\n");
        return -1;
    } else if (SPD_in->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0\n");
        return -1;
    } else if (SPD_in->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL\n");
        return -1;
    } else if (SPD_in->wavelength[0] > 780 || SPD_in->wavelength[SPD_in->arr_len - 1] < 380) {
    	fprintf(stderr, "ERROR: wavelength range must overlap 380nm - 780nm\n");
        return -1;
    }
    
    //ajust range of data to 380nm - 780nm to match melanopic lux weightings data
    //cannot be done on normalised data as normalisation will no longer make sense
    struct SPD_data *correct_range_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    correct_range_data->wav_interval = 5;
    correct_range_data->arr_len = (780 - 380) / 5 + 1;
    correct_range_data->wavelength = (long double*)calloc(correct_range_data->arr_len, sizeof(long double));
    switch (vector_flag) {
        case 1:
            correct_range_data->spec_irradiance = (long double*)malloc(sizeof(long double) * correct_range_data->arr_len);
            correct_range_data->spec_radiance = NULL;
            break;
        case 3:
            correct_range_data->spec_radiance = (long double*)malloc(sizeof(long double) * correct_range_data->arr_len);
            correct_range_data->spec_irradiance = NULL;
            break;
        default:
            free(correct_range_data->wavelength);
            free(correct_range_data);
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            return -1;
    }
    ajust_wavlength_range(SPD_in, correct_range_data, 380, 780, vector_flag);
	
	//adjust interval of data to 5nm
	struct SPD_data *correct_interval_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
	correct_interval_data->wav_interval = 5;
    correct_interval_data->arr_len = (780 - 380) / 5 + 1;
    correct_interval_data->wavelength = (long double*)calloc(correct_interval_data->arr_len, sizeof(long double));
    correct_interval_data->spec_radiance = (long double*)calloc(correct_interval_data->arr_len, sizeof(long double));
    
	ajust_wavelength_interval(correct_range_data, correct_interval_data, 5, vector_flag);
	
	//read in melanpoic curve data
	FILE *fpmel = fopen("./data/melanopic_weightings.csv", "r");
    if (fpmel == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/melanopic_weightings.csv file.\n");
        exit(EXIT_FAILURE);
    }
	
	//skip headder of 7 lines
	int templinelen = 9999;
    char templine[templinelen];
    for (int skip = 0; skip < 7; skip++)
    {
        fgets(&templine, templinelen, fpmel);
    }
	
	long double *weight = (long double*)malloc(sizeof(long double) * 81);
	long double *weight2 = (long double*)malloc(sizeof(long double) * 81);
	long double *results = (long double*)malloc(sizeof(long double) * 81);
	long double *results2 = (long double*)malloc(sizeof(long double) * 81);
	for (int i = 0; i < 81; i++) {
        fscanf(fpmel, "%*Lf,%*Lf,%*Lf,%*Lf,%*Lf,%Lf,%Lf,%*Lf,%*Lf\n", &weight[i], &weight2[i]);
        switch (vector_flag) {			//multiply
        	case 1:
        		results[i] = correct_interval_data->spec_irradiance[i] * weight[i] * 5;
        		results2[i] = correct_interval_data->spec_irradiance[i] * weight2[i] * 5;
        		break;
        	case 3:
        		results[i] = correct_interval_data->spec_radiance[i] * weight[i] * 5;
        		results2[i] = correct_interval_data->spec_radiance[i] * weight2[i] * 5;
        		break;
        	default:
        		break;
        }
    }
    fclose(fpmel);
	
	//sum
	long double aopic_luminance = 0;
	long double lux = 0;
	//printf("\nmelanopic weightings = results:\n");
    for(int i = 0; i < 81; i++)
    {
    	aopic_luminance += results[i];
    	lux += results2[i];
    	//printf("%.0Lf %Lf %Lf\n", correct_interval_data->wavelength[i], weight[i], results[i]);
    }
	
	aopic_luminance *= 72983.25;
	//aopic_luminance = aopic_luminance / lux * 1.218;	//melanopic ratio
	
	//free() temporary data
	free(results);
	free(results2);
	free(weight);
	free(weight2);
    free(correct_range_data->wavelength);
    free(correct_range_data->spec_radiance);
    free(correct_range_data);
    free(correct_interval_data->wavelength);
    free(correct_interval_data->spec_radiance);
    free(correct_interval_data);
    
    return aopic_luminance;
}

//
//	Calculates circadian luminance
//	A combination of all visual and non-visual photoreceptors in the eye
//	http://lucasgroup.lab.manchester.ac.uk/research/measuringmelanopicilluminance/
//
//
//  EXTERNAL REQUIREMENTS: none
//  INPUTS: - struct SPD_data* with wavelength, arr_len and a spectral value filled in
//          - int flag: 1 = spec_irradiance, 3 = spec_radiance
//  OUTPUT: - equivalent melanopic luminance, long double
//			- -1 on error
//
long double circadian_luminance(struct SPD_data *SPD_in, int vector_flag)
{	
	long double aopic_luminance = 0;
	
	//check arguments
    if (SPD_in == NULL) {
        fprintf(stderr, "ERROR: pointer to SPD data to normalise cannot be NULL\n");
        return -1;
    } else if (SPD_in->arr_len <= 0) {
        fprintf(stderr, "ERROR: array length cannot be 0\n");
        return -1;
    } else if (SPD_in->wavelength == NULL) {
        fprintf(stderr, "ERROR: pointer to wavelength vector cannot be NULL\n");
        return -1;
    } else if (SPD_in->wavelength[0] > 780 || SPD_in->wavelength[SPD_in->arr_len - 1] < 380) {
    	fprintf(stderr, "ERROR: wavelength range must overlap 380nm - 780nm\n");
        return -1;
    }
    
    //ajust range of data to 380nm - 780nm to match melanopic lux weightings data
    //cannot be done on normalised data as normalisation will no longer make sense
    struct SPD_data *correct_range_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
    correct_range_data->wav_interval = 5;
    correct_range_data->arr_len = (780 - 380) / 5 + 1;
    correct_range_data->wavelength = (long double*)calloc(correct_range_data->arr_len, sizeof(long double));
    switch (vector_flag) {
        case 1:
            correct_range_data->spec_irradiance = (long double*)malloc(sizeof(long double) * correct_range_data->arr_len);
            correct_range_data->spec_radiance = NULL;
            break;
        case 3:
            correct_range_data->spec_radiance = (long double*)malloc(sizeof(long double) * correct_range_data->arr_len);
            correct_range_data->spec_irradiance = NULL;
            break;
        default:
            free(correct_range_data->wavelength);
            free(correct_range_data);
            fprintf(stderr, "ERROR: vector flag supplied is not valid\n");
            return -1;
    }
    ajust_wavlength_range(SPD_in, correct_range_data, 380, 780, vector_flag);
	
	//adjust interval of data to 5nm
	struct SPD_data *correct_interval_data = (struct SPD_data*)calloc(1, sizeof(struct SPD_data));
	correct_interval_data->wav_interval = 5;
    correct_interval_data->arr_len = (780 - 380) / 5 + 1;
    correct_interval_data->wavelength = (long double*)calloc(correct_interval_data->arr_len, sizeof(long double));
    correct_interval_data->spec_radiance = (long double*)calloc(correct_interval_data->arr_len, sizeof(long double));
    
	ajust_wavelength_interval(correct_range_data, correct_interval_data, 5, vector_flag);
	
	//read in melanpoic curve data
	FILE *fpmel = fopen("./data/melanopic_weightings.csv", "r");
    if (fpmel == NULL) {
        fprintf(stderr, "ERROR: could not open ./data/melanopic_weightings.csv file.\n");
        exit(EXIT_FAILURE);
    }
	
	//skip headder of 7 lines
	int templinelen = 9999;
    char templine[templinelen];
    for (int skip = 0; skip < 7; skip++)
    {
        fgets(&templine, templinelen, fpmel);
    }
	
	long double *weight = (long double*)malloc(sizeof(long double) * 81);
	long double *results = (long double*)malloc(sizeof(long double) * 81);
	for (int i = 0; i < 81; i++) {
        fscanf(fpmel, "%*Lf,%*Lf,%*Lf,%*Lf,%*Lf,%*Lf,%*Lf,%Lf,%*Lf\n", &weight[i]);
        switch (vector_flag) {			//multiply
        	case 1:
        		results[i] = correct_interval_data->spec_irradiance[i] * weight[i] * 5;
        		break;
        	case 3:
        		results[i] = correct_interval_data->spec_radiance[i] * weight[i] * 5;
        		break;
        	default:
        		break;
        }
    }
    fclose(fpmel);
	
	//sum
	//printf("\nmelanopic weightings = results:\n");
    for(int i = 0; i < 81; i++)
    {
    	aopic_luminance += results[i];
    	//printf("%.0Lf %Lf %Lf\n", correct_interval_data->wavelength[i], melweight[i], results[i]);
    }
	
	aopic_luminance *= 72983.25;
	
	//free() temporary data
	free(results);
	free(weight);
    free(correct_range_data->wavelength);
    free(correct_range_data->spec_radiance);
    free(correct_range_data);
    free(correct_interval_data->wavelength);
    free(correct_interval_data->spec_radiance);
    free(correct_interval_data);
    
    return aopic_luminance;
}
