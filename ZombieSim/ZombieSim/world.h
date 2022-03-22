#ifndef WORLD_H
#define WORLD_H
// VS2010 Compiler wants fopen_s instead of fopen - we've found fopen_s to be occasionally unreliable, incorrectly giving PERMISSION_DENIED errors for some reason.
#define _CRT_SECURE_NO_WARNINGS
// So that removes the warnings.

#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include "lodepng.h"
#include "maths.h"
#include <vector>

#define SUSCEPTIBLE 1
#define CONTACTED   2
#define RECOVERED   3

class world { // The world as this node sees it.
  public:
    float* people_lat;
    float* people_lon;
    unsigned short* people_patch;
    unsigned char* people_status;
    unsigned int n_people;

    unsigned short n_patches;
    float* patch_top_lat;
    float* patch_left_lon;
    unsigned int* patch_pop;
    unsigned int* patch_first_index;
    
    
    float** qmatrix;

    int thread_count;

    std::vector<int>** contactQ;    // Contacts[thread][timestep in window]
    std::vector<int>** recoveryQ;   // Contacts[thread][timestep in window]
    #define WINDOW_SIZE 32
    unsigned char WINDOW_START;

    
    // Bits for the output images

    #define PNG_WIDTH 260  
    #define PNG_HEIGHT 248
    #define LEFT_LON -11
    #define TOP_LAT 61
    #define PNG_SCALE 20.0

    unsigned char image[PNG_WIDTH*PNG_HEIGHT];
    std::vector<unsigned char> png_buffer;
    int*** infected_grid;
    int*** immune_grid;
    

    // Parameters
    
    DTYPE k_cut,k_a,k_b,r0,s_amp,s_cen,s_per,s_off;
	DTYPE vac_lon,vac_lat,vac_prop,vac_rad;
	DTYPE i_isolate,i_limit_travel;
	DTYPE h_latent,h_inf,h_sympt;

	DTYPE n_seeds,seed1,seed2;
	DTYPE seed_radius,seed_lon,seed_lat;
      

    // Stats

    int current_infected; 
    int new_cases;       
    int* new_cases_threaded;
   
};

#endif