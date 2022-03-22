#include "initialise.h"

void initialise(world* w, const char* patch_file, const char* people_file,const char* ini_file) {
  w->thread_count = omp_get_max_threads();
  if (w->thread_count > 1) w->thread_count -= 1; // Leave a core for java. 
                                                 // This actually makes things faster, due to less thread competition...

  omp_set_num_threads(w->thread_count);

  printf("Zombie ");
#pragma omp parallel
{
  printf("%d ", omp_get_thread_num());
}
  printf("\n");
  printf("Read patch file %s\n", patch_file);
  FILE* f_pat = fopen(&patch_file[0], "rb");
  FILE* f_peo = fopen(&people_file[0], "rb");
  unsigned int i,j;
  int k;
  size_t dummy;

  dummy = fread(&w->n_patches, 4, 1, f_pat);
  w->patch_left_lon = new float[w->n_patches];
  w->patch_top_lat = new float[w->n_patches];
  w->patch_pop = new unsigned int[w->n_patches];
  w->patch_first_index = new unsigned int[w->n_patches];
  for (i = 0; i < w->n_patches; i++) {
    dummy = fread(&w->patch_left_lon[i], 4, 1, f_pat);
    dummy = fread(&w->patch_top_lat[i], 4, 1, f_pat);
    dummy = fread(&w->patch_first_index[i], 4, 1, f_pat);
    dummy = fread(&w->patch_pop[i], 4, 1, f_pat);
  }

  w->infected_grid = new int**[PNG_WIDTH];
  w->immune_grid = new int**[PNG_WIDTH];
  for (i = 0; i < PNG_WIDTH; i++) {
    w->infected_grid[i] = new int*[PNG_HEIGHT];
    w->immune_grid[i] = new int*[PNG_HEIGHT];
    for (j = 0; j < PNG_HEIGHT; j++) {
      w->infected_grid[i][j] = new int[w->thread_count];
      w->immune_grid[i][j] = new int[w->thread_count];
      for (k = 0; k < w->thread_count; k++) {
        w->infected_grid[i][j][k] = 0;
        w->immune_grid[i][j][k] = 0;
      }
    }
  }
  

  printf("Read people file %s\n", people_file);
  dummy=fread(&w->n_people, 4, 1, f_peo);
  w->people_lat = new float[w->n_people];
  w->people_lon = new float[w->n_people];
  w->people_status = new unsigned char[w->n_people];
  w->people_patch = new unsigned short[w->n_people];
  for (i = 0; i < w->n_people; i++) {
    dummy = fread(&w->people_lon[i], 4, 1, f_peo);
    dummy = fread(&w->people_lat[i], 4, 1, f_peo); 
    w->people_status[i] = SUSCEPTIBLE;
  }
  
  int p = 0;
  for (i = 0; i < w->n_patches; i++) {
    for (j = w->patch_first_index[i]; j < w->patch_first_index[i] + w->patch_pop[i]; j++) {
      w->people_patch[p++] = (unsigned short) i;
    }
  }

  fclose(f_pat);
  fclose(f_peo);
  w->qmatrix = new float*[w->n_patches];
  for (i = 0; i < w->n_patches; i++) {
    w->qmatrix[i] = new float[w->n_patches];
  }
  
  FILE* iniFile = fopen(&ini_file[0], "r");
  dummy = fscanf(iniFile, DSCAN, &w->k_cut);
  dummy = fscanf(iniFile, DSCAN, &w->k_a);
  dummy = fscanf(iniFile, DSCAN, &w->k_b);
  dummy = fscanf(iniFile, DSCAN, &w->r0);
  dummy = fscanf(iniFile, DSCAN, &w->s_amp);
  dummy = fscanf(iniFile, DSCAN, &w->s_per);
  dummy = fscanf(iniFile, DSCAN, &w->s_off);
  dummy = fscanf(iniFile, DSCAN, &w->s_cen);
  dummy = fscanf(iniFile, DSCAN, &w->vac_lon);
  dummy = fscanf(iniFile, DSCAN, &w->vac_lat);
  dummy = fscanf(iniFile, DSCAN, &w->vac_rad);
  dummy = fscanf(iniFile, DSCAN, &w->vac_prop);
  dummy = fscanf(iniFile, DSCAN, &w->i_isolate);
  dummy = fscanf(iniFile, DSCAN, &w->i_limit_travel);
  dummy = fscanf(iniFile, DSCAN, &w->h_latent);
  dummy = fscanf(iniFile, DSCAN, &w->h_inf);
  dummy = fscanf(iniFile, DSCAN, &w->h_sympt);
  dummy = fscanf(iniFile, DSCAN, &w->n_seeds);
  dummy = fscanf(iniFile, DSCAN, &w->seed1);
  dummy = fscanf(iniFile, DSCAN, &w->seed2);
  dummy = fscanf(iniFile, DSCAN, &w->seed_radius);
  dummy = fscanf(iniFile, DSCAN, &w->seed_lon);
  dummy = fscanf(iniFile, DSCAN, &w->seed_lat);

  fclose(iniFile);
  w->h_latent -= 0.5;
  w->h_inf -= 0.5;
    
  w->contactQ = new std::vector<int>*[w->thread_count];
  w->recoveryQ = new std::vector<int>*[w->thread_count];
  w->new_cases = 0; 
  w->current_infected = 0;
  for (k = 0; k < w->thread_count; k++) {
    w->contactQ[k] = new std::vector<int>[WINDOW_SIZE];
    w->recoveryQ[k] = new std::vector<int>[WINDOW_SIZE];
    for (j = 0; j < WINDOW_SIZE; j++) {
      w->contactQ[k][j].reserve(100000);
      w->recoveryQ[k][j].reserve(100000);
    }
  }
  w->png_buffer.reserve(20000);
  w->WINDOW_START = 0;
  w->new_cases_threaded = new int[w->thread_count];

}



void calculateQ(world* w) {
  // This loop parallelises embarrassingly, since qk1,k' is independent from qk2,k'
  
  #pragma omp parallel for shared(w) default(none) schedule(static, 1)
  for (int thread_no = 0; thread_no < w->thread_count; thread_no++) {
    for (int k = thread_no; k < w->n_patches; k += w->thread_count) {
      double Z_k = 0;
      for (int k_prime = 0; k_prime < w->n_patches; k_prime++) {
        Z_k += (kernel_F(w, patch_distance(w, k, k_prime)) * w->patch_pop[k_prime]);
      }
      Z_k = 1.0 / Z_k;

      // Now generate cumulative distribution

      double cumulative = 0;
      for (int k_prime = 0; k_prime < w->n_patches; k_prime++) {
        cumulative += (kernel_F(w, patch_distance(w, k, k_prime)) * w->patch_pop[k_prime] * Z_k);
        w->qmatrix[k][k_prime] = (float)cumulative;
      }
      w->qmatrix[k][w->n_patches - 1] = 2.0;
    }
  }
}

