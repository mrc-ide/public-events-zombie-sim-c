#include "initialise.h"

void doContactsAndRecoveries(world* w,int day) {

  // This handles all the contacts and recoveries - I fused the two together to increase throughput, since contacts and recoveries never
  // happen in the same timestep.

  // Thread safety works by writing to 2-D vectors [t][time] - for thread t.
  // Efficient reading happens by ALL threads visiting ALL [t][time], but mod arithmetic to only process 1/thread_count of them.
  


  for (int i = 0; i < w->thread_count; i++) w->new_cases_threaded[i] = 0;
  double p_season = w->s_cen + (w->s_amp * cos((TWO_PI * (day + w->s_off)) / w->s_per));
  if (p_season < 0) p_season = 0;

 #pragma omp parallel for schedule(static,1) 
  for (int thread = 0; thread < w->thread_count; thread++) {
    unsigned int p = 0;
    for (int i = 0; i < w->thread_count; i++) {
      for (unsigned int j = 0; j < w->contactQ[i][w->WINDOW_START].size(); j++) {
        if (p%w->thread_count == thread) { // Parallelisation basically happens on this line.
          
          bool symptomatic = (rnd(thread) < w->h_sympt);
          int n_contacts = 0;
          if (w->i_isolate > 0) {
            if (symptomatic) {
              if (rnd(thread) < w->i_isolate) {
                n_contacts = -1;
              }
            }
          }
          if (n_contacts == 0) {
            n_contacts = poi(w->r0 * p_season, thread);
            int inf_id = w->contactQ[i][w->WINDOW_START].at(j);
            int inf_period = (int)(1.0 + (w->h_inf * rnd(thread)));   // Infectious period of host always >0 - so recovery always happens later.
            int j = 0;
            while (j < n_contacts) {
              int candidate_patch = getCommunityContactPatch(w, rnd(thread), w->people_patch[inf_id]);
              int candidate_susc = w->patch_first_index[candidate_patch] + (int)(rnd(thread)*w->patch_pop[candidate_patch]);
              DTYPE D_kk2 = patch_distance(w, candidate_patch, w->people_patch[inf_id]);
              DTYPE r_ij = haversine(w->people_lon[inf_id], w->people_lat[inf_id], w->people_lon[candidate_susc], w->people_lat[candidate_susc]);
              DTYPE S_ij = kernel_F(w, r_ij) / kernel_F(w, D_kk2);
              if (symptomatic) {
                if (w->i_limit_travel > 0) {
                  if (r_ij > w->i_limit_travel) S_ij = -1;
                }
              }
              if (S_ij >= 0) {
                if (rnd(thread) < S_ij) {
                  j++;
                  bool got_susc = false;

                  // Below is the one section that makes simulations "non-repeatable" (ie, no guaranteed same result)
                  //
                  // It IS thread-safe - we won't get nasty numerical messups, but there IS a race condition, and it
                  // DOES come up often, since you have zombies near each other competing for their victims with a local travel kernel.
                  //
                  // The race is for which thread (zombie) enters this critical section first, to bite the same victim.
                  // And that changes the random number sequence, so we diverge from the next timestep.
                  //
                  // But we want this to be fast, rather than repeatable, so hey ho.
                  //
                  // The real solution in the global sim is to store the time of biting, and after the parallel loop, revisit
                  // and pick the earliest bite, which then re-enforces determinism.

#pragma omp critical
                  {
                    got_susc = (w->people_status[candidate_susc] == SUSCEPTIBLE);   // Protected read
                    if (got_susc) w->people_status[candidate_susc] = CONTACTED;     // and write. Race!
                  }

                  if (got_susc) {
                    w->new_cases_threaded[thread]++;  // Thread-private write.
                    int latent = (int)(1.0 + (w->h_latent * rnd(thread)));  // At least one timestep, so we'll always deal with new contact in later step.
                    int t_contact = (int)(rnd(thread) * inf_period);        // Assuming uniformly distributed infectious and latent period.
                    int daysAway = (w->WINDOW_START + t_contact + latent) % WINDOW_SIZE;
                    w->contactQ[thread][daysAway].push_back(candidate_susc); // Thread-private write.
                  }
                }
              }
            }

            int daysAway = (w->WINDOW_START + inf_period) % WINDOW_SIZE; // Schedule recovery.
            w->infected_grid[(int)(PNG_SCALE*(w->people_lon[inf_id] - LEFT_LON))][(int)(PNG_SCALE * (TOP_LAT - w->people_lat[inf_id]))][thread]++;
            w->recoveryQ[thread][daysAway].push_back(inf_id);  // Thread-private write
          }
        }
        p++;
      }
    }
    if (thread == 0) {
      w->current_infected = p; // All threads visit "p" people - they only do something interesting with (1/thread_count) of them though. 
                               // Hence as long as one of the threads sets w->current_infected, we're good.
    }
    
    // Recoveries below - kept in the same openmp loop to save fork overhead.

    p = 0;
    for (int i = 0; i < w->thread_count; i++) {
      for (unsigned int j = 0; j < w->recoveryQ[i][w->WINDOW_START].size(); j++) {
        if (p%w->thread_count == thread) {   // Parallelisation step again.
          int inf_id = w->recoveryQ[i][w->WINDOW_START].at(j);
          w->people_status[inf_id] = RECOVERED;
          w->immune_grid[(int)(PNG_SCALE * (w->people_lon[inf_id] - LEFT_LON))][(int)(PNG_SCALE * (TOP_LAT - w->people_lat[inf_id]))][thread]++;
          w->infected_grid[(int)(PNG_SCALE * (w->people_lon[inf_id] - LEFT_LON))][(int)(PNG_SCALE * (TOP_LAT - w->people_lat[inf_id]))][thread]--;
        }
        p++;
      }
    }
  }

  for (int t = 0; t < w->thread_count; t++) {
    w->new_cases += w->new_cases_threaded[t];
    w->contactQ[t][w->WINDOW_START].clear();
    w->recoveryQ[t][w->WINDOW_START].clear();
  }
}

void randomSeed(world* w, unsigned int no_seeds, DTYPE seed_lon, DTYPE seed_lat, DTYPE seed_radius) {
  unsigned int i = 0;
  if (seed_lon < -50) {
    while (i < no_seeds) {
      int p = (int) (rnd(0) * w->n_people);
      if (w->people_status[p] == SUSCEPTIBLE) {
        w->contactQ[0][0].push_back(p);
        w->people_status[p] = CONTACTED;
        i++;
      }
    }
  } else {
    unsigned int candidates = 0;
    for (i = 0; i < w->n_people; i++) {
      if (haversine(seed_lon, seed_lat, w->people_lon[i], w->people_lat[i]) < seed_radius) candidates++;
    }
    if (no_seeds > candidates) no_seeds = candidates;
    int remember_candidates = candidates;
    if (no_seeds > 0) {
      double p_match = (double) no_seeds / (double) candidates;
      while (no_seeds > 0) {
        for (i = 0; i < w->n_people; i++) {
          if (haversine(seed_lon, seed_lat, w->people_lon[i], w->people_lat[i]) < seed_radius) {
            if (rnd(0) < p_match) {
              if (w->people_status[i] == SUSCEPTIBLE) {
                w->contactQ[0][0].push_back(i);
                w->people_status[i] = CONTACTED;
                no_seeds--;
                if (no_seeds == 0) i = w->n_people;
              }
            }
            candidates--;
            if (candidates == 0) {
              i = w->n_people;
              candidates = remember_candidates;
            }
          }
        }
      }
    }
  }
}

void updateImage(world *w) {

  #pragma omp parallel for shared(w) default(none) schedule(static,1)
  for (int thread_no = 0; thread_no < w->thread_count; thread_no++) {
    for (int x = thread_no; x < PNG_WIDTH; x += w->thread_count) {
      for (int y = 0; y < PNG_HEIGHT; y++) {
        int total_inf = 0;
        int total_imm = 0;
        for (int z = 0; z < w->thread_count; z++) {
          total_inf += w->infected_grid[x][y][z];
          total_imm += w->immune_grid[x][y][z];
        }
          
        if (total_inf + total_imm !=0) {
          unsigned char val = (unsigned char) (1 + (254.0 * ((float) total_imm / (total_imm + total_inf))));
          w->image[(PNG_WIDTH * y) + x] = val;
        }
      }
    }
  }
}

void saveImage(world *w,int no) {
  std::string filename;
  filename.append("mov");
  std::stringstream noConverter;
  if (no<1000) filename.append("0");
  if (no<100) filename.append("0");
  if (no<10) filename.append("0");
  noConverter << (no);
  filename.append(noConverter.str());
  filename.append(".png");
  char* fpointer = &filename[0];

  //create encoder and set settings and info (optional)
  
  LodePNG::Encoder encoder;
  
  
  // Black is black...
  // This makes the zombie palette: start at 0,255,0 (green), progress to 127,255,0   (sort of greeny red)
  // Then, from 127,255,0 to 255,0,0. 
  // Then add one more to make sure we have 256 entries.
  
  // encoder.addPalette((unsigned char)0, (unsigned char)0, (unsigned char)0, (unsigned char)0);
  // for (int i=0; i<127; i++) encoder.addPalette((unsigned char)i,(unsigned char)255,(unsigned char) 0,(unsigned char)255);
  // for (int i=0; i<127; i++) encoder.addPalette((unsigned char)127+i,(unsigned char)255-(2*i),(unsigned char)255-(2*i),(unsigned char)255);
  // encoder.addPalette((unsigned char)255,(unsigned char)0,(unsigned char)0,(unsigned char)255);
  
  // More conventionally, this one will start at red, and go to green.
 
  /*
  encoder.addPalette((unsigned char)0, (unsigned char)0, (unsigned char)0, (unsigned char)0);
  for (int i=0; i<127; i++) encoder.addPalette((unsigned char)255,(unsigned char)(i),(unsigned char)0,(unsigned char)255);
  for (int i=0; i<127; i++) encoder.addPalette((unsigned char)(255-(2*i)),(unsigned char)(i+127),(unsigned char)(i/4),(unsigned char)255);
  encoder.addPalette((unsigned char)0,(unsigned char)255,(unsigned char)0,(unsigned char)255);
  */

  // And this one will start at red, and go to blue.
 
  encoder.addPalette((unsigned char) 0, (unsigned char) 0, (unsigned char) 0, (unsigned char) 0);
  for (int i = 0; i < 127; i++) encoder.addPalette((unsigned char) 255,(unsigned char)(i),(unsigned char) 0,(unsigned char) 255);
  for (int i = 0; i < 127; i++) encoder.addPalette((unsigned char) (255 - (2 * i)),(unsigned char) (127 - i),(unsigned char) (i * 2),(unsigned char) 255);
  encoder.addPalette((unsigned char) 0, (unsigned char) 0, (unsigned char) 255, (unsigned char) 255);
 
  
   //both the raw image and the encoded image must get colorType 3 (palette)
  encoder.getInfoPng().color.colorType = 3; //if you comment this line, and store the palette in InfoRaw instead (use getInfoRaw() in the previous lines), then you get the same image in a RGBA PNG.
  encoder.getInfoRaw().color.colorType = 3;

  /* end */
  encoder.getSettings().zlibsettings.windowSize = 2048;

  //encode and save
 
  encoder.encode(w->png_buffer, &(w->image[0]), PNG_WIDTH, PNG_HEIGHT);
  LodePNG::saveFile(w->png_buffer, fpointer);
  w->png_buffer.clear();
}

void saveData(world *w, int no) {
  std::string filename;
  filename.append("data");
  std::stringstream noConverter;
  if (no<1000) filename.append("0");
  if (no<100) filename.append("0");
  if (no<10) filename.append("0");
  noConverter << (no);
  filename.append(noConverter.str());
  std::string filename2 = filename;
  filename.append(".txt");
  char* fpointer = &filename[0];
  FILE* f = fopen(&fpointer[0],"w");
  fprintf(f,"%d\t%d\t%d\n", no, w->current_infected, w->new_cases);
  printf("%d\t%d\t%d\n", no, w->current_infected, w->new_cases);
  fflush(f);
  fclose(f);
  filename.append("2");
  f = fopen(&fpointer[0],"w");
  fprintf(f, "\n");
  fclose(f);
}

void vaccinate(world* w) {
  for (unsigned int i = 0; i < w->n_people; i++) {
    if (rnd(0) < w->vac_prop) {
      if (haversine(w->vac_lon, w->vac_lat, w->people_lon[i], w->people_lat[i]) < w->vac_rad) {
        w->people_status[i] = RECOVERED;
      }
    }
  }
}

void mainLoop(world* w) {
  clock_t time = clock();
  int day = 0;
  bool activity = true;
  vaccinate(w);
  while ((w->current_infected > 0) || (activity)) {
    w->new_cases = 0;
    w->current_infected = 0;
    doContactsAndRecoveries(w, day);
    
    w->WINDOW_START = (w->WINDOW_START + 1) % WINDOW_SIZE;

    day++;
    if (w->current_infected == 0) {
      activity = false;
      for (int i = 0; i < w->thread_count; i++) {
        for (int j = 0; j < WINDOW_SIZE; j++) {
          if ((w->contactQ[i][j].size() > 0) || (w->recoveryQ[i][j].size() > 0)) {
            j = WINDOW_SIZE;
            i = w->thread_count;
            activity = true;
          }
        }
      }
    }
    updateImage(w);
    saveImage(w, day);
    saveData(w, day);
  }
  time = clock() - time;
  double dt = (double) time / CLOCKS_PER_SEC;
  printf("time = %lf\n", dt);

}

int main()  {
  world* w = new world();
  
#ifdef _WIN32
  const char* patch_file = "..\\data\\patches.bin";
  const char* person_file = "..\\data\\people.bin";
  const char* ini_file = ".\\zom.ini";
#else
  const char* patch_file = "../data/patches.bin";
  const char* person_file = "../data/people.bin";
  const char* ini_file = "./zom.ini";
#endif
  initialise(w, patch_file, person_file, ini_file);
  calculateQ(w);
  initSeeds((int) w->seed1, (int) w->seed2);
  printf("Zombie Server\n");
  printf("People: %d, Patches: %d\n", w->n_people, w->n_patches);
  randomSeed(w,(int) w->n_seeds, w->seed_lon, w->seed_lat, w->seed_radius);
  mainLoop(w);
  fflush(stdout);
  FILE* f = fopen("end.zom", "w");
  fprintf(f, "The End\n");
  fclose(f);
  int s = 0;
  for (unsigned int i = 0; i < w->n_people; i++) {
    if (w->people_status[i] != SUSCEPTIBLE) s++;

  }
  printf("tot people = %d,tot=%d\n", w->n_people, s);
  return 0;
}
