#ifndef INIT_H
#define INIT_H

#include "omp.h"
#include "world.h"
#include "maths.h"
#include "lodepng.h"
#include <ctime>


void initialise(world* w, const char* patch_file, const char* people_file,const char* iniFile);
void calculateQ(world* w);

#endif
