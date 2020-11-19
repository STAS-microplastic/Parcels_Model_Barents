/*
 *
 * @date: December 1 2018
 * @author: Philippe Delandmeter
 *
 * Useful functions using open for loop for MP particle postprocessing
*/


#include <stdio.h>
#include <stdlib.h>

void histogram(int **pxi, int **pyi, int **histogram_map, int npart, int ntimes, int nlat, int nlon)
{
  int (*xi)[ntimes] = (int (*)[ntimes]) pxi;
  int (*yi)[ntimes] = (int (*)[ntimes]) pyi;
  int (*hist_map)[nlon] = (int (*)[nlon]) histogram_map;

  int pi,ti;
  for (pi=0; pi<npart; pi++){
    for (ti=0; ti<ntimes; ti++){
      if ((yi[pi][ti] >= 0) && (yi[pi][ti] <nlat) &&
          (xi[pi][ti] >= 0) && (xi[pi][ti] <nlon))
        hist_map[yi[pi][ti]][xi[pi][ti]] += 1;
    }
  }
}

void cumulative_age(int **pxi, int **pyi, float **page, float **cumulative_age_map, int npart, int ntimes, int nlat, int nlon)
{
  int (*xi)[ntimes] = (int (*)[ntimes]) pxi;
  int (*yi)[ntimes] = (int (*)[ntimes]) pyi;
  float (*age)[ntimes] = (float (*)[ntimes]) page;
  float (*age_map)[nlon] = (float (*)[nlon]) cumulative_age_map;

  int pi,ti;
  for (pi=0; pi<npart; pi++){
    for (ti=0; ti<ntimes; ti++){
      if ((yi[pi][ti] >= 0) && (yi[pi][ti] <nlat) &&
          (xi[pi][ti] >= 0) && (xi[pi][ti] <nlon))
        age_map[yi[pi][ti]][xi[pi][ti]] += age[pi][ti];
    }
  }
}

void cumulative_depth_lim(int **pxi, int **pyi, float **pdepth, float **cumulative_depth_count_map, int npart, int ntimes, int nlat, int nlon, float depth_lim)
{
  int (*xi)[ntimes] = (int (*)[ntimes]) pxi;
  int (*yi)[ntimes] = (int (*)[ntimes]) pyi;
  float (*depth)[ntimes] = (float (*)[ntimes]) pdepth;
  float (*depth_count_map)[nlon] = (float (*)[nlon]) cumulative_depth_count_map;

  int pi,ti;
  for (pi=0; pi<npart; pi++){
    for (ti=0; ti<ntimes; ti++){
      if ((yi[pi][ti] >= 0) && (yi[pi][ti] <nlat) &&
          (xi[pi][ti] >= 0) && (xi[pi][ti] <nlon)){
        if (depth[pi][ti] >= depth_lim)
          depth_count_map[yi[pi][ti]][xi[pi][ti]] += 1.;
      }
    }
  }
}


void touched(int **pxi, int **pyi, float **touched_map, int npart, int ntimes, int nlat, int nlon)
{
  int (*xi)[ntimes] = (int (*)[ntimes]) pxi;
  int (*yi)[ntimes] = (int (*)[ntimes]) pyi;
  float (*touch_map)[nlon] = (float (*)[nlon]) touched_map;
  int *touch_tmp = malloc(nlon*nlat*sizeof(int));

  int pi,ti,i;
  for (pi=0; pi<npart; pi++){
    for (i=0; i<nlon*nlat; i++)
       touch_tmp[i] = 0;
    for (ti=0; ti<ntimes; ti++){
      if ((yi[pi][ti] >= 0) && (yi[pi][ti] <nlat) &&
          (xi[pi][ti] >= 0) && (xi[pi][ti] <nlon)){
        if (touch_tmp[ xi[pi][ti]+nlon*yi[pi][ti] ] == 0){
          touch_map[yi[pi][ti]][xi[pi][ti]] += 1.;
          touch_tmp[ xi[pi][ti]+nlon*yi[pi][ti] ] = 1;
        }
      }
    }
  }

  free(touch_tmp);
}

void zone_concentration(int **pxi, int **pyi, int **page, int **zones, float **zone_concentration,
                        int npart, int ntimes, int nlat, int nlon, int nzone, int nstep)
{
  int (*xi)[ntimes] = (int (*)[ntimes]) pxi;
  int (*yi)[ntimes] = (int (*)[ntimes]) pyi;
  int (*age)[ntimes] = (int (*)[ntimes]) page;
  int (*zone_map)[nlon] = (int (*)[nlon]) zones;
  float (*zone_count)[nstep] = (float (*)[nstep]) zone_concentration;

  int pi,ti, z, t;
  for (pi=0; pi<npart; pi++){
    for (ti=0; ti<ntimes; ti++){
      t = age[pi][ti];
      if ((yi[pi][ti] >= 0) && (yi[pi][ti] <nlat) &&
          (xi[pi][ti] >= 0) && (xi[pi][ti] <nlon)){
        z = zone_map[yi[pi][ti]][xi[pi][ti]] - 1;
        if ((t >= nstep) || (z >= nzone) || (z < 0) || (t<0) )
          printf("error here %d %d (%d %d) [%d %d]\n", z, t, nzone, nstep,yi[pi][ti],xi[pi][ti]);
        else
          zone_count[z][t] +=1;
      }
      else{
        zone_count[3][t] +=1; // Arctic particle
      }
    }
  }
}

