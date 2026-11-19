/*
   causticfinder calculates caustics based on Hans Witt's Method
   (PhD thesis, Hamburg, 1991)

   Copyright 2007, 2008  Robert W. Schmidt <rschmidt@ari.uni-heidelberg.de>
  
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.


   **********************************************************************

   USAGE NOTICE

   star positions in input file given in lens plane Einstein radii
                                      ER = sqrt(4GM/c^2 Dd Dds/Ds)

   First version 01 05 2007
   last update 06 02 2008
*/

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <omp.h>

#define M 200

gsl_complex diffp (gsl_complex, double, double, double, int, gsl_complex*);

int main(int argc, char **argv){
  int i,k,n,m,nstars;
  float x,y,mass;
  double gamma_o, norm, phi, sigmac, normalize;
  gsl_complex z, zp, t, sum, dum, dum2;
  gsl_complex *solns, *stars;

  
  if( argc != 2 ){
      fprintf(stderr,"Usage: causticfinder path/to/map/\n");
      exit(-1);
  }

  std::string path = argv[1];
  
  
  
  /* read stars */
  std::vector<float> x_vec,y_vec;
  std::ifstream input(path+"/lens_pos.dat");
  while( true ){
    input >> x >> y >> mass;
    if( input.eof() ) break;
    x_vec.push_back(x);
    y_vec.push_back(y);
  }
  input.close();
  nstars = x_vec.size();

  /* all objects have the same mass (yet) */
  mass=1.0;




  // GERLUMPH code: read mapmeta
  std::ifstream mapmeta(path+"mapmeta.dat");
  float avgmu,avgN,ss;
  int res;
  float kappa,gamma,smf;
  mapmeta >> avgmu >> avgN;
  mapmeta >> res;
  mapmeta >> ss;
  mapmeta >> kappa >> gamma >> smf;
  mapmeta.close();
  // GERLUMPH code

  std::cout << avgmu << "avgN\n";
  std::cout << res << "res\n";
  
  /* ML parameters */
  sigmac  = kappa*smf;
  gamma_o = -gamma/(1.0-sigmac);
  
  /* The normalized lens equation is written with the Einstein radius including sigmac, set to unity if you're using normalized units */
  //normalize = 1.0;
  normalize=sqrt(fabs(1.0-sigmac));


  
  
  stars = (gsl_complex *) calloc(nstars,sizeof(gsl_complex));
  solns = (gsl_complex *) calloc(2*(nstars),sizeof(gsl_complex));

  /* clear arrays */
  for (i=1;i<=nstars;i++) {
    stars[i-1]=gsl_complex_rect(0.0,0.0);
    solns[2*(i-1)]=gsl_complex_rect(0.0,0.0);
    solns[2*i-1]=gsl_complex_rect(0.0,0.0);
  }
  
  for(i=1;i<=nstars;i++) {
    stars[i-1] = gsl_complex_rect(x_vec[i-1]*normalize,y_vec[i-1]*normalize);
  }

  
  /* determine starting values for all critical curve loops */
  zp=gsl_complex_rect(0.5,0.0);

  for (m=1;m<=(2*nstars);m++) {
    if( m%300 == 0 ){
      printf("setting up %5.1f %% \n",100*(float)m/(2*nstars));
    }

    zp=gsl_complex_mul(zp,gsl_complex_exp(gsl_complex_rect(0.0,0.1)));

    dum=gsl_complex_rect(1.0,0.0);
    norm=gsl_complex_abs(zp);
    int aux = 0;
    while (norm>0.000001) {
      if ( aux > 10000 ) break;
      t=zp;

      sum=gsl_complex_rect(0.0,0.0);
      
      for (k=1;k<=m-1;k++) {
	sum=gsl_complex_add(sum,gsl_complex_div(dum,gsl_complex_sub(t,solns[k-1])));
      }
      dum2=diffp(t,0.0,gamma_o,mass,nstars,stars);
      zp=gsl_complex_sub(t,gsl_complex_div(dum,gsl_complex_sub(dum2,sum)));

      norm=gsl_complex_abs(gsl_complex_sub(t,zp));
      aux = aux + 1;
    }
    solns[m-1]=zp;
  }

  printf("setup done\n");

  FILE *output = fopen((path+"solns.dat").c_str(),"w");
  for (m=1;m<=(2*nstars);m++) {
  fprintf(output,"%f %f\n", GSL_REAL(solns[m-1]),GSL_IMAG(solns[m-1]));
  }
  fclose(output);


  return 0;
}


gsl_complex diffp (gsl_complex z, double phi, double gamma, double mass, int nstars, gsl_complex *stars){
  int i;
  gsl_complex sum, sum1, sum2, sum3, dum, dum2, dum3;

  sum1=gsl_complex_rect(0.0,0.0);
  dum=gsl_complex_rect(-2.0,0.0);
  for (i=1;i<=nstars;i++) {
    dum2=gsl_complex_div(dum,gsl_complex_sub(stars[i-1],z));
    sum1=gsl_complex_add(sum1,dum2);
  }

  sum2=gsl_complex_rect(0.0,0.0);
  dum=gsl_complex_rect(2.0*mass,0.0);
  for (i=1;i<=nstars;i++) {
    dum2=gsl_complex_pow_real(gsl_complex_sub(stars[i-1],z),3.0);
    sum2=gsl_complex_add(sum2,gsl_complex_div(dum,dum2));    
  }

  sum3=gsl_complex_rect(0.0,0.0);
  dum=gsl_complex_rect(mass,0.0);
  for (i=1;i<=nstars;i++) {
    dum2=gsl_complex_pow_real(gsl_complex_sub(stars[i-1],z),2.0);
    sum3=gsl_complex_add(sum3,gsl_complex_div(dum,dum2));    
  }
  dum=gsl_complex_rect(gamma,0.0);
  dum2=gsl_complex_exp(gsl_complex_rect(0.0,phi));
  dum3=gsl_complex_sub(dum,dum2);

  sum=gsl_complex_add(sum1,gsl_complex_div(sum2,gsl_complex_add(dum3,sum3)));

  return sum;
}
