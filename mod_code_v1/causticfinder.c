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

#define M 200

void set_source_pixel(int i,int j,float** source_plane);
gsl_complex diffp (gsl_complex, double, double, double, int, gsl_complex*);
gsl_complex diffp2 (gsl_complex, double, double, double, int, gsl_complex*);
gsl_complex defl (gsl_complex, double, double, int, gsl_complex*);

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
  std::cout << "stars read" << std::endl;




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


/*
  output=fopen("solns.dat","w");
  for (m=1;m<=(2*nstars);m++) {
  fprintf(output,"%f %f\n", GSL_REAL(solns[m-1]),GSL_IMAG(solns[m-1]));
  }
  fclose(output);
*/

  // GERLUMPH code: read map properties and initialize
  double ss2 = ss/2.0;
  
  // define source plane grid, identical to a map
  float** source_plane = (float**) calloc(res,sizeof(float*));
  for(int i=0;i<res;i++){
    source_plane[i] = (float*) calloc(res,sizeof(float));
    for(int j=0;j<res;j++){
      source_plane[i][j] = (float) 0.0;
    }
  }

  float scalingFy = ((float)(res/2)/(float)(ss2));
  float pixsize = (float)(ss)/(float)(res);
  int i_bin,j_bin;
  // GERLUMPH code  
  


  /* determine critical curves and map to caustics */
  gsl_complex* curve = (gsl_complex*) malloc(M*sizeof(gsl_complex));
  FILE* output = fopen((path+"caustics.txt").c_str(),"w");
  for(n=1;n<=(2*nstars);n++){
    if( n%200 == 0 ){
      printf("criticals+caustics %5.1f %% \n",100*(float)n/(2*nstars));
    }

    zp = solns[n-1];
    
    for(m=1;m<=M;m++){
      t    = gsl_complex_rect(0.0,0.0);      
      phi  = m*2.9*M_PI/M;
      /*norm = gsl_complex_abs(zp);*/
      norm = 10;
      while( norm>0.000001 ){
	t    = zp;
	dum  = diffp2(t,phi,gamma_o,mass,nstars,stars);
	zp   = gsl_complex_sub(t,dum);
	norm = gsl_complex_abs(gsl_complex_sub(t,zp));
      }
      z = defl(zp,gamma_o,mass,nstars,stars);
      fprintf(output,"%f %f %f %f\n", GSL_REAL(zp)/normalize,GSL_IMAG(zp)/normalize,GSL_REAL(z)*normalize,GSL_IMAG(z)*normalize);
      curve[m-1] = z;
      
      if( m > 1 ){
	float x0 = GSL_REAL(curve[m-2])*normalize;
	float y0 = GSL_IMAG(curve[m-2])*normalize;
	float x1 = GSL_REAL(curve[m-1])*normalize;
	float y1 = GSL_IMAG(curve[m-1])*normalize;
	
	if( (fabs(x0) < ss2) && (fabs(y0) < ss2) ){
	  float len = hypot(x0-x1,y0-y1);
	  int steps = (int)(len/(0.5*pixsize));
	  float xx,yy;
	  for(int p=0;p<steps;p++){
	    xx = x0+p*(x1-x0)/steps;
	    yy = y0+p*(y1-y0)/steps;
	    if( (fabs(xx) < ss2) && (fabs(yy) < ss2) ){
	      i_bin = (int)((xx+ss2)*scalingFy);
	      j_bin = (int)((yy+ss2)*scalingFy);
	      set_source_pixel(i_bin,j_bin,source_plane);
	    }
	  }
	}	
      }

      
      /*
	curve[(n-1)*M+m-1]=zp;
	printf("%f %f\n",GSL_REAL(zp),GSL_IMAG(zp));
      */
    }
  }
  fclose(output);
  free(curve);

  
  // GERLUMPH code: Write caustic map
  int pix;
  std::ofstream out_bin(path+"caustics.bin",std::ios::out|std::ios::binary);
  for(int i=0;i<res;i++){
    for(int j=0;j<res;j++){
      pix = int(source_plane[i][j]);
      out_bin.write((const char*) (&pix),sizeof(int));
    }
  }
  out_bin.close();
  // GERLUMPH code

  
  
  return 0;
}


void set_source_pixel(int i,int j,float** source_plane){
  std::vector<int> ii{i-1,i-1,i-1,i,i,i,i+1,i+1,i+1};
  std::vector<int> jj{j-1,j,j+1,j-1,j,j+1,j-1,j,j+1};

  for(int k=0;k<ii.size();k++){
    if( 0 < ii[k] && ii[k] < 10000 && 0 < jj[k] && jj[k] < 10000 ){
      source_plane[ii[k]][jj[k]] = 1;
    }
  } 
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

gsl_complex diffp2 (gsl_complex z, double phi, double gamma, double mass, int nstars, gsl_complex *stars){
  int i;
  gsl_complex sum, sum1, sum2, sum3, dum, dum2, dum3;

  sum1=gsl_complex_rect(0.0,0.0);
  dum=gsl_complex_rect(mass,0.0);
  for (i=1;i<=nstars;i++) {
    dum2=gsl_complex_pow_real(gsl_complex_sub(stars[i-1],z),2.0);
    sum1=gsl_complex_add(sum1,gsl_complex_div(dum,dum2));    
  }

  sum2=gsl_complex_rect(0.0,0.0);
  dum=gsl_complex_rect(2.0*mass,0.0);
  for (i=1;i<=nstars;i++) {
    dum2=gsl_complex_pow_real(gsl_complex_sub(stars[i-1],z),3.0);
    sum2=gsl_complex_add(sum2,gsl_complex_div(dum,dum2));    
  }

  dum=gsl_complex_rect(gamma,0.0);
  dum2=gsl_complex_exp(gsl_complex_rect(0.0,phi));
  dum3=gsl_complex_sub(dum,dum2);

  sum=gsl_complex_div(gsl_complex_add(dum3,sum1),sum2);

  return sum;
}

gsl_complex defl (gsl_complex z, double gamma, double mass, int nstars, gsl_complex *stars){
  int i;
  gsl_complex sum, sum1, dum, dum2;

  sum1=gsl_complex_rect(0.0,0.0);
  dum=gsl_complex_rect(mass,0.0);
  for (i=1;i<=nstars;i++) {
    dum2=gsl_complex_conjugate(gsl_complex_sub(stars[i-1],z));
    sum1=gsl_complex_add(sum1,gsl_complex_div(dum,dum2));    
  }

  dum=gsl_complex_mul(gsl_complex_conjugate(z),gsl_complex_rect(gamma,0.0));
  dum2=gsl_complex_add(z,dum);

  sum=gsl_complex_add(dum2,sum1);

  return sum;
}
