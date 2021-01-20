#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#include<stdio.h>
#include<math.h>
#include<mex.h>

double surfacetension(double ori1, double ori2);

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
 /* 
  *  updatelevelsetdata(presence,grains,ID,ORI);
  *
  *  CAUTION: MODIFIES MATLAB INPUT *IN PLACE*.
  */

  mxArray *indices, *grainlevsetvals, *grainconvvals, *locs, *sharedd;
  double *pindices, *plocs, *pgrainlevsetvals, *pgrainconvvals, *id, *ori, *psharedd;
  double sum, mink, st;
  int N,i,j,k,p,ell,b,dims,nograins,gind,gind2,idk,idk1,idell,nosharedd,Nin;
  double temp[100][100], phi[100], minphi[100];
  
  dims = mxGetM(prhs[0]); /* Number of pixels. */
  N = mxGetM(prhs[1]);    /* Number of grains. */  

  id = (double *) mxGetData(prhs[2]);  /* List of grain IDs. */
  ori = (double *) mxGetData(prhs[3]); /* Grain orientations. */
  
 for (j=0;j<dims;j++){ /* Loop over pixels. */
//     printf("%d\n", j);
    indices = mxGetCell(prhs[0],j); /* Grains near this pixel. */
    pindices = (double *) mxGetData(indices);
    locs = mxGetCell(prhs[0],2*dims+j); /* Location of pixel in grain's data. */
    plocs = (double *) mxGetData(locs);
    nograins = mxGetN(indices); /* Number of grains near this pixel. */
    
    for (k=0;k<nograins;k++){ /* Loop over grains. */
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainconvvals = mxGetCell(prhs[1],2*N+gind-1);
      pgrainconvvals = (double *) mxGetData(grainconvvals);
//       idk = (int) id[gind-1];
      sharedd = mxGetCell(prhs[4],(gind-1)); /* Grains in nhd of grain k. */
      nosharedd = mxGetN(sharedd);
      Nin = mxGetM(grainconvvals);
      for (p=0;p<nosharedd;p++){
        temp[k][p] = pgrainconvvals[p*Nin  + i];
      }
    }

    /* These lines implement the redistribution step in Esedoglu-Otto algorithm: */
    /* Form the "phi" functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
//       printf("%d\n", gind);
      idk1 = (int) id[gind-1];
//       idk1 = (gind);
      sum = 0.0;
     // idk = (int) id[gind-1]; /* id of the k-th grain in the local list. */
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          gind = (int) pindices[ell]; /* Index of grain in list of all grains. */
          //idell = (int) id[gind-1]; /* id of the ell-th grain in the local list. */
          //st = surfacetension(ori[idk-1],ori[idell-1]);
//           idk = (int) id[gind-1];
          sharedd = mxGetCell(prhs[4],(gind-1));/* Grains in nhd of grain k. */
          psharedd = (double *) mxGetData(sharedd);
          nosharedd = mxGetN(sharedd);
          idk1 = (int) id[gind-1];
          
          for (p=0;p<nosharedd;p++){
              if (psharedd[p] == (idk1)) b=p;
          }
//           printf("%f\n", temp[ell][b]);
          st=1.0;
          sum = sum + st*temp[ell][b];
        }
      }
      phi[k] = sum;
//       printf("%f\n", phi[k]);
    }
    
    /* Minimization over the "phi" functions involved in forming level set functions: */
    for (k=0;k<nograins;k++){
      mink = 1e100;
      for (ell=0;ell<nograins;ell++){
        if (ell != k) {
          mink = min( mink , phi[ell] );
        }
      }
      minphi[k] = mink;
    }    

    /* Form the level set functions: */
    for (k=0;k<nograins;k++){
      gind = (int) pindices[k]; /* Index of grain in list of all grains. */
      i = (int) plocs[k]-1; /* Location of pixel within grain's data. */
      grainlevsetvals = mxGetCell(prhs[1],N+gind-1);
      pgrainlevsetvals = (double *) mxGetData(grainlevsetvals);
      
      pgrainlevsetvals[i] = (minphi[k] - phi[k]);
      
//       if (nograins==1) pgrainlevsetvals[i] = temp[k][0];
//       if (nograins==1) pgrainlevsetvals[i] = 1;
      
      idk1 = (int) id[gind-1];
      idk = gind;
      sharedd = mxGetCell(prhs[4],(idk-1)); /* Grains in nhd of grain k. */
      psharedd = (double *) mxGetData(sharedd);
      nosharedd = mxGetN(sharedd);
      
      if (nograins==1){
          for (p=0;p<nosharedd;p++){
              if (psharedd[p] == (idk1)) b=p;
          }
          pgrainlevsetvals[i] = temp[k][b];
      }
    }
    
   } /* (for j). Loop over pixels ends. */
    
}

double surfacetension(double ori1, double ori2)
{
  double ang1, minang, st;
  if (ori2>ori1) ang1 = 6.28 - ori2 + ori1;
  else ang1 = 6.28 - ori1 + ori2;
  minang = min( ang1, fabs(ori1-ori2) );
  
  /* Read-Shockley with Brandon angle */
  /*
  if (minang>1.5) st = 1;
  else st = minang / 1.5 * ( 1-log(minang/1.5) );
  */

  /* Equal surface tensions: */  
  st = 1.0;
  
  return st;
}
