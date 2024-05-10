/*
*
* Copyright [1993-2018] [Charles Kooperberg]
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*   Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*
*/
#include <math.h>
#include <stdio.h>
#define NC 50
#include "R.h"

void F77_NAME(xdsifa)(double[][NC], int *, int *, int *, int *);
void F77_NAME(xdsisl)(double[][NC], int *, int *, int *, double *);
void F77_NAME(xdsidi)(double[][NC], int *, int *, int *, double *, int *, double *, int *);
void F77_NAME(xssort)(double *, double *, int *, int *);

static double knots[NC],coef[NC][4][NC],zheta[NC],czheta,xg[NC],dfunpar[6];
static int nknots,ng3[NC],ng4[NC];

void logcensor(int *idelete, int *iknotauto, double sample[], int nsample[], double bound[], int SorC[], int *ynknots, double yknots[], double ycoef[], double *alpha, double wk[], double wk2[], double logl[]);
static int removeknot(double info[][NC], double coef2[][NC]);
static void fits(double xcoef2[][NC], double xzheta[], double xczheta, double ycoef[], int xiknots[], int xnknots);
static double liter(double info[][NC], double sufficient[][2], double bound[], int SorC[], int nsample[], double sample[], int accuracy, int *itrouble);
static void setbounds(double bound[], double cbound[], int nsample[]);
static double erroradjust(double shift[]);
static double error2(double shift[], double rvr[]);
static double likeli(double candidate[], int nsample[], double sample[], double bound[], int accuracy);
static double linsearch(double shift[], double oldll, double bound[], int nsample[], double sample[], int accuracy);
static double onesearch(double rt, double shift[], int accuracy, double bound[], int *err, int nsample[], double sample[]);
static int numbertester(double aa);
static double expin(int version, double t1, double t2, double a[]);
static double expin2(int version, double t1, double t2, double aa[], double b1, double b0);
static double dens3(double x);
static double dens33(double x);
static double numint(double k1, double k2, int accuracy);
static double numint2(double k1, double k2, int accuracy);
static double numints(double vv[], double k1, double k2, int accuracy, int ip);
static double fun2(double x);
static double fun48(double w, double x, double vv[], int ip);
static void coeff(double coef2[][NC]);
static void start1(double crossprods[][NC], double derivatives[], double sample[], int nsample[]);
static void start2(double crossprods[][NC], double derivatives[], double coef2[][NC], int nkstart, int iremove);
static void suffstat1(double suffcombine[][2], double sample[], int nsample[]);
static void suffstat2(double suffcombine[][2], double coef2[][NC], double sufficient[][2]);
static void knotplace(int iknots[], double rknots[], int iknotauto, double bound[], double sample[], int nsample[], int SorC[], double smp2[], double smp3[], double qt[]);
static int knotnumber(int idelete, int nsample[], int nknots, int SorC[]);
static int piecedens(double sample[], double smp2[], double smp3[], int nsample[]);
static double middle(double info[][NC], double shift[], double sufficient[][2], double bound[], int accuracy, int nsample[], double sample[], double zheta[], int what);
static void intnum2(double x1, double x2, double qolint[][NC +1], double shift[], double info[][NC], int n1, int n2, int what);
static void intnum3(double x, double qolint[][NC +1], double d1, double e1, int vs, double bd, double shift[], double info[][NC], int n1, int n2, int what);
static void intnum4(double x, double qolint[][NC +1], double d1, double e1, int vs, double bd, double shift[], double info[][NC], int n1, int n2, int what);
static double tails(double info[][NC], double shift[], double coef[][4][NC], double bound[], double knots[], double zheta[], int nknots, int what);

/******************************************************************************/
/* this is the main program                                                   */
/* remove follows at the end                                                  */

void logcensor(int *idelete, int *iknotauto, double sample[], int nsample[], double bound[], int SorC[], int *ynknots, double yknots[], double ycoef[], double *alpha, double wk[], double wk2[], double logl[])

/* these quantities are defined in the file where they originated and lhead.h */

{
   int itrouble,accuracy=0,xnknots=0.;
   double qt[2];
/* functions - see the functions themselves */

double info[NC][NC],loglikelihood,sufficient[NC][2],coef2[NC][NC],xczheta=0.;
/* stuf used in iter
   info - the information matrix
   loglikelihood                                  
   sufficient - sufficient statistics
   coef,coef2 - 2 matrices defining the splines as a function of the knots, see
                lcoef.c for the exact definitions.                            */

double derivatives[NC],crossprods[NC][NC],suffcombine[NC][2],xzheta[NC];
/* derivative crossprods used to compute the starting values - see lstart.c    
   suffcombine used to compute the sufficient statistics - see lsuff.c        */

double aic,aicmin,r1,rknots[NC],xcoef2[NC][NC];
/* local double stuff:
   r1 - utility
   aic - akaike information criterion
   aicmin - minimum aic encountered
   alpha - alpha value in aic               
   rknots - copy of knots                                                     */

int i,j,nkstart,iremove=0,iknots[NC],xiknots[NC];
/* local integers
   i,j,k - counter, utility                                                   
   nkstart - number of knots at the beginning of the algorithm
   iremove - number of the knot that is removed                               */



/******************************************************************************/
/* compute the number of knots (to start)                                     */
   for(i=0;i<NC;i++) ng3[i]=0;
   for(i=0;i<NC;i++) ng4[i]=0;
   nknots = *ynknots;
   for(i=0;i<NC;i++)knots[i]=yknots[i];

/* 14 is the minimum sample size                                              */
   if(nsample[0] < 14){
      if(SorC[0] == 0)
         (void)Rprintf("sample is too small\n");
      else
         SorC[0] = 2;
      return;
   }

/* determine the number of starting knots                                     */
   nknots = knotnumber(*idelete,nsample,nknots,SorC);

/* place the knots                                                            */
   knotplace(iknots,rknots,*iknotauto,bound,sample,nsample,SorC,wk,wk2,qt);
   bound[2]=qt[0]+qt[1]*bound[2];
   bound[4]=qt[0]+qt[1]*bound[4];

/* initialize aicmin to a very large value                                    */
   aicmin = pow(10.,100.);

   nkstart = nknots;
   for(i=0;i<nknots;i++)xg[i]=knots[i];
   for(i=nsample[1]+2*nsample[2];i<nsample[0]+nsample[2]-nsample[4];i++){
      for(j=0;knots[j]<sample[i]&&j<nknots;j++);
      if(j!=0&&j!=nknots)
         if(sample[i]-knots[j-1]<knots[j]-sample[i])j=j-1;
      if(j==nknots)j=j-1;
      ng3[j]++;
   }
   for(i=nsample[1]+2*nsample[2]+nsample[3];i<nsample[0]+nsample[2];i++){
      for(j=0;knots[j]<sample[i]&&j<nknots;j++);
      if(j!=0&&j!=nknots)
         if(sample[i]-knots[j-1]<knots[j]-sample[i])j=j-1;
      if(j==nknots)j=j-1;
      ng4[j]++;
   }
   for(i=0;i<nknots;i++)yknots[i]=knots[i];

/* some possible errors                                                       */
   if(SorC[0] == -2)return;
   if(SorC[0] == 647)return;
   if(SorC[0] == 23)return;

/* Compute stuff later used for computing sufficient statistics.              */
   suffstat1(suffcombine,sample,nsample);

/* Compute stuff later used for computing the starting values.                */
   start1(crossprods,derivatives,sample,nsample);

/* Compute coefficient matrix.                                                */
   itrouble = 0;
   do{   
      coeff(coef2); 

/* Compute sufficient statistics.                                             */
      suffstat2(suffcombine,coef2,sufficient);

/* Compute starting values                                                    */
      start2(crossprods,derivatives,coef2,nkstart,iremove);

/* Compute zheta-hat (main iteration loop).                                   */
      if(itrouble == 1){
         for(i=0;i<nknots-1;i++)zheta[i]=0.;
         loglikelihood = liter(info,sufficient,bound,SorC,nsample,sample,5,
         &itrouble);
      }
      else{
      loglikelihood = liter(info,sufficient,bound,SorC,nsample,sample,accuracy,
         &itrouble);
      if(itrouble == 2){
         for(i=0;i<nknots-1;i++)zheta[i]=0.;
         itrouble = 1;
         loglikelihood = liter(info,sufficient,bound,SorC,nsample,sample,5,
         &itrouble);
      }}

/* this means that there were problems (like no convergence)                 */
      SorC[27]=0;
      if(SorC[0] < 0 || itrouble == 5||itrouble==17){
         if(nknots==nkstart){
             SorC[0] = -1;
             SorC[27]=0;
             return;
         }
         if(SorC[0]==-647)
            (void)Rprintf("Smallest number of knots tried is %d\n",nknots+1);
         SorC[27]=nknots;
         nknots=0;
      }
      
/* If we have a fixed number of knots we are done now.           
   Copy the solution in the x-arrays                                         */
      if(nknots>0){
      if(*idelete == 0){
         xczheta = czheta;
         xnknots = nknots;
         for(i=0; i<NC; i++){
            xzheta[i] = zheta[i];
            xiknots[i] = iknots[i];
            for(j=0; j<NC; j++)
               xcoef2[i][j] = coef2[i][j];
         }
         logl[0]=loglikelihood*nsample[0]+log(qt[1])*nsample[1]; /* &&&&& */
      }
   
/* We are deleting knots. Compute aic. Is it an improvement?                  */
      else {
         aic = -2. * loglikelihood * nsample[0] + *alpha * (nknots-1);      
         logl[nknots-1]=loglikelihood*nsample[0]+log(qt[1])*nsample[1]; /* &&&&& */
         if(aic <= aicmin){

/* Then we mark the solution   */
            aicmin = aic;
            xczheta = czheta;
            xnknots = nknots;
            for(i=0; i<NC; i++){
               xzheta[i] = zheta[i];
               xiknots[i] = iknots[i];
               for(j=0;j<NC;j++)
                  xcoef2[i][j]=coef2[i][j];
            }

/* that is we fill it here, iknots contains the rank numbers of the remaining
   knots                                                                      */
         }
         else{
            r1 = -2. * loglikelihood * nsample[0] + *alpha * 2;
            if(r1 >= aicmin) nknots = 0;
         }

      
/* If there were more than 4 knots, we should remove one.                     */
         if(nknots >= 4){

/* Select the one to remove.                                                  */
            iremove = removeknot(info,coef2);

/* Remove its remainders from all sort of arrays. That is, shift the ends of the
   array 1 closer to 0, crossprods gets shifted in 2 directions               */
            for(i=iremove; i<nknots; i++){
               knots[i] = knots[i+1];
               iknots[i] = iknots[i+1];
               derivatives[i+1] = derivatives[i+2];
               suffcombine[i+2][0] = suffcombine[i+3][0];
               suffcombine[i+2][1] = suffcombine[i+3][1];
               crossprods[0][i+1] = crossprods[0][i+2];
               crossprods[i+1][0] = crossprods[i+2][0];
            }
            for(j=0; j<nknots; j++){
               for(i=iremove; (i<nknots && j>=iremove);i++)
                  crossprods[i+1][j+1] = crossprods[i+2][j+2];
               for(i=iremove;(i<nknots && j<iremove);i++)
                  crossprods[i+1][j+1] = crossprods[i+2][j+1];
               for(i=0;(i<iremove && j>=iremove);i++)
                  crossprods[i+1][j+1] = crossprods[i+1][j+2];
            }
         }

/* this takes care that we stop                                               */
         else
           if(nknots == 3)nknots=2;
      }}
   } while(nknots>=3 && *idelete >0);

/* Now write the solution down.   This prints the density                     */
   nknots = nkstart;
   fits(xcoef2,xzheta,xczheta,ycoef,xiknots,xnknots);
   *ynknots = nkstart;
   bound[2]=(bound[2]-qt[0])/qt[1];
   bound[4]=(bound[4]-qt[0])/qt[1];
   for(i=0;i<nkstart;i++)
       yknots[i]=(yknots[i]-qt[0])/qt[1];
   for(i=0;i<nkstart;i++)
       ycoef[i+2]=ycoef[i+2]*qt[1]*qt[1]*qt[1];
   ycoef[0]=ycoef[0]+ycoef[1]*qt[0]+log(qt[1]);
   ycoef[1]=ycoef[1]*qt[1];
}
/******************************************************************************/

static int removeknot(double info[][NC],double coef2[][NC])

/* all described in originating files and lhead.h                             */

{
   double ratmax=0.,wk1[NC],se[NC],phi[NC],x1[2];
   int i,j,infol,k,x2[3],x3=1,kpvt[NC],irmax=1;
/* local:

	   wk1[]      - workarray for linpack
	   kpvt       - workarray for linpack
	   x1,x2,x3   - for linpack
	   infol      - for linpack
	   i j        - counters
	   phi        - linear combination of zhetas
	   se         - standard errors of phi
	   ratmax     - maximum ratio se/phi

	   output:

	   irmax      - index of maximum ratio                                        */

	/* Take linear combinations of zheta such that phi is zheta(phi) for
	   the truncated power basis. (Which is not a basis.)                         */
	   for(i=0;i< nknots-1;i++){
	      phi[i] = 0.;
	      for(j=0;j< nknots-1;j++)
		 phi[i] = phi[i] + zheta[j] * coef2[j][i+2];
	      phi[i] = fabs(phi[i]);
	   }

	/* Invert the information matrix, giving the approximate covariance 
	   matrix for zheta                                                           */
	   i= nknots-1;
	   j=NC;
	   F77_CALL(xdsifa)(info,&j,&i,kpvt,&infol);
	   F77_CALL(xdsidi)(info,&j,&i,kpvt,x1,x2,wk1,&x3);

/* For convenience: symmetarize                                               */
   for(i=1;i< nknots-1;i++)
      for(j=0;j<i;j++)
         info[j][i] = info[i][j];

/* Take linear combinations, to get the standard errors of phi      
   we don't use the tail ones                                                 */
   for(i=1;i< nknots-1;i++){
      se[i] = 0.;
      for(j=0;j< nknots-1;j++)
         for(k=0;k< nknots-1;k++)
            se[i] = se[i] + coef2[j][i+2] * coef2[k][i+2] * info[j][k];
      if(se[i]>0){
      se[i] = sqrt(se[i]);

/* Select for which knot se/phi takes it maximal value */
      if(se[i] > phi[i] * ratmax){
         ratmax = se[i] / phi[i];
         irmax = i;
      }}
   }
   nknots = nknots-1;
   return irmax;
}
/*****************************************************************************/
static void fits(double xcoef2[][50], double xzheta[], double xczheta, double ycoef[], int xiknots[], int xnknots)
{
   int i,j;
   for(i=0;i< NC;i++)ycoef[i]=0.;
   ycoef[0]=-log(xczheta);
   for(i=0;i< xnknots;i++){
      ycoef[xiknots[i]+2] = 0.;
      for(j=0;j< xnknots-1;j++)
         ycoef[xiknots[i]+2] = ycoef[xiknots[i]+2] + xzheta[j] * xcoef2[j][i+2];
   }
   for(j=0;j< xnknots-1;j++){
      ycoef[0]=ycoef[0] + xzheta[j] * xcoef2[j][0];
      ycoef[1]=ycoef[1] + xzheta[j] * xcoef2[j][1];
   }
}
/******************************************************************************/
/* this is the main iteration loop                                            */
/* setbounds follows at the end                                               */


static double liter(double info[][50], double sufficient[][2], double bound[], int SorC[], int nsample[], double sample[], int accuracy, int *itrouble)

{
   int counter=0,infol,i1,i2,i3,kpvt[NC],jaja1,i4,iii[4],ithere,i7,nrc=0,nrc2=0;
   double zerror=0.,oldlikelihood,cbound[7],shift[NC],dd[2];
   double candidate[NC],newlikelihood=0.,work[NC][NC];
   double one,rvr[100],zerrorx;

/* local
   i1 i2 i3    - counters
   oldlikelihood - loglikelihood previous iteration
   newlikelihood - loglikelihood present iteration
   infol,kpvt  - for linpack
   work        - for linpack
   shift       - used for the shift and score
   zerror       - stop criterion
   cbound      - see below
   candidate   - candidate for new zheta
   counter     - number of iterations since last reset of boundaries          */

   one=0.99999;
/* sets integration bounds                                                    */
   setbounds(bound,cbound,nsample);

/* accuracy = 0 means that we approximate censoring                           */
   if(accuracy == 0 && nsample[0]==nsample[1])accuracy=1;

/* start of the iteration                                                     */
   if(SorC[0] == 0) (void)Rprintf("%d \n",nknots);
   for(i1=1; i1<500; i1++){

/* if we go extremely far out in the tails, and zheta[first] or zheta[last]
   still doesn't have the right sign, we fix the bounds at the max value      */
      ithere = 0.;
      if(bound[3]<0.5||bound[3]>1.5)ithere=5;
      else{
         if(bound[4]>=25*knots[nknots-1]-24*knots[0]){
            if(SorC[0]==0&&SorC[24]==0)(void)Rprintf("a very long right tail");
            SorC[24]=nknots;
            *itrouble=17;
            for(i2=1;i2<5;i2++)bound[i2]=cbound[i2];
            return 0.;
         }
      }
      if(bound[1]<0.5||bound[1]>1.5)ithere=ithere+5;
      else{
         if(bound[2]<=-24*knots[nknots-1]+25*knots[0]){
            if(SorC[0]==0&&SorC[23]==0)(void)Rprintf("a very long left tail");
            SorC[23]=nknots;
            *itrouble=17;
            for(i2=1;i2<5;i2++)bound[i2]=cbound[i2];
            return 0.;
         }
      }

/* after 250 iterations we do the censoring exact                             */
      if(accuracy == 0 && i1==250){
         accuracy=1;
         counter = 0;
      }
      jaja1=0;
      counter++;

/* if something changed in the bounds since last time counter = 0, and we have 
   to recompute czheta and the loglikelihood                                  */
      i7=0;
      do{
      if(counter==1){
         do{
            czheta = middle(info,shift,sufficient,bound,accuracy,
                                                        nsample,sample,zheta,0);
/* check against overflow and NAs                                             */
            i2=-1;
            if(czheta>0.0)i2=i2+1;
            if(czheta<2.0)i2=i2+1;
            if(numbertester(czheta)==1)i2=i2-4;
            if(i2<0){
               for(i2=0;i2<nknots-1;i2++)zheta[i2]=zheta[i2]/10.;
               if(i7>7 && *itrouble == 0){
                  *itrouble = 2;
                  for(i2=1;i2<5;i2++)bound[i2]=cbound[i2];
                  return 0.;
               }
               i2=-6;
               i7++;
               if(i7==200){
                  *itrouble=5;  
                  return 1.;
               }
            }
         }while(i2<0);
         newlikelihood=likeli(zheta,nsample,sample,bound,accuracy);
      }
      i2=1000*nknots+i1;
      
/* compute score (stored in shift), info and czheta                           */
      czheta = middle(info,shift,sufficient,bound,accuracy,
                                                        nsample,sample,zheta,2);
           
/* store the loglikelihood of the previous iteration                          */
      oldlikelihood = newlikelihood;

/* copy info                                                                  */
      for(i2=0;i2<nknots-1;i2++)
         for(i3=0;i3<nknots-1;i3++) work[i2][i3] = info[i2][i3];
/* solve the system                                                           */
      i2=NC;
      i3=nknots-1;
      i4 = 100;
      F77_CALL(xdsifa)(work,&i2,&i3,kpvt,&infol);
/* check for underflows overflows and NA                                      */
      if(infol!=0) {for(i2=0;i2<nknots-1;i2++) zheta[i2] = zheta[i2]/100.;
         if(i7>7 && *itrouble == 0){
            *itrouble = 2;
            for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
            return 0.;
         }
         i7++;
         if(i7==200){
            *itrouble = 5;
            return 1.;
         }
      }
      }while(infol!=0);
      iii[1]=0;
      iii[2]=0;
      if(nsample[2]+nsample[3]+nsample[4] >0)
         F77_CALL(xdsidi)(work,&i2,&i3,kpvt,dd,iii,rvr,&i4); 
      for(i4=0;i4<nknots-1;i4++)rvr[i4]=shift[i4];
      iii[3]=iii[1]+iii[2];
      if(iii[2]+iii[1]==0)nrc=0;
      if(iii[2]+iii[1]==0)nrc2=0;
      if(iii[1]+iii[2]>0&&nrc+nrc2>=5&&nrc>=2){
         nrc2++;
         if(nrc2==15 && *itrouble == 0){
              *itrouble = 2;
              for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
              return 0.;
          }
         if(nrc2==48 && *itrouble != 0){
              *itrouble = 5;
              for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
              return 1.;
          }
          iii[2]=0;
          iii[1]=0;
          nrc=0;
      }
      if(iii[2]+iii[1]==0){
         F77_CALL(xdsisl)(work,&i2,&i3,kpvt,shift);  

/* compute the stop criterion and adjust stepsize, if too large               */
         i7=0;
         do{
            zerrorx=error2(shift,rvr);
            zerrorx=zerrorx*100.;
            zerror = erroradjust(shift);
            if(iii[3]==0)zerror=zerrorx;
            i2=-1;
            if(zerror>0.0)i2=i2+2;
            if(zerror<2.0)i2=i2+2;
            if(numbertester(zerror)==1)i2=i2-4;
            if(i2<0){
               for(i2=0;i2<nknots-1;i2++)shift[i2]=shift[i2]/10.;
               if(i7>7 && *itrouble == 0){
                  *itrouble = 2;
                  for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
                  return 0.;
               }
               i7++;
               if(i7==200){
                  *itrouble = 5;
                  return 1.;
               }
               i2=-7;
            }
         }while(i2<0);

/* counter is the number of iterations since the last adjustment of integration
   boundaries. Thus if we find an zerror, we put it to 0 and start all over.
   Temporarily we deduct 10000 from it if we have to half the step size first */
         i7=0;
         do{
            if(counter<-1000)counter=counter+10000;

/* We should check whether zheta's that should be negative (the ones for tail
   basis functions if we are integrating to plus/minus infinity) stay negative.
   If counter become 0 we essentially leave the loop here.....................
   We then go from case A to case B...........................................*/
            if(zheta[0] + shift[0] >= 0. && bound[1] < 0.5){
               cbound[5] = -cbound[5];
               bound[2] = cbound[5] * (knots[0] - knots[1]) + knots[0];
               bound[1] = 1;
               counter = 0;
            }
            if(zheta[nknots-2] + shift[nknots-2] >= 0. && bound[3] < 0.5){
               cbound[6] = -cbound[6];
               bound[4] = cbound[6] * (knots[nknots-1] - knots[nknots-2]) +
                              knots[nknots-1];
               bound[3] = 1.;
               counter = 0;
            }

/* Compute the C(zheta) for the candidate zheta-hat. We can
   then compute the new log-likelihood.                                       */
            if(counter >0){
               for(i2=0;i2<nknots-1;i2++) candidate[i2] = shift[i2] + zheta[i2];
               czheta = middle(info,shift,sufficient,bound,accuracy,
                                                    nsample,sample,candidate,0);
               i2=0;
               if(czheta>0.0)i2=i2+1;
               if(czheta<2.0)i2=i2+1;
               if(numbertester(czheta)==1)i2=i2-4;
               if(i2>0){
               newlikelihood=likeli(candidate,nsample,sample,bound,accuracy);

/* If the loglikelihood really decreases, we step size and go back a bit      */
               if(newlikelihood < oldlikelihood &&
                  ((newlikelihood/oldlikelihood < one) ||
                  (oldlikelihood/newlikelihood < one)) && zerror > 0.00001){
		  jaja1++;
                  if(jaja1<12 || i1<3){
                     for(i2=0;i2<nknots-1;i2++) shift[i2] = shift[i2]/2.;
                     counter = counter - 10000;
                  }
               }
               }
               else{
                  if(i7 > 7 && *itrouble == 0){
                     *itrouble = 2;
                     for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
                     return 0.;
                  }
                  i7++;
                  if(i7==200){
                     *itrouble = 5;
                     return 1.;
                  }
                  for(i2=0;i2<nknots-1;i2++) shift[i2] = shift[i2]/10.;
                  counter = counter - 10000;
               }
            }
         }while(counter <0);
/* Only if the counter is negative we have to go back.                        */
      }
      if(jaja1==25)for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
      if(jaja1==25)return newlikelihood;
      if(iii[2]+iii[1] >0){
         nrc++;
         for(i2=0;i2<nknots-1;i2++) shift[i2] = rvr[i2];
         for(i2=0;i2<nknots-1;i2++) rvr[i2] = zheta[i2];
         (void)erroradjust(shift);
         if(cbound[5]<0. && bound[1] < 0.5 &&
            (cbound[5]>-6.||zheta[0]+shift[0]>0)){
            cbound[5] = -cbound[5];
            bound[2] = cbound[5] * (knots[0] - knots[1]) + knots[0];
            bound[1] = 1;
            counter = 0;
         }
         if(cbound[6]<0. && bound[3] < 0.5 && 
            (cbound[6]>-6.||zheta[nknots-2]+shift[nknots-2]>0)){
            cbound[6] = -cbound[6];
            bound[4] = cbound[6] * (knots[nknots-1] - knots[nknots-2]) +
                           knots[nknots-1];
            bound[3] = 1.;
            counter = 0;
         }
         if(counter ==0){
            czheta = middle(
                   info,shift,sufficient,bound,accuracy,nsample,sample,zheta,0);
            (void)
            middle(info,shift,sufficient,bound,accuracy,nsample,sample,zheta,1);
            oldlikelihood=likeli(zheta,nsample,sample,bound,accuracy);
            for(i2=0;i2<nknots-1;i2++) shift[i2] = rvr[i2];
            for(i2=0;i2<nknots-1;i2++) rvr[i2] = zheta[i2];
            (void)erroradjust(shift);
         }
         oldlikelihood = linsearch(
                             shift,oldlikelihood,bound,nsample,sample,accuracy);
         czheta = middle(
                   info,shift,sufficient,bound,accuracy,nsample,sample,zheta,0);
         (void)
            middle(info,shift,sufficient,bound,accuracy,nsample,sample,zheta,1);
         (void)erroradjust(shift);
         oldlikelihood = linsearch(
                             shift,oldlikelihood,bound,nsample,sample,accuracy);
         for(i2=0;i2<nknots-1;i2++) shift[i2]= zheta[i2] - rvr[i2];
         oldlikelihood = linsearch(
                             shift,oldlikelihood,bound,nsample,sample,accuracy);
         czheta = middle(
                   info,shift,sufficient,bound,accuracy,nsample,sample,zheta,0);
         for(i2=0;i2<nknots-1;i2++) shift[i2]= 0.;
         for(i2=0;i2<nknots-1;i2++)rvr[i2]=zheta[i2]-rvr[i2];
         zerror = 100.*erroradjust(rvr);
         newlikelihood=oldlikelihood;
      }

/* O.K. the new zheta is good, let's store it.                                */
      if(counter > 0){
         for(i2=0;i2<nknots-1;i2++) zheta[i2] = zheta[i2] + shift[i2];

/* Is there convergence? If this is not the final iteration, we're not picky. */
         if(ithere<10) zerror=zerror/100.;
         if(zerror < 0.000001){

/* If we were integrating to +/-infinity or if we exceeded the bounds, done.  */
            if(ithere >9){
                if(accuracy == 5){
                   for(i3=1;i3<5;i3++)bound[i3]=cbound[i3];
                   rvr[99] = middle(
                          info,shift,sufficient,bound,5,nsample,sample,zheta,0);
                   if(iii[1]+iii[2]!=0)
                   rvr[99] = middle(
                          info,shift,sufficient,bound,5,nsample,sample,zheta,0);
                   i2=0;
                   if(rvr[99]>0.0)i2=i2+1;
                   if(rvr[99]<2.0)i2=i2+1;
                   if(numbertester(rvr[99])==1)i2=i2-4;
                   if(i2>0){
                      czheta=rvr[99];
                      newlikelihood=likeli(zheta,nsample,sample,bound,5);
                   }
                   return newlikelihood;
                } 
                else {
                   accuracy = 5;
                   counter = 0;
                }
            }
            else {
          
/* If we were not integrating to +/- infinity, we either double the
   integration tails (if a tail zheta was possitive) or we go for it:
   +/- infinity (That means cbound[5] and cbound[6] become negative.)         */
               counter=0;
/* Left tail                                                                  */
               if(cbound[5] < 0)                                /* case A     */
                  cbound[5] = 2. * cbound[5];
               else{        
                  if(cbound[1] < 0.5 && zheta[0] >=0.){         /* case B     */
                     cbound[5] = 2. * cbound[5];
                     bound[2] = cbound[5] * (knots[0] - knots[1]) + knots[0];
                  }
                  if(cbound[1] < 0.5 && zheta[0] <0. ){         /* case B to A*/
                     cbound[5] = -2. * cbound[5];
                     bound[1] = 0;
                  }
                  if(cbound[1] > 0.5 && bound[1] > 0.5 && bound[1] < 1.5){ /*E*/
                     cbound[5] = 2. * cbound[5];
                     bound[2] = cbound[5] * (knots[0] - knots[1]) + knots[0];
                     if(bound[2] < cbound[2]){                  /* goto C/D   */
                        bound[2] = cbound[2];
                        bound[1] = 2.;   
                     }
                  }
               }
/* Right tail                                                                 */
               if(cbound[6] < 0)                                /* case A     */
                  cbound[6] = 2. * cbound[6];
               else{        
                  if(cbound[3] < 0.5 && zheta[nknots-2] >=0.){  /* case B     */
                     cbound[6] = 2. * cbound[6];
                     bound[4] = cbound[6] * (knots[nknots-1] - knots[nknots-2])
                                                             + knots[nknots-1];
                  }
                  if(cbound[3] < 0.5 && zheta[nknots-2] <0. ){  /* case B to A*/
                     cbound[6] = -2. * cbound[6];
                     bound[3] = 0.;
                  }
                  if(cbound[3] > 0.5 && bound[3] > 0.5 &&
                                                    bound[3] < 1.5){ /* case E*/
                     cbound[6] = 2. * cbound[6];
                     bound[4] = cbound[6] * (knots[nknots-1] - knots[nknots-2])
                                                             + knots[nknots-1];
                     if(bound[4] > cbound[4]){                  /* goto C/D   */
                        bound[4] = cbound[4];
                        bound[3] = 2;   
                     }
                  }
               }
            }
         }
      }
   }

/* If we ended up here, there was no convergence in 300 iterations.           */
   if(SorC[0] == 0){
      (void)Rprintf("no convergence was achieved with %d knots\n",nknots);
         SorC[0] = -647;
   }
   else
      SorC[0] = -SorC[0];
   return newlikelihood;
}
static void setbounds(double bound[], double cbound[], int nsample[])
{
/* set the integration boundaries. O.k., this is quite complicated. There are 2
   arrays that determine how far we are integrating in each tail.
   For the lower tail:

   bound[1] - is there a lower bound to which we should integrate (on entry)
          are we right now integrating to a lower bound (during)
          is there a lower bound to which we should integrate (on exit)
          (0 = -infinity; 1 = lower bound);
   cbound[1]- copies the begin value of bound[1] for during the iterations.
   bound[2] - if bound[1] == 1/2 the lower bound of integration (on entry)
          if bound[1] == 1/2 the lower bound of integration (during)
          if bound[1] == 1/2 the lower bound of integration (on exit)
   cbound[2]- copies the begin value of bound[2] for during the integration.
   cbound[5] - technical limit on the integration. If cbound[5] < 0 we can
          integrate to -infinity (if we want to), if cbound[5] > 0 the
          furthest we want to integrate is to
          cbound[5]*(knots[1]-knots[2])+knots[1];
          if cbound[5] is smaller than 0 it is twice its last possitive value.
   For the upper tail: bound[3], bound[4], cbound[3], cbound[4], cbound[6]    */

/* the lower tail                                                             */
   cbound[2] = bound[2];
   cbound[1] = bound[1];
   cbound[5] = 1.;
   if(nsample[3]*3>nsample[1])cbound[5]=0.5;
   bound[2] = cbound[5] * (knots[0] - knots[1]) + knots[0];

/* there are 5 possible situations:
A  cbound[1]=0, zheta[0]<0    need: cbound[5]=-1, bound[1]=0
B               zheta[0]>=0   need: cbound[5]=1,  bound[1]=1
C  cbound[1]=1, zheta[0]<0    need: cbound[5]=1,  bound[1]=2, bound[2]=cbound[2]
D               zheta[0]>=0, bound[2]<cbound[2]:
                                     cbound[5]=1, bound[1]=2, bound[2]=cbound[2]
E                            bound[2]>=cbound[2]:
                                     cbound[5]=1,  bound[1]=1              */
   if(zheta[0]<0. && cbound[1] < 0.5 && nsample[0] == nsample[1])
      cbound[5] = -1.;        
   if(cbound[1] < 0.5 && (zheta[0] >= 0. || nsample[0]!=nsample[1]))
      bound[1] = 1;
   if(cbound[1] > 0.5) bound[1] = 2;
   if(cbound[1] > 0.5 && zheta[1] >= 0. && bound[2] >= cbound[2]) bound[1] = 1;
   if(cbound[1] > 0.5 && (zheta[0]<0. || (zheta[0] >= 0. && 
                   bound[2] < cbound[2]))) bound[2] = cbound[2];

/* the upper tail almost similar                                              */
   cbound[4] = bound[4];
   cbound[3] = bound[3];
   cbound[6] = 1;
   if(nsample[4]*3>nsample[1])cbound[6]=0.5;
   bound[4] = cbound[6] * (knots[nknots-1] - knots[nknots-2]) + knots[nknots-1];
   if(zheta[nknots-2]<0. && cbound[3] < 0.5 && nsample[0] == nsample[1])
      cbound[6] = -cbound[6]; 
   if(cbound[3] < 0.5 && (zheta[nknots-2] >=0. || nsample[0]!=nsample[1]))
      bound[3] = 1;
   if(cbound[3] > 0.5) bound[3] = 2;
   if(cbound[3] > 0.5 && zheta[nknots-2] >= 0. && cbound[4] >= bound[4])
      bound[3] = 1;
   if(cbound[3] > 0.5 && (zheta[nknots-2]<0. || (zheta[nknots-2] >= 0. &&
              bound[4] > cbound[4]))) bound[4] = cbound[4];
      
}
/******************************************************************************/
/* this file contains a few miscelaneous routines

   erroradjust - computes the error criterion
   likeli      - computes the loglikelihood
   linsearch   - computes the winning stepsize for steepest decent
   onesearch   - computes czheta and the loglikelihood for one zheta+shift    */

/******************************************************************************/

/* this function computes the stopcriterion (zerror) and adjust the stepsize 
   (shift) if this is too large.                                              */

static double erroradjust(double shift[])

{
   double r1,r2;
   int i;
/* all utility numbers */

/* the zerror is the sum of (shift/zheta)^2, except where zheta is very small  */
   r1 = 0.;
   for(i=0; i<nknots-1; i++){
      r2 = zheta[i] * zheta[i];
      if(r2 < 0.0000001) r2 = 0.0000001;
      r1 = r1 + shift[i] * shift[i] / r2;
   }
   if(numbertester(r1)!=1)r2 = sqrt(r1);

/* if the zerror is infinity we compute an alternate zerror criterion to prevent
   trouble                                                                    */
   else{
      r1 = 0.;
      for(i=0; i<nknots-1; i++){
         r2 = fabs(zheta[i]);
         if(r2<0.00001) r2 = 0.00001;
         r1 = r1 + fabs(shift[i]) / r2; 
      }
      r2=r1;
   }

/* If the zerror would be too big, we only take a smaller step: this
   is to prevent overshoot.                                                   */
   for(i=0;i<nknots-1 && r1 > 1000.;i++) shift[i]=shift[i] * 3. / r2;

   return r1;
}
static double error2(double shift[], double rvr[])
{
   int i,j;
   double r=0.;
   for(i=0;i<nknots-1;i++)r=r+shift[i]*rvr[i]/2.;
   j=-1;
   if(r>=0.0)j++;
   if(numbertester(r)==1)j=j-2;
   if(j<0)r=1000;
   return r;
}
/******************************************************************************/

/* this routine computes the loglikelihood                                    */

static double likeli(double candidate[], int nsample[], double sample[], double bound[], int accuracy)
{
   double r0,r1,likl,r3[NC+1],aa[6],bb[6];
   int i1,i2,i3,i4,iv,iw;
   r0=exp((double)(-740));

/* the function uses numint and dens33 - which always use zheta, but the
   routine doesn't always compute the likelihood in zheta, sometimes at
   another place - so we have to swap zheta and candidate                     */
   for(i1=0;i1<(nknots-1);i1++){
      r1=zheta[i1];
      zheta[i1]=candidate[i1];
      candidate[i1]=r1;
   }

/* the stuff for the exact data is easy                                       */
   likl=0.;
   for(i1=0;i1<nsample[1];i1++){
         likl=likl+dens33(sample[i1]);
   }

/* compute the integral of the density between each two knots in r3           */
   if(bound[3]>0.5)iv=4;
   else iv=3;
   if(bound[1]>0.5)iw=2;
   else iw=1;
   aa[1]=0.;
   aa[2]=0.;
   aa[3]=1.;
   aa[4]=zheta[nknots-2]*coef[nknots-2][1][nknots];
   aa[5]=zheta[nknots-2]*coef[nknots-2][0][nknots]
     +zheta[nknots-3]*coef[nknots-3][0][nknots]-log(czheta);
   if(nknots<4)
   aa[5]=zheta[nknots-2]*coef[nknots-2][0][nknots]-log(czheta);
   r3[nknots]=expin(iv,knots[nknots-1],bound[4],aa);
   bb[1]=0.;
   bb[2]=0.;
   bb[3]=1.;
   bb[4]=zheta[0] * coef[0][1][0];
   bb[5]=zheta[0] * coef[0][0][0]-log(czheta);
   r3[0]=expin(iw,knots[0],bound[2],bb);
   for(i1=1;i1<nknots;i1++)
      r3[i1]=numint(knots[i1-1],knots[i1],accuracy);

/* depending on the accuracy we either deal with the censored data approximate
   or exact                                                                   */
   if(accuracy!=0){
/* exact first the rounded data                                               */
      for(i1=0;i1<nsample[2];i1++){
         i2=i1*2+nsample[1];
/* between which knots are the two endpoints                                  */
         for(i3=0;knots[i3]<sample[i2] && i3 < nknots;i3++);
         for(i4=0;knots[i4]<sample[i2+1] && i4 < nknots;i4++);
/* is it the same interval                                                    */
         if(i3==i4)
            likl=likl+log(numint(sample[i2],sample[i2+1],0));
/* or a different interval                                                    */
         else{
            r1=numint(sample[i2],knots[i3],0);
            r1=r1+numint(knots[i4-1],sample[i2+1],0);
            for(i3=i3+1;i3<i4;i3++)r1=r1+r3[i3];
            likl=likl+log(r1);
         }
      }
/* the right censored data                                                    */
      for(i1=0;i1<nsample[3];i1++){
         i2=i1+nsample[1]+2*nsample[2];
/* in which interval is the datapoint                                         */
         for(i3=0;knots[i3]<sample[i2] && i3 < nknots;i3++);
/* somewhere in the middle                                                    */
         if(i3<nknots){
            r1=numint(sample[i2],knots[i3],0);
            for(i3=i3+1;i3<=nknots;i3++) r1=r1+r3[i3];
/* if r1<=r0 the likl is effectively -infinity                                */
            if(r1>r0) likl=likl+log(r1);
            else likl=likl-1000;
         }
/* in the tail                                                                */
         else likl=likl+log(expin(iv,sample[i2],bound[4],aa));
      }
/* the left censored data                                                    */
      for(i1=0;i1<nsample[4];i1++){
         i2=i1+nsample[1]+2*nsample[2]+nsample[3];
/* in which interval is the datapoint                                         */
         for(i3=0;knots[i3]<sample[i2] && i3 < nknots;i3++);
/* somewhere in the middle                                                    */
         if(i3>0){
            r1=numint(knots[i3-1],sample[i2],0);
            for(i4=0;i4<i3;i4++) r1=r1+r3[i4];
/* if r1<=r0 the likl is effectively -infinity                                */
            if(r1>r0) likl=likl+log(r1);
            else likl=likl-1000;
         }
/* in the tail                                                                */
         else likl=likl+log(expin(iw,sample[i2],bound[2],bb));
      }
   }
   else{
/* approximate, first right censored - essentially as above   
   there are ng3 points at xg                                                 */
      if(nsample[3]>0){
         i1=0;
         for(i2=0;i2<NC;i2++){
            if(ng3[i2]>0){
               i1=i1+ng3[i2];
               for(i3=0;knots[i3]<xg[i2] && i3<nknots;i3++);
               if(i3<nknots){
                  r1=numint(xg[i2],knots[i3],0);
                  for(i3=i3+1;i3<=nknots;i3++) r1=r1+r3[i3];
                  if(r1>r0) likl=likl+ng3[i2]*log(r1);
                  else likl=likl-1000.*ng3[i2];
               }
               else likl=likl+ng3[i2]*log(expin(iv,xg[i2],bound[4],aa));
               if(i1==nsample[3])i2=NC+3;
            }   
         }   
      }   
/* approximate, now left censored - essentially as above   
   there are ng4 points at xg                                                 */
      if(nsample[4]>0){
         i1=0;
         for(i2=0;i2<NC;i2++){
            if(ng4[i2]>0){
               i1=i1+ng4[i2];
               for(i3=0;knots[i3]<xg[i2] && i3<nknots;i3++);
               if(i3>0){
                  r1=numint(knots[i3-1],xg[i2],0);
                  for(i4=0;i4<i3;i4++) r1=r1+r3[i4];
                  if(r1>r0) likl=likl+ng4[i2]*log(r1);
                  else likl=likl-1000.*ng4[i2];
               }
               else likl=likl+ng4[i2]*log(expin(iw,xg[i2],bound[2],bb));
               if(i1==nsample[4])i2=NC+3;
            }   
         }   
      }   
/* the interval censored stuff is easy                                        */
      if(nsample[2]>0){
         for(i1=0;i1<nsample[2];i1++){
            r1=(sample[nsample[1]+i1*2]+sample[nsample[1]+i1*2+1])/2.;
            likl=likl+dens33(r1);
         }
      }
   }
/* swap zheta and candidate back                                              */
   for(i1=0;i1<(nknots-1);i1++){
      r1=zheta[i1];
      zheta[i1]=candidate[i1];
      candidate[i1]=r1;
   }
   return likl/(double)nsample[0];
}
/******************************************************************************/

/* this routine does a steepest decent search in the direction shift          */

static double linsearch(double shift[], double oldll, double bound[], int nsample[], double sample[], int accuracy)
{
   double rt,rl=1000.,rr=1000.,maxll,maxrt,ll;
   int i,err;
   if(bound[1]<0.5 && shift[0]*zheta[0]<0)rl=-zheta[0]/shift[0];
   if(bound[3]<0.5 && shift[nknots-2]*zheta[nknots-2]<0)
                           rr=-zheta[nknots-2]/shift[nknots-2];
/* rl and rr are the maximum stepsize, given the sign of zheta and integration
   boundaries                                                                 */
   if(rl<rr)rr=rl;
   maxll=oldll;
   maxrt=-100;
/* rt is the first stepsize guess, on a 2log basis                            */
   if(rr>0.5)rt=-2.;
   else rt=(double)floor(log(rr)/log((double)2))-2.;
/* onesearch computes the loglikelihood for one rt (stepsize)                 */
   ll=onesearch(rt,shift,accuracy,bound,&err,nsample,sample);
/* if for the basis stepsize the likelihood goes up we try more               */
   if(err==0 && ll>maxll){
      do{
         maxrt=rt;
         maxll=ll;
         rt=rt+2.;
         ll=onesearch(rt,shift,accuracy,bound,&err,nsample,sample);
      }while(rt<8.5 && ll>maxll && err==0);
   }
/* if for the basis stepsize the likelihood goes down we try less             */
   else {
      do{
         rt=rt-2.;
         ll=onesearch(rt,shift,accuracy,bound,&err,nsample,sample);
      }while(rt>-14.1 && (ll < maxll || err==1));
   }
/* write down the winning combination                                         */
   if(err==0 && ll>maxll){
      maxrt=rt;
      maxll=ll;
   }
   if(maxrt > -50.){
      rt = maxrt + 1;
      ll=onesearch(rt,shift,accuracy,bound,&err,nsample,sample);
      if(err==0 && ll>maxll){
         maxrt=rt;
         maxll=ll;
      }
      maxrt=pow(2.,maxrt);
      for(i=0;i<nknots-1;i++)zheta[i]=zheta[i]+maxrt*shift[i];
   }
   return maxll;
}
   
/******************************************************************************/

/* this routine computes the loglikelihood in zheta+2^rt shift                */

static double onesearch(double rt, double shift[], int accuracy, double bound[], int *err, int nsample[], double sample[])
{
   int i,j;
   double rv[NC],ll,dm1[NC][NC],dm2[NC],dm3[NC][2];
   *err=0;
   rt=pow(2.,rt);
   for(i=0;i<nknots-1;i++)rv[i]=zheta[i]+rt*shift[i];
/* compute czheta                                                             */
   czheta=middle(dm1,dm2,dm3,bound,accuracy,nsample,sample,rv,0);
/* zerror conditions                                                           */
   if(czheta<=0){
      *err=1;
      return 0.;
   }
   if(log(czheta)>200){
      *err=1;
      return 0.;
   }
   j=0;
   if(czheta > -5)j++; 
   if(czheta < 5)j++; 
   if(j==0){
      *err=1;
      return 0.;
   }
/* compute likeli                                                             */
   ll=likeli(rv,nsample,sample,bound,accuracy);
   return ll;
}
static int numbertester(double aa)
/* if aa = -Inf: 0
      aa = +Inf: 1
      aa =  NaN: 2
      otherwise: 3 */
{
   int i1=0,i2=0,i3=0,i4=0;
   if(aa< 2.){
      i1=1;
   }
   if(aa> 0.){
      i2=1;
   }
   if(aa< pow(10.,200.)){
      i3=1;
   }
   if(aa> -pow(10.,200.)){
      i4=1;
   }
   if(i1+i2+i3+i4>=3){
      return 3;
   }
   if(i2==1 && i4==1){
      return 1;
   }
   if(i1==1 && i3==1){
      return 0;
   }
   return 2;
}
/******************************************************************************/
/* this file contains the following functions:
expin  - computes an exponential integral
expin2 - computes another exponential integral
dens3  - computes a logspline-density in one point
dens33 - computes the log of a logspline-density in one point
numint - computes a numerical integral
numints- computes a vector of numerical integrals                             */

/******************************************************************************/
/*   this function computes (analytical) the integral:   
                                                          
      t1 t1 inf t2                                         
      /   /  /  /     2.          a4*x+a5                       
      |   |  |  |(a1*x +a2*x+a3)*e       dx                        
      /   /  /  /                                            
    -inf t2 t1 t1                                              
      1   2.  3  4     <==== version                              */
                                                                 

static double expin(int version, double t1, double t2, double a[])

/* input: a,b,c,d,e,t1,t2,version: see figure above
   local: a1,b1,c1: as a,b and c, but for the primitive.                      
          f1,f2 half-products                                    */
/* this version does not contain much information - the best way to figure out
   what is happening is to compute the integrals above, and then check below  */

{
   double a1,b1,c1,f1,f2;
   int i1,i2;
/* if d=0 and version is 1 or 3 well...........                               */
   if(a[4]!=0 || version == 1 || version ==3){
      a1 = a[1]/a[4];
      b1 = (a[2]-2*a1)/a[4];
      c1 = (a[3]-b1)/a[4];
      f1 = a1*t1*t1+b1*t1+c1; 
      i1 = 1;
      if(f1<0) i1 = -1;
      f1 = log(fabs(f1)) + a[4]*t1+a[5];
      if(f1>2000.) f1=2000.;
      if(version==2 || version == 4){
         f2 = a1*t2*t2+b1*t2+c1;
         i2 = 1;
         if(f2<0) i2= -1;
         f2 = log(fabs(f2)) + a[4]*t2+a[5];
         if(f2 > 2000.) f2=2000.;
         if(version == 2)return i1*exp(f1)-i2*exp(f2);
         return i2*exp(f2)-i1*exp(f1);
      } 
      if(version==1) return i1*exp(f1);
      return -i1*exp(f1);
   }
   a1 = (((a[1]/3)*(t2-t1)+a[2]/2)*(t2-t1)+a[3])*(t2-t1)*exp(a[5]);
   if(version==4)return a1;
   return -a1;
}

/******************************************************************************/
/* This function computes a similar integral, but with a higher order leading
polinomial                                                                    */

static double expin2(int version, double t1, double t2, double aa[], double b1, double b0)

/* input: aa,b1,b0,t1,t2,version: see figure above
   local: u6,u5,u4,u3,u2,u1,u0: as a4,a3,a2,a1,a0, but for the primitive.            
          f1,f2 half-products                                    */
/* this version does not contain much information - the best way to figure out
   what is happening is to compute the integrals above, and then check below  */

{
   double u6,u5,u4,u3,u2,u1,u0,f1,f2;
   int i1,i2;
/* if b1=0 and version is 1 or 3 well...........                              */
   if(b1!=0 || version ==1 || version == 3){
      u6 = aa[6]/b1;
      u5 = (aa[5]-6*u6)/b1;
      u4 = (aa[4]-5*u5)/b1;
      u3 = (aa[3]-4*u4)/b1;
      u2 = (aa[2]-3*u3)/b1;
      u1 = (aa[1]-2*u2)/b1;
      u0 = (aa[0]-u1)/b1;
      f1 = (((((u6*t1+u5)*t1+u4)*t1+u3)*t1+u2)*t1+u1)*t1+u0;
      i1 = 1;
      if(f1<0) i1 = -1;
      f1 = log(fabs(f1)) + b1*t1+b0;
      if(f1>2000.) f1=2000.;
      if(version==2 || version == 4){
         f2 = (((((u6*t2+u5)*t2+u4)*t2+u3)*t2+u2)*t2+u1)*t2+u0;
         i2 = 1;
         if(f2<0) i2= -1;
         f2 = log(fabs(f2)) + b1*t2+b0;
         if(f2 > 2000.) f2=2000.;
         if(version == 2)return i1*exp(f1)-i2*exp(f2);
         return i2*exp(f2)-i1*exp(f1);
      } 
      if(version==1) return i1*exp(f1);
      return -i1*exp(f1);
   }
   u6 = (((aa[6]/7*(t2-t1)+aa[5]/6)*(t2-t1)+aa[4]/5)*(t2-t1)+aa[3]/4)*(t2-t1);
   u6 = (((u6 + aa[2]/3)*(t2-t1)+aa[1]/2)*(t2-t1)+aa[0])*(t2-t1)*exp(b0);
   if(version==4)return u6;
   return -u6;
}
/******************************************************************************/
/* this function return a value from a logspline density                      */

static double dens3(double x)
/* point of interest                                                          */

{
      return exp(dens33(x));
}

/******************************************************************************/
/* this function return the log of a value from a logspline density           */
static double dens33(double x)
/* point of interest                                                          */

{
   int j,k;
   double f;

/* circle through the knots                                                   */
      for(j=0; (j<nknots && x > knots[j]); j++); 
      f = -log(czheta);

/* find in between which knots the point is, and compute the spline           */
      for(k=j-3; k<j+1; k++){
         if(k>=0 && k<=nknots){
            f = f + zheta[k] * (coef[k][0][j] + x * (coef[k][1][j] 
                  + x * (coef[k][2][j] + x * coef[k][3][j])));
         }
      }
      return f;
}
/******************************************************************************/
/* Numerical integration using gaussian quadrature.
   See Abromowitz and Stegun.

   k1 and k2: lower and upper integration bounds
   fun      : function to be integrated          */


static double numint(double k1, double k2, int accuracy)

/* Intgerals using Gauss-Legendre quadrature with 12 points    
   y1,y2,... - abisces
   w1,w2,... - weight
   accuracy  - accuracy
   r1 and r2 - from (k1,k2) to (-1,1)         */

{
   double r1,r2,w[33],y[33];
   int i;
   r1 = ((k2 - k1) / 2);
   r2 = ((k2 + k1) / 2);

   if(accuracy == 0 || accuracy == 1){
      w[1] = 0.467913934572691 * r1; y[1] = 0.238619186083197 * r1;
      w[2] = 0.360761573048139 * r1; y[2] = 0.661209386466265 * r1;
      w[3] = 0.171324429379170 * r1; y[3] = 0.932469514203152 * r1;
      r1 = 0.;
      for(i=1;i<4;i++)r1 = r1 + w[i]*(dens3(r2-y[i])+dens3(r2+y[i]));
      return r1;
   }

   w[ 1]=  0.00178328072169643 * r1; y[1 ]=  0.99930504173577217 * r1;
   w[ 2]=  0.00414703326056247 * r1; y[2 ]=  0.99634011677195533 * r1;
   w[ 3]=  0.00650445796897836 * r1; y[3 ]=  0.99101337147674429 * r1;
   w[ 4]=  0.00884675982636395 * r1; y[4 ]=  0.98333625388462598 * r1;
   w[ 5]=  0.01116813946013113 * r1; y[5 ]=  0.97332682778991098 * r1;
   w[ 6]=  0.01346304789671864 * r1; y[6 ]=  0.96100879965205377 * r1;
   w[ 7]=  0.01572603047602472 * r1; y[7 ]=  0.94641137485840277 * r1;
   w[ 8]=  0.01795171577569734 * r1; y[8 ]=  0.92956917213193957 * r1;
   w[ 9]=  0.02013482315353021 * r1; y[9 ]=  0.91052213707850282 * r1;
   w[10]=  0.02227017380838325 * r1; y[10]=  0.88931544599511414 * r1;
   w[11]=  0.02435270256871087 * r1; y[11]=  0.86599939815409277 * r1;
   w[12]=  0.02637746971505466 * r1; y[12]=  0.84062929625258032 * r1;
   w[13]=  0.02833967261425948 * r1; y[13]=  0.81326531512279754 * r1;
   w[14]=  0.03023465707240248 * r1; y[14]=  0.78397235894334139 * r1;
   w[15]=  0.03205792835485155 * r1; y[15]=  0.75281990726053194 * r1;
   w[16]=  0.03380516183714161 * r1; y[16]=  0.71988185017161088 * r1;
   w[17]=  0.03547221325688239 * r1; y[17]=  0.68523631305423327 * r1;
   w[18]=  0.03705512854024005 * r1; y[18]=  0.64896547125465731 * r1;
   w[19]=  0.03855015317861563 * r1; y[19]=  0.61115535517239328 * r1;
   w[20]=  0.03995374113272034 * r1; y[20]=  0.57189564620263400 * r1;
   w[21]=  0.04126256324262353 * r1; y[21]=  0.53127946401989457 * r1;
   w[22]=  0.04247351512365359 * r1; y[22]=  0.48940314570705296 * r1;
   w[23]=  0.04358372452932345 * r1; y[23]=  0.44636601725346409 * r1;
   w[24]=  0.04459055816375657 * r1; y[24]=  0.40227015796399163 * r1;
   w[25]=  0.04549162792741814 * r1; y[25]=  0.35722015833766813 * r1;
   w[26]=  0.04628479658131442 * r1; y[26]=  0.31132287199021097 * r1;
   w[27]=  0.04696818281621002 * r1; y[27]=  0.26468716220876742 * r1;
   w[28]=  0.04754016571483031 * r1; y[28]=  0.21742364374000708 * r1;
   w[29]=  0.04799938859645831 * r1; y[29]=  0.16964442042399283 * r1;
   w[30]=  0.04834476223480295 * r1; y[30]=  0.12146281929612056 * r1;
   w[31]=  0.04857546744150343 * r1; y[31]=  0.07299312178779904 * r1;
   w[32]=  0.04869095700913972 * r1; y[32]=  0.02435029266342443 * r1;
   r1 = 0.;
   for(i=1;i<33;i++)r1 = r1 + w[i]*(dens3(r2-y[i])+dens3(r2+y[i]));
   return r1;
}
/******************************************************************************/
/* Numerical integration using gaussian quadrature.
   See Abromowitz and Stegun.

   k1 and k2: lower and upper integration bounds
   fun      : function to be integrated          */


static double numint2(double k1,double k2,int accuracy)

/* Intgerals using Gauss-Legendre quadrature with 12 points    
   y1,y2,... - abisces
   w1,w2,... - weight
   accuracy  - accuracy
   r1 and r2 - from (k1,k2) to (-1,1)         */

{
   double r1,r2,w[33],y[33];
   int i;
   r1 = ((k2 - k1) / 2);
   r2 = ((k2 + k1) / 2);

   if(accuracy == 0 || accuracy == 1){
      w[1] = 0.467913934572691 * r1; y[1] = 0.238619186083197 * r1;
      w[2] = 0.360761573048139 * r1; y[2] = 0.661209386466265 * r1;
      w[3] = 0.171324429379170 * r1; y[3] = 0.932469514203152 * r1;
      r1 = 0.;
      for(i=1;i<4;i++)r1 = r1 + w[i]*(fun2(r2-y[i])+fun2(r2+y[i]));
      return r1;
   }

   w[ 1]=  0.00178328072169643 * r1; y[1 ]=  0.99930504173577217 * r1;
   w[ 2]=  0.00414703326056247 * r1; y[2 ]=  0.99634011677195533 * r1;
   w[ 3]=  0.00650445796897836 * r1; y[3 ]=  0.99101337147674429 * r1;
   w[ 4]=  0.00884675982636395 * r1; y[4 ]=  0.98333625388462598 * r1;
   w[ 5]=  0.01116813946013113 * r1; y[5 ]=  0.97332682778991098 * r1;
   w[ 6]=  0.01346304789671864 * r1; y[6 ]=  0.96100879965205377 * r1;
   w[ 7]=  0.01572603047602472 * r1; y[7 ]=  0.94641137485840277 * r1;
   w[ 8]=  0.01795171577569734 * r1; y[8 ]=  0.92956917213193957 * r1;
   w[ 9]=  0.02013482315353021 * r1; y[9 ]=  0.91052213707850282 * r1;
   w[10]=  0.02227017380838325 * r1; y[10]=  0.88931544599511414 * r1;
   w[11]=  0.02435270256871087 * r1; y[11]=  0.86599939815409277 * r1;
   w[12]=  0.02637746971505466 * r1; y[12]=  0.84062929625258032 * r1;
   w[13]=  0.02833967261425948 * r1; y[13]=  0.81326531512279754 * r1;
   w[14]=  0.03023465707240248 * r1; y[14]=  0.78397235894334139 * r1;
   w[15]=  0.03205792835485155 * r1; y[15]=  0.75281990726053194 * r1;
   w[16]=  0.03380516183714161 * r1; y[16]=  0.71988185017161088 * r1;
   w[17]=  0.03547221325688239 * r1; y[17]=  0.68523631305423327 * r1;
   w[18]=  0.03705512854024005 * r1; y[18]=  0.64896547125465731 * r1;
   w[19]=  0.03855015317861563 * r1; y[19]=  0.61115535517239328 * r1;
   w[20]=  0.03995374113272034 * r1; y[20]=  0.57189564620263400 * r1;
   w[21]=  0.04126256324262353 * r1; y[21]=  0.53127946401989457 * r1;
   w[22]=  0.04247351512365359 * r1; y[22]=  0.48940314570705296 * r1;
   w[23]=  0.04358372452932345 * r1; y[23]=  0.44636601725346409 * r1;
   w[24]=  0.04459055816375657 * r1; y[24]=  0.40227015796399163 * r1;
   w[25]=  0.04549162792741814 * r1; y[25]=  0.35722015833766813 * r1;
   w[26]=  0.04628479658131442 * r1; y[26]=  0.31132287199021097 * r1;
   w[27]=  0.04696818281621002 * r1; y[27]=  0.26468716220876742 * r1;
   w[28]=  0.04754016571483031 * r1; y[28]=  0.21742364374000708 * r1;
   w[29]=  0.04799938859645831 * r1; y[29]=  0.16964442042399283 * r1;
   w[30]=  0.04834476223480295 * r1; y[30]=  0.12146281929612056 * r1;
   w[31]=  0.04857546744150343 * r1; y[31]=  0.07299312178779904 * r1;
   w[32]=  0.04869095700913972 * r1; y[32]=  0.02435029266342443 * r1;
   r1 = 0.;
   for(i=1;i<33;i++)r1 = r1 + w[i]*(fun2(r2-y[i])+fun2(r2+y[i]));
   return r1;
}
/***************************************/
static double numints(double vv[],double k1,double k2,int accuracy,int ip)

/* Intgerals using Gauss-Legendre quadrature with 12 points    
   y1,y2,... - abisces
   w1,w2,... - weight
   accuracy  - accuracy
   r1 and r2 - from (k1,k2) to (-1,1)         */

{
   double y[33],w[33],r1,r2;
   int i1;
   r1 = ((k2 - k1) / 2);
   r2 = ((k2 + k1) / 2);

   if(accuracy == 0 || accuracy == 1){
      w[1 ]= 0.467913934572691 * r1; y[1 ]= 0.238619186083197 * r1;
      w[2 ]= 0.360761573048139 * r1; y[2 ]= 0.661209386466265 * r1;
      w[3 ]= 0.171324429379170 * r1; y[3 ]= 0.932469514203152 * r1;
      for(i1=0;i1<ip;i1++)vv[i1]=0.;
      for(i1=1;i1<4;i1++){ 
         fun48(w[i1],r2-y[i1],vv,ip);
         fun48(w[i1],r2+y[i1],vv,ip);
      }
      return 0.;
   }
   w[1 ]=  0.00178328072169643 * r1; y[1 ]=  0.99930504173577217 * r1;
   w[2 ]=  0.00414703326056247 * r1; y[2 ]=  0.99634011677195533 * r1;
   w[3 ]=  0.00650445796897836 * r1; y[3 ]=  0.99101337147674429 * r1;
   w[4 ]=  0.00884675982636395 * r1; y[4 ]=  0.98333625388462598 * r1;
   w[5 ]=  0.01116813946013113 * r1; y[5 ]=  0.97332682778991098 * r1;
   w[6 ]=  0.01346304789671864 * r1; y[6 ]=  0.96100879965205377 * r1;
   w[7 ]=  0.01572603047602472 * r1; y[7 ]=  0.94641137485840277 * r1;
   w[8 ]=  0.01795171577569734 * r1; y[8 ]=  0.92956917213193957 * r1;
   w[9 ]=  0.02013482315353021 * r1; y[9 ]=  0.91052213707850282 * r1;
   w[10]=  0.02227017380838325 * r1; y[10]=  0.88931544599511414 * r1;
   w[11]=  0.02435270256871087 * r1; y[11]=  0.86599939815409277 * r1;
   w[12]=  0.02637746971505466 * r1; y[12]=  0.84062929625258032 * r1;
   w[13]=  0.02833967261425948 * r1; y[13]=  0.81326531512279754 * r1;
   w[14]=  0.03023465707240248 * r1; y[14]=  0.78397235894334139 * r1;
   w[15]=  0.03205792835485155 * r1; y[15]=  0.75281990726053194 * r1;
   w[16]=  0.03380516183714161 * r1; y[16]=  0.71988185017161088 * r1;
   w[17]=  0.03547221325688239 * r1; y[17]=  0.68523631305423327 * r1;
   w[18]=  0.03705512854024005 * r1; y[18]=  0.64896547125465731 * r1;
   w[19]=  0.03855015317861563 * r1; y[19]=  0.61115535517239328 * r1;
   w[20]=  0.03995374113272034 * r1; y[20]=  0.57189564620263400 * r1;
   w[21]=  0.04126256324262353 * r1; y[21]=  0.53127946401989457 * r1;
   w[22]=  0.04247351512365359 * r1; y[22]=  0.48940314570705296 * r1;
   w[23]=  0.04358372452932345 * r1; y[23]=  0.44636601725346409 * r1;
   w[24]=  0.04459055816375657 * r1; y[24]=  0.40227015796399163 * r1;
   w[25]=  0.04549162792741814 * r1; y[25]=  0.35722015833766813 * r1;
   w[26]=  0.04628479658131442 * r1; y[26]=  0.31132287199021097 * r1;
   w[27]=  0.04696818281621002 * r1; y[27]=  0.26468716220876742 * r1;
   w[28]=  0.04754016571483031 * r1; y[28]=  0.21742364374000708 * r1;
   w[29]=  0.04799938859645831 * r1; y[29]=  0.16964442042399283 * r1;
   w[30]=  0.04834476223480295 * r1; y[30]=  0.12146281929612056 * r1;
   w[31]=  0.04857546744150343 * r1; y[31]=  0.07299312178779904 * r1;
   w[32]=  0.04869095700913972 * r1; y[32]=  0.02435029266342443 * r1;
   for(i1=0;i1<ip;i1++)vv[i1]=0.;
   for(i1=1;i1<33;i1++){
      fun48(w[i1],r2-y[i1],vv,ip);
      fun48(w[i1],r2+y[i1],vv,ip);
   }
   return 0.;
}
/******************************************************************************/

/* These are 5 functions that are used in integration routines.     

   dfunpar[] - parameters used by the integrayion routines.                   */

static double fun2(double x)
/* This function is:
                            2      3
               p1 + x p2 + x p3 + x p4
              e
                                                                              */
{
   return exp(dfunpar[1] + x * (dfunpar[2] + x * (dfunpar[3] + x *dfunpar[4])));
}

/******************************************************************************/
static double fun48(double w,double x,double vv[],int ip)
/* This function is:
           p0
          x  f(x)
                 +
                      where f(x) is a density provided by dens3               */
{
   double vx;
   int i1;
   vx=w*dens3(x);
   vv[0]=vv[0]+vx;
   for(i1=1;i1<ip;i1++){
      vx=vx*x;
      vv[i1]=vv[i1]+vx;
   }
   return 0.;
}
/******************************************************************************/
/* this file contains functions that are executed - before the iterations start-
   by the main program.

   coeff  - computes the coefficient matrices for the basis functions          
   start1 - gathers statistics used to compute the starting values
   start2 - computes the starting values
   sufficient1 - gathers statistics used to compute the sufficientt statistics
   sufficient2 - computes the sufficientt statistics
   knotplace   - places the knots
   knotnumber  - determines the number of starting knots
   piecedens   - determines a piecewise constant density for initial knot placem
*******************************************************************************/

static void coeff(double coef2[][NC])
/* This function computes the coefficients of the basis functions from the knots
   coef2:
   first index: basis function number-1,
   second index: 0:1, 1:x, 2:(x-t1)+^3,  3:(x-t2)+^3, 4:(x-t3)+^3,........
   coef: 
   between knot(i) and knot(i+1) the coefficient of x^power of basisfunction(j)
   first index: basis function number-1 (j-1)
   second index: power of x
   third index: interval (i)                                                  

   Basis functions B(2)-B(nknots-3) are multiples of B-splines.
   Further the coefficients are choosen such that the quadratic and  
   cubic terms in both tails are 0;  this leads to differnt basis
   functions for B(1), B(nknots-2) and B(nknots-1). B(1) is linear left
   of the first knot. B(nknots-2) is constant to the right of the last
   knot and B(nknots-1) is linear to the right of the last knot.              */


{
   int i,j,k;
   double z0,z1;
/* local
   i j k       - counter
   z0,z1       - value of constants of two succesive basisfunctions           */

/* Initializations                                                            */
   for(i=0; i<nknots; i++){
      for(j=0; j<nknots+3; j++){
         coef2[i][j]=0.;
         if(j<nknots+2)
            for(k=0; k<4; k++)
               coef[i][k][j]=0.;
      }
   }

/* The coefficients for the two tail basis functions are easy to compute      */
   coef2[0][2] = 1.;
   coef2[0][3] = (knots[0] - knots[2]) / (knots[2] - knots[1]);
   coef2[0][4] = (knots[1] - knots[0]) / (knots[2] - knots[1]);
   coef2[0][1] = -3. * (pow(knots[0], 2.) 
                     + coef2[0][3] * pow(knots[1], 2.)
                     + coef2[0][4] * pow(knots[2], 2.));
   coef2[0][0] = - knots[nknots-1] * coef2[0][1]
                 - coef2[0][2] * pow((knots[nknots-1] - knots[0]), 3.)
                 - coef2[0][3] * pow((knots[nknots-1] - knots[1]), 3.)
                 - coef2[0][4] * pow((knots[nknots-1] - knots[2]), 3.);
   coef2[0][5] = 0.;
   coef2[nknots-2][nknots-1] = 1.;
   coef2[nknots-2][nknots]   = (knots[nknots-3] - knots[nknots-1]) / 
                               (knots[nknots-1] - knots[nknots-2]);
   coef2[nknots-2][nknots+1] = (knots[nknots-2] - knots[nknots-3]) /
                               (knots[nknots-1] - knots[nknots-2]);    

/* we first create basis functions that are 0 before knot[i] and constant
   after knot [i+3]                                                           */

   if(nknots >= 4){
      for(i=1; i<nknots-2; i++){
         coef2[i][i+1] = 1.;
         coef2[i][i+4] = (knots[i+1]-knots[i-1]) * (knots[i-1]-knots[i])
                     / ((knots[i+1]-knots[i+2]) * (knots[i]-knots[i+2]));
         coef2[i][i+3] = (coef2[i][i+4] * (knots[i]-knots[i+2])
                     + knots[i] - knots[i-1]) / (knots[i+1]-knots[i]);
         coef2[i][i+2] = -1. - coef2[i][i+3] - coef2[i][i+4];
      }

/* In the following part we subtract a number of times one basis
   function from another - so that basis function i becomes 0 after knot[i+4] */
      if(nknots >= 5){
         for(i=1; i<nknots-3; i++){
            z0 = 0.;
            z1 = 0.;
            for(j=2; j<nknots+1; j++){
               z0 = z0 + coef2[i][j]   * pow((knots[nknots-1]- knots[j-2]), 3.);
               z1 = z1 + coef2[i+1][j] * pow((knots[nknots-1]- knots[j-2]), 3.);
            }
            for(j=2; j<nknots+2; j++)
                  coef2[i][j] = coef2[i][j] - (z0 / z1) * coef2[i+1][j];
         }
      }
   }

/* Now the coef matrix. First basis function 1.                               */

   for(k=0; k<3; k++){
      coef[0][1][k] = coef2[0][1];
      coef[0][0][k] = coef2[0][0];
   }

/* The rest is a bit tricking with the correct indices                        */
   for(i=0; i<nknots-1; i++){
      for(j=i; j<i+4; j++){ 
         for(k=i+1; k<j+2; k++){ 
            if(j > 0 && j < nknots+1 && (i != 0 || j != 3)){
               if(k != 1){
                  coef[i][0][j] = coef[i][0][j] - 
                                  coef2[i][k] * pow(knots[k-2], 3.);
                  coef[i][1][j] = coef[i][1][j] + 
                                  3. * coef2[i][k] * pow(knots[k-2], 2.);
                  coef[i][2][j] = coef[i][2][j] - 3. * coef2[i][k] * knots[k-2];
                  coef[i][3][j] = coef[i][3][j] + coef2[i][k];
               }
            }
         }  
      }  
   }    
}
/******************************************************************************/

/* These are two functions to compute the starting values.
   The first one gathers statistics from all the datapoints and knots.
   The second one computes the starting values from these statistics          */

static void start1(double crossprods[][NC],double derivatives[],double sample[],int nsample[])

/* the objective of this routine is the computation of derivatives and
   crossproducts                                                              */


/* sample,nsample,nknots, accuracy and knots see lhead.h and originating file
   derivatives[0]   = sum((sample))'' which is 0....
   derivatives[i]   = sum(([(sample-knot[i-1])+]^3)'')
   crossprods[0][0] = sum((sample)'*(sample)')=nsample[0]
   crossprods[i][0] = sum(([(sample-knot[i-1])+]^3)'*(sample)')
   crossprods[0][i] = sum(([(sample-knot[i-1])+]^3)'*(sample)')
   crossprods[i][j] = sum(([(sample-knot[i-1])+]^3)'*([(sample-knot[j-1])+]^3)')
   they are used in start2 to compute starting values                         */

{
   int i,j,k,i2,i3;
   double xs,uuu,rr,rs;
/* counters, utility                                                          */

/* Initializations                                                            */
   for(i=0; i<nknots+1; i++){
      derivatives[i] = 0.;
      for(j=0; j<nknots+1; j++)
         crossprods[i][j] = 0.;
   }
   crossprods[0][0] = nsample[0];

   rr=sample[0];
   rs=sample[0];
/* First the uncensored part                                                  */
   for(i=0;i<nsample[0]+nsample[2];i++) if(sample[i]>rr) rr=sample[i];
   for(i=0;i<nsample[0]+nsample[2];i++) if(sample[i]<rs) rs=sample[i];
   if(nsample[1]>0)for(j=0; j<nknots; j++){
      for(i=0;i<nsample[1];i++){
         if(sample[i]>knots[j]){
            derivatives[j+1] = derivatives[j+1] + 6 * (sample[i] - knots[j]);
            crossprods[j+1][0]=crossprods[j+1][0]+3*pow(sample[i]-knots[j],2.);
            for(k=0; k<=j; k++)
               crossprods[j+1][k+1] = crossprods[j+1][k+1] + 
                  9 * pow(((sample[i] - knots[j]) * (sample[i] - knots[k])),2.);
         }
      }
   }
/* The interval censored part - we take the midpoints of the intervals        */
   if(nsample[2]>0)for(j=0; j<nknots; j++){
      for(i=0;i<nsample[2];i++){
         xs=(sample[nsample[1]+2*i]+sample[nsample[1]+2*i+1])/2.;
         if(xs>knots[j]){
            derivatives[j+1] = derivatives[j+1] + 6 * (xs - knots[j]);
            crossprods[j+1][0]=crossprods[j+1][0]+3*pow(xs-knots[j],2.);
            for(k=0; k<=j; k++)
               crossprods[j+1][k+1] = crossprods[j+1][k+1] + 
                  9 * pow(((xs - knots[j]) * (xs - knots[k])),2.);
         }
      }
   }
/* The right censored part                                                    */
   if(nsample[3]>0)for(j=0; j<nknots; j++){
      for(i2=0;i2<nsample[3];i2++){
         i=nsample[1]+2*nsample[2]+i2;
         for(i3=0;i3<4;i3++)uuu=(double)(i3)*rr/4.+(4-i3)*sample[i]/4.;
         if(uuu>knots[j]){
            derivatives[j+1] = derivatives[j+1] + 1.5 * (uuu - knots[j]);
            crossprods[j+1][0]=crossprods[j+1][0]+0.75*pow(uuu-knots[j],2.);
            for(k=0; k<=j; k++)
               crossprods[j+1][k+1] = crossprods[j+1][k+1] + 
                  2.25 * pow(((uuu - knots[j]) * (uuu - knots[k])),2.);
         }
      }
   }

/* The left  censored part                                                    */
   if(nsample[4]>0)for(j=0; j<nknots; j++){
      for(i2=0;i2<nsample[4];i2++){
         i=nsample[1]+2*nsample[2]+nsample[3]+i2;
         for(i3=0;i3<4;i3++)uuu=(double)(i3)*rs/4.+(4-i3)*sample[i]/4.;
         if(uuu>knots[j]){
            derivatives[j+1] = derivatives[j+1] + 1.5 * (uuu - knots[j]);
            crossprods[j+1][0]=crossprods[j+1][0]+0.75*pow(uuu-knots[j],2.);
            for(k=0; k<=j; k++)
               crossprods[j+1][k+1] = crossprods[j+1][k+1] + 
                  2.25 * pow(((uuu - knots[j]) * (uuu - knots[k])),2.);
         }
      }
   }

/* symmetarize crossprods                                                     */
   for(j=0; j<nknots; j++)
      for(k=j+1; k<nknots+1; k++)
         crossprods[j][k] = crossprods[k][j];
}

/******************************************************************************/
static void start2(double crossprods[][NC],double derivatives[],
    double coef2[][NC],int nkstart,int iremove)

/* this function combines derivatives and crossprods, and then inverts a
   system to compute starting values                                          */


{
   int i,j,k,l,j2,l2,infox,kpvt[200];
   double work[NC][NC],r1,r2,r3,r4;

/* local:
   i j j2 k l l2  - counters
   work           - workmatrix, 
   infox kpvt     - reguired by inversion routine                             */

/* Initializations */
   r1=zheta[0];
   r2=zheta[1];
   r3=zheta[nknots-2];
   r4=zheta[nknots-1];
   if(nknots>4){
      for(i=0; i<nknots; i++){
         zheta[i] = 0.;
         for(j=0; j<nknots; j++)
            work[i][j] = 0.;
      }

/* The left hand side is a matrix with sums of crossproducts of derivatives of
   basis functions.  We make them by combining coefficients (coef2) and
   crossproducts (crossprods) of the powerbasis.  The right hand side is a
   vector with second derivatives of the basis functions.  Again combining
   coefficients (coef2) now with derivatives.                                 */

      for(i=0; i<nknots-1; i++){
/* these are the only combinations for which coef2 is not equal to 0          */
         for(j2=i+1; (j2<i+6 && j2<nknots+2); j2++){
            j = j2;
            if(i == 0 && j == 5)j = 1;
            zheta[i] = zheta[i] - coef2[i][j] * derivatives[j-1];
            for(k=i; k<nknots-1; k++){

/* these are the only combinations for which coef2 is not equal to 0          */
               for(l2=k+1; (l2<k+6 && l2<nknots+2); l2++){
                  l = l2;
                  if(k == 0 && l == 5)l = 1;
                  work[i][k] = work[i][k] +
                     coef2[i][j] * coef2[k][l] * crossprods[j-1][l-1];
               }
            }
         }
      }
      
/* Make symmetric                                                             */
      for(i=1; i<nknots-1; i++)
         for(j=0; j<i; j++)
            work[i][j]=work[j][i];

/* Solving the system finishes the job.                                       */
      i=nknots-1;
      j=NC;
      F77_CALL(xdsifa)(work,&j,&i,kpvt,&infox);
      F77_CALL(xdsisl)(work,&j,&i,kpvt,zheta);
/* If we had previous iterations - we're sometimes better of using that       */
      if(nkstart!=nknots){
         if(iremove>2){
            zheta[0]=r1;
            zheta[1]=r2;
         }
         if(iremove<3){
            zheta[0]=r1/25.;
            zheta[1]=r2/25.;
         }
         if(iremove<nknots-2){
            zheta[nknots-3]=r3;
            zheta[nknots-2]=r4;
         }
         if(iremove>nknots-3){
            zheta[nknots-3]=r3/25.;
            if(nknots==4&&iremove==3)zheta[1]=(r2+r3)/50.;
            zheta[nknots-2]=r4/25.;

         }
      }
/* If the zhetas are too large, we are better of if we make them smaller      */
      else{
         r1=0.;
         for(i=0;i<nknots-2;i++)r1=r1+fabs(zheta[i]);
         if(r1>10000)
            for(i=0;i<nknots-2;i++)zheta[i]=zheta[i]*10000/r1;
      }
   }
   else{
      for(i=0;i<3;i++) zheta[i]=0.;
   }
}
/******************************************************************************/

/* The first routine collects statistics about the data. The second one
   uses these to compute the sufficiemt statistics                            */

static void suffstat1(double suffcombine[][2],double sample[],int nsample[])
/* these quantities are defined in lhead.h and the file where they originate  
   suffcombine[0][a] = 1
   suffcombine[1][a] = sum(sample)/nsample[0]
   suffcombine[i][a] = sum(((sample-knot[i-1])+)^3)/nsample[0]                    
   where a is 0 for the exact sample and 1 for the interval censored stuff.
   they are later (in suffstat2) combined to compute the sufficient statistics*/

{
   int i,j;
   double rr;
/* counters, utility                                                          */

/* initialization */
   for(i=1;i<nknots+2;i++){
      suffcombine[i][0] = 0.;
      suffcombine[i][1] = 0.;
   }

   suffcombine[0][0] = (double)nsample[1]/(double)nsample[0];
   suffcombine[0][1] = (double)nsample[2]/(double)nsample[0];

/* for the complete observations                                              */
   for(j=0; j<nsample[1]; j++){
      suffcombine[1][0] = suffcombine[1][0] + sample[j]/(double)nsample[0];
      for(i=2; i<nknots+2; i++)
         if(sample[j] > knots[i-2])
           suffcombine[i][0] = suffcombine[i][0] + 
                           pow(sample[j]  - knots[i-2], 3.)/(double)nsample[0] ;
   }
/* for the interval censored observations                                     */
   for(j=0; j<nsample[2]; j++){
      rr=(sample[nsample[1]+2*j]+sample[nsample[1]+2*j+1])/2.;
      suffcombine[1][1] = suffcombine[1][1] + rr/(double)nsample[0];
      for(i=2; i<nknots+2; i++)
         if(rr > knots[i-2])
           suffcombine[i][1] = suffcombine[i][1] + 
                          pow(rr  - knots[i-2], 3.)/(double)nsample[0] ;
   }
}

/******************************************************************************/

static void suffstat2(double suffcombine[][2],double coef2[][NC],double sufficient[][2])
/* all defined in lhead.h and the file where they originate.
   suffcombine defined in suffstat1                                           */

{
   int i,j;
/* counters                                                                   */

/* over the basisfunctions                                                    */
   for(i=0;i<nknots;i++){
      sufficient[i][0] = 0.;
      sufficient[i][1] = 0.;
         
/* over the intervals                                                         */
      for(j=0;j<nknots+2;j++){
         sufficient[i][0] = sufficient[i][0] + suffcombine[j][0] * coef2[i][j];
         sufficient[i][1] = sufficient[i][1] + suffcombine[j][1] * coef2[i][j];
      }
   }
}
/******************************************************************************/
/* These routines determine the position of the knots and the number of knots.*/

static void knotplace(int iknots[], double rknots[], int iknotauto, double bound[],
	double sample[], int nsample[], int SorC[], double smp2[], double smp3[], double qt[])
                                    
/* these quantities are defined in lhead.h and the files where they originate */

{
   int i,j=0,j2,k,kk,ll,ia=0,il,kx,w1,w2,w3,w4,w5,w6;
/* local integers
   i k     - counters
   j j2    - is there an odd or an even number of knots?                      */

   double s,eps,eps1,eps2,v,w,five=4,u1,u2,u3=0.,u4=0.;
/* local   - double precision
   five    - parameter for knot rule
   s       - this parameter has to become (nsample[0]+1)/2
   eps     - current estimate for five rule
   eps1    - lower bound for eps
   eps2    - upper bound for eps
   v w     - combinations of counters                                         */

/* the major part of this function is computing the knots for the
   five rule. The essence is to find eps such that s=(nsample[0]-1)/2
   The usual make it a bit larger, smaller, till we have one just
   too large and one just too small.                                          */

   w1=nsample[0];
   w2=nsample[1];
   w3=nsample[2];
   w4=nsample[3];
   w5=nsample[4];
   w6=nsample[5];
   if(SorC[16]==1){
      nsample[0]=nsample[1]=w2;
      nsample[2]=nsample[3]=nsample[4]=0;
      nsample[5]=w2-1;
   }

   qt[0]=-1.;
   qt[1]=1.;
   if(iknotauto == 0&&nknots>12){
      qt[0]=knots[4];
      qt[1]=knots[nknots-5];
   }
   if(iknotauto == 0&&nknots>9 && nknots <13){
      qt[0]=knots[3];
      qt[1]=knots[nknots-4];
   }
   if(iknotauto == 0&&nknots>6 && nknots <10){
      qt[0]=knots[2];
      qt[1]=knots[nknots-3];
   }
   if(iknotauto == 0&&nknots==3){
      qt[0]=(knots[0]+knots[1])/2.;
      qt[1]=(knots[2]+knots[1])/2.;
   }
   if(iknotauto == 0&&nknots>3 &&nknots<7){
      qt[0]=knots[1];
      qt[1]=knots[nknots-2];
   }
   if(iknotauto == 1){

/* compute the piecewise density                                              */
      il = piecedens(sample,smp2,smp3,nsample);
      kk = floor((double)il/2+2.1);
      if(nknots>kk && nsample[5]==0){
          nknots=kk;
          if(SorC[0] == 0)
               (void)
                  Rprintf("running with maximum number of degrees of freedom");
          else
             SorC[20] = 1;
      }
      if(nknots==nsample[5]+1)SorC[20]=1;

/* Check whether there are not too many knots                                 */
      kk=0;
      ll=0;
      do{
         ll=kk+ll;
         kk=0;
         if(bound[1] < 0.5 && bound[3] < 0.5){

/* all knots are a minimum of "five" apart                                    */

            if((nknots-1) * five >= nsample[0] - 1){
               i = floor(1. + (nsample[0] - 1.) / five);
               if(SorC[0] == 0){
                  (void)Rprintf("too many knots, at most %d knots possible\n",i);
                  SorC[0] = -647;
               }
               else {
                  SorC[0] = -2;
                  SorC[1] = i;
               }
               return;
            }

/*  place the 2 extreme knots                                                 */

            rknots[0] = 1.;
            rknots[nknots-1] = nsample[0];

/*  j and j2 are this way to deal both with odd and even situations           */

            j = ceil((nknots-1)/2.);
            j2 = floor((nknots-1)/2.);

/*  eps1 and eps2 are lower and upper bound on eps, eps is our first guess.   */

            eps1 = five - pow(((nsample[0] - 1) / five),(1. / (j - 1)));
            if(eps1 > 0.)
               eps1 = 0.;
            eps2 = five - 1.;

/* s should become exactly (nsample[0]-1)/2: in that case the knots
   are symmetric and cover exactly the whole range.
   We here compute what s would be for the present value of eps 
   s is the span of all the knots: i.e. the location of the middle knot
   if the first knot is at 1, nknots as specified and eps as guessed.         */

            do{
               eps = (eps1 + eps2) / 2.;
               s = 1.;
               w = five;
               for(i=1; i<=j2; i++){
                  v = i;
/* s is the location after adding another knot                                */
                  s = s + w;   
/* we store the rknots - in case they are good                                */
                  rknots[i] = s;
                  rknots[nknots-i-1] = nsample[0] + 1 - s;
                  v = five - v * eps;
/* w is what is going to be the next gap                                      */
                  if(v < 1.) v = 1.;
                  w = w * v;
               }
/* Are there an odd or even number of gaps?                                   */
               if(j * 2 == nknots) 
/* even - no more knots to place.                                             */
                  s = s + w/2.;     
               else
/* odd - the last knot                                                        */
                  rknots[j] = (nsample[0] + 1.) / 2.;
   
/* Is eps too large or too small?                                             */
               if(2. * s >= nsample[0]+1) 
                  eps1 = eps;
               else
                  eps2 = eps;
   
/* Are eps1 and eps2 close together                                           */
            } while(eps2 - eps1 > .001);   
         }
         else{
/* if both sides there are finite limits we put them equidistant              */
            if(bound[1] > 0.5 && bound[3] > 0.5){
               for(i=0;i<nknots;i++)rknots[i]=
                   1.+(double)(nsample[0]-1)*((double)i)/(double)(nknots-1);
            }
            else{
/* one side bounded, one side unbounded half of the five-rule              
   which side is bounded?                                                     */
               rknots[0]=1.;
               rknots[nknots-1]=nsample[0];
               eps1 = 0.;
               eps2 = five;
               do {
                  eps = (eps1 + eps2)/2.;
                  s = 1.;
                  w = five;
                  for(i=1;i<nknots-1;i++){
                     v = i;
/* s is the location after adding another knot                                */
                     s = s+w;
/* we store the rknots - in case they are good                                */
                     if(bound[3] > 0.5)
                        rknots[i] = s;
                     else
                        rknots[nknots-i-1]=nsample[0]+1.-s;
                     v = five - v *eps;
/* w is what is going to be the next gap                                      */
                     if(v < 1.) v = 1.;
                     w = w * v;
                 }
/* Is eps too large or too small?                                             */
                 if(s + w >(double)nsample[0])
                    eps1 = eps;
                 else
                    eps2 = eps;
/* Are eps1 and eps2 close together                                           */
              } while(eps2 - eps1 > .001);
            }
         }
         if(nknots==3)rknots[1]=nsample[0]/2.;

/* Translate rknots in knots                                                  */
         if(ll==0){
/* the first two knots are easy                                               */
            knots[0] = smp2[0];     
            knots[nknots-1] = smp2[il];

 /* cycle through the endpoints, average                                      */
            k = 0;
            for(i=1;i<nknots-1;i++){
               rknots[i]=(rknots[i]-0.5);
               kx=k;
               for(k=kx;k<=il;k++){
                  if(smp3[k]>=rknots[i]){
                     knots[i]=((rknots[i]-smp3[k-1])*smp2[k]+
                             (smp3[k]-rknots[i])*smp2[k-1])/(smp3[k]-smp3[k-1]);
                     ia=k;
                     k=k+4*il;
                  }
               }
               k=ia;
            }
         }

/* find the first and the third quartile                                      */
         for(i=1;i<il;i++){
            if(smp3[i]>0.25*nsample[0]){
               u3=smp2[i];  
               i=i+2*il;
               qt[0]=u3;
            }
         }
         for(i=1;i<il+1;i++){
            if(smp3[i]>0.75*nsample[0]){
               u4=smp2[i];  
               i=i+2*il;
               qt[1]=u4;
            }
         }

/* knots that are close together at the end, are an indicator of a discontinuity
   we replace the knots if this happens. kk is an indicator that somthing like
   this is the case - we fix this by, temporarily, cheating and saying that
   there is a discontinuity (changing bound), ll reminds us to change back, by
   the way, u2 is the IQR of the data                                         */

         if(bound[1]<0.5){
            u1=knots[2]-knots[0];
            u2=u4-u3;
            if(u1<u2/30.){
               if(SorC[0]==0){
                  (void)
                  Rprintf("possible discontinuity of the density at left end\n");
                  (void)
                  Rprintf("considder rerunning it with lbound=min(sample)\n");
               }
               else SorC[21]=1;
               kk=kk+1;
               bound[1]=1;
               bound[2]=knots[0];
            }
         }
         if(bound[3]<0.5){
            u1=knots[nknots-1]-knots[nknots-3];
            u2=u4-u3;
            if(u1<u2/30.){
               if(SorC[0]==0){
                  (void)
                 Rprintf("possible discontinuity of the density at right end\n");
                  (void)
                  Rprintf("considder rerunning it with ubound=max(sample)\n");
               }
               else SorC[22]=1;
               kk=kk+2;
               bound[3]=1;
               bound[4]=knots[nknots-1];
            }
         }
      }while(kk>0);
      if(ll==1||ll==3)bound[1]=0.;
      if(ll==2||ll==3)bound[3]=0.;
   
/* if after previous changes these knots are still to close together, the
   density might be real high at the end - we toss out knots in this case     */
   u1=knots[nknots-1]-knots[nknots-3];
   u2=u4-u3;
   if(u1<u2/30){
      if(SorC[0]==0){
      (void)Rprintf("potential very high value of the density at right end\n");
      (void)Rprintf("running program with fewer knots than usual\n");
      }
      else SorC[22]=3;
      knots[nknots-3]=(knots[nknots-3]+knots[nknots-2])/2.;
      knots[nknots-2]=knots[nknots-1];
      nknots=nknots-1;
   }
   do{
      u1=knots[nknots-1]-knots[nknots-3];
      u2=u4-u3;
      if(u1<u2/30){
         knots[nknots-3]=(knots[nknots-3]+knots[nknots-2])/2.;
         nknots=nknots-1;
      }
   }while(u1<u2/30);
   u1=knots[2]-knots[0];
   u2=u4-u3;
   if(u1<u2/30){
      if(SorC[0]==0){
         (void)Rprintf("potential very high value of the density at left end\n");
         (void)Rprintf("running program with fewer knots than usual\n");
      }
      else SorC[21]=3;
      knots[1]=(knots[1]+knots[2])/2.;
      for(i=2;i<nknots-1;i++)knots[i]=knots[i+1];
      nknots=nknots-1;
   }
   do{
      u1=knots[2]-knots[0];
      u2=u4-u3;
      if(u1<u2/30){
         knots[1]=(knots[1]+knots[2])/2.;
         for(i=1;i<nknots-1;i++)knots[i]=knots[i+1];
         nknots=nknots-1;
      }
   }while(u1<u2/30);
   }

/* Remove double knots,either store zerror messages, or print them.            */

   j = 0;                          
   k = 0;
   for(i=1; i<nknots; i++){
      if(knots[i] > knots[j]){
         j++;
         knots[j]=knots[i];
      }
      else {
         k++;
         if(SorC[0] == 8){
            (void)Rprintf("===> warning: knot %d removed - double knot\n",i+1);
            if(k == 2){
               (void)
               Rprintf("* several double knots suggests that your data is *\n");
               (void)
               Rprintf("* strongly rounded: attention might be required.  *\n");
            }
         }
         else {
            SorC[0] = 2 + k;
            if(SorC[0] == 23)return;
            SorC[k] = i+1;
         }
      }
   }
   nknots = j+1;
   nsample[0]=w1;
   nsample[1]=w2;
   nsample[2]=w3;
   nsample[3]=w4;
   nsample[4]=w5;
   nsample[5]=w6;
   u3=2./(qt[1]-qt[0]);
   u4=1.-2.*qt[1]/(qt[1]-qt[0]);
   qt[0]=u4;
   qt[1]=u3;
   for(i=0;i<nknots;i++){
      rknots[i] = u4 + u3 * knots[i];
      knots[i] = u4 + u3 * knots[i];
      iknots[i] = i;
   }           
   for(i=0;i<nsample[0]+nsample[2];i++)sample[i]=u4+u3*sample[i];
}

/******************************************************************************/
/* determines the number of (starting) knots                                  */

static int knotnumber(int idelete,int nsample[],int nknots,int SorC[])

/* all defined in lhead.c and originating file                                */
{

   int nn;
   double r;

   nn=nsample[1]+nsample[2]+(nsample[3]+nsample[4])/2;
   if(nknots>=NC || (nknots > nsample[5]+1 && nsample[5]>1)){ 
      if(SorC[0]==0)
         (void)Rprintf("can not run with that many knots\n");
      else
         SorC[20]=1;
      nknots=NC-1;
      if(nknots > nsample[5]+1)nknots=nsample[5]+1;
   }
   if(nknots > 0)
    return nknots;
   r = 2.5*pow((double)nn,0.2);
   if(idelete > 0)r=1.6*r;
   if(idelete > 0 && nn < 51)r=(double)nn/6.;
   if(r>NC-1.1)r=NC-1.1;
   if(r<2.5)r=2.5;
   if(nsample[5]!=0 && r>nsample[5]+1)r=nsample[5]+0.1;
   return ceil(r);
}
/******************************************************************************/
static int piecedens(double sample[],double smp2[],double smp3[],int nsample[])
{
   int i,j,il,k,m,n;
   if(nsample[1]>0){
      i=1;
      j=nsample[1];
      F77_CALL(xssort)(sample,smp3,&j,&i);
   }
   if(nsample[2]>0){
      for(i=0;i<nsample[2];i++){
         smp2[i]=sample[nsample[1]+2*i];
         smp3[i]=sample[nsample[1]+2*i+1];
      }
      i=2;
      j=nsample[2];
      F77_CALL(xssort)(smp2,smp3,&j,&i);
      for(i=0;i<nsample[2];i++){
         sample[nsample[1]+2*i]=smp2[i];
         sample[nsample[1]+2*i+1]=smp3[i];
      }
   }
   if(nsample[3]>0){
      for(i=0;i<nsample[3];i++) smp2[i]=sample[i+nsample[1]+2*nsample[2]];
      i=1;
      j=nsample[3];
      F77_CALL(xssort)(smp2,smp3,&j,&i);
      for(i=0;i<nsample[3];i++) sample[i+nsample[1]+2*nsample[2]]=smp2[i];
   }
   if(nsample[4]>0){
      for(i=0;i<nsample[4];i++)
         smp2[i]=sample[i+nsample[1]+2*nsample[2]+nsample[3]];
      i=1;
      j=nsample[4];
      F77_CALL(xssort)(smp2,smp3,&j,&i);
      for(i=0;i<nsample[4];i++)
         sample[i+nsample[1]+2*nsample[2]+nsample[3]]=smp2[i];
   } 
   for(i=0;i<nsample[1];i++){
      smp2[i]=sample[i];
      smp3[i]=1.;
   }
   for(i=0;i<nsample[3]+nsample[4];i++){
      j=nsample[1]+2*nsample[2]+i;
      smp2[j]=sample[j];
      smp3[j]=1.;
   }
   if(nsample[2]>0){
      for(i=0;i<nsample[2];i++){
         j=nsample[1]+i*2;
         smp2[j]=sample[j];
         smp3[j]=0.;
         smp2[j+1]=sample[j+1];
         smp3[j+1]=0.;
      }
   }
   i=2;
   j=nsample[0]+nsample[2];
   F77_CALL(xssort)(smp2,smp3,&j,&i);
   i=0;
   for(j=1;j<nsample[0]+nsample[2];j++){
      if(smp2[j]==smp2[i])
         smp3[i]=smp3[i]+smp3[j];
      else{
         i++;
         smp2[i]=smp2[j];
         smp3[i]=smp3[j];
      }
   }
   il=i;
      
   for(i=il-1;i>=0;i=i-1){
      if(smp3[i]>0){
         smp3[i+1]=smp3[i+1]+0.5*smp3[i];
         smp3[i]=0.5*smp3[i];
      }
   }
   if(nsample[2]>0){
      k=0;
      for(i=0;i<nsample[2];i++){
         j=nsample[1]+i*2;
         for(m=k;m<=il;m++){
            if(sample[j]<smp2[m]){
               k=m;
               m=2*il+1;
               for(n=k;n<=il;n++){
                  if(sample[j+1]<smp2[n])n=2*il+1;
                  else
                  smp3[n]=smp3[n]+(smp2[n]-smp2[n-1])/(sample[j+1]-sample[j]);
               }
            }
         }
      }
   }
   for(i=1;i<=il;i++) smp3[i]=smp3[i]+smp3[i-1];
   return il;
}
/******************************************************************************/
/* These 2 routines comoute the informationmatrix (info) the score
   vector (score) and C(zheta) (czheta).  The routine consists
   of 2 parts. The middle1 routine, that does the numerical integrations
   in the middle, and tails that does the exact computations for
   the tails.                                                                 */

/* computes the numerical integrals for shift(score) info and czheta for the
   middle.  calls tails for the tails. */


static double middle(double info[][NC],double shift[],double sufficient[][2],
   double bound[],int accuracy,int nsample[],double sample[],double zheta[],int what)

{
   double qolint[7][NC+1],aa[7],d1,e1,cth,d2,e2;
   int i1,i2,i3,i4,i5,version,version2,k0;
/* local: i1,i2,i3,i4,i5 - counter
          polint[]       - integrals of the form x**i.exp(polynomial)         */

/* The tails can be computed analytically                                     */
   cth = tails(info,shift,coef,bound,knots,zheta,nknots,what);

/* what == 0: we are only interested in czheta                                */
   if(what==0){
      for(i1=1; i1<nknots; i1++){               
         for(i2=1; i2<5; i2++){
            dfunpar[i2]=0.;
            for(i3=i1-3; i3<i1+1; i3++){
               if(i3>=0 && i3 < nknots-1)
                  dfunpar[i2]=dfunpar[i2]+coef[i3][i2-1][i1]*zheta[i3];
            }
         }
         cth = cth + numint2(knots[i1-1],knots[i1],accuracy);
      }
      return cth;
   }

/* The integrals are computed numerically per interval between knots          */
   for(i1=1;i1<nknots;i1++){               

/* polint gets the result of a numerical integration of the density
   times x^i2                                                                 */
      (void)numints(aa,knots[i1-1],knots[i1],accuracy,7);
      for(i2=0;i2<7;i2++)qolint[i2][i1]=aa[i2];

/* For the scorefunction we need integrals of a basisfunction times
   the density, this is a combination of qolint[]. Basisfunctions
   do not exists in all intervals.                                            */
      for(i2=i1-3;i2<=i1;i2++){
         if(i2>=0 && i2 < nknots-1 && (i2!=0 || i1!=3)){
            shift[i2] = shift[i2] + qolint[0][i1] * coef[i2][0][i1]
                                  + qolint[1][i1] * coef[i2][1][i1]
                                  + qolint[2][i1] * coef[i2][2][i1]
                                  + qolint[3][i1] * coef[i2][3][i1];

/* For the information matrix we need integrals of a basisfunction
   times a basisfunction times the density, this is a combination of
   qolints. Many combinations do not exist.                                   */
            if(what==2){
               for(i3=i1-3;i3<=i2;i3++)
                  if(i3>=0 && i3 < nknots-1 && (i3!=0 || i1!=3))
                     for(i5=0;i5<4;i5++)
                        for(i4=0;i4<4;i4++)
                           info[i2][i3] = info[i2][i3] + qolint[i5+i4][i1]
                                     * coef[i2][i5][i1] * coef[i3][i4][i1];     
            }
         }
      }
   }
/* The following lines make the informationmatrix and score function
   from what is stored into score and info up to now.                         */
   if(what==2){
      for(i1=0;i1<nknots-1;i1++){
         for(i2=0;i2<=i1;i2++){
            info[i1][i2] = info[i1][i2] - shift[i1] * shift[i2];
            info[i2][i1] = info[i1][i2];
         }
      }
   }
/* if there is no censoring, we are done                                      */
   if(nsample[0]!=nsample[1]){

      version=3;
      if(bound[3] > 0.5) version=4;
      version2=1;
      if(bound[1] > 0.5) version2=2;
      d1 = zheta[nknots-2] * coef[nknots-2][1][nknots];
      d2 = zheta[0] * coef[0][1][0];
      e2 = zheta[0] * coef[0][0][0] -log(czheta);
      if(nknots>3)
         e1 = zheta[nknots-2] * coef[nknots-2][0][nknots]
            + zheta[nknots-3] * coef[nknots-3][0][nknots]-log(czheta);
      else e1 = zheta[nknots-2] * coef[nknots-2][0][nknots]-log(czheta);
      for(i2=0;i2<7;i2++)aa[i2]=0.;
      for(i2=0;i2<7;i2++){
         aa[i2]=1.;
         if(i2!=0)aa[i2-1]=0.;
         qolint[i2][nknots]=expin2(version,knots[nknots-1],bound[4],aa,d1,e1);
      }
      for(i2=0;i2<7;i2++)aa[i2]=0.;
      for(i2=0;i2<7;i2++){
         aa[i2]=1.;
         if(i2!=0)aa[i2-1]=0.;
         qolint[i2][0]=expin2(version2,knots[0],bound[2],aa,d2,e2);
      }
   
      if(accuracy!=0){
         if(nsample[3]>0){
            k0=1;
            for(i2=0;i2<nsample[3];i2++){
            if(i2==nsample[3]-1||sample[i2+nsample[1]+2*nsample[2]]!=
            sample[i2+nsample[1]+2*nsample[2]+1]){
            intnum3(sample[i2+nsample[1]+2*nsample[2]],
                   qolint,d1,e1,version,bound[4],shift,info,k0,nsample[0],what);
            k0=0;
            }
            k0++;
            }
         }
         if(nsample[4]>0){
            k0=1;
            for(i2=0;i2<nsample[4];i2++){
            if(i2==nsample[4]-1||sample[i2+nsample[1]+2*nsample[2]+nsample[3]]!=
            sample[i2+nsample[1]+2*nsample[2]+nsample[3]+1]){
            intnum4(sample[i2+nsample[1]+2*nsample[2]+nsample[3]],
                  qolint,d2,e2,version2,bound[2],shift,info,k0,nsample[0],what);
            k0=0;
            }
            k0++;
            }
         }
         if(nsample[2]>0){
            k0=1;
            for(i2=0;i2<nsample[2];i2++){
               i1=nsample[1]+2*i2;
               if(i2==nsample[2]-1 || sample[i1]!=sample[i1+2] ||
                  sample[i1+1]!=sample[i1+3]){
                  intnum2(sample[i1],sample[i1+1],qolint,shift,info,
                                                            k0,nsample[0],what);
                  k0=0;
               }
               k0++;
            }
         }
      }
      else{
         if(nsample[3]>0){
            i1=0;
            for(i2=0;i2<NC;i2++){
               if(ng3[i2]>0){
                  i1=i1+ng3[i2];
                  intnum3(xg[i2],qolint,d1,e1,version,bound[4],shift,info,
                          ng3[i2],nsample[0],what);
                  if(i1==nsample[3])i2=NC+3;
               }
            }
         }
         if(nsample[4]>0){
            i1=0;
            for(i2=0;i2<NC;i2++){
               if(ng4[i2]>0){
                  i1=i1+ng4[i2];
                  intnum4(xg[i2],qolint,d2,e2,version2,bound[2],shift,info,
                          ng4[i2],nsample[0],what);
                  if(i1==nsample[4])i2=NC+3;
               }
            }
         }
      }
   }
   for(i1=0;i1<nknots-1;i1++)
      shift[i1] = sufficient[i1][0] - shift[i1];
   if(accuracy==0 && nsample[2]>0)
      for(i1=0;i1<nknots-1;i1++)
         shift[i1] = sufficient[i1][1] + shift[i1];
   return czheta;
}

/******************************************************************************/
static void intnum2(double x1,double x2,double qolint[][NC+1],double shift[],
   double info[][NC],int n1,int n2,int what)
{
   int i1,i2,i3,i4,i5,jl,jr,im;
   double z0,z1[NC],z2[NC][NC],y1[7],y2[7];

/* jl - in which interval is the left endpoint
   jr - in which interval is the right endpoint                               */

   im = 3*what+1;

/* we are dealing here with one interval censored point - the lower and the
   upper bound can either be in the same interval (i2==i3) or not. In the first
   case we have to compute one set of extra integrals, in the second case two.*/
   for(jl=0; knots[jl]<x1 && jl < nknots; jl++);
   for(jr=0; knots[jr]<x2 && jr < nknots; jr++);
   for(i4=0;i4<7;i4++)y1[i4]=0;
   for(i4=0;i4<7;i4++)y2[i4]=0;
   if(jl==jr){
      (void)numints(y1,x1,x2,0,im);
   }
   else{
      (void)numints(y1,x1,knots[jl],0,im);
      (void)numints(y2,knots[jr-1],x2,0,im);
   }

/* z1 is going to be the first order contribution, z2 the second order, z0
   the normalizing constant. Here we initialize them                          */
   z0=0.;
   for(i1=0;i1<nknots;i1++){
      z1[i1]=0.;
      if(what==2)for(i2=0;i2<nknots;i2++)z2[i1][i2]=0.;
   }

/* for each interval strictly between jl and jr we now compute them           */
   for(i1=1; i1<nknots; i1++){
     if(i1>jl && i1<jr) z0 = z0 + qolint[0][i1];
     for(i3=i1-3; i3<=i1 && i3<nknots-1; i3++) if(i3>=0){
	if(i1>jl && i1<jr) 
           for(i2=0; i2<4; i2++)
	      z1[i3] = z1[i3] + coef[i3][i2][i1] * qolint[i2][i1];
        if(what==2){
           for(i4=i1-3; i4<=i1 && i4<nknots-1;i4++) 
              if(i4>=0 && i1>jl && i1<jr)
                 for(i2=0; i2<4; i2++)
                    for(i5=0;i5<4;i5++)
	               z2[i3][i4] = z2[i3][i4] + coef[i3][i2][i1]
	                             * qolint[i2+i5][i1] * coef[i4][i5][i1];
        }
     }
   }
   z0 = z0 + y1[0] + y2[0];

/* we now compute them for the interval in which the left endpoint is         */
   for(i3=jl-3; i3<=jl && i3<nknots-1; i3++) if(i3>=0){
      for(i2=0; i2<4; i2++)
         z1[i3] = z1[i3] + coef[i3][i2][jl] * y1[i2];
      if(what  == 2){
         for(i4=jl-3; i4<=jl && i4<nknots-1; i4++) 
            if(i4>=0)
               for(i2=0;i2<4;i2++)
                  for(i5=0;i5<4;i5++) 
                     z2[i3][i4] = z2[i3][i4] + coef[i3][i2][jl]
                                   * y1[i2+i5] * coef[i4][i5][jl];
      } 
   } 
/* we now compute them for the interval in which the right endpoint is - since
   y2 = 0 if this interval is the same as the left one, we do not count double*/
   for(i3=jr-3; i3<=jr && i3<nknots-1; i3++) if(i3>=0){
      for(i2=0; i2<4; i2++)
         z1[i3] = z1[i3] + coef[i3][i2][jr] * y2[i2];
      if(what  == 2){
         for(i4=jr-3; i4<=jr && i4<nknots-1; i4++)
            if(i4>=0)
               for(i2=0;i2<4;i2++)
                  for(i5=0;i5<4;i5++) 
                     z2[i3][i4] = z2[i3][i4] + coef[i3][i2][jr]
                                   * y2[i2+i5] * coef[i4][i5][jr];
      } 
   } 

/* now we update shift and info                                               */
   if(z0>0.){
      for(i3=0; i3<nknots; i3++){
	 shift[i3] = shift[i3] - n1 * z1[i3] / (z0 * n2);
         if(what == 2){
	    for(i4=0; i4<nknots; i4++)
               info[i3][i4] = info[i3][i4]+
	          n1 * (z1[i3]*z1[i4] / z0 - z2[i3][i4]) / (n2*z0);
         }
      }
   }
}
/******************************************************************************/
static void intnum3(double x,double qolint[][NC+1],double d1,double e1,int vs,
   double bd,double shift[],double info[][NC],int n1,int n2,int what)
{
   int i1,i2,i3,i4,i5,jin,im;
   double aa[7],z0,yy[7],z1[NC],z2[NC][NC];

   im = 3*what+1;

/* in which interval is the point                                             */
   for(jin=0; jin<7; jin++) aa[jin]=0.;
   for(jin=0; jin<7; jin++) yy[jin]=0.;
   for(jin=0; knots[jin]<x && jin < nknots; jin++);
   if(jin<nknots)
      (void)numints(yy,x,knots[jin],0,im);
   else{
      for(i1=0; i1<7; i1++){
         aa[i1]=1.;
         if(i1!=0) aa[i1-1]=0.;
         yy[i1]=expin2(vs,x,bd,aa,d1,e1);
      }
   }

/* initialize                                                                 */
   z0=0.;
   for(i1=0; i1<nknots; i1++){
      z1[i1]=0.;
      if(what==2) for(i2=0; i2<nknots; i2++) z2[i1][i2]=0.;
   }

/* combine for the intervals which are strictly right of the observation      */
   for(i1=1; i1<=nknots; i1++){
      if(i1>jin) z0=z0+qolint[0][i1];
      for(i3=i1-3; i3<=i1 && i3<nknots-1; i3++) if(i3>=0){
         if(i1>jin)
            for(i2=0; i2<4; i2++)
               z1[i3]=z1[i3]+coef[i3][i2][i1]*qolint[i2][i1];
         if(what==2){
            for(i4=i1-3; i4<=i1 && i4<nknots-1; i4++)
               if(i4>=0 && i1>jin)
                  for(i2=0; i2<4; i2++)
                     for(i5=0; i5<4; i5++)
                        z2[i3][i4] = z2[i3][i4] + coef[i3][i2][i1]
                                   * qolint[i2+i5][i1] * coef[i4][i5][i1];
         }
      }
   }

/* add the one in the interval in which the point is                          */
   z0=z0+yy[0];
   for(i3=jin-3; i3<=jin && i3<nknots-1; i3++) if(i3>=0){
      for(i2=0; i2<4; i2++) z1[i3] = z1[i3] + coef[i3][i2][jin] * yy[i2];
      if(what == 2){
         for(i4=jin-3; i4<=jin && i4<nknots-1; i4++)
            if(i4>=0)
               for(i2=0; i2<4; i2++)
                  for(i5=0; i5<4; i5++)
                     z2[i3][i4] = z2[i3][i4] + coef[i3][i2][jin] * yy[i2+i5]
                                  * coef[i4][i5][jin];
      }
   } 

/* update shift and info                                                      */
   if(z0>0.){
      for(i3=0; i3<nknots-1; i3++){
         shift[i3] = shift[i3] - n1 * z1[i3] / (z0 * n2);
         if(what == 2){
            for(i4=0; i4<nknots-1; i4++) info[i3][i4] = info[i3][i4] 
                         + n1 * (z1[i3] * z1[i4] / z0 - z2[i3][i4]) / (n2 * z0);
         }
      }
   }
}
/******************************************************************************/
static void intnum4(double x,double qolint[][NC+1],double d1,double e1,int vs,double bd,
       double shift[],double info[][NC],int n1,int n2,int what)
{
   int i1,i2,i3,i4,i5,jin,im;
   double aa[7],z0,yy[7],z1[NC],z2[NC][NC];

   im = 3*what+1;

/* in which interval is the point                                             */
   for(jin=0; jin<7; jin++) aa[jin]=0.;
   for(jin=0; jin<7; jin++) yy[jin]=0.;
   for(jin=0; knots[jin]<x && jin < nknots; jin++);
   if(jin>0)
      (void)numints(yy,knots[jin-1],x,0,im);
   else{
      for(i1=0; i1<7; i1++){
         aa[i1]=1.;
         if(i1!=0) aa[i1-1]=0.;
         yy[i1]=expin2(vs,x,bd,aa,d1,e1);
      }
   }

/* initialize                                                                 */
   z0=0.;
   for(i1=0; i1<nknots; i1++){
      z1[i1]=0.;
      if(what==2) for(i2=0; i2<nknots; i2++) z2[i1][i2]=0.;
   }

/* combine for the intervals which are strictly right of the observation      */
   for(i1=0; i1<nknots; i1++){
      if(i1<jin) z0=z0+qolint[0][i1];
      for(i3=i1-3; i3<=i1 && i3<nknots-1; i3++) if(i3>=0){
         if(i1<jin)
            for(i2=0; i2<4; i2++)
               z1[i3]=z1[i3]+coef[i3][i2][i1]*qolint[i2][i1];
         if(what==2){
            for(i4=i1-3; i4<=i1 && i4<nknots-1; i4++)
               if(i4>=0 && i1<jin)
                  for(i2=0; i2<4; i2++)
                     for(i5=0; i5<4; i5++)
                        z2[i3][i4] = z2[i3][i4] + coef[i3][i2][i1]
                                   * qolint[i2+i5][i1] * coef[i4][i5][i1];
         }
      }
   }

/* add the one in the interval in which the point is                          */
   z0=z0+yy[0];
   for(i3=jin-3; i3<=jin && i3<nknots-1; i3++) if(i3>=0){
      for(i2=0; i2<4; i2++) z1[i3] = z1[i3] + coef[i3][i2][jin] * yy[i2];
      if(what == 2){
         for(i4=jin-3; i4<=jin && i4<nknots-1; i4++)
            if(i4>=0)
               for(i2=0; i2<4; i2++)
                  for(i5=0; i5<4; i5++)
                     z2[i3][i4] = z2[i3][i4] + coef[i3][i2][jin] * yy[i2+i5]
                                  * coef[i4][i5][jin];
      }
   } 

/* update shift and info                                                      */
   if(z0>0.){
      for(i3=0; i3<nknots-1; i3++){
         shift[i3] = shift[i3] - n1 * z1[i3] / (z0 * n2);
         if(what == 2){
            for(i4=0; i4<nknots-1; i4++) info[i3][i4] = info[i3][i4] 
                         + n1 * (z1[i3] * z1[i4] / z0 - z2[i3][i4]) / (n2 * z0);
         }
      }
   }
}
/******************************************************************************/
static double tails(double info[][NC],double shift[],double coef[][4][NC],
    double bound[],double knots[],double zheta[],int nknots,int what)

/* tails1 computes the integrals in the tails.                                */



{
   double a[6],cth;
   int i,j=0,version;
/* local:
   a,version:    for expin
   i,j      :    counters                                                     */
      
/* Initialize cth, score and info                                             */
   cth = 0.;
   if(what>=1){
      for(i=0;i<nknots;i++){
         shift[i] = 0.;
         if(what==2) for(j=0;j<nknots;j++) info[i][j]=0.;
      }
   }

/* Now we are going to compute the integrals, first the ones  
   for the left tail, then the right tail, a large number of   
   integrals are zero, and they do not have to be computed  

   First the left-tail, only the first basisfunction is not
   equal to 0 here. It only has a linear and a constant term
   in the tails.                                                              */

   version=1;
   if(bound[1] > 0.5) version=2;
      
/* The numbers a1,b1,c1,d1,e1 are constants in the integrals                  */
   a[1] = 0.;
   a[4] = zheta[0] * coef[0][1][0];
   a[5] = zheta[0] * coef[0][0][0];
   if(what==0){
      a[2] = 0.;
      a[3] = 1.;
      cth = expin(version,knots[0],bound[2],a);
   }
   else{
      a[2] = coef[0][1][0];
      a[3] = coef[0][0][0];
      shift[0] = expin(version,knots[0],bound[2],a) /czheta;
      if(what==2){
         a[1]=coef[0][1][0]*coef[0][1][0];
         a[2]=coef[0][1][0]*coef[0][0][0]*2.;
         a[3]=coef[0][0][0]*coef[0][0][0];
         info[0][0] =expin(version,knots[0],bound[2],a)/czheta;
      }
   } 
      
/* Now the right tail. Only the last 2 basisfunctions are not
   equal to 0 here. The last one has a linear and a constant
   term, the one but last one only a constant term.                           */
   version=3;
   if(bound[3] > 0.5) version=4;
   a[1] = 0.;
   a[4] = zheta[nknots-2] * coef[nknots-2][1][nknots];
   if(nknots>3)
      a[5] = zheta[nknots-2] * coef[nknots-2][0][nknots]
           + zheta[nknots-3] * coef[nknots-3][0][nknots];
   else a[5] = zheta[nknots-2] * coef[nknots-2][0][nknots];
   if(what == 0){
      a[2] = 0.;
      a[3] = 1.;
      cth = cth + expin(version,knots[nknots-1],bound[4],a);
   }
   else {
      if(nknots>3){
         a[2] = 0.;
         a[3] = coef[nknots-3][0][nknots];
         shift[nknots-3] = expin(version,knots[nknots-1],bound[4],a)/czheta;
      }
      a[2] = coef[nknots-2][1][nknots];
      a[3] = coef[nknots-2][0][nknots];
      shift[nknots-2] = expin(version,knots[nknots-1],bound[4],a)/czheta;
   }
   if(what == 2 && nknots>3){
      a[2] = 0.;
      a[3] = coef[nknots-3][0][nknots] * coef[nknots-3][0][nknots];
      info[nknots-3][nknots-3]=expin(version,knots[nknots-1],bound[4],a)/czheta;
      a[2] = coef[nknots-3][0][nknots] * coef[nknots-2][1][nknots];
      a[3] = coef[nknots-3][0][nknots] * coef[nknots-2][0][nknots];
      info[nknots-2][nknots-3]=expin(version,knots[nknots-1],bound[4],a)/czheta;
      info[nknots-3][nknots-2] = info[nknots-2][nknots-3];
   }
   if(what == 2){
      a[1] = coef[nknots-2][1][nknots] * coef[nknots-2][1][nknots];
      a[2] = coef[nknots-2][0][nknots] * coef[nknots-2][1][nknots] *2.;
      a[3] = coef[nknots-2][0][nknots] * coef[nknots-2][0][nknots];
      info[nknots-2][nknots-2]=expin(version,knots[nknots-1],bound[4],a)/czheta;
   }
   return cth;
}
/******************************************************************************/
/* static int where(x,knots,lk)
double x,knots[];
int lk;
{
   int i;
   if(x<knots[0])return 0;
   if(x>=knots[lk-1])return lk;
   for(i=1;i<lk-1;i++)
      if(x<knots[i])return i;
   return lk-1;
} */
/******************************************************************************/
