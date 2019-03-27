 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #define MAX(x, y) (x > y ? x : y)
 #define MIN(x, y) (x < y ? x : y)
 #define NTAB 32
     
    // Global Variable Declarations
    int idum;
  
    // ran11 save statement
    int iv[NTAB], iy;

    // ran1 function declaration
    double ran11();

 int main(void) 
 {
    //  Variables Declarations
     int thet[104], gclust[104], nqi[104], nsi[104];
     double mclust[104], nqu[104], naq[104], nst[104], nas[104];
     double mstr, dpt;
     double xmax[104];
     char reffram[8];
     double ran1;

     double rapgap, bimp, RAu, drtbin, dtbin;
     int id_enh1, id_enh2, id_enh3, id_enh4, id_enh5, Nevents;

     const int i2maxbin = 38;
     const int imaxbin = 20;
     const int nID = 10;
     struct dNdpt
     {
         int nID;
         int i2maxbin;
         int imaxbin;
         int imaxbin1;
         int imaxbin2;
     };
     struct xnorm 
     {
         int nID;
     };

     const double PI = acos(-1.0);
     const double PI2 = PI * PI;
     const double HC3 = pow( 197.32858 , 3 );

    //  IMPLICIT REAL*8 (A-H,O-Z) Variables Declarations
     double dphi, xx;
     

     FILE * fptr;
     fptr = fopen("evgen_par.dat","r");
     char temp[100];

    //  id-nr. to be enhanced
     fscanf(fptr, "%d", &id_enh1);
     fgets(temp, 100, fptr);

     fscanf(fptr, "%d", &id_enh2);
     fgets(temp, 100, fptr);

     fscanf(fptr, "%d", &id_enh3);
     fgets(temp, 100, fptr);

     fscanf(fptr, "%d", &id_enh4);
     fgets(temp, 100, fptr);

     fscanf(fptr, "%d", &id_enh5);
     fgets(temp, 100, fptr);

    //  total rapidity gap
     fscanf(fptr, "%lf", &rapgap);
     fgets(temp, 100, fptr);

    //  impact parameter (irrelevant; appears only in output)
     fscanf(fptr, "%lf", &bimp);
     fgets(temp, 100, fptr);

    //  initial Radius
     fscanf(fptr, "%lf", &RAu);
     fgets(temp, 100, fptr);

    //  # of events
     fscanf(fptr, "%d", &Nevents);
     fgets(temp, 100, fptr);

    //  pt bin width, in MeV
     fscanf(fptr, "%lf", &dpt);
     fgets(temp, 100, fptr);

    //  rt bin width, in units of RAu
     fscanf(fptr, "%lf", &drtbin);
     fgets(temp, 100, fptr);

    //  t bin width, in units of RAu
     fscanf(fptr, "%lf", &dtbin);
     fgets(temp, 100, fptr);

    //  random number seed (<0)
     fscanf(fptr, "%d", &idum);
     fgets(temp, 100, fptr);

     fclose(fptr);

     dphi = 2.0*PI/imaxbin;
     printf("%f\n", dphi);

     xx = ran11 (); 
    /* ran1 Works!! */ 
    //  printf("%f\n",xx);
    //  printf("%d\n",idum); 
     
 }

 double ran11()
 {

     double ran1;
     const int IA = 16807;
     const int IM = 2147483647;
     const double AM = (1.0)/IM;
     const int IQ = 127773;
     const int IR = 2836;
     const int NDIV = 1 + (IM - 1)/NTAB;
     const double EPS = 1.2e-7;
     const double RNMX = 1.0 - EPS;

    //  printf("%d, %d, %e, %d, %d, %d, %d, %e, %f", IA, IM, AM, IQ, IR, NTAB, NDIV, EPS, RNMX); tamam

     int j, k; 

     if ( (idum <= 0) || (iy == 0) )
     {
         idum = MAX(-idum , 1);
        //  printf("%d\n", idum); tamam
         for( j = NTAB + 8; j >= 1; j--)
         {
             k = idum/IQ;
             idum = IA * (idum - k*IQ) - IR*k;
            //  printf("%d\n", idum); tamam
             if(idum < 0)
             {
                  idum = idum + IM;
             }
             if(j <= NTAB) 
             {
                 iv[j-1] = idum;
                //  printf("iv[%d] = %d\n", j, iv[j-1]); tamam
             }
         }

         iy = iv[0]; 
        //  printf ("%d\n", iy); tamam
     }

     k = idum/IQ;
     idum = IA * (idum - k*IQ) - IR*k;
     if (idum < 0)
     {
         idum = idum + IM;
     }
     j=1+iy/NDIV;
     iy=iv[j-1];
     iv[j-1] = idum;
     ran1 = MIN( (AM*iy) , RNMX );

    //  printf("%f\n",ran1); // Zalfol

     return ran1;
 }