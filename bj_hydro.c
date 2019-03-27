#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)
#define MAX3(X, Y, Z) (MAX(X, Y) > Z ? MAX(X, Y) : Z)
#define MIN3(X, Y, Z) (MIN(X, Y) < Z ? MIN(X, Y) : Z)

/* ** BJ_HYDRO Finished ** */

/* Global Variables Declarations */
double gm, B, dx, t, eta, alp, de;
float dt;
int eos;
double press, temp, velo, sound;
double total0[2];
int kp[1000], km[1000];
double pars[10];
double presstab[5001], temptab[5001];

/* File pointers declarations for write files (rhllew00 --> rhllew21) */
FILE *fptr00, *fptr01, *fptr02, *fptr03, *fptr04, *fptr05, *fptr06, *fptr07, *fptr08, *fptr09, *fptr10, *fptr11, *fptr12, *fptr13, *fptr14, *fptr15, *fptr16, *fptr17, *fptr18, *fptr19, *fptr20, *fptr21;

/* File pointer for 'fort.30' file */
FILE *fptr30;

/* File pointer for 'hyper.dat' file */
FILE *fptr40;

/* Helper Functions Declarations */
double *getParDatData();
int fortNumberOfLines();

/* Function Declaration */
void sinit(double e[1000], double gamma[1000], double *ks);
void init(double e[1000], double p[1000], double v[1000], double gamma[1000], double elab[1000], double m[1000]);
void untang(double elab[1000], double m[1000], double e[1000], double p[1000], double v[1000], double gamma[1000]);
double velo1(double el, double ml);
double press1(double a);
double press_table(double e);
double temp1(double e);
double temp_table(double e);
double sound1(double a);
double minmod(double a, double b);
double slope1(double a, double b, double ks);
void prop(double elab[1000], double m[1000], double ks);
void fileo(double elab[1000], double m[1000], double e[1000], double p[1000], double v[1000], double gamma[1000], int it, int ind);
void PressInit();
void wrhyps(double e[1000], double v[1000], double tc);
void sort(double t0);

/* Functions and Subroutine Implementations */
// subroutine sinit(e,gamma,ks) |Done: Data Checked|
void sinit(double e[1000], double gamma[1000], double *ks)
{
    double r, epc, ra, lam;
    int i;

    /* Getting data from par.dat file as an Array of integers */
    double *parFileData = getParDatData();

    /* Assigning necessary subroutine data*/
    dx = parFileData[2];
    // printf("dx: %f\n",dx);

    /* Assigning necessary subroutine data*/
    lam = parFileData[3];
    // printf("lam: %f\n",lam);
    dt = lam * dx;
    //  printf("dt: %f\n",dt);

    /* Assigning necessary subroutine data*/
    eta = parFileData[4];
    // printf("eta: %f\n",eta);

    /* Initialize ratio of degrees of freedom in QGP to hadron matter,
    adiabatic index of the EoS, epc, the initial energy density of 
    the slab in units of the critical pressure (resp. in units of the
    e0=147.7MeV/fm^3 for tab.-EoS), and the radius ra of the system.
    Then initialize Bag constant B (in units of the critical pressure), 
    and the fields e and gamma. */

    r = 37.0 / 3.0;
    gm = 1.0 / 3.0;
    // printf("%f\n",gm);

    /* Assigning necessary subroutine data*/
    epc = parFileData[5];
    // printf("epc: %f\n",epc);

    /* Assigning necessary subroutine data*/
    ra = parFileData[6];
    // printf("ra: %f\n",ra);

    B = r - 1.0;
    // printf("B: %f",B); -. Tamam

    for (int i = 1; i < 1000; i++)
    {
        /* code */
        if ((fabs((float)(i)-500.5e0) * dx) <= ra)
        {
            e[i - 1] = epc;
        }
        else
        {
            e[i - 1] = 0.0;
        }

        gamma[i] = 1.0;
    }

    // Initialize slope approximation
    *ks = 1.e0;

    // Initialize flags:
    // eos=1         equation of state is ultrarelativistic ideal gas
    // eos=0         ur-idgas + Bag Model QGP, 1st order phase transition
    // eos=2         all hadrons + Bag Model QGP, 1st order phase transition
    eos = parFileData[7];
    // printf("eos: %d\n",eos);

    //     Initialize geometry: alp = 0d0     longitudinal expansion
    //                          alp = 1d0     cylindrical expansion
    //                          alp = 2d0     spherical expansion
    //                          alp = 3d0     cylindrical plus long. scaling
    alp = parFileData[8];
    // printf("alp: %f\n",alp);
}

// subroutine init(e,p,v,gamma,elab,m) |Done Data Checked|
void init(double e[1000], double p[1000], double v[1000], double gamma[1000], double elab[1000], double m[1000])
{
    int i;

    // Start calculation. Null vector total0
    total0[0] = 0.0;
    total0[1] = 0.0;

    // Start calculation. Each cell is treated separately.
    for (i = 0; i < 1000; i++)
    {
        // Calculate index pointer (used in function prop).
        kp[i] = i + 1;
        if (i == 999)
        {
            kp[i] = 999;
        }

        km[i] = i - 1;
        if (i == 0)
        {
            km[i] = 0;
        }

        if (e[i] <= 0.0)
        {
            // Zero or negative temperatures correspond to vacuum.
            // The vacuum is supposed to have velocity -1 in (-x)-direction
            // and velocity +1 in (+x)-direction.
            e[i] = 0.0;
            p[i] = 0.0;

            if (i <= 500)
                v[i] = -1.0;

            if (i > 500)
                v[i] = 1.0;

            elab[i] = 0.0;
            m[i] = 0.0;
            gamma[i] = 1.e+8;
        }
        else
        {
            if ((fabs(gamma[i] - 1.0)) < (1.e-16))
                v[i] = 0.0;

            else if (gamma[i] < 0.0)
                v[i] = -sqrt(1.0 - 1.0 / (gamma[i] * gamma[i]));

            else
                v[i] = sqrt(1.0 - 1.0 / (gamma[i] * gamma[i]));

            gamma[i] = fabs(gamma[i]);

            p[i] = press1(e[i]);
            m[i] = gamma[i] * gamma[i] * (e[i] + p[i]) * v[i];
            elab[i] = gamma[i] * gamma[i] * (e[i] + p[i]) - p[i];
        }

        if (alp != 0.e0)
        {
            total0[0] = total0[0] + elab[i] * pow(dx * fabs((double)(i)-500.5e0), alp) * dx;
            total0[1] = total0[1] + m[i] * pow(dx * fabs((double)(i)-500.5e0), alp) * dx;
        }
        else
        {
            total0[0] = total0[0] + elab[i] * dx;
            total0[1] = total0[1] + m[i] * dx;
        }
        // printf("%d  %d  %d\n", i, km[i], kp[i]);
    }

    return;
}

// subroutine untang(elab,m,e,p,v,gamma) |Done: Data Checked|
void untang(double elab[1000], double m[1000], double e[1000], double p[1000], double v[1000], double gamma[1000])
{
    int i;

    for (i = 0; i < 1000; i++)
    {
        if ((fabs(m[i]) < 1.e-16) && (fabs(elab[i]) < 1.e-16))
        {
            e[i] = 0.0;
            if (i < 500)
                v[i] = -1.0;
            if (i >= 500)
                v[i] = 1.0;

            gamma[i] = 1.e8;
        }
        else if (fabs(m[i]) < 1.e-16)
        {
            e[i] = elab[i];
            v[i] = 0.e0;
            gamma[i] = 1.e0;
        }
        else
        {
            v[i] = velo1(elab[i], m[i]);
            // printf("%f\n", v[i]);
            if (fabs(v[i]) < 1.e0 - 1.e-16)
            {
                e[i] = elab[i] - v[i] * m[i];
                gamma[i] = (1.0) / sqrt(1.0 - v[i] * v[i]);
            }
            else
            {
                e[i] = 0.e0;
                if (i < 500)
                    v[i] = -1.e0;
                if (i >= 500)
                    v[i] = 1.e0;

                gamma[i] = 1.e8;
            }
        }
        p[i] = press1(e[i]);

        // printf("%d   %f    %f    %f    %f\n", i+1, e[i], p[i], v[i], gamma[i]);
    }

    return;
}

// function velo(el,ml) |Done: but missing the write on terminate statements|
double velo1(double el, double ml)
{
    // printf("%f  %f\n", el, ml);
    double v, f, root, temp;
    int i;

    i = 0;
    f = 0.0;

   while(fabs((v-f)/f) > 1.e-16)
   {
       v = f;
       temp = press1(el - v*ml);
       root = sqrt(1.0 - v*v);
       f = ml/(el + temp);
    //    printf("%d %f\n", i+1, f);
       i = i + 1;
       if (i > 100)
        {
            fprintf(fptr00, "No root found for Elab=  %f , M= %f", el, ml);
            fprintf(fptr00, "Program terminated.");
            exit(0);
        }

   }

    velo = f;
    return velo;
}

// function press(a) |Done|
double press1(double a)
{

    /* This function-subprogram determines the pressure of the
       underlying EoS.
       It is used in subroutine velo, untang.
       a is the energy density in the local rest frame. */

    if (eos == 1)
    {
        press = gm * a;
    }
    else if (eos == 0)
    {
        if (a >= (3.0 + 4.0 * B))
        {
            press = gm * (a - 4.0 * B);
        }
        else if (a >= 3.0)
        {
            press = 1.0;
        }
        else
        {
            press = gm * a;
        }
    }
    else if (eos == 2)
    {
        press = press_table(a);
    }
    return press;
}

// function press_table(e) |Done|
double press_table(double e)
{
    double press_table, p1, p2;
    int i = (int)(e / de);
    double offe = e - ((double)(i)) * de;

    if (i >= 0)
    {
        if (i > 4999)
        {
            press_table = (e - 10.5e0) / 3.e0;
            return press_table;
        }
        p1 = presstab[i];
        p2 = presstab[i + 1];
        press_table = (p2 - p1) / de * offe + p1;
    }

    /* Ask the doctor where did he get the press from, and what does the fortran function return */
    if (press < 0.0)
    {
        press = 0.0;
    }

    return press_table;
}

// function temp(e) |Done: Data Checked|
double temp1(double e)
{
    double ttc;

    if(eos == 1)
    {
        double t = gm/(1.0 + gm);
        ttc = pow((e/3.0) , t);
    } else if (eos == 0) {
        if(e >= (3.0 + 4.0*B)) 
        {
            ttc = pow(((e-B)/3.0/(1.0+B)) , t);
        } else if(e >= 3.0) {
            ttc = 1.0;
        } else {
            ttc = pow((e/3.0) , t);
        }
    } else if(eos == 2) {
        ttc = temp_table(e);
    }
    temp = ttc;
    // printf("%f\n", temp);
    return temp;
}

// function temp_table(e) |Done|
double temp_table(double e)
{
    int i;
    double offe, temp_table, p1, p2;

    i = (int)(e / de);
    // printf("%d\n",i);
    offe = e - (double)(i)*de;

    if (i >= 0)
    {
        if (i > 4999)
        {
            temp_table = 1.e+3;
            return temp_table;
        }
        p1 = temptab[i];
        p2 = temptab[i + 1];
        temp_table = (p2 - p1) / de * offe + p1;
    }
    else
    {
        temp_table = 0.e0;
    }

    return temp_table;
}

// function sound(a) |Done|
double sound1(double a)
{
    if (eos == 1)
    {
        sound = sqrt(gm);
    }
    else
    {
        sound = sqrt((1.e0) / (3.e0));
    }

    /* Ask Doctor: commented code */
    // if(a >= ((3.e0)+(4.e0)*B))
    // {
    //     sound = sqrt(gm);
    // } else if(a >= 3.e0)
    // {
    //     sound = 1.e-6;
    //     sound = sqrt(gm);
    // } else {
    //     sound = sqrt(gm);
    // }

    return sound;
}

// function minmod(a,b) |Done|
double minmod(double a, double b)
{
    double minmod;

    if ((a * b) <= 0.0)
    {
        minmod = 0.0;
    }
    else if (a < 0.0)
    {
        minmod = MAX(a, b);
    }
    else
    {
        minmod = MIN(a, b);
    }

    return minmod;
}

// function slope(a,b,ks) |Done|
double slope1(double a, double b, double ks)
{
    double slope;
    slope = minmod(a, b) * ks;

    return slope;
}

// subroutine prop(elab,m,ks) |Done: - Data Checked|
void prop(double elab[1000], double m[1000], double ks)
{
    /* This subroutine-subprogram propagates the laboratory frame
    quantities Elab,M. It is used in the main program.
    Since a second order scheme is used, laboratory frame quantities
    are computed at cell interfaces. */

    // Variable decarations
    double mt[1000], elabt[1000], pos;
    double dd, sslope, bp, bm, flux, csmean, vmean, sqrl, sqrr, slope;
    double denomi, sound;
    double csp[1000], csm[1000];
    double a1[1000], a2[1000], a3[1000];
    double elm[1000], mm[1000], em[1000], pm[1000], vm[1000], gamm[1000];
    double elp[1000], mp[1000], ep[1000], pp[1000], vp[1000], gamp[1000];
    double dpv[1000], gv[1000];
    double fem[1000], fmm[1000];
    double fep[1000], fmp[1000];
    int i;

    // Start calculation. dd gets its half-step value
    dd = 0.5 * dt / dx;
    // printf("%f\n",dd); // Checked

    /* Step 1: Slope and interface values are calculated. */

    // For Elab:
    for (i = 0; i < 1000; i++)
    {
        dpv[i] = (elab[kp[i]] - elab[i]) / dx;
        // printf("%d  %f\n",i+1, dpv[i]);
    }

    for (i = 0; i < 1000; i++)
    {
        sslope = 0.5 * dx * slope1(dpv[km[i]], dpv[i], ks);
        elp[i] = elab[i] + sslope;
        elm[i] = elab[i] - sslope;
        // printf("%d    %f   %f\n", i+1, elp[i], elm[i]);
    }

    // For M:
    for (i = 0; i < 1000; i++)
    {
        dpv[i] = (m[kp[i]] - m[i]) / dx;
        // printf("%d  %f\n", i+1, dpv[i]);
    }

    for (i = 0; i < 1000; i++)
    {
        sslope = 0.5 * dx * slope1(dpv[km[i]], dpv[i], ks);
        if ((m[i] + sslope) > 0.0)
        {
            mp[i] = MIN((m[i] + sslope), elp[i]);
        }
        else
        {
            mp[i] = MAX((m[i] + sslope), -elp[i]);
            if (mp[i] == -0.0)
            {
                mp[i] = 0.0;
            }
        }

        if ((2.0 * m[i] - mp[i]) > 0.0)
        {
            mm[i] = MIN((2.0 * m[i] - mp[i]), elm[i]);
            if (mm[i] == -0.0)
            {
                mm[i] = 0.0;
            }
        }
        else
        {
            mm[i] = MAX((2.0 * m[i] - mp[i]), -elm[i]);
            if (mm[i] == -0.0)
            {
                mm[i] = 0.0;
            }
        }

        // printf("%d  %f  %f\n", i+1, mp[i], mm[i]);
    }

    /* Step 2: Calculate rest frame variables via untang. */

    untang(elp, mp, ep, pp, vp, gamp);
    untang(elm, mm, em, pm, vm, gamm);

    for (i = 0; i < 1000; i++)
    {
        // printf ("%d   %f   %f   %f   %f   %f   %f\n", i+1, elp[i], mp[i], ep[i], pp[i], vp[i], gamp[i]);
        // printf ("%d   %f   %f   %f   %f   %f   %f\n", i+1, elm[i], mm[i], em[i], pm[i], vm[i], gamm[i]);
    }

    /* Step 3: Calculate interface values in half step. */

    for (i = 0; i < 1000; i++)
    {
        flux = dd * (mp[i] * vp[i] + pp[i] - mm[i] * vm[i] - pm[i]);
        mp[i] = mp[i] - flux;
        mm[i] = mm[i] - flux;
        flux = dd * ((pp[i] + elp[i]) * vp[i] - (pm[i] + elm[i]) * vm[i]);
        elp[i] = elp[i] - flux;
        elm[i] = elm[i] - flux;

        //    Check physical consistence of half step values.
        if ((elp[i] - fabs(mp[i])) < 0.0)
        {
            if (fabs(mp[i]) > 1.e-16)
            {
                mp[i] = elp[i] * mp[i] / fabs(mp[i]);
                //    fprintf(fptr00,"   M+(',%d,') changed to Elab+ at Step 3",i); Commented in the source code

                // After correction, linear slope is enforced (but not if it
                // produces acausalities in elm(i), mm(i) ).

                if (mm[i] > 0.0)
                {
                    mm[i] = MIN((2.0 * m[i] - mp[i]), elm[i]);
                }
                else
                {
                    mm[i] = MAX((2.0 * m[i] - mp[i]), -elm[i]);
                }
            }
            else
            {
                mp[i] = 0.0;
            }
        }

        if ((elm[i] - fabs(mm[i])) < 0.0)
        {
            elm[i] = fabs(mm[i]);
            //    fprintf(fptr00,"   Elab-(',%d,') changed to ,M- at Step 3",i); Commented in the source code

            //    After correction, linear slope is enforced (but not if it
            //    produces acausalities in elp(i), mp(i) ).

            elp[i] = MAX((2.0 * elab[i] - elm[i]), fabs(mp[i]));
        }

        // printf("%d  %f   %f   %f   %f\n", i+1, mp[i], mm[i], elm[i], elp[i]);
    }

    /* Step 4: Update local rest frame variables at cell interfaces via untang */

    untang(elp, mp, ep, pp, vp, gamp);
    untang(elm, mm, em, pm, vm, gamm);

    for (i = 0; i < 1000; i++)
    {
        // printf ("%d   %f   %f   %f   %f   %f   %f\n", i+1, elp[i], mp[i], ep[i], pp[i], vp[i], gamp[i]);
        // printf ("%d   %f   %f   %f   %f   %f   %f\n", i+1, elm[i], mm[i], em[i], pm[i], vm[i], gamm[i]);
    }

    /* Step 5: Calculate sound velocity at the interface. */

    for (i = 0; i < 1000; i++)
    {
        csp[i] = sound1(ep[i]);
        csm[i] = sound1(em[kp[i]]);
        // printf("%d   %f   %f\n", i+1, csp[i], csm[i]);
    }

    /* Step 6: Estimate signal velocities. */

    for (i = 0; i < 1000; i++)
    {
        // Roe-mean for the velocities

        sqrl = sqrt(elp[i]);
        sqrr = sqrt(elm[kp[i]]);
        denomi = sqrl + sqrr;
        if (denomi > 1.e-16)
        {
            vmean = (vp[i] * sqrl + vm[kp[i]] * sqrr) / denomi;
            csmean = sqrt((pow(csp[i], 2) * sqrl + pow(csm[kp[i]], 2) * sqrr) / denomi + eta * sqrl * sqrr / denomi / denomi * pow((vp[i] - vm[kp[i]]), 2));
            bp = MAX3(0.0, (vmean + csmean) / (1.0 + vmean * csmean), (vm[kp[i]] + csm[kp[i]]) / (1.0 + vm[kp[i]] * csm[kp[i]]));
            bm = MIN3(0.0, (vmean - csmean) / (1.0 - vmean * csmean), (vp[i] - csp[i]) / (1.0 - vp[i] * csp[i]));
        }
        else
        {
            bp = 1.0;
            bm = -1.0;
        }

        /* Step 7, first part: calculate coefficients a. */

        a1[i] = bp / (bp - bm);
        a2[i] = bm / (bp - bm);
        a3[i] = bp * bm / (bp - bm);

        // printf("%d   %f   %f   %f\n", i+1, a1[i], a2[i], a3[i]);
    }

    /* Step 7, second part: calculate fluxes. */
    for (i = 0; i < 1000; i++)
    {
        fep[i] = (pp[i] + elp[i]) * vp[i];
        fmp[i] = mp[i] * vp[i] + pp[i];
        fem[i] = (pm[i] + elm[i]) * vm[i];
        fmm[i] = mm[i] * vm[i] + pm[i];
        // printf("%d   %f   %f   %f   %f\n", i+1, fep[i], fmp[i], fem[i], fmm[i]);
    }

    /* Step 8 and 9 : Calculation of the total numerical fluxes and the transported quantities at full time step. */

    dd = dd + dd;

    for (i = 0; i < 1000; i++)
    {
        gv[i] = a1[i] * fep[i] - a2[i] * fem[kp[i]] + a3[i] * (em[kp[i]] - ep[i]);
        // printf("%d   %f\n", i+1, gv[i]);
    }

    for (i = 0; i < 1000; i++)
    {
        elabt[i] = elab[i] + dd * (gv[km[i]] - gv[i]);
        // printf("%d   %f\n", i+1, elabt[i]);
    }

    for (i = 0; i < 1000; i++)
    {
        gv[i] = a1[i] * fmp[i] - a2[i] * fmm[kp[i]] + a3[i] * (mm[kp[i]] - mp[i]);
        // printf("%d   %f\n", i+1, gv[i]);
    }

    for (i = 0; i < 1000; i++)
    {
        mt[i] = m[i] + dd * (gv[km[i]] - gv[i]); // Update: The loop does not stuck here
        // printf("%d  %f\n", i+1, mt[i]);
    }

    // printf("%d  %f\n", km[400] - 1, gv[km[400] - 1] - gv[400]);

    /* Step 9a: remove acausalities */

    for (i = 0; i < 1000; i++)
    {

        if ((elabt[i] - fabs(mt[i])) < 0.0) //Update: This is nit where the loop stucks
        {
            if (fabs(mt[i]) > 1.e-16)
            {
                mt[i] = elabt[i] * mt[i] / fabs(mt[i]) * (1.0 - 1.e-10);
                // fprintf(fptr00,"MT(',%d,') changed to Elabt at Step 9a",i); Commented in the source code
            }
            else
            {
                mt[i] = 0.0;
            }
        }

        // printf("%d   %f\n", i+1, mt[i]);
    }

    /* Step 9b: Sod's method to correct for geometry
             First call untangle, since the source terms
             contain variables in the local rest frame.
             These variables are temporarily stored in
             the fields ep,pp,vp,gamp. */

    untang(elabt, mt, ep, pp, vp, gamp); //Update the code stucks here
    for (i = 0; i < 1000; i++)
    {
        pos = ((double)(i)-500.5e0) * dx;

        if (alp == 3.e0)
        {
            elab[i] = elabt[i] - dt * (vp[i] / pos + 1.0 / t) * (elabt[i] + pp[i]);
            m[i] = mt[i] - dt * (vp[i] / pos + 1.0 / t) * mt[i];
        }
        else
        {
            elab[i] = elabt[i] - dt * alp * vp[i] / pos * (elabt[i] + pp[i]);
            m[i] = mt[i] - dt * alp * vp[i] / pos * mt[i];
        }

        // printf("%d  %f  %f\n", i+1, elab[i], m[i]);
    }

    /* Step 10: The last consistency check is to remove acausalities produced in the full step propagation. */

    for (i = 0; i < 1000; i++)
    {
        if ((elab[i] - fabs(m[i])) < 0.0)
        {
            if (fabs(m[i]) > 1.e-16)
            {
                m[i] = elab[i] * m[i] / fabs(m[i]) * (1.e0 - 1.e-10);
                // fprintf(fptr00,"M(',%d,') changed to Elab at Step 10",i); Commented in the source code
            }
            else
            {
                m[i] = 0.0;
            }
        }

        // printf("%d   %f\n", i+1, m[i]);
    }

    /* Step 11: Cells with smaller lab energy density than a cut-off are regarded as vacuum. */
    for (i = 0; i < 1000; i++)
    {
        if (elab[i] < 1.e-16)
        {
            elab[i] = 0.0;
            m[i] = 0.0;
        }
        // printf("%d  %f  %f\n", i + 1, elab[i], m[i]);
    }

    return;
}

// subroutine fileo(elab,m,e,p,v,gamma,it,ind) |Done: Data Checked|
void fileo(double elab[1000], double m[1000], double e[1000], double p[1000], double v[1000], double gamma[1000], int it, int ind)
{
    // Variables declarations
    double total[2], pos[1000];
    double ttc, temp;
    int i, index;

    // Total is a vector which contains the sum of the primary variables
    // over all cells. It is set zero in the beginning.

    total[0] = 0.e0;
    total[1] = 0.e0;

    // Pos is a vector whose element i contains the centre position of
    // cell no. i.

    pos[0] = (0.5e0) * dx;
    for (i = 1; i < 1000; i++)
    {
        /* code */
        pos[i] = ((double)(i + 1) - (0.5e0)) * dx;
    }

    // Start writing.
    index = 7 + it / ind;
    fprintf(fptr00, "time=  %f  file index = %d\n", t, (index - (7 - 1)));

    // printf("%d\n",index);
    // Write quantities for all cells with pos(i)>0.
    switch (index)
    {
    case 7:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr01, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 8:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr02, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 9:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr03, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 10:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr04, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 11:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr05, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 12:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr06, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 13:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr07, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 14:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr08, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 15:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr09, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 16:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr10, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 17:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr11, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 18:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr12, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 19:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr13, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 20:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr14, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 21:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr15, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 22:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr16, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 23:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr17, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 24:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr18, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 25:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr19, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 26:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr20, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    case 27:
        for (i = 500; i < 1000; i++)
        {
            ttc = temp1(e[i]);

            double elem1 = (5.e2) * dx;
            elem1 = pos[i] - elem1;
            fprintf(fptr21, "%e  %e  %e  %e  %e  %e\n", elem1, ttc, elab[i], v[i], e[i], m[i]);
        }
        break;
    }

    for (i = 500; i < 1000; i++)
    {
        if (alp != 0.e0)
        {
            total[0] = total[0] + elab[i] * pow(dx * fabs((double)(i)-500.5e0), alp) * dx;
            total[1] = total[1] + m[i] * pow(dx * fabs((double)(i)-500.5e0), alp) * dx;
        }
        else
        {
            total[0] = total[0] + elab[i] * dx;
            total[0] = total[0] + m[i] * dx;
        }
    }
}

// subroutine PressInit |Done: Data Checked|
void PressInit()
{
    double d1, d2;

    FILE *fptr81;
    fptr81 = fopen("eos.dat", "r");

    double temp;

    fscanf(fptr81, "%lf\n", &de);
    for (int i = 0; i <= 5000; i++)
    {
        fscanf(fptr81, "%lf %lf %lf %lf %lf \n", &d1, &presstab[i], &d2, &temptab[i], &temp);
    }
}

// subroutine wrhyps(e,v,tc) implicit real*8 (a-h,o-z) |Done: - Data Checked|
void wrhyps(double e[1000], double v[1000], double tc)
{
    // Variable Declarations
    double pos, posi, ttc, ttcm, ei, vi, del_x;
    int i, j, index;

    pos = -(1000.0 * 0.5 + 0.5) * dx;
    index = 1;
    pos = pos + dx;

    // printf("%f\n",pos);

    for (i = 1; i < 1000; i++)
    {
        pos = pos + dx;

        // Find points which have given energy density
        ttc = e[i];
        ttcm = e[i - 1];

        // Search routine uses the fact that for given time T (or e)
        // first increases for increasing cell index.

        if ((ttc >= tc) && (index == 1))
        {
            j = abs(i - 500);
            fprintf(fptr30, "%f        %f        %f       %f        %f\n", pos, t, e[i], v[i], 0.0);
            index = 2;
            // printf("%d  %f  %f  %f  %f\n", i ,pos, t, e[i], v[i]);
        }
        else if ((ttc < tc) && (index == 2))
        {
            j = abs(i - 500);
            del_x = (tc - ttc) / (ttc - ttcm); // interpolate to position between
            posi = pos + del_x * dx;           // the cells
            ei = e[i] + del_x * (e[i] - e[i - 1]);
            vi = v[i] + del_x * (v[i] - v[i - 1]);
            fprintf(fptr30, "%f         %f        %f        %f       %f\n", posi, t, ei, vi, 0.0);
            index = 1;
            // printf("%d  %f  %f  %f  %f\n", i ,posi, t, ei, vi);
        }
    }


    return;
}

// subroutine sort(t0) |Incomplete|
void sort(double t0)
{
    // Variable decarations
    double x[10000], t[10000], v[10000], abs[10000];
    double e[10000], rho[10000];
    double pos, time, vel, eps;
    double xn[10000], tn[10000], vn[10000], min, xx;
    double en[10000], rhon[10000];
    double xn2[10000], tn2[10000], vn2[10000];
    double en2[10000], rhon2[10000];
    double rhob, fo_temp, fo_muq, fo_mus;
    int i, j, k, l, n, ind, count;

    /* Read in cells on hypersurface. */

    fptr40 = fopen("hyper.dat", "w");

    n = -1;
    int nol = fortNumberOfLines();
    for (i = 0; i < nol - 1; i++)
    {
        fscanf(fptr30,"%lf  %lf  %lf  %lf  %lf\n", &pos, &time, &eps, &vel, &rhob);
        // if(i >= nol - 1)
        // {
        //     pos = 0;
        //     time = 0;
        //     eps = 0;
        //     vel = 0;
        //     rhob = 0;
        // }
        if (pos >= 0.0)
        {
            n = n + 1;
            x[n] = pos;
            t[n] = time;
            v[n] = vel;
            e[n] = eps;
            rho[n] = rhob;
        }
        // printf("%d  %f\n", i+1, x[n]);
        // printf("%d  %f\n", i+1, t[n]);
        // printf("%d  %f\n", i+1, v[n]);
        // printf("%d  %f\n", i+1, e[n]);
        // printf("%d  %f\n", i+1, rho[n]);
    }
    // printf("%d\n",n);

    if (n == -1)
    {
        return;
    }

    /* Sort (x,t) combinations. */
    xn[0] = x[0];
    tn[0] = t[0];
    vn[0] = v[0];
    en[0] = e[0];
    rhon[0] = rho[0];

    // printf("%d  %f  %f  %f  %f  %f\n", i+1, xn[0], tn[0], vn[0], en[0], rhon[0]);

    for (i = 0; i < n ; i++)
    {
        /* Calculate vector of distances between remaining points
        and point i. */
        for (j = i + 1; j <= n; j++)
        {
            abs[j] = pow((tn[i] - t[j]), 2.0) + pow((xn[i] - x[j]), 2.0);
            // printf("%d %f\n", j+1, abs[j]);
        }
        // Error Here Stopped

        min = abs[i + 1];
        // printf("%d %f\n",i+1, min);
        l = i + 1;
        // printf("%d\n",l);
        

        for (k = i + 2; k <= n; k++)
        {
            if (abs[k] < min)
            {
                min = abs[k];
                l = k - 1;
            }
        }

        xn[i + 1] = x[l];
        tn[i + 1] = t[l];
        vn[i + 1] = v[l];
        en[i + 1] = e[l];
        rhon[i + 1] = rho[l];

        x[l] = x[i + 1];
        t[l] = t[i + 1];
        v[l] = v[i + 1];
        e[l] = e[i + 1];
        rho[l] = rho[i + 1];

        // printf("%d  %f\n", i+1, x[i]);
        // printf("%d  %f\n", i+1, t[i]);
        // printf("%d  %f\n", i+1, v[i]);
        // printf("%d  %f\n", i+1, e[i]);
        // printf("%d  %f\n", i+1, rho[i]);
    }

    /* If positions are the same for different times, take
    only the arithmetic mean time at that position. */
    // printf("%f\n", t0);

    fprintf(fptr40, "%f  %f  %f  %f  %f  %f  %f\n", 1.0, t0, 0.0, 0.0, 0.0, 0.0, 0.0);
    i = 0;
    l = 1;

T31:
    xn2[i] = xn[l - 1];
    tn2[i] = tn[l - 1];
    en2[i] = en[l - 1];
    vn2[i] = vn[l - 1];
    rhon2[i] = rhon[l - 1];

    count = 1;

    /* Search next 100 entries if position is equal. */

    for (j = 0; j < 100; j++)
    {
        if (xn[j + l] == xn[l - 1])
        {
            tn2[i] = tn2[i] + tn[j + l];
            en2[i] = en2[i] + en[j + l];
            vn2[i] = vn2[i] + vn[j + l];
            rhon2[i] = rhon2[i] + rhon[j + l];
            count = count + 1;
        }
        else
        {
            goto T33;
        }
    }

T33:
    tn2[i] = tn2[i] / count;
    en2[i] = en2[i] / count;
    vn2[i] = vn2[i] / count;
    rhon2[i] = rhon2[i] / count;

    //  fo_temp=temp(en2[i],rhon2[i])
    fo_temp = temp1(en2[i]);
    // printf("%f\n", fo_temp);
    fo_muq = 0.0;
    fo_mus = 0.0;

    fprintf(fptr40, "%f  %f  %f  %f  %f  %f  %f\n", xn2[i], tn2[i], en2[i], vn2[i], fo_temp, fo_muq, fo_mus);

    i = i + 1;
    l = l + count;

    if (l < n + 1)
    {
        goto T31;
    }

    if (xn2[i - 1] > dx)
    {
        xx = xn2[i - 1] - dx;

    T34:
        fprintf(fptr40, "%f  %f  %f  %f  %f  %f  %f\n", xx, tn2[i - 1], en2[i - 1], 0.0, fo_temp, fo_muq, fo_mus);

        if (xx > dx)
        {
            xx = xx - dx;
            goto T34;
        }
    }

    fprintf(fptr40, "%f  %f  %f  %f  %f  %f  %f\n", 0.0, tn2[i - 1], 0.0, 0.0, 0.0, 0.0, 0.0);

    fclose(fptr30);
    fclose(fptr40);

    return;
}

int main(void)
{
    double e[1000], p[1000], v[1000], gamma[1000], elab[1000], m[1000], t0, eps_hyper, ks;
    int itplot, it, ind;

    /* First part of initialization (continued in sinit) */
    double *parFileData = getParDatData();

    /* Initialize t0 and number of time steps */

    t0 = parFileData[0];
    itplot = parFileData[1];
    sinit(e, gamma, &ks);
    if (eos == 2)
        PressInit();
    init(e, p, v, gamma, elab, m);
    eps_hyper = parFileData[9];

    t = t0;
    fptr00 = fopen("rhllew00.dat", "w");
    fptr01 = fopen("rhllew01.dat", "w");
    fptr02 = fopen("rhllew02.dat", "w");
    fptr03 = fopen("rhllew03.dat", "w");
    fptr04 = fopen("rhllew04.dat", "w");
    fptr05 = fopen("rhllew05.dat", "w");
    fptr06 = fopen("rhllew06.dat", "w");
    fptr07 = fopen("rhllew07.dat", "w");
    fptr08 = fopen("rhllew08.dat", "w");
    fptr09 = fopen("rhllew09.dat", "w");
    fptr10 = fopen("rhllew10.dat", "w");
    fptr11 = fopen("rhllew11.dat", "w");
    fptr12 = fopen("rhllew12.dat", "w");
    fptr13 = fopen("rhllew13.dat", "w");
    fptr14 = fopen("rhllew14.dat", "w");
    fptr15 = fopen("rhllew15.dat", "w");
    fptr16 = fopen("rhllew16.dat", "w");
    fptr17 = fopen("rhllew17.dat", "w");
    fptr18 = fopen("rhllew18.dat", "w");
    fptr19 = fopen("rhllew19.dat", "w");
    fptr20 = fopen("rhllew20.dat", "w");
    fptr21 = fopen("rhllew21.dat", "w");

    /* intermediate hypersurface file */
    fptr30 = fopen("fort.30", "w");

    fprintf(fptr00, "et=  %f Gamma-1=  %e\nmax. no. of steps= %d dx= %f, dt= %f\n", eta, gm, itplot, dx, dt);
    ind = itplot / 50;
    fileo(elab, m, e, p, v, gamma, 0, ind);

    for (it = 0; it < itplot; it++) /* main loop (time-steps) */
    {
        t = t + dt;
        prop(elab, m, ks);
        untang(elab, m, e, p, v, gamma);
        if ((it == ind) || (it == (ind * 2)) || (it == (ind * 3)) ||
            (it == (ind * 4)) || (it == (ind * 5)) || (it == (ind * 6)) ||
            (it == (ind * 7)) || (it == (ind * 8)) || (it == (ind * 9)) ||
            (it == (ind * 10)) || (it == (ind * 11)) || (it == (ind * 12)) ||
            (it == (ind * 13)) || (it == (ind * 14)) || (it == (ind * 15)) ||
            (it == (ind * 16)) || (it == (ind * 17)) || (it == (ind * 18)) ||
            (it == (ind * 19)) || (it == (ind * 20)))
        {
            fileo(elab, m, e, p, v, gamma, it, ind);
            // printf("fileo called at it = %d\n", it);
        }

        /* write out cells on the hypersurface */
        wrhyps(e, v, eps_hyper);
    } /* main loop (time-steps) */

    fclose(fptr30);

    for (int i = 0; i < 1000; i++)
    {
        // printf ("%d   %f   %f   %f   %f   %f   %f\n", i+1, elab[i], m[i], e[i], p[i], v[i], gamma[i]);
    }

    fclose(fptr00);
    fclose(fptr01);
    fclose(fptr02);
    fclose(fptr03);
    fclose(fptr04);
    fclose(fptr05);
    fclose(fptr06);
    fclose(fptr07);
    fclose(fptr08);
    fclose(fptr09);
    fclose(fptr10);
    fclose(fptr11);
    fclose(fptr12);
    fclose(fptr13);
    fclose(fptr14);
    fclose(fptr15);
    fclose(fptr16);
    fclose(fptr17);
    fclose(fptr18);
    fclose(fptr19);
    fclose(fptr20);
    fclose(fptr21);

    fptr30 = fopen("fort.30" , "r");
    sort(t0);
    fclose(fptr30);

    

    return 0;
}

// Function to get data from par.dat file into an array
double *getParDatData()
{
    FILE *fpointer;
    fpointer = fopen("par.dat", "r");

    char line[100];
    char delimt[] = " ";
    int i = -1;

    while (!feof(fpointer))
    {
        fgets(line, 100, fpointer);
        char *ptr = strtok(line, delimt);

        while (ptr != NULL)
        {
            i++;
            pars[i] = strtod(ptr, &ptr);
            ptr = strtok(NULL, delimt);
        }
    }

    fclose(fpointer);

    return pars;
}

int fortNumberOfLines() {
    FILE *fpointer;
    fpointer = fopen("fort.30", "r");
    char line[100];
    int count = 0;

    while (!feof(fpointer))
    {
        fgets(line, 100, fpointer);
        count ++;
    }

    fclose(fpointer);

    return count;
}
