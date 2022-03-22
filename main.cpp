
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

#define N 3
#define ND 50
#define PI 3.14159

int nm = N - 1;
int ndm = ND - 1;

int nstep = 1001;
int pstep = 100;

double dx = 1.0e-7;
double dtime = 4.0e-10;

double gamma0 = 0.1;
double delta = 5.0 * dx;

double mobi = 2.0e-5;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double F10 = 1.5e6;
double F20 = 1.5e6;

double Ds = 0.1e-6;
double Dl = 0.2e-5;

double temp = 600.0; // K
double Te = 800.0;
double ce = 0.4;
double cl = 0.4;

double ml1 = -2000.0;
double kap1 = 0.25;
double c01e = ce + (temp - Te) / ml1;
double c1e = c01e * kap1;

double ml2 = 2000.0;
double kap2 = 2.5;
double c02e = ce + (temp - Te) / ml2;
double c2e = c02e * kap2;

double c0, dc0;

double cddtt, dev1_s, dev2_s, dev1_l, dev2_l;

double con[ND][ND][ND], con_new[ND][ND][ND], con1[ND][ND][ND], con2[ND][ND][ND], con0[ND][ND][ND];

double mij[N][N], aij[N][N], wij[N][N], fij[N][N];

int anij[N][N];
double thij[N][N];
double vpij[N][N];
double etaij[N][N];

double phi[N][ND][ND][ND], phi_new[N][ND][ND][ND];

double phis_ijk, phis_ipjk, phis_imjk, phis_ijpk, phis_ijmk, phis_ijkp, phis_ijkm;
double cons_ijk, cons_ipjk, cons_imjk, cons_ijpk, cons_ijmk, cons_ijkp, cons_ijkm;

int phinum;
int phiNum[ND][ND][ND];
int phiIdx[N + 1][ND][ND][ND];

int i, j, im, ip, jm, jp, k, km, kp, l;
int ii, jj, kk;
int n1, n2, n3;
int istep;

double astre = 0.05;
double th, vp, eta;
double thetax, thetay;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;
double phidxy, phidxz, phidyz;
double phiabs;

double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;

double phidxp, phidyp, phidzp;
double phidxpx, phidypx, phidzpx;
double phidxpy, phidypy, phidzpy;
double phidxpz, phidypz, phidzpz;
double ep, epdx, epdy, epdz;
double term0;
double termx, termx0, termx1, termx0dx, termx1dx;
double termy, termy0, termy1, termy0dy, termy1dy;
double termz, termz0, termz1, termz0dz, termz1dz;

double pddtt, sum1;
// double termiikk, termjjkk;

double dF;

void datasave(int step);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << dtime * Dl / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * mobi * gamma0 << endl;
    cout << "unpper limit of driving force is: " << delta / dtime / mobi / PI << endl;
    cout << "The dimension is" << endl;
    cout << "***************   " << ND * dx << " m"
         << "    ***************" << endl;
    cout << "**                                          **" << endl;
    cout << "**********************************************" << endl;

    for (i = 0; i <= nm; i++)
    {
        for (j = 0; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            anij[i][j] = 0;
            thij[i][j] = 0.0;
            vpij[i][j] = 0.0;
            etaij[i][j] = 0.0;
            if ((i == 0) || (j == 0))
            {
                anij[i][j] = 1;
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                anij[i][j] = 0;
            }
        }
    }

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {

                if (i <= ND / 8 && j < ND / 2)
                // if ((i - ND / 2) * (i - ND / 2) + (j - ND / 2) * (j - ND / 2) + (k - ND / 2) * (k - ND / 2) < 49)
                {
                    phi[1][i][j][k] = 1.0;
                    con1[i][j][k] = c1e;
                    phi[2][i][j][k] = 0.0;
                    con2[i][j][k] = c2e;
                    phi[0][i][j][k] = 0.0;
                    con0[i][j][k] = c01e;
                }
                else if (i <= ND / 8 && j >= ND / 2)
                {
                    phi[2][i][j][k] = 1.0;
                    con2[i][j][k] = c2e;
                    phi[1][i][j][k] = 0.0;
                    con1[i][j][k] = c1e;
                    phi[0][i][j][k] = 0.0;
                    con0[i][j][k] = c02e;
                }
                else
                {
                    phi[1][i][j][k] = 0.0;
                    con1[i][j][k] = c1e;
                    phi[2][i][j][k] = 0.0;
                    con2[i][j][k] = c2e;
                    phi[0][i][j][k] = 1.0;
                    con0[i][j][k] = cl;
                }
                con[i][j][k] = phi[1][i][j][k] * con1[i][j][k] + phi[2][i][j][k] * con2[i][j][k] + phi[0][i][j][k] * con0[i][j][k];
                con_new[i][j][k] = con[i][j][k];
                sum1 += con[i][j][k];
            }
        }
    }
    c0 = sum1 / ND / ND / ND;
    cout << "nominal concentration is: " << c0 << endl;

start:;

    if ((((int)(istep) % pstep) == 0))
    {
        datasave(istep);
        cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;

        sum1 = 0.0;
        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                for (k = 0; k <= ndm; k++)
                {
                    sum1 += con[i][j][k];
                }
            }
        }
        cout << "nominal concentration is: " << sum1 / ND / ND / ND << endl;
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndm)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndm;
                }
                if (j == ndm)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndm;
                }
                if (k == ndm)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndm;
                }

                phinum = 0;
                for (ii = 0; ii <= nm; ii++)
                {
                    if ((phi[ii][i][j][k] > 0.0) ||
                        ((phi[ii][i][j][k] == 0.0) && (phi[ii][ip][j][k] > 0.0) ||
                         (phi[ii][im][j][k] > 0.0) ||
                         (phi[ii][i][jp][k] > 0.0) ||
                         (phi[ii][i][jm][k] > 0.0) ||
                         (phi[ii][i][j][kp] > 0.0) ||
                         (phi[ii][i][j][km] > 0.0)))
                    {
                        phinum++;
                        phiIdx[phinum][i][j][k] = ii;
                    }
                }
                phiNum[i][j][k] = phinum;
            }
        }
    }

    // Calculate the concentration field in solid and liqud phase based on local concentration and equilibrium rule.
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                con1[i][j][k] = c1e;
                con2[i][j][k] = c2e;
                // liquid in pure phase 1 or phase 2 or their mixer
                if (phi[0][i][j][k] == 0.0)
                {
                    con0[i][j][k] = c01e * phi[1][i][j][k] + c02e * phi[2][i][j][k];
                }
                // liquid in phase 2 or phase 1 or their mixer
                else if (phi[0][i][j][k] > 0.0 || ((phi[1][i][j][k] >= 0.0) || (phi[2][i][j][k] >= 0.0)))
                {
                    con0[i][j][k] = (con[i][j][k] - phi[2][i][j][k] * con2[i][j][k] - phi[1][i][j][k] * con1[i][j][k]) / phi[0][i][j][k];
                    if (con0[i][j][k] >= c01e)
                    {
                        con0[i][j][k] = c01e;
                    }
                    else if (con0[i][j][k] <= c02e)
                    {
                        con0[i][j][k] = c02e;
                    }
                }
                // The local concentation should be re-assigned after setting cons and conl (very important)
                con[i][j][k] = con1[i][j][k] * phi[1][i][j][k] + phi[2][i][j][k] * con2[i][j][k] + con0[i][j][k] * phi[0][i][j][k];
            }
        }
    }

    // Evolution Equation of Concentration field
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndm)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndm;
                }
                if (j == ndm)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndm;
                }
                if (k == ndm)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndm;
                }

                phis_ijk = phi[1][i][j][k] + phi[2][i][j][k];
                phis_ipjk = phi[1][ip][j][k] + phi[2][ip][j][k];
                phis_imjk = phi[1][im][j][k] + phi[2][im][j][k];
                phis_ijpk = phi[1][i][jp][k] + phi[2][i][jp][k];
                phis_ijmk = phi[1][i][jm][k] + phi[2][i][jm][k];
                phis_ijkp = phi[1][i][j][kp] + phi[2][i][j][kp];
                phis_ijkm = phi[1][i][j][km] + phi[2][i][j][km];

                cons_ijk = con1[i][j][k] + con2[i][j][k];
                cons_ipjk = con1[ip][j][k] + con2[ip][j][k];
                cons_imjk = con1[im][j][k] + con2[im][j][k];
                cons_ijpk = con1[i][jp][k] + con2[i][jp][k];
                cons_ijmk = con1[i][jm][k] + con2[i][jm][k];
                cons_ijkp = con1[i][j][kp] + con2[i][j][kp];
                cons_ijkm = con1[i][j][km] + con2[i][j][km];

                dev1_s = 0.25 * ((phis_ipjk - phis_imjk) * (cons_ipjk - cons_imjk) + (phis_ijpk - phis_ijmk) * (cons_ijpk - cons_ijmk) + (phis_ijkp - phis_ijkm) * (cons_ijkp - cons_ijkm)) / dx / dx;
                dev1_l = 0.25 * ((phi[0][ip][j][k] - phi[0][im][j][k]) * (con0[ip][j][k] - con0[im][j][k]) + (phi[0][i][jp][k] - phi[0][i][jm][k]) * (con0[i][jp][k] - con0[i][jm][k]) + (phi[0][i][j][kp] - phi[0][i][j][km]) * (con0[i][j][kp] - con0[i][j][km])) / dx / dx;
                dev2_s = (phis_ijk) * (cons_ipjk + cons_imjk + cons_ijpk + cons_ijmk + cons_ijkp + cons_ijkm - 6.0 * cons_ijk) / dx / dx;
                dev2_l = phi[0][i][j][k] * (con0[ip][j][k] + con0[im][j][k] + con0[i][jp][k] + con0[i][jm][k] + con0[i][j][kp] + con0[i][j][km] - 6.0 * con0[i][j][k]) / dx / dx;

                cddtt = Ds * (dev1_s + dev2_s) + Dl * (dev1_l + dev2_l);
                con_new[i][j][k] = con[i][j][k] + cddtt * dtime;
                // ch2[i][j] = ch[i][j] + cddtt * dtime + (2. * DRND(1.) - 1.) * 0.001;
            }
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                con[i][j][k] = con_new[i][j][k];
            }
        }
    }

    sum1 = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                sum1 += con[i][j][k];
            }
        }
    }
    dc0 = sum1 / ND / ND / ND - c0;

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                con[i][j][k] = con[i][j][k] - dc0;
                if (con[i][j][k] > 1.0)
                {
                    con[i][j][k] = 1.0;
                }
                if (con[i][j][k] < 0.0)
                {
                    con[i][j][k] = 0.0;
                }
            }
        }
    }

    // Evolution Equation of Phase Fields
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                ip = i + 1;
                im = i - 1;
                jp = j + 1;
                jm = j - 1;
                kp = k + 1;
                km = k - 1;
                if (i == ndm)
                {
                    ip = 0;
                }
                if (i == 0)
                {
                    im = ndm;
                }
                if (j == ndm)
                {
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndm;
                }
                if (k == ndm)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndm;
                }

                for (n1 = 1; n1 <= phiNum[i][j][k]; n1++)
                {
                    ii = phiIdx[n1][i][j][k];
                    pddtt = 0.0;
                    for (n2 = 1; n2 <= phiNum[i][j][k]; n2++)
                    {
                        jj = phiIdx[n2][i][j][k];
                        sum1 = 0.0;
                        for (n3 = 1; n3 <= phiNum[i][j][k]; n3++)
                        {
                            kk = phiIdx[n3][i][j][k];

                            // calculate the interface normal and deirivatives of the phase field
                            phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
                            phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
                            phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

                            phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                            phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                            phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

                            phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0 / dx / dx;
                            phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0 / dx / dx;
                            phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0 / dx / dx;

                            phiabs = phidx * phidx + phidy * phidy + phidz * phidz;

                            if (anij[ii][kk] == 1 && phiabs != 0.0)
                            {
                                epsilon0 = sqrt(aij[ii][kk]);

                                th = thij[ii][kk];
                                vp = vpij[ii][kk];
                                eta = etaij[ii][kk];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                termiikk = term0 + termx + termy + termz;
                            }
                            else
                            {
                                termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                            }

                            if (anij[jj][kk] == 1 && phiabs != 0.0)
                            {
                                epsilon0 = sqrt(aij[jj][kk]);

                                th = thij[jj][kk];
                                vp = vpij[jj][kk];
                                eta = etaij[jj][kk];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 2.0));

                                epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));
                                epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0));

                                term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs;
                                termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs;
                                termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs;

                                termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs, 2.0);
                                termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs, 2.0);
                                termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs, 2.0);

                                termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs, 2.0);
                                termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs, 2.0);
                                termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs, 2.0);

                                termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);
                                termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs, 3.0);

                                termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                termjjkk = term0 + termx + termy + termz;
                            }
                            else
                            {
                                termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                            }

                            sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];

                            // termiikk = aij[ii][kk] * (phi[kk][ip][j][k] + phi[kk][im][j][k] + phi[kk][i][jp][k] + phi[kk][i][jm][k] + phi[kk][i][j][kp] + phi[kk][i][j][km] - 6.0 * phi[kk][i][j][k]) / (dx * dx);

                            // termjjkk = aij[jj][kk] * (phi[kk][ip][j][k] + phi[kk][im][j][k] + phi[kk][i][jp][k] + phi[kk][i][jm][k] + phi[kk][i][j][kp] + phi[kk][i][j][km] - 6.0 * phi[kk][i][j][k]) / (dx * dx);

                            // sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
                        }
                        if (ii == 1 && jj == 0)
                        {
                            dF = F10 * (c01e - con0[i][j][k]);
                        }
                        else if (ii == 0 && jj == 1)
                        {
                            dF = -F10 * (c01e - con0[i][j][k]);
                        }
                        else if (ii == 2 && jj == 0)
                        {
                            dF = -F20 * (c02e - con0[i][j][k]);
                        }
                        else if (ii == 0 && jj == 2)
                        {
                            dF = F20 * (c02e - con0[i][j][k]);
                        }
                        else
                        {
                            dF = 0.0;
                        }
                        pddtt += -2.0 * mij[ii][jj] / double(phiNum[i][j][k]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
                    }
                    phi_new[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime;
                    if (phi_new[ii][i][j][k] >= 1.0)
                    {
                        phi_new[ii][i][j][k] = 1.0;
                    }
                    if (phi_new[ii][i][j][k] <= 0.0)
                    {
                        phi_new[ii][i][j][k] = 0.0;
                    }
                }
            } // i
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                for (l = 0; l <= nm; l++)
                {
                    phi[l][i][j][k] = phi_new[l][i][j][k];
                }
            }
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                sum1 = 0.0;
                for (l = 0; l <= nm; l++)
                {
                    sum1 += phi[l][i][j][k];
                }
                for (l = 0; l <= nm; l++)
                {
                    phi[l][i][j][k] = phi[l][i][j][k] / sum1;
                }
            }
        }
    }

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
    return 0;
}

void datasave(int step)
{
    FILE *stream;
    char buffer[30];
    sprintf(buffer, "data/phi/3d%d.vtk", step);
    stream = fopen(buffer, "a");

    fprintf(stream, "# vtk DataFile Version 1.0\n");
    fprintf(stream, "phi_%e.vtk\n", step);
    fprintf(stream, "ASCII\n");
    fprintf(stream, "DATASET STRUCTURED_POINTS\n");
    fprintf(stream, "DIMENSIONS %d %d %d\n", ND, ND, ND);
    fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
    fprintf(stream, "\n");
    fprintf(stream, "POINT_DATA %d\n", ND * ND * ND);
    fprintf(stream, "SCALARS scalars float\n");
    fprintf(stream, "LOOKUP_TABLE default\n");

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                fprintf(stream, "%e\n", phi[0][i][j][k]);
            }
        }
    }
    fclose(stream);

    FILE *streamc;
    char bufferc[30];
    sprintf(bufferc, "data/con/3d%d.vtk", step);
    streamc = fopen(bufferc, "a");

    fprintf(streamc, "# vtk DataFile Version 1.0\n");
    fprintf(streamc, "con_%e.vtk\n", step);
    fprintf(streamc, "ASCII\n");
    fprintf(streamc, "DATASET STRUCTURED_POINTS\n");
    fprintf(streamc, "DIMENSIONS %d %d %d\n", ND, ND, ND);
    fprintf(streamc, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(streamc, "ASPECT_RATIO 1.0 1.0 1.0\n");
    fprintf(streamc, "\n");
    fprintf(streamc, "POINT_DATA %d\n", ND * ND * ND);
    fprintf(streamc, "SCALARS scalars float\n");
    fprintf(streamc, "LOOKUP_TABLE default\n");

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (k = 0; k <= ndm; k++)
            {
                fprintf(streamc, "%e\n", con[i][j][k]);
            }
        }
    }
    fclose(streamc);
}