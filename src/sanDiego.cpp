#include <cmath>
#include <vector>
#include <iostream>

namespace mech {

  int    mm = 3;
  int    kk = 9;
  int    ii = 24;
  int    ff = 2;
  double OneAtm      = 1.01325e5;
  double OneThird    = 1.0 / 3.0;
  double GasConstant = 8314.4621;

  std::vector<double> mw = { 2.015880e+00, 1.007940e+00, 3.199880e+01, 1.599940e+01, 1.700734e+01, 3.300674e+01, 3.401468e+01, 1.801528e+01, 2.801348e+01 };

  void getSpecificHeats_R(double& T, std::vector<double>& cp0_R) {

    double tt0 = T;
    double tt1 = T * tt0;
    double tt2 = T * tt1;
    double tt3 = T * tt2;
    double tt4 = 1.0 / T;
    double tt5 = tt4 / T;

    if(tt0 > 1.000000e+03) {
      cp0_R[0] = 3.337279e+00 - 4.940247e-05 * tt0 + 4.994568e-07 * tt1 - 1.795664e-10 * tt2 + 2.002554e-14 * tt3;
    } else {
      cp0_R[0] = 2.344331e+00 + 7.980521e-03 * tt0 - 1.947815e-05 * tt1 + 2.015721e-08 * tt2 - 7.376118e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[1] = 2.500000e+00 - 2.308430e-11 * tt0 + 1.615619e-14 * tt1 - 4.735152e-18 * tt2 + 4.981974e-22 * tt3;
    } else {
      cp0_R[1] = 2.500000e+00 + 7.053328e-13 * tt0 - 1.995920e-15 * tt1 + 2.300816e-18 * tt2 - 9.277323e-22 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[2] = 3.282538e+00 + 1.483088e-03 * tt0 - 7.579667e-07 * tt1 + 2.094706e-10 * tt2 - 2.167178e-14 * tt3;
    } else {
      cp0_R[2] = 3.782456e+00 - 2.996734e-03 * tt0 + 9.847302e-06 * tt1 - 9.681295e-09 * tt2 + 3.243728e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[3] = 2.569421e+00 - 8.597411e-05 * tt0 + 4.194846e-08 * tt1 - 1.001778e-11 * tt2 + 1.228337e-15 * tt3;
    } else {
      cp0_R[3] = 3.168267e+00 - 3.279319e-03 * tt0 + 6.643064e-06 * tt1 - 6.128066e-09 * tt2 + 2.112660e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[4] = 2.864729e+00 + 1.056504e-03 * tt0 - 2.590828e-07 * tt1 + 3.052187e-11 * tt2 - 1.331959e-15 * tt3;
    } else {
      cp0_R[4] = 4.125306e+00 - 3.225449e-03 * tt0 + 6.527647e-06 * tt1 - 5.798536e-09 * tt2 + 2.062374e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[5] = 4.017211e+00 + 2.239820e-03 * tt0 - 6.336581e-07 * tt1 + 1.142464e-10 * tt2 - 1.079085e-14 * tt3;
    } else {
      cp0_R[5] = 4.301798e+00 - 4.749121e-03 * tt0 + 2.115829e-05 * tt1 - 2.427639e-08 * tt2 + 9.292251e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[6] = 4.165003e+00 + 4.908317e-03 * tt0 - 1.901392e-06 * tt1 + 3.711860e-10 * tt2 - 2.879083e-14 * tt3;
    } else {
      cp0_R[6] = 4.276113e+00 - 5.428224e-04 * tt0 + 1.673357e-05 * tt1 - 2.157708e-08 * tt2 + 8.624544e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[7] = 3.033992e+00 + 2.176918e-03 * tt0 - 1.640725e-07 * tt1 - 9.704199e-11 * tt2 + 1.682010e-14 * tt3;
    } else {
      cp0_R[7] = 4.198641e+00 - 2.036434e-03 * tt0 + 6.520402e-06 * tt1 - 5.487971e-09 * tt2 + 1.771978e-12 * tt3;
    };

    if(tt0 > 1.000000e+03) {
      cp0_R[8] = 2.926640e+00 + 1.487977e-03 * tt0 - 5.684760e-07 * tt1 + 1.009704e-10 * tt2 - 6.753351e-15 * tt3;
    } else {
      cp0_R[8] = 3.298677e+00 + 1.408240e-03 * tt0 - 3.963222e-06 * tt1 + 5.641515e-09 * tt2 - 2.444854e-12 * tt3;
    };

  };

  void getEnthalpies_RT(double& T, std::vector<double>& h0_RT) {

    double tt0 = T;
    double tt1 = T * tt0;
    double tt2 = T * tt1;
    double tt3 = T * tt2;
    double tt4 = 1.0 / T;
    double tt5 = tt4 / T;
    double tt6 = std::log(tt0) * tt4;

    if(tt0 > 1.000000e+03) {
      h0_RT[0] = 3.337279e+00 - 4.940247e-05 * 0.50 * tt0 + 4.994568e-07 * OneThird * tt1 - 1.795664e-10 * 0.25 * tt2 + 2.002554e-14 * 0.20 * tt3 - 9.501589e+02 * tt4;
    } else {
      h0_RT[0] = 2.344331e+00 + 7.980521e-03 * 0.50 * tt0 - 1.947815e-05 * OneThird * tt1 + 2.015721e-08 * 0.25 * tt2 - 7.376118e-12 * 0.20 * tt3 - 9.179352e+02 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[1] = 2.500000e+00 - 2.308430e-11 * 0.50 * tt0 + 1.615619e-14 * OneThird * tt1 - 4.735152e-18 * 0.25 * tt2 + 4.981974e-22 * 0.20 * tt3 + 2.547366e+04 * tt4;
    } else {
      h0_RT[1] = 2.500000e+00 + 7.053328e-13 * 0.50 * tt0 - 1.995920e-15 * OneThird * tt1 + 2.300816e-18 * 0.25 * tt2 - 9.277323e-22 * 0.20 * tt3 + 2.547366e+04 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[2] = 3.282538e+00 + 1.483088e-03 * 0.50 * tt0 - 7.579667e-07 * OneThird * tt1 + 2.094706e-10 * 0.25 * tt2 - 2.167178e-14 * 0.20 * tt3 - 1.088458e+03 * tt4;
    } else {
      h0_RT[2] = 3.782456e+00 - 2.996734e-03 * 0.50 * tt0 + 9.847302e-06 * OneThird * tt1 - 9.681295e-09 * 0.25 * tt2 + 3.243728e-12 * 0.20 * tt3 - 1.063944e+03 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[3] = 2.569421e+00 - 8.597411e-05 * 0.50 * tt0 + 4.194846e-08 * OneThird * tt1 - 1.001778e-11 * 0.25 * tt2 + 1.228337e-15 * 0.20 * tt3 + 2.921758e+04 * tt4;
    } else {
      h0_RT[3] = 3.168267e+00 - 3.279319e-03 * 0.50 * tt0 + 6.643064e-06 * OneThird * tt1 - 6.128066e-09 * 0.25 * tt2 + 2.112660e-12 * 0.20 * tt3 + 2.912226e+04 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[4] = 2.864729e+00 + 1.056504e-03 * 0.50 * tt0 - 2.590828e-07 * OneThird * tt1 + 3.052187e-11 * 0.25 * tt2 - 1.331959e-15 * 0.20 * tt3 + 3.718858e+03 * tt4;
    } else {
      h0_RT[4] = 4.125306e+00 - 3.225449e-03 * 0.50 * tt0 + 6.527647e-06 * OneThird * tt1 - 5.798536e-09 * 0.25 * tt2 + 2.062374e-12 * 0.20 * tt3 + 3.381538e+03 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[5] = 4.017211e+00 + 2.239820e-03 * 0.50 * tt0 - 6.336581e-07 * OneThird * tt1 + 1.142464e-10 * 0.25 * tt2 - 1.079085e-14 * 0.20 * tt3 + 1.118567e+02 * tt4;
    } else {
      h0_RT[5] = 4.301798e+00 - 4.749121e-03 * 0.50 * tt0 + 2.115829e-05 * OneThird * tt1 - 2.427639e-08 * 0.25 * tt2 + 9.292251e-12 * 0.20 * tt3 + 2.948080e+02 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[6] = 4.165003e+00 + 4.908317e-03 * 0.50 * tt0 - 1.901392e-06 * OneThird * tt1 + 3.711860e-10 * 0.25 * tt2 - 2.879083e-14 * 0.20 * tt3 - 1.786179e+04 * tt4;
    } else {
      h0_RT[6] = 4.276113e+00 - 5.428224e-04 * 0.50 * tt0 + 1.673357e-05 * OneThird * tt1 - 2.157708e-08 * 0.25 * tt2 + 8.624544e-12 * 0.20 * tt3 - 1.770258e+04 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[7] = 3.033992e+00 + 2.176918e-03 * 0.50 * tt0 - 1.640725e-07 * OneThird * tt1 - 9.704199e-11 * 0.25 * tt2 + 1.682010e-14 * 0.20 * tt3 - 3.000430e+04 * tt4;
    } else {
      h0_RT[7] = 4.198641e+00 - 2.036434e-03 * 0.50 * tt0 + 6.520402e-06 * OneThird * tt1 - 5.487971e-09 * 0.25 * tt2 + 1.771978e-12 * 0.20 * tt3 - 3.029373e+04 * tt4;
    };

    if(tt0 > 1.000000e+03) {
      h0_RT[8] = 2.926640e+00 + 1.487977e-03 * 0.50 * tt0 - 5.684760e-07 * OneThird * tt1 + 1.009704e-10 * 0.25 * tt2 - 6.753351e-15 * 0.20 * tt3 - 9.227977e+02 * tt4;
    } else {
      h0_RT[8] = 3.298677e+00 + 1.408240e-03 * 0.50 * tt0 - 3.963222e-06 * OneThird * tt1 + 5.641515e-09 * 0.25 * tt2 - 2.444854e-12 * 0.20 * tt3 - 1.020900e+03 * tt4;
    };

  };

  void getEntropies_R(double& T, std::vector<double>& s0_R) {

    double tt0 = T;
    double tt1 = T * tt0;
    double tt2 = T * tt1;
    double tt3 = T * tt2;
    double tt4 = 1.0 / T;
    double tt5 = tt4 / T;
    double tt6 = std::log(T);

    if(tt0 > 1.000000e+03) {
      s0_R[0] = 3.337279e+00 * tt6 - 4.940247e-05 * tt0 + 4.994568e-07 * 0.50 * tt1 - 1.795664e-10 * OneThird * tt2 + 2.002554e-14 * 0.25 * tt3 - 3.205023e+00;
    } else {
      s0_R[0] = 2.344331e+00 * tt6 + 7.980521e-03 * tt0 - 1.947815e-05 * 0.50 * tt1 + 2.015721e-08 * OneThird * tt2 - 7.376118e-12 * 0.25 * tt3 + 6.830102e-01;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[1] = 2.500000e+00 * tt6 - 2.308430e-11 * tt0 + 1.615619e-14 * 0.50 * tt1 - 4.735152e-18 * OneThird * tt2 + 4.981974e-22 * 0.25 * tt3 - 4.466829e-01;
    } else {
      s0_R[1] = 2.500000e+00 * tt6 + 7.053328e-13 * tt0 - 1.995920e-15 * 0.50 * tt1 + 2.300816e-18 * OneThird * tt2 - 9.277323e-22 * 0.25 * tt3 - 4.466829e-01;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[2] = 3.282538e+00 * tt6 + 1.483088e-03 * tt0 - 7.579667e-07 * 0.50 * tt1 + 2.094706e-10 * OneThird * tt2 - 2.167178e-14 * 0.25 * tt3 + 5.453231e+00;
    } else {
      s0_R[2] = 3.782456e+00 * tt6 - 2.996734e-03 * tt0 + 9.847302e-06 * 0.50 * tt1 - 9.681295e-09 * OneThird * tt2 + 3.243728e-12 * 0.25 * tt3 + 3.657676e+00;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[3] = 2.569421e+00 * tt6 - 8.597411e-05 * tt0 + 4.194846e-08 * 0.50 * tt1 - 1.001778e-11 * OneThird * tt2 + 1.228337e-15 * 0.25 * tt3 + 4.784339e+00;
    } else {
      s0_R[3] = 3.168267e+00 * tt6 - 3.279319e-03 * tt0 + 6.643064e-06 * 0.50 * tt1 - 6.128066e-09 * OneThird * tt2 + 2.112660e-12 * 0.25 * tt3 + 2.051933e+00;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[4] = 2.864729e+00 * tt6 + 1.056504e-03 * tt0 - 2.590828e-07 * 0.50 * tt1 + 3.052187e-11 * OneThird * tt2 - 1.331959e-15 * 0.25 * tt3 + 5.701641e+00;
    } else {
      s0_R[4] = 4.125306e+00 * tt6 - 3.225449e-03 * tt0 + 6.527647e-06 * 0.50 * tt1 - 5.798536e-09 * OneThird * tt2 + 2.062374e-12 * 0.25 * tt3 - 6.904330e-01;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[5] = 4.017211e+00 * tt6 + 2.239820e-03 * tt0 - 6.336581e-07 * 0.50 * tt1 + 1.142464e-10 * OneThird * tt2 - 1.079085e-14 * 0.25 * tt3 + 3.785102e+00;
    } else {
      s0_R[5] = 4.301798e+00 * tt6 - 4.749121e-03 * tt0 + 2.115829e-05 * 0.50 * tt1 - 2.427639e-08 * OneThird * tt2 + 9.292251e-12 * 0.25 * tt3 + 3.716662e+00;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[6] = 4.165003e+00 * tt6 + 4.908317e-03 * tt0 - 1.901392e-06 * 0.50 * tt1 + 3.711860e-10 * OneThird * tt2 - 2.879083e-14 * 0.25 * tt3 + 2.916157e+00;
    } else {
      s0_R[6] = 4.276113e+00 * tt6 - 5.428224e-04 * tt0 + 1.673357e-05 * 0.50 * tt1 - 2.157708e-08 * OneThird * tt2 + 8.624544e-12 * 0.25 * tt3 + 3.435051e+00;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[7] = 3.033992e+00 * tt6 + 2.176918e-03 * tt0 - 1.640725e-07 * 0.50 * tt1 - 9.704199e-11 * OneThird * tt2 + 1.682010e-14 * 0.25 * tt3 + 4.966770e+00;
    } else {
      s0_R[7] = 4.198641e+00 * tt6 - 2.036434e-03 * tt0 + 6.520402e-06 * 0.50 * tt1 - 5.487971e-09 * OneThird * tt2 + 1.771978e-12 * 0.25 * tt3 - 8.490322e-01;
    };

    if(tt0 > 1.000000e+03) {
      s0_R[8] = 2.926640e+00 * tt6 + 1.487977e-03 * tt0 - 5.684760e-07 * 0.50 * tt1 + 1.009704e-10 * OneThird * tt2 - 6.753351e-15 * 0.25 * tt3 + 5.980528e+00;
    } else {
      s0_R[8] = 3.298677e+00 * tt6 + 1.408240e-03 * tt0 - 3.963222e-06 * 0.50 * tt1 + 5.641515e-09 * OneThird * tt2 - 2.444854e-12 * 0.25 * tt3 + 3.950372e+00;
    };

  };

  void getGibbsFunctions_RT(double& T, std::vector<double>& g0_RT) {

    std::vector<double> h0_RT(kk, 0.0);
    std::vector<double> s0_R(kk, 0.0);

    getEnthalpies_RT(T, h0_RT);
    getEntropies_R(T, s0_R);
    for(int k = 0; k < kk; ++k) { g0_RT[k] = h0_RT[k] - s0_R[k]; }

  };

  void getEquilibriumConstants(double& T, std::vector<double>& keq) {

    double            p0 = OneAtm;
    double              RT = GasConstant * T;
    double              C0 = p0 / RT;
    std::vector<double> g0_RT(kk, 0.0);

    getGibbsFunctions_RT(T, g0_RT);
    for(int k = 0; k < kk; ++k) { g0_RT[k] = exp(g0_RT[k]); }

    keq[0] = ( g0_RT[3] * g0_RT[4] ) / ( g0_RT[1] * g0_RT[2] );
    keq[1] = ( g0_RT[1] * g0_RT[4] ) / ( g0_RT[0] * g0_RT[3] );
    keq[2] = ( g0_RT[1] * g0_RT[7] ) / ( g0_RT[0] * g0_RT[4] );
    keq[3] = ( g0_RT[4] * g0_RT[4] ) / ( g0_RT[7] * g0_RT[3] );
    keq[4] =  C0 * ( g0_RT[0] ) / ( g0_RT[1] * g0_RT[1] );
    keq[5] =  C0 * ( g0_RT[7] ) / ( g0_RT[1] * g0_RT[4] );
    keq[6] =  C0 * ( g0_RT[2] ) / ( g0_RT[3] * g0_RT[3] );
    keq[7] =  C0 * ( g0_RT[4] ) / ( g0_RT[1] * g0_RT[3] );
    keq[8] =  C0 * ( g0_RT[5] ) / ( g0_RT[3] * g0_RT[4] );
    keq[9] =  C0 * ( g0_RT[5] ) / ( g0_RT[1] * g0_RT[2] );
    keq[10] = ( g0_RT[4] * g0_RT[4] ) / ( g0_RT[1] * g0_RT[5] );
    keq[11] = ( g0_RT[0] * g0_RT[2] ) / ( g0_RT[1] * g0_RT[5] );
    keq[12] = ( g0_RT[7] * g0_RT[3] ) / ( g0_RT[1] * g0_RT[5] );
    keq[13] = ( g0_RT[2] * g0_RT[4] ) / ( g0_RT[5] * g0_RT[3] );
    keq[14] = ( g0_RT[7] * g0_RT[2] ) / ( g0_RT[5] * g0_RT[4] );
    keq[15] = ( g0_RT[7] * g0_RT[2] ) / ( g0_RT[5] * g0_RT[4] );
    keq[16] =  C0 * ( g0_RT[6] ) / ( g0_RT[4] * g0_RT[4] );
    keq[17] = ( g0_RT[6] * g0_RT[2] ) / ( g0_RT[5] * g0_RT[5] );
    keq[18] = ( g0_RT[6] * g0_RT[2] ) / ( g0_RT[5] * g0_RT[5] );
    keq[19] = ( g0_RT[0] * g0_RT[5] ) / ( g0_RT[1] * g0_RT[6] );
    keq[20] = ( g0_RT[7] * g0_RT[4] ) / ( g0_RT[1] * g0_RT[6] );
    keq[21] = ( g0_RT[7] * g0_RT[5] ) / ( g0_RT[6] * g0_RT[4] );
    keq[22] = ( g0_RT[7] * g0_RT[5] ) / ( g0_RT[6] * g0_RT[4] );
    keq[23] = ( g0_RT[5] * g0_RT[4] ) / ( g0_RT[6] * g0_RT[3] );

  };

  void getTemperature(double& h, double& Told, std::vector<double>& y, double& T) {

    double            tol   = 1.0e-06;
    int               niter = 500;
    double              RT;
    double              To;
    double              Tp;
    double              dT = 1.0;
    double              fy = h;
    double              jy = 0.0;
    std::vector<double> hk(kk, 0.0);
    std::vector<double> cpk(kk, 0.0);

    To = Told;
    Tp = Told;

    for(int it = 0; it < niter; ++it) { 

      RT = GasConstant * To;
      getSpecificHeats_R(To, cpk);
      getEnthalpies_RT(To, hk);
      for(int k = 0; k < kk; ++k) { hk[k]  = RT * hk[k] / mw[k]; }
      for(int k = 0; k < kk; ++k) { cpk[k] = GasConstant * cpk[k] / mw[k]; }
      for(int k = 0; k < kk; ++k) { fy -=  hk[k] * y[k]; }
      for(int k = 0; k < kk; ++k) { jy -= cpk[k] * y[k]; }
      dT  = -fy / jy;
      Tp  =  To + dT;
      To  =  Tp;

      if( (std::fabs(dT) < tol)) {
	T = Tp;
	return;
      }

      fy = h;
      jy = 0.0;

    }

    T = Tp;

  };

  void getFalloffRates(double& T, std::vector<double>& C, std::vector<double>& kfwd) {

    int TROE = 110;
    double lpr;
    double cc;
    double nn;
    double f1;
    double logT = std::log(T);
    double invT = 1.0/T;
    std::vector<double> khi(ff,0.0);
    std::vector<double> klo(ff,0.0);
    std::vector<double> pr(ff,0.0);
    std::vector<double> work(ff,0.0);
    std::vector<int> falloffType(ff,100);

    khi[0] = std::exp(2.226013e+01 + 4.400000e-01 * logT);
    klo[0] = std::exp(3.168281e+01 - 1.400000e+00 * logT);

    khi[1] = std::exp(2.528239e+01 - 2.700000e-01 * logT);
    klo[1] = std::exp(4.476435e+01 - 3.200000e+00 * logT);

    pr[0] = 2.50e+00 * C[0] + 1.00e+00 * C[1] + 1.00e+00 * C[2] + 1.00e+00 * C[3] + 1.00e+00 * C[4] + 1.00e+00 * C[5] + 1.00e+00 * C[6] + 1.60e+01 * C[7] + 1.00e+00 * C[8]; 
    pr[1] = 2.50e+00 * C[0] + 1.00e+00 * C[1] + 1.00e+00 * C[2] + 1.00e+00 * C[3] + 1.00e+00 * C[4] + 1.00e+00 * C[5] + 1.00e+00 * C[6] + 6.00e+00 * C[7] + 1.00e+00 * C[8]; 

    for(int i = 0; i < ff; ++i) { pr[i] *= (klo[i] / khi[i]); }

    falloffType[0] = 110;
    falloffType[1] = 110;

    work[0] = (1.0 - 5.000000e-01) * std::exp( - 1.000000e+30 * T) + 5.000000e-01 * std::exp( - 1.000000e-30 * T);
    work[1] = (1.0 - 4.300000e-01) * std::exp( - 1.000000e-30 * T) + 4.300000e-01 * std::exp( - 1.000000e+30 * T);

    for(int i = 0; i < ff; ++i) {
      lpr =  std::log10(pr[i]);
      if(falloffType[i] == TROE) {
        cc      = -0.40 - 0.67 * std::log10(work[i]);
        nn      =  0.75 - 1.27 * std::log10(work[i]);
        f1      =  (lpr + cc)/(nn - 0.14 * (lpr + cc));
        work[i] =  std::log10(work[i])/(1 + f1 * f1);
        work[i] =  std::pow(10.0, work[i]);
      }
      work[i] =  (pr[i] * work[i])/(1 + pr[i]);
    }

    kfwd[9] = khi[0] * work[0];
    kfwd[16] = khi[1] * work[1];

  };

  void getRateCoefficients(double& T, std::vector<double>& C, std::vector<double>& kfwd, std::vector<double>& krev) {

    double logT = std::log(T);
    double invT = 1.0/T;
    std::vector<double> keq(ii,0.0);

    getEquilibriumConstants(T, keq);

    kfwd[0] = std::exp(3.119207e+01 - 7.000000e-01 * logT - 8.589852e+03 * invT);
    kfwd[1] = std::exp(3.923952e+00 + 2.670000e+00 * logT - 3.165569e+03 * invT);
    kfwd[2] = std::exp(1.397251e+01 + 1.300000e+00 * logT - 1.829343e+03 * invT);
    kfwd[3] = std::exp(6.551080e+00 + 2.330000e+00 * logT - 7.320979e+03 * invT);
    kfwd[4] = std::exp(2.789339e+01 - 1.000000e+00 * logT);
    kfwd[5] = std::exp(3.822766e+01 - 2.000000e+00 * logT);
    kfwd[6] = std::exp(2.254296e+01 - 5.000000e-01 * logT);
    kfwd[7] = std::exp(2.918071e+01 - 1.000000e+00 * logT);
    kfwd[8] = std::exp(2.280271e+01);
    kfwd[10] = std::exp(2.498312e+01 - 1.484161e+02 * invT);
    kfwd[11] = std::exp(2.353267e+01 - 4.140977e+02 * invT);
    kfwd[12] = std::exp(2.415725e+01 - 8.659610e+02 * invT);
    kfwd[13] = std::exp(2.371900e+01);
    kfwd[14] = std::exp(2.683251e+01 - 5.500055e+03 * invT);
    kfwd[15] = std::exp(2.411777e+01 + 2.501665e+02 * invT);
    kfwd[17] = std::exp(1.908337e+01 + 7.090056e+02 * invT);
    kfwd[18] = std::exp(2.535799e+01 - 5.556583e+03 * invT);
    kfwd[19] = std::exp(2.385876e+01 - 4.000620e+03 * invT);
    kfwd[20] = std::exp(2.302585e+01 - 1.804085e+03 * invT);
    kfwd[21] = std::exp(2.505268e+01 - 3.659888e+03 * invT);
    kfwd[22] = std::exp(2.127715e+01 - 1.599622e+02 * invT);
    kfwd[23] = std::exp(9.172639e+00 + 2.000000e+00 * logT - 2.008548e+03 * invT);

    kfwd[4] *= ( 2.500e+00 * C[0] + 1.000e+00 * C[1] + 1.000e+00 * C[2] + 1.000e+00 * C[3] + 1.000e+00 * C[4] + 1.000e+00 * C[5] + 1.000e+00 * C[6] + 1.200e+01 * C[7] + 1.000e+00 * C[8] ); 
    kfwd[5] *= ( 2.500e+00 * C[0] + 1.000e+00 * C[1] + 1.000e+00 * C[2] + 1.000e+00 * C[3] + 1.000e+00 * C[4] + 1.000e+00 * C[5] + 1.000e+00 * C[6] + 1.200e+01 * C[7] + 1.000e+00 * C[8] ); 
    kfwd[6] *= ( 2.500e+00 * C[0] + 1.000e+00 * C[1] + 1.000e+00 * C[2] + 1.000e+00 * C[3] + 1.000e+00 * C[4] + 1.000e+00 * C[5] + 1.000e+00 * C[6] + 1.200e+01 * C[7] + 1.000e+00 * C[8] ); 
    kfwd[7] *= ( 2.500e+00 * C[0] + 1.000e+00 * C[1] + 1.000e+00 * C[2] + 1.000e+00 * C[3] + 1.000e+00 * C[4] + 1.000e+00 * C[5] + 1.000e+00 * C[6] + 1.200e+01 * C[7] + 1.000e+00 * C[8] ); 
    kfwd[8] *= ( 2.500e+00 * C[0] + 1.000e+00 * C[1] + 1.000e+00 * C[2] + 1.000e+00 * C[3] + 1.000e+00 * C[4] + 1.000e+00 * C[5] + 1.000e+00 * C[6] + 1.200e+01 * C[7] + 1.000e+00 * C[8] ); 

    getFalloffRates(T, C, kfwd);

    for(int i = 0; i < ii; ++i) { krev[i] = kfwd[i] * keq[i]; }

  };

  void getNetRatesOfProgress(double& T, std::vector<double>& C, std::vector<double>& Rnet) {

    std::vector<double> kfwd(ii,0.0);
    std::vector<double> krev(ii,0.0);
    std::vector<double> Rfwd(ii,0.0);
    std::vector<double> Rrev(ii,0.0);

    getRateCoefficients(T, C, kfwd, krev);

    Rfwd[0] = kfwd[0] * C[1] * C[2];
    Rrev[0] = krev[0] * C[3] * C[4];

    Rfwd[1] = kfwd[1] * C[0] * C[3];
    Rrev[1] = krev[1] * C[1] * C[4];

    Rfwd[2] = kfwd[2] * C[0] * C[4];
    Rrev[2] = krev[2] * C[1] * C[7];

    Rfwd[3] = kfwd[3] * C[7] * C[3];
    Rrev[3] = krev[3] * C[4] * C[4];

    Rfwd[4] = kfwd[4] * C[1] * C[1];
    Rrev[4] = krev[4] * C[0];

    Rfwd[5] = kfwd[5] * C[1] * C[4];
    Rrev[5] = krev[5] * C[7];

    Rfwd[6] = kfwd[6] * C[3] * C[3];
    Rrev[6] = krev[6] * C[2];

    Rfwd[7] = kfwd[7] * C[1] * C[3];
    Rrev[7] = krev[7] * C[4];

    Rfwd[8] = kfwd[8] * C[3] * C[4];
    Rrev[8] = krev[8] * C[5];

    Rfwd[9] = kfwd[9] * C[1] * C[2];
    Rrev[9] = krev[9] * C[5];

    Rfwd[10] = kfwd[10] * C[1] * C[5];
    Rrev[10] = krev[10] * C[4] * C[4];

    Rfwd[11] = kfwd[11] * C[1] * C[5];
    Rrev[11] = krev[11] * C[0] * C[2];

    Rfwd[12] = kfwd[12] * C[1] * C[5];
    Rrev[12] = krev[12] * C[7] * C[3];

    Rfwd[13] = kfwd[13] * C[5] * C[3];
    Rrev[13] = krev[13] * C[2] * C[4];

    Rfwd[14] = kfwd[14] * C[5] * C[4];
    Rrev[14] = krev[14] * C[7] * C[2];

    Rfwd[15] = kfwd[15] * C[5] * C[4];
    Rrev[15] = krev[15] * C[7] * C[2];

    Rfwd[16] = kfwd[16] * C[4] * C[4];
    Rrev[16] = krev[16] * C[6];

    Rfwd[17] = kfwd[17] * C[5] * C[5];
    Rrev[17] = krev[17] * C[6] * C[2];

    Rfwd[18] = kfwd[18] * C[5] * C[5];
    Rrev[18] = krev[18] * C[6] * C[2];

    Rfwd[19] = kfwd[19] * C[1] * C[6];
    Rrev[19] = krev[19] * C[0] * C[5];

    Rfwd[20] = kfwd[20] * C[1] * C[6];
    Rrev[20] = krev[20] * C[7] * C[4];

    Rfwd[21] = kfwd[21] * C[6] * C[4];
    Rrev[21] = krev[21] * C[7] * C[5];

    Rfwd[22] = kfwd[22] * C[6] * C[4];
    Rrev[22] = krev[22] * C[7] * C[5];

    Rfwd[23] = kfwd[23] * C[6] * C[3];
    Rrev[23] = krev[23] * C[5] * C[4];

    for(int i = 0; i < ii; ++i) { Rnet[i] = Rfwd[i] - Rrev[i]; }

  };

  void getNetProductionRates(double& p, double& T, std::vector<double>& y, std::vector<double>& omega) {

    double W;
    double rho;
    std::vector<double> C(kk,0.0);
    std::vector<double> Rnet(ii,0.0);

    W   = 0.0;
    for(int k = 0; k < kk; ++k) { W += y[k] / mw[k]; }
    W   = 1.0 / W;
    rho = (p * W) / (GasConstant * T);
    for(int k = 0; k < kk; ++k) { C[k] = rho * y[k] / mw[k]; }

    getNetRatesOfProgress(T, C, Rnet);

    omega[0] =  - Rnet[1] - Rnet[2] + Rnet[4] + Rnet[11] + Rnet[19];
    omega[1] =  - Rnet[0] + Rnet[1] + Rnet[2] - Rnet[4] - Rnet[4] - Rnet[5] - Rnet[7] - Rnet[9] - Rnet[10] - Rnet[11] - Rnet[12] - Rnet[19] - Rnet[20];
    omega[2] =  - Rnet[0] + Rnet[6] - Rnet[9] + Rnet[11] + Rnet[13] + Rnet[14] + Rnet[15] + Rnet[17] + Rnet[18];
    omega[3] =  + Rnet[0] - Rnet[1] - Rnet[3] - Rnet[6] - Rnet[6] - Rnet[7] - Rnet[8] + Rnet[12] - Rnet[13] - Rnet[23];
    omega[4] =  + Rnet[0] + Rnet[1] - Rnet[2] + Rnet[3] + Rnet[3] - Rnet[5] + Rnet[7] - Rnet[8] + Rnet[10] + Rnet[10] + Rnet[13] - Rnet[14] - Rnet[15] - Rnet[16] - Rnet[16] + Rnet[20] - Rnet[21] - Rnet[22] + Rnet[23];
    omega[5] =  + Rnet[8] + Rnet[9] - Rnet[10] - Rnet[11] - Rnet[12] - Rnet[13] - Rnet[14] - Rnet[15] - Rnet[17] - Rnet[17] - Rnet[18] - Rnet[18] + Rnet[19] + Rnet[21] + Rnet[22] + Rnet[23];
    omega[6] =  + Rnet[16] + Rnet[17] + Rnet[18] - Rnet[19] - Rnet[20] - Rnet[21] - Rnet[22] - Rnet[23];
    omega[7] =  + Rnet[2] - Rnet[3] + Rnet[5] + Rnet[12] + Rnet[14] + Rnet[15] + Rnet[20] + Rnet[21] + Rnet[22];

  };

}
