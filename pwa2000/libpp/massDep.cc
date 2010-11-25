#include "massDep.h"


using namespace std;
	

extern particleDataTable PDGtable;


complex<double> flat::val(const particle& p) 
{ return complex<double>(1,0); }


complex<double> breitWigner::val(const particle& p) 
{
  const double m0     = p.Mass ();
  const double Gamma0 = p.Width ();
  const double q      = p.q ();
  const double q0     = p.q0 ();   
  const double m      = ~(p.get4P ());
  const int    l      = p.Decay ()->L ();

  const double    GammaV = Gamma0 * (m0 / m) * (q / q0) * (pow(F(l, q), 2) / pow(F(l, q0), 2));
  complex<double> ret    = (m0 * Gamma0) / (m0 * m0 - m * m - complex<double>(0, 1) * m0 * GammaV);
  return ret;
}


// rho(1450)/rho(1700)
// see DOI: 10.1007/BF01552547
complex<double> rhoPrime::val(const particle& p) 
{
  
  const double m1     = 1.465;
  const double Gamma1 = 0.235;
  const double m2     = 1.720;
  const double Gamma2 = 0.220;

  const double q      = p.q ();
  const double q0     = p.q0 ();   
  const double m      = ~(p.get4P ());
  const int    l      = p.Decay ()->L ();

  const double    GammaV1 = Gamma1 * (m1 / m) * (q / q0) * (pow(F(l, q), 2) / pow(F(l, q0), 2));
  const double    GammaV2 = Gamma2 * (m2 / m) * (q / q0) * (pow(F(l, q), 2) / pow(F(l, q0), 2));
  
  cout << m << "     " << GammaV1 << "    " << GammaV2 << endl;

  complex<double> amp1    = (m1 * Gamma1) / (m1 * m1 - m * m - complex<double>(0, 1) * m1 * GammaV1);
  complex<double> amp2    = (m2 * Gamma2) / (m2 * m2 - m * m - complex<double>(0, 1) * m2 * GammaV2);
  

  complex<double> ret= 1./7 * ( 4. * amp1 - 3. * amp2 );
  return ret;

}


AMP_M::AMP_M() {
  ves_sheet = 0;
  _Pmax     = 2;
  _Nmax     = 5;

  _rho = matrix<complex<double> >(2, 2);
  _M   = matrix<complex<double> >(2, 2);
  _T   = matrix<complex<double> >(2, 2);

  _f = matrix<complex<double> >(1, 2);
  _f.el(0, 0) = 0.1968;
  _f.el(0, 1) = -0.0154;

  _a = vector<matrix<complex<double> > >(2, matrix<complex<double> >(2, 2));
  _a[0].el(0, 0) = 0.1131;
  _a[0].el(0, 1) = 0.0150;
  _a[0].el(1, 0) = 0.0150;
  _a[0].el(1, 1) = -0.3216;
  _a[1].el(0, 0) = _f.el(0, 0) * _f.el(0, 0);
  _a[1].el(0, 1) = _f.el(0, 0) * _f.el(0, 1);
  _a[1].el(1, 0) = _f.el(0, 1) * _f.el(0, 0);
  _a[1].el(1, 1) = _f.el(0, 1) * _f.el(0, 1);

	
  _c = vector<matrix<complex<double> > >(5, matrix<complex<double> >(2, 2));
  _c[0].el(0, 0) = 0.0337;
  _c[1].el(0, 0) = -0.3185;
  _c[2].el(0, 0) = -0.0942;
  _c[3].el(0, 0) = -0.5927;
  _c[4].el(0, 0) = 0.1957; 
  _c[0].el(0, 1) = _c[0].el(1, 0) = -0.2826;
  _c[1].el(0, 1) = _c[1].el(1, 0) = 0.0918;
  _c[2].el(0, 1) = _c[2].el(1, 0) = 0.1669;
  _c[3].el(0, 1) = _c[3].el(1, 0) = -0.2082;
  _c[4].el(0, 1) = _c[4].el(1, 0) = -0.1386;
  _c[0].el(1, 1) = 0.3010;
  _c[1].el(1, 1) = -0.5140;
  _c[2].el(1, 1) = 0.1176;
  _c[3].el(1, 1) = 0.5204;
  _c[4].el(1, 1) = -0.3977; 

  _sP = matrix<double>(1, 2);
  _sP.el(0, 0) = -0.0074;
  _sP.el(0, 1) = 0.9828;
}


complex<double> AMP_M::val(const particle& p) 
{
  extern particleDataTable PDGtable;

  _M.el(0, 0) = 0;
  _M.el(0, 1) = 0;
  _M.el(1, 0) = 0;
  _M.el(0, 1) = 0;

  complex <double> ret;
  complex <double> i(0, 1);

  double pi_mass    = PDGtable.get("pi").Mass();
  double pi0_mass   = PDGtable.get("pi0").Mass();
  double K_mass     = PDGtable.get("K").Mass();
  double K0_mass    = PDGtable.get("K0").Mass();
  double Kmean_mass = 0.5 * (K_mass + K0_mass);
  double My_mass    = ~(p.get4P());
  double s          = My_mass * My_mass;
  if (fabs(s - _sP.el(0, 1)) < 1e-6) {
    My_mass += 1e-6;
    s = My_mass * My_mass;
  }

  complex<double> q_pipi   = q(My_mass, pi_mass,    pi_mass);
  complex<double> q_pi0pi0 = q(My_mass, pi0_mass,   pi0_mass);
  complex<double> q_KK     = q(My_mass, K_mass,     K_mass);
  complex<double> q_K0K0   = q(My_mass, K0_mass,    K0_mass);
  complex<double> q_KmKm   = q(My_mass, Kmean_mass, Kmean_mass);
  _rho.el(0, 0) = 0.5 * ((2.0 * q_pipi) / My_mass + (2.0 * q_pi0pi0) / My_mass);
  _rho.el(1, 1) = 0.5 * ((2.0 * q_KK)   / My_mass + (2.0 * q_K0K0)   / My_mass);

  if (ves_sheet) {
    if (q_KmKm.imag() > 0.0)
      q_KmKm *= -1;
    _rho.el(0, 0) = (2.0 * q_pipi) / My_mass;
    _rho.el(1, 1) = (2.0 * q_KmKm) / My_mass;
  }
  //_rho.print();

  double scale = (s / (4.0 * Kmean_mass * Kmean_mass)) - 1;
  //cout << "scale=" << scale << endl;

  for (int p = 0; p < _Pmax; ++p) {
    complex<double> fa = 1. / (s - _sP.el(0, p));
    //cout << "fa" << p <<"="<< fa << endl;
    _M += fa * _a[p];
  }
  for (int n = 0; n < _Nmax; ++n) {
    complex<double> sc = pow(scale, n);
    //cout << "sc" << n <<"="<< sc << endl;
    _M += sc *_c[n];
  }
	
  // Modification: Here we set off-diagonal terms to 0
  _M.el(0, 1) = 0;
  _M.el(1, 0) = 0;
  //_M.print();

  _T = (_M - i * _rho).inv();
  //_T.print();

  ret = _T.el(0, 0);
  return ret;
}


complex<double> AMP_ves::val(const particle& p) {
  complex<double>        bw, amp_m;
  static complex<double> coupling(-0.3743, 0.3197);
  static double          massf0 = 0.9837, widthf0 = 0.0376;

  double pi_mass = PDGtable.get("pi").Mass();
  double My_mass = ~(p.get4P());

  ves_sheet = 1;
  amp_m     = AMP_M::val(p);

  if (My_mass > 2.0 * pi_mass) {
    double p, p0, gam, denom, A, B, C;
    p     = q(My_mass, pi_mass, pi_mass).real();
    p0    = q(massf0,  pi_mass, pi_mass).real();
    gam   = widthf0 * (p / p0);
    A     = massf0 * massf0 - My_mass * My_mass;
    B     = massf0 * gam;
    C     = B * (My_mass / p);
    denom = C / (A * A + B * B);
    bw    = denom * complex<double>(A, B);
  }

  return amp_m - coupling * bw;
}


AMP_kach::AMP_kach(): AMP_M()
{
  // Change parameters according to Kachaev's prescription
  _c[4].el(0, 0) = 0; // was 0.1957;
  _c[4].el(1, 1) = 0; // was -0.3977;

  _a[0].el(0, 1) = 0; // setting off-diagonal values to 0
  _a[0].el(1, 0) = 0; // 
 
  // _a[1] are the f's from the AMP paper! Setting to 0!
  _a[1].el(0, 0) = 0;
  _a[1].el(0, 1) = 0;
  _a[1].el(1, 0) = 0;
  _a[1].el(1, 1) = 0;
}
