void converter_U(double A, double *U, double angle);

double calculation_Id(double L, double R, double p,
                      double U, double *I, double w);


double calculation_Iq(double L, double R, double p, double F,
                      double U, double *I, double w);

double calculation_Mn(double Mn_0, double Mn_1, double Mn_2,
                      double w, double kr, double nr);

double PID(double *error, double dt);