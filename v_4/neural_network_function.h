double* read_matrix(char file_name[], int N, int M);

void write_matrix(char file_name[], double* A, int N, int M);  

double* straight(double* A, double* B, double *b, int N, int M);

double* nnout(double* A, double* B, double *b, int N, int M);

double* delta_out(double out, double error);

double* back(double* W, double* delta, double *out, int N, int M);

void correction_w(double E, double A, double* w, double* dw, double* d, double* out, int N, int M);

void correction_b(double E, double A, double* b, double* db, double* d, int N);
