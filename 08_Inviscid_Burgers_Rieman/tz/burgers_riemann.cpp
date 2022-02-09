#include "../../inc/_.hpp"

double wcL(double v1, double v2, double v3, double v4, double v5) {
  double eps = 1.0e-6;

  double s1 = (13.0 / 12.0) * SQR(v1 - 2.0 * v2 + v3) +
              0.25 * SQR(v1 - 4.0 * v2 + 3.0 * v3);
  double s2 = (13.0 / 12.0) * SQR(v2 - 2.0 * v3 + v4) + 0.25 * SQR(v2 - v4);
  double s3 = (13.0 / 12.0) * SQR(v3 - 2.0 * v4 + v5) +
              0.25 * SQR(3.0 * v3 - 4.0 * v4 + v5);

  double c1 = 1.0e-1 / (SQR(eps + s1));
  double c2 = 6.0e-1 / (SQR(eps + s2));
  double c3 = 3.0e-1 / (SQR(eps + s3));

  double w1 = c1 / (c1 + c2 + c3);
  double w2 = c2 / (c1 + c2 + c3);
  double w3 = c3 / (c1 + c2 + c3);

  double q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3;
  double q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0;
  double q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0;

  return (w1 * q1 + w2 * q2 + w3 * q3);
}

double wcR(double v1, double v2, double v3, double v4, double v5) {
  double eps = 1.0e-6;

  double s1 = (13.0 / 12.0) * SQR(v1 - 2.0 * v2 + v3) +
              0.25 * SQR(v1 - 4.0 * v2 + 3.0 * v3);
  double s2 = (13.0 / 12.0) * SQR(v2 - 2.0 * v3 + v4) + 0.25 * SQR(v2 - v4);
  double s3 = (13.0 / 12.0) * SQR(v3 - 2.0 * v4 + v5) +
              0.25 * SQR(3.0 * v3 - 4.0 * v4 + v5);

  double c1 = 3.0e-1 / (SQR(eps + s1));
  double c2 = 6.0e-1 / (SQR(eps + s2));
  double c3 = 1.0e-1 / (SQR(eps + s3));

  double w1 = c1 / (c1 + c2 + c3);
  double w2 = c2 / (c1 + c2 + c3);
  double w3 = c3 / (c1 + c2 + c3);

  double q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + v3 / 3.0;
  double q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0;
  double q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0;

  return (w1 * q1 + w2 * q2 + w3 * q3);
}

void wenoL(int n, VectorD const& u, VectorD& f) {
  for (int i = 0; i < n; ++i) {
    int ir = (i + 1) % n;
    int irr = (i + 2) % n;
    int il = (i - 1 + n) % n;
    int ill = (i - 2 + n) % n;

    double v1, v2, v3, v4, v5;
    v1 = u(ill);
    v2 = u(il);
    v3 = u(i);
    v4 = u(ir);
    v5 = u(irr);

    f(ir) = wcL(v1, v2, v3, v4, v5);
  }
}

void wenoR(int n, VectorD const& u, VectorD& f) {
  for (int i = 0; i < n; ++i) {
    int ir = (i + 1) % n;
    int irr = (i + 2) % n;
    int il = (i - 1 + n) % n;
    int ill = (i - 2 + n) % n;

    double v1, v2, v3, v4, v5;
    v1 = u(ill);
    v2 = u(il);
    v3 = u(i);
    v4 = u(ir);
    v5 = u(irr);

    f(i) = wcR(v1, v2, v3, v4, v5);
  }
}

void rusanov(int nx, VectorD const& u, VectorD const& uL, VectorD const& uR,
             VectorD& f, VectorD const& fL, VectorD const& fR) {
  auto ps = MakeVectorD(nx);
  ps(0) = max(abs(u(0)), abs(u(nx - 1)));
  for (int i = 1; i < nx; ++i) {
    ps(i) = max(abs(u(i)), abs(u(i - 1)));
  }

  for (int i = 0; i < nx; ++i) {
    f(i) = 0.5 * (fR(i) + fL(i)) - 0.5 * ps(i) * (uR(i) - uL(i));
  }
}

void fluxes(int nx, VectorD const& u, VectorD& f) {
  for (int i = 0; i < nx; ++i) {
    f(i) = 0.5 * u(i) * u(i);
  }
}

void rhs(int nx, double dx, VectorD const& u, VectorD& r) {
  auto f = MakeVectorD(nx);
  auto uL = MakeVectorD(nx);
  auto uR = MakeVectorD(nx);
  auto fL = MakeVectorD(nx);
  auto fR = MakeVectorD(nx);

  wenoL(nx, u, uL);
  wenoR(nx, u, uR);

  fluxes(nx, uL, fL);
  fluxes(nx, uR, fR);

  rusanov(nx, u, uL, uR, f, fL, fR);

  for (int i = 0; i < nx - 1; ++i) {
    r(i) = -(f(i + 1) - f(i)) / dx;
  }
  r(nx - 1) = -(f(0) - f(nx - 1)) / dx;
}

void numerical(int nx, int ns, int nt, double dx, double dt, MatrixD& u,
               VectorD const& x) {
  auto un = MakeVectorD(nx);
  auto ut = MakeVectorD(nx);
  auto r = MakeVectorD(nx);

  int freq = nt / ns;

  for (int i = 0; i < nx; ++i) {
    un(i) = sin(2.0 * PI * x(i));
    u(0, i) = un(i);
  }

  for (int j = 0; j < nt; ++j) {
    rhs(nx, dx, un, r);

    for (int i = 0; i < nx; ++i) {
      ut(i) = un(i) + dt * r(i);
    }

    rhs(nx, dx, ut, r);

    for (int i = 0; i < nx; ++i) {
      ut(i) = 0.75 * un(i) + 0.25 * ut(i) + 0.25 * dt * r(i);
    }

    rhs(nx, dx, ut, r);

    for (int i = 0; i < nx; ++i) {
      un(i) =
          (1.0 / 3.0) * un(i) + (2.0 / 3.0) * ut(i) + (2.0 / 3.0) * dt * r(i);
    }

    if ((j + 1) % freq == 0) {
      // printf("t = %g\n", (j + 1) * dt);
      int k = j / freq + 1;
      for (int i = 0; i < nx; ++i) {
        u(k, i) = un(i);
      }
    }
  }
}

int main() {
  int nx = 200;
  double dx = 1.0 / nx;

  double tm = 0.25;
  int nt = 2500;
  double dt = tm / nt;
  int ns = 10;

  auto u = MakeMatrixD(ns + 1, nx);

  auto x = LinspaceD(0.5 * dx, 1.0 - 0.5 * dx, nx);

  numerical(nx, ns, nt, dx, dt, u, x);

  auto solution = fopen("solution_riemann.csv", "w");
  for (int i = 0; i < nx; ++i) {
    fprintf(solution, "%g ", x(i));
  }
  fprintf(solution, "\n");

  for (int j = 0; j < ns + 1; ++j) {
    for (int i = 0; i < nx; ++i) {
      fprintf(solution, "%g ", u(j, i));
    }
    fprintf(solution, "\n");
  }

  fclose(solution);
}
