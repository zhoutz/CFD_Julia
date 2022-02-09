#include "../../../inc/_.hpp"

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
  double v1, v2, v3, v4, v5;
  int i;

  i = 0;
  v1 = u(n - 2);
  v2 = u(n - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(i + 2);
  f(i) = wcL(v1, v2, v3, v4, v5);

  i = 1;
  v1 = u(n - 1);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(i + 2);
  f(i) = wcL(v1, v2, v3, v4, v5);

  for (int i = 2; i < n - 1; ++i) {
    v1 = u(i - 2);
    v2 = u(i - 1);
    v3 = u(i);
    v4 = u(i + 1);
    v5 = u(i + 2);
    f(i) = wcL(v1, v2, v3, v4, v5);
  }

  i = n - 1;
  v1 = u(i - 2);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(1);
  f(i) = wcL(v1, v2, v3, v4, v5);
}

void wenoR(int n, VectorD const& u, VectorD& f) {
  double v1, v2, v3, v4, v5;
  int i;

  i = 1;
  v1 = u(n - 1);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(i + 2);
  f(i) = wcR(v1, v2, v3, v4, v5);

  for (int i = 2; i < n - 1; ++i) {
    v1 = u(i - 2);
    v2 = u(i - 1);
    v3 = u(i);
    v4 = u(i + 1);
    v5 = u(i + 2);
    f(i) = wcR(v1, v2, v3, v4, v5);
  }

  i = n - 1;
  v1 = u(i - 2);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(1);
  f(i) = wcR(v1, v2, v3, v4, v5);

  i = n;
  v1 = u(i - 2);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(1);
  v5 = u(2);
  f(i) = wcR(v1, v2, v3, v4, v5);
}

void rhs(int nx, double dx, VectorD const& u, VectorD& r) {
  auto uL = MakeVectorD(nx);
  auto uR = MakeVectorD(nx + 1);

  wenoL(nx, u, uL);
  wenoR(nx, u, uR);

  for (int i = 1; i < nx; ++i) {
    if (u(i) >= 0.0) {
      r(i) = -u(i) * (uL(i) - uL(i - 1)) / dx;
    } else {
      r(i) = -u(i) * (uR(i + 1) - uR(i)) / dx;
    }
  }

  int i = 0;
  if (u(i) >= 0.0) {
    r(i) = -u(i) * (uL(i) - uL(nx - 1)) / dx;
  } else {
    r(i) = -u(i) * (uR(i + 1) - uR(nx)) / dx;
  }
}

void numerical(int nx, int ns, int nt, double dx, double dt, MatrixD& u,
               VectorD const& x) {
  auto un = MakeVectorD(nx + 1);
  auto ut = MakeVectorD(nx + 1);
  auto r = MakeVectorD(nx);

  int freq = nt / ns;

  for (int i = 0; i < nx + 1; ++i) {
    un(i) = sin(2.0 * PI * x(i));
  }

  for (int j = 0; j < nt; ++j) {
    rhs(nx, dx, un, r);

    for (int i = 0; i < nx; ++i) {
      ut(i) = un(i) + dt * r(i);
    }
    ut(nx) = ut(0);

    rhs(nx, dx, ut, r);

    for (int i = 0; i < nx; ++i) {
      ut(i) = 0.75 * un(i) + 0.25 * ut(i) + 0.25 * dt * r(i);
    }
    ut(nx) = ut(0);

    rhs(nx, dx, ut, r);

    for (int i = 0; i < nx; ++i) {
      un(i) =
          (1.0 / 3.0) * un(i) + (2.0 / 3.0) * ut(i) + (2.0 / 3.0) * dt * r(i);
    }
    ut(nx) = ut(0);

    if ((j + 1) % freq == 0) {
      printf("t = %g\n", (j + 1) * dt);
      int k = j / freq;
      for (int i = 0; i < nx + 1; ++i) {
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

  auto u = MakeMatrixD(ns, nx + 1);

  auto x = MakeVectorD(nx + 1);
  for (int i = 0; i < nx + 1; ++i) {
    x(i) = i * dx;
  }

  numerical(nx, ns, nt, dx, dt, u, x);

  auto solution = fopen("solution_p.csv", "w");
  for (int i = 0; i < nx + 1; ++i) {
    fprintf(solution, "%g ", x(i));
  }
  fprintf(solution, "\n");

  for (int j = 0; j < ns; ++j) {
    for (int i = 0; i < nx + 1; ++i) {
      fprintf(solution, "%g ", u(j, i));
    }
    fprintf(solution, "\n");
  }

  fclose(solution);
}
