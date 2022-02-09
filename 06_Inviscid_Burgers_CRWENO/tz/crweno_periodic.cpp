#include "../../inc/_.hpp"

void tdms(VectorD const& a, VectorD const& b, VectorD const& c,
          VectorD const& d, VectorD& u, int s, int e) {
  auto q = MakeVectorD(e);
  u(s) = d(s) / b(s);
  double bb = b(s);
  for (int i = s + 1; i <= e; ++i) {
    q(i) = c(i - 1) / bb;
    bb = b(i) - q(i) * a(i);
    u(i) = (d(i) - a(i) * u(i - 1)) / bb;
  }
  for (int i = e - 1; i >= s; --i) {
    u(i) = u(i) - q(i + 1) * u(i + 1);
  }
}

void ctdms(VectorD const& a, VectorD const& b, VectorD const& c, double alpha,
           double beta, VectorD const& r, VectorD& x, int s, int e) {
  auto bb = MakeVectorD(e);
  auto u = MakeVectorD(e);
  auto z = MakeVectorD(e);

  double gamma = -b(s);
  bb(s) = b(s) - gamma;
  bb(e) = b(e) - alpha * beta / gamma;

  for (int i = s + 1; i <= e - 1; ++i) {
    bb(i) = b(i);
  }

  tdms(a, bb, c, r, x, s, e);

  u(s) = gamma;
  u(e) = alpha;
  for (int i = s + 1; i <= e - 1; ++i) {
    u(i) = 0.0;
  }

  tdms(a, bb, c, u, z, s, e);

  double fact =
      (x(s) + beta * x(e) / gamma) / (1.0 + z(s) + beta * z(e) / gamma);

  for (int i = s; i <= e; ++i) {
    x(i) = x(i) - fact * z(i);
  }
}

void crwcL(double v1, double v2, double v3, double v4, double v5, double& a1,
           double& a2, double& a3, double& b1, double& b2, double& b3) {
  double eps = 1.0e-6;

  double s1 = (13.0 / 12.0) * SQR(v1 - 2.0 * v2 + v3) +
              0.25 * SQR(v1 - 4.0 * v2 + 3.0 * v3);
  double s2 = (13.0 / 12.0) * SQR(v2 - 2.0 * v3 + v4) + 0.25 * SQR(v2 - v4);
  double s3 = (13.0 / 12.0) * SQR(v3 - 2.0 * v4 + v5) +
              0.25 * SQR(3.0 * v3 - 4.0 * v4 + v5);

  double c1 = 2.0e-1 / (SQR(eps + s1));
  double c2 = 5.0e-1 / (SQR(eps + s2));
  double c3 = 3.0e-1 / (SQR(eps + s3));

  double w1 = c1 / (c1 + c2 + c3);
  double w2 = c2 / (c1 + c2 + c3);
  double w3 = c3 / (c1 + c2 + c3);

  a1 = (2.0 * w1 + w2) / 3.0;
  a2 = (w1 + 2.0 * w2 + 2.0 * w3) / 3.0;
  a3 = w3 / 3.0;

  b1 = w1 / 6.0;
  b2 = (5.0 * w1 + 5.0 * w2 + w3) / 6.0;
  b3 = (w2 + 5.0 * w3) / 6.0;
}

void crwcR(double v1, double v2, double v3, double v4, double v5, double& a1,
           double& a2, double& a3, double& b1, double& b2, double& b3) {
  double eps = 1.0e-6;

  double s1 = (13.0 / 12.0) * SQR(v1 - 2.0 * v2 + v3) +
              0.25 * SQR(v1 - 4.0 * v2 + 3.0 * v3);
  double s2 = (13.0 / 12.0) * SQR(v2 - 2.0 * v3 + v4) + 0.25 * SQR(v2 - v4);
  double s3 = (13.0 / 12.0) * SQR(v3 - 2.0 * v4 + v5) +
              0.25 * SQR(3.0 * v3 - 4.0 * v4 + v5);

  double c1 = 3.0e-1 / (SQR(eps + s1));
  double c2 = 5.0e-1 / (SQR(eps + s2));
  double c3 = 2.0e-1 / (SQR(eps + s3));

  double w1 = c1 / (c1 + c2 + c3);
  double w2 = c2 / (c1 + c2 + c3);
  double w3 = c3 / (c1 + c2 + c3);

  a1 = w1 / 3.0;
  a2 = (w3 + 2.0 * w2 + 2.0 * w1) / 3.0;
  a3 = (2.0 * w3 + w2) / 3.0;

  b1 = (w2 + 5.0 * w1) / 6.0;
  b2 = (5.0 * w3 + 5.0 * w2 + w1) / 6.0;
  b3 = w3 / 6.0;
}

void crwenoL(int n, VectorD const& u, VectorD& f) {
  auto a = MakeVectorD(n);
  auto b = MakeVectorD(n);
  auto c = MakeVectorD(n);
  auto r = MakeVectorD(n);

  double v1, v2, v3, v4, v5;
  double a1, a2, a3, b1, b2, b3;
  int i;

  i = 0;
  v1 = u(n - 2);
  v2 = u(n - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(i + 2);

  crwcL(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
  a(i) = a1;
  b(i) = a2;
  c(i) = a3;
  r(i) = b1 * u(n - 1) + b2 * u(i) + b3 * u(i + 1);

  i = 1;
  v1 = u(n - 1);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(i + 2);

  crwcL(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
  a(i) = a1;
  b(i) = a2;
  c(i) = a3;
  r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(i + 1);

  for (int i = 2; i < n - 1; ++i) {
    v1 = u(i - 2);
    v2 = u(i - 1);
    v3 = u(i);
    v4 = u(i + 1);
    v5 = u(i + 2);

    crwcL(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
    a(i) = a1;
    b(i) = a2;
    c(i) = a3;
    r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(i + 1);
  }

  i = n - 1;
  v1 = u(i - 2);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(1);

  crwcL(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
  a(i) = a1;
  b(i) = a2;
  c(i) = a3;
  r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(i + 1);

  double alpha = c(n - 1);
  double beta = a(0);

  ctdms(a, b, c, alpha, beta, r, f, 0, n - 1);
}

void crwenoR(int n, VectorD const& u, VectorD& f) {
  auto a = MakeVectorD(n + 1);
  auto b = MakeVectorD(n + 1);
  auto c = MakeVectorD(n + 1);
  auto r = MakeVectorD(n + 1);

  double v1, v2, v3, v4, v5;
  double a1, a2, a3, b1, b2, b3;
  int i;

  i = 1;
  v1 = u(n - 1);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(i + 2);

  crwcR(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
  a(i) = a1;
  b(i) = a2;
  c(i) = a3;
  r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(i + 1);

  for (int i = 2; i < n - 1; ++i) {
    v1 = u(i - 2);
    v2 = u(i - 1);
    v3 = u(i);
    v4 = u(i + 1);
    v5 = u(i + 2);

    crwcR(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
    a(i) = a1;
    b(i) = a2;
    c(i) = a3;
    r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(i + 1);
  }

  i = n - 1;
  v1 = u(i - 2);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(i + 1);
  v5 = u(1);

  crwcR(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
  a(i) = a1;
  b(i) = a2;
  c(i) = a3;
  r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(i + 1);

  i = n;
  v1 = u(i - 2);
  v2 = u(i - 1);
  v3 = u(i);
  v4 = u(1);
  v5 = u(2);

  crwcR(v1, v2, v3, v4, v5, a1, a2, a3, b1, b2, b3);
  a(i) = a1;
  b(i) = a2;
  c(i) = a3;
  r(i) = b1 * u(i - 1) + b2 * u(i) + b3 * u(1);

  double alpha = c(n);
  double beta = a(1);

  ctdms(a, b, c, alpha, beta, r, f, 1, n);
}

void rhs(int nx, double dx, VectorD const& u, VectorD& r) {
  auto uL = MakeVectorD(nx);
  auto uR = MakeVectorD(nx + 1);

  crwenoL(nx, u, uL);
  crwenoR(nx, u, uR);

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
      // printf("t = %g\n", (j + 1) * dt);
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
