#include "../../../inc/_.hpp"

void rhs(int nx, double dx, VectorD const& u, VectorD& r) {
  for (int i = 1; i < nx; ++i) {
    r(i) = -u(i) * (u(i + 1) - u(i - 1)) / (2.0 * dx);
  }
}

void numerical(int nx, int ns, int nt, double dx, double dt, MatrixD& u,
               VectorD const& x) {
  auto un = MakeVectorD(nx + 1);
  auto ut = MakeVectorD(nx + 1);
  auto r = MakeVectorD(nx);

  int freq = nt / ns;

  un(0) = 0.0;
  for (int i = 1; i < nx; ++i) {
    un(i) = sin(2.0 * PI * x(i));
  }
  un(nx) = 0.0;

  ut(0) = 0.0;
  ut(nx) = 0.0;

  for (int j = 0; j < nt; ++j) {
    rhs(nx, dx, un, r);
    for (int i = 1; i < nx; ++i) {
      ut(i) = un(i) + dt * r(i);
    }
    rhs(nx, dx, ut, r);
    for (int i = 1; i < nx; ++i) {
      ut(i) = 0.75 * un(i) + 0.25 * ut(i) + 0.25 * dt * r(i);
    }
    rhs(nx, dx, ut, r);
    for (int i = 1; i < nx; ++i) {
      un(i) =
          (1.0 / 3.0) * un(i) + (2.0 / 3.0) * ut(i) + (2.0 / 3.0) * dt * r(i);
    }

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

  auto solution = fopen("solution.csv", "w");
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
