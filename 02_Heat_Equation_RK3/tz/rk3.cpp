#include "../../inc/_.hpp"

double compute_l2norm(VectorD const& r) {
  double rms = 0.0;
  for (int i = 1; i < r.size() - 1; ++i) {
    rms += r(i) * r(i);
  }
  return sqrt(rms / ((r.size() - 2)));
}

void rhs(int nx, double dx, VectorD const& u, VectorD& r, double alpha) {
  for (int i = 1; i < nx; ++i) {
    r(i) = alpha * (u(i + 1) - 2.0 * u(i) + u(i - 1)) / (dx * dx);
  }
}

void numerical(int nx, int nt, double dx, double dt, VectorD const& x,
               MatrixD& u, double alpha) {
  auto un = MakeVectorD(nx + 1);
  auto ut = MakeVectorD(nx + 1);
  auto r = MakeVectorD(nx);

  for (int i = 0; i < nx + 1; ++i) {
    un(i) = -sin(PI * x(i));
    u(0, i) = un(i);
  }

  un(0) = 0.0;
  un(nx) = 0.0;

  ut(0) = 0.0;
  ut(nx) = 0.0;

  for (int j = 1; j < nt + 1; ++j) {
    rhs(nx, dx, un, r, alpha);
    for (int i = 1; i < nx; ++i) {
      ut(i) = un(i) + dt * r(i);
    }
    rhs(nx, dx, ut, r, alpha);
    for (int i = 1; i < nx; ++i) {
      ut(i) = 0.75 * un(i) + 0.25 * ut(i) + 0.25 * dt * r(i);
    }
    rhs(nx, dx, ut, r, alpha);
    for (int i = 1; i < nx; ++i) {
      un(i) =
          (1.0 / 3.0) * un(i) + (2.0 / 3.0) * ut(i) + (2.0 / 3.0) * dt * r(i);
      u(j, i) = un(i);
    }
  }
}

int main() {
  double x_l = -1.0;
  double x_r = 1.0;
  double dx = 0.025;
  int nx = (x_r - x_l) / dx;

  double dt = 0.0025;
  double t = 1.0;
  int nt = t / dt;

  double alpha = 1 / (PI * PI);

  auto u = MakeMatrixD(nt + 1, nx + 1);
  auto x = MakeVectorD(nx + 1);
  auto u_e = MakeVectorD(nx + 1);
  auto error = MakeVectorD(nx + 1);

  for (int i = 0; i < nx + 1; ++i) {
    x(i) = x_l + dx * i;
    u_e(i) = -exp(-t) * sin(PI * x(i));
  }

  numerical(nx, nt, dx, dt, x, u, alpha);

  auto u_error = MakeVectorD(nx + 1);
  for (int i = 0; i < nx + 1; ++i) {
    u_error(i) = u(nt, i) - u_e(i);
  }
  double rms_error = compute_l2norm(u_error);
  double max_error = 0;
  for (auto x : u_error) {
    max_error = max(max_error, abs(x));
  }

  auto output = fopen("output.txt", "w");
  fprintf(output, "Error details: \n");
  fprintf(output, "L-2 Norm = %.20f\n", rms_error);
  fprintf(output, "Maximum Norm = %.20f\n", max_error);

  auto field_final = fopen("field_final.csv", "w");
  fprintf(field_final, "x ue un uerror \n");
  for (int i = 0; i < nx + 1; ++i) {
    fprintf(field_final, "%.16f %.16f %.16f %.16f \n", x(i), u_e(i), u(nt, i),
            u_error(i));
  }

  fclose(output);
  fclose(field_final);
}
