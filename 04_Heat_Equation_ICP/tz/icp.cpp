#include "../../inc/_.hpp"

double compute_l2norm(VectorD const& r) {
  double rms = 0.0;
  for (int i = 1; i < r.size() - 1; ++i) {
    rms += r(i) * r(i);
  }
  return sqrt(rms / ((r.size() - 2)));
}

void tdms(VectorD const& a, VectorD const& b, VectorD const& c,
          VectorD const& d, VectorD& u, int N) {
  auto q = MakeVectorD(N);
  u(0) = d(0) / b(0);
  double bb = b(0);
  for (int i = 1; i < N; ++i) {
    q(i) = c(i - 1) / bb;
    bb = b(i) - q(i) * a(i);
    u(i) = (d(i) - a(i) * u(i - 1)) / bb;
  }
  for (int i = N - 2; i >= 0; --i) {
    u(i) = u(i) - q(i + 1) * u(i + 1);
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

  auto x = MakeVectorD(nx + 1);
  auto u_e = MakeVectorD(nx + 1);
  auto u_n = MakeMatrixD(nt + 1, nx + 1);
  auto a = MakeVectorD(nx + 1);
  auto b = MakeVectorD(nx + 1);
  auto c = MakeVectorD(nx + 1);
  auto r = MakeVectorD(nx + 1);
  auto p = MakeVectorD(nx + 1);

  for (int i = 0; i < nx + 1; ++i) {
    x(i) = x_l + dx * i;
    u_n(0, i) = -sin(PI * x(i));
    u_e(i) = -exp(-t) * sin(PI * x(i));
  }

  u_n(0, 0) = 0.0;
  u_n(0, nx) = 0.0;

  double alpha1 = alpha * dt / (2.0 * dx * dx);

  for (int k = 1; k < nt + 1; ++k) {
    a(0) = 0.0;
    b(0) = 1.0;
    c(0) = 0.0;
    r(0) = 0.0;
    for (int i = 1; i < nx; ++i) {
      a(i) = 12.0 / (dx * dx) - 2.0 / (alpha * dt);
      b(i) = -24.0 / (dx * dx) - 20.0 / (alpha * dt);
      c(i) = 12.0 / (dx * dx) - 2.0 / (alpha * dt);
      r(i) =
          -2.0 / (alpha * dt) *
              (u_n(k - 1, i + 1) + 10.0 * u_n(k - 1, i) + u_n(k - 1, i - 1)) -
          12.0 / (dx * dx) *
              (u_n(k - 1, i + 1) - 2.0 * u_n(k - 1, i) + u_n(k - 1, i - 1));
    }
    a(nx) = 0.0;
    b(nx) = 1.0;
    c(nx) = 0.0;
    r(nx) = 0.0;

    tdms(a, b, c, r, p, nx + 1);

    for (int i = 0; i < nx + 1; ++i) {
      u_n(k, i) = p(i);
    }
  }

  auto u_error = MakeVectorD(nx + 1);
  for (int i = 0; i < nx + 1; ++i) {
    u_error(i) = u_n(nt, i) - u_e(i);
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
    fprintf(field_final, "%.16f %.16f %.16f %.16f \n", x(i), u_e(i), u_n(nt, i),
            u_error(i));
  }

  fclose(output);
  fclose(field_final);
}
