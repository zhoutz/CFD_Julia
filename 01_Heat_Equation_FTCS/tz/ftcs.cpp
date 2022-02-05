#include "../../inc/_.hpp"

double compute_l2norm(Tensor<double, 1> const& r) {
  double rms = 0.0;
  for (int i = 1; i < r.size() - 1; ++i) {
    rms += r(i) * r(i);
  }
  return sqrt(rms / ((r.size() - 2)));
}

double x_l = -1.0;
double x_r = 1.0;
double dx = 0.025;
int nx = (x_r - x_l) / dx;

double dt = 0.0025;
double t = 1.0;
int nt = t / dt;

double alpha = 1 / (PI * PI);

auto x = MakeTensor<double>(nx + 1);
auto u_e = MakeTensor<double>(nx + 1);
auto error = MakeTensor<double>(nx + 1);
auto un = MakeTensor<double>(nt + 1, nx + 1);

int main() {
  for (int i = 0; i < nx + 1; ++i) {
    x(i) = x_l + dx * i;
    un(0, i) = -sin(PI * x(i));
    u_e(i) = -exp(-t) * sin(PI * x(i));
  }
  un(0, 0) = 0.0;
  un(0, nx) = 0.0;

  double beta = alpha * dt / (dx * dx);

  for (int k = 1; k < nt + 1; ++k) {
    for (int i = 1; i < nx; ++i) {
      un(k, i) = un(k - 1, i) + beta * (un(k - 1, i + 1) - 2.0 * un(k - 1, i) +
                                        un(k - 1, i - 1));
    }
    un(k, 0) = 0.0;
    un(k, nx) = 0.0;
  }

  auto u_error = MakeTensor<double>(nx + 1);
  for (int i = 0; i < nx + 1; ++i) {
    u_error(i) = un(nt, i) - u_e(i);
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
    fprintf(field_final, "%.16f %.16f %.16f %.16f \n", x(i), u_e(i), un(nt, i),
            u_error(i));
  }

  fclose(output);
  fclose(field_final);
}
