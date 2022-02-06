#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <memory>
#include <numeric>
#include <vector>

inline namespace tz {
inline namespace ts {

template <class T, int N>
struct Tensor {
  std::array<int, N> shape;
  std::array<int, N> strides;
  std::unique_ptr<T[]> data;

  int size() const {
    return std::accumulate(shape.begin(), shape.end(), 1,
                           std::multiplies<int>());
  }

  void UpdateStrides() {
    std::exclusive_scan(shape.rbegin(), shape.rend(), strides.rbegin(), 1,
                        std::multiplies<int>());
  }

  template <class... Ints>
  Tensor(Ints... shape_) : shape{shape_...}, data(new T[this->size()]) {
    static_assert(sizeof...(shape_) == N);
    this->UpdateStrides();
  }

  template <class... Ints>
  T& operator()(Ints... indices) {
    static_assert(sizeof...(indices) <= N);
    std::array<int, sizeof...(Ints)> indices_arr{indices...};
    int index = std::inner_product(indices_arr.rbegin(), indices_arr.rend(),
                                   strides.rbegin(), 0);
    return data[index];
  }

  template <class... Ints>
  T operator()(Ints... indices) const {
    static_assert(sizeof...(indices) <= N);
    std::array<int, sizeof...(Ints)> indices_arr{indices...};
    int index = std::inner_product(indices_arr.rbegin(), indices_arr.rend(),
                                   strides.rbegin(), 0);
    return data[index];
  }

  T* begin() { return data.get(); }
  T* end() { return data.get() + this->size(); }

  T const* begin() const { return data.get(); }
  T const* end() const { return data.get() + this->size(); }
};

template <class T, class... Ints>
auto MakeTensor(Ints... ints) {
  return Tensor<T, sizeof...(Ints)>(ints...);
}

template <class... Ints>
auto MakeTensorD(Ints... ints) {
  return Tensor<double, sizeof...(Ints)>(ints...);
}

using VectorD = Tensor<double, 1>;
using MatrixD = Tensor<double, 2>;

auto MakeVectorD(int i0) { return VectorD(i0); }

auto MakeMatrixD(int i0, int i1) { return MatrixD(i0, i1); }

}  // namespace ts
}  // namespace tz
