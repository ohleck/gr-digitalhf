// -*- C++ -*-

#ifndef _VITERBI_HPP_
#define _VITERBI_HPP_

#include <cassert>
#include <algorithm>
#include <array>
#include <bitset>
#include <vector>

// soft-decision viterbi decoder
// based on Phil Karn's libfec
template<std::size_t N, std::size_t K>
class viterbi {
public:

  enum { M = 1<<(K-1) };
  typedef std::vector<int> vec_type;
  typedef std::array<int, M> arr_type;

  viterbi(std::array<std::uint32_t, N> const& polys) // ={0x6d,0x4f}
    : _decisions() // len<<(K-1))
    , _metric()
    , _bits()
    , _prev()
    , _last_max_metric(0) {
    make_tables(polys);
  }

  void reset() {
    std::fill_n(_metric.begin(), _metric.size(), 0);
    _last_max_metric = 0;
  }

  void resize(size_t len) {
    _decisions.resize(len<<(K-1));
  }

  void update(int j, std::array<std::uint8_t,N>const& sym) {
    int s[N];
    for (int l=0; l<N; ++l)
      s[l] = sym[l] ^ 255;

    arr_type new_metric;
    auto jdec = decision(j);

    int mmin[2] = {65535, 65535};
    for (int i=0; i<M; i+=2) {
      int const p0 = _prev[i][0];
      int const p1 = _prev[i][1];
      int m0[2] = { _metric[p0], _metric[p1] };
      int m1[2] = { _metric[p0], _metric[p1] };
      for (int l=0; l<N; ++l) {
	m0[0] += _bits[p0][0][l] ^ s[l];
	m0[1] += _bits[p1][0][l] ^ s[l];
	m1[0] += _bits[p0][1][l] ^ s[l];
	m1[1] += _bits[p1][1][l] ^ s[l];
      }
      jdec[i  ] = m0[0] < m0[1];
      jdec[i+1] = m1[0] < m1[1];

      new_metric[i  ] = m0[jdec[i  ]];
      new_metric[i+1] = m1[jdec[i+1]];
      mmin[0] = std::min(mmin[0], new_metric[i  ]);
      mmin[1] = std::min(mmin[1], new_metric[i+1]);
    }
    // avoid path metric overflow
    int const imin = std::min(mmin[0], mmin[1]);
    if (imin > (1<<15)) {
      _last_max_metric -= imin;
      for (int i=0; i<M; ++i)
	_metric[i] = new_metric[i] - imin;
    } else {
      std::copy(new_metric.begin(), new_metric.end(), _metric.begin());
    }
  }

  float chainback(std::vector<uint8_t>& v) {
    return chainback(v.begin(), v.end());
  }

  float chainback(std::vector<uint8_t>::iterator begin,
                  std::vector<uint8_t>::iterator end) {
    assert(std::distance(begin, end) == ssize_t((_decisions.size()>>(K-1))));

    auto const imax = std::max_element(_metric.begin(), _metric.end());
    int idx_max = std::distance(_metric.begin(), imax);
    for (int k=_decisions.size()>>(K-1); k!=0; --k) {
      begin[k-1] = idx_max&1;
      //idx_max    = _prev[idx_max][decision(k-1)[idx_max]];
      idx_max    = (idx_max>>1) + (decision(k-1)[idx_max] != 0)*M/2;
    }
    int const max_metric = *imax;
    float const quality  = float(max_metric - _last_max_metric)/255.0;
    _last_max_metric = max_metric;
    return quality;
  }

protected:
  vec_type::iterator decision(int i) {
    return _decisions.begin() + (i<<(K-1));
  }
  void make_tables(std::array<std::uint32_t, N> const& polys) {
    for (int i=0, n=1<<K; i<n; ++i) {
      for (int l=0; l<N; ++l) {
	std::bitset<K> const b(polys[l]&i);
	_bits[i>>1][i%2][l] = 255*(b.count()%2);
      }
    }
    for (int i=0; i<M; ++i) {
      _prev[i][0] = (i>>1);
      _prev[i][1] = _prev[i][0] + M/2;
    }
  }
private:
  vec_type _decisions;
  arr_type _metric;
  int      _bits[M][2][N];
  int      _prev[M][2];
  int      _last_max_metric;
} ;

#endif // _VITERBI2_HPP_
