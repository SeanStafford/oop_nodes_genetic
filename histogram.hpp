#include <iostream>
#include <cmath>
#include <map>
#include <functional>

class Histogram {
  public:
  Histogram() {};
  void Add( double v ) {
    WeightedAdd( v, 1.0 );
  }
  void WeightedAdd( double v, double w ) {
    int key = val_to_binidx(v);
    if( histo.find(key) == histo.end() ) { histo[key] = 0; }
    histo[key] += w;
  }
  std::map<double,double> Frequency() const {
    std::map<double,double> result;
    if( histo.begin() == histo.end() ) { return result; }
    int key_min = histo.begin()->first;
    int key_max = histo.rbegin()->first;
    for( int i = key_min; i <= key_max; i++ ) {
      double val = binidx_to_val( i );
      double n = ( histo.find(i) == histo.end() ) ? 0 : histo.at(i);
      double freq = n / binidx_to_binsize( i );
      result[val] = freq;
    }
    return result;
  }
  std::map<double,double> PDF() const {
    double sum = 0;
    for( auto keyval : histo ) {
      sum += keyval.second;
    }
    std::map<double,double> result;
    for( const auto & keyval : Frequency() ) {
      result[ keyval.first ] = keyval.second / sum;
    }
    return result;
  }
  void Clear() { histo.clear(); }

  protected:
  std::function<int(double)> val_to_binidx;
  std::function<double(int)> binidx_to_val;
  std::function<double(int)> binidx_to_binsize;
  std::map<int, double> histo;
};

class HistoNormalBin : public Histogram {
public:
  HistoNormalBin(double bin_size) : Histogram(), bin(bin_size) {
    val_to_binidx = [=](double v)->int {
      return static_cast<int>( std::floor(v/bin) );
    };
    binidx_to_val = [=](int i)->double {
      return i * bin;
    };
    binidx_to_binsize = [=](int)->double {
      return bin;
    };
  }
private:
  const double bin;
};

class HistoLogBin : public Histogram {
public:
  HistoLogBin(double bin_base = 2.0, double min_bin_idx = -13) : Histogram(), base(bin_base), min(pow(base,min_bin_idx)), min_idx(min_bin_idx) {
    val_to_binidx = [=](double v)->int {
      if( v >= min ) {
        return static_cast<int>( std::floor( log(v)/log(base) ) );
      }
      else {
        return min_idx;
      }
    };
    binidx_to_val = [=](int i)->double {
      return pow(base, i);
    };
    binidx_to_binsize = [=](int i)->double {
      if( i != min_idx ) {
        return pow(base, i+1) - pow(base, i);
      }
      else {
        return pow(base, min_idx+1);
      }
    };
  }
private:
  const double base;
  const double min;
  const int min_idx;
};

