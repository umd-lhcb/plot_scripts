#ifndef H_AXIS
#define H_AXIS

#include <cstddef>

#include <vector>
#include <set>
#include <string>

#include "core/named_func.hpp"

class Axis{
public:
  Axis(const std::vector<double> &bins,
       const NamedFunc &var,
       const std::string &title = "unset",
       const std::set<double> &cut_vals = {});
  Axis(std::size_t nbins,
       double xmin, double xmax,
       const NamedFunc &var,
       const std::string &title = "unset",
       const std::set<double> &cut_vals = {});
  Axis(const std::vector<double> &bins,
       const std::vector<NamedFunc> &vars,
       const std::string &title = "unset",
       const std::set<double> &cut_vals = {});
  Axis(std::size_t nbins,
       double xmin, double xmax,
       const std::vector<NamedFunc> &vars,
       const std::string &title = "unset",
       const std::set<double> &cut_vals = {});

  Axis(const Axis &) = default;
  Axis& operator=(const Axis &) = default;
  Axis(Axis &&) = default;
  Axis& operator=(Axis &&) = default;
  ~Axis() = default;

  std::size_t Nbins() const;
  Axis & Bins(const std::vector<double> &bins);
  Axis & Bins(std::size_t nbins, double xmin, double xmax);
  const std::vector<double> & Bins() const;
  double AvgBinWidth() const;

  std::string Title() const;

  NamedFunc var_;
  std::vector<NamedFunc> vars_;//!< Variables to be plotted
  std::string title_;//!< Axis title without units
  std::string units_;//!< Units of Axis::var_
  std::set<double> cut_vals_;//!< Values of HistoDef::var_ for which to plot a line

private:
  std::vector<double> bins_;//!<List of bin edges
  
  static std::vector<double> GetEdges(std::size_t nbins, double xmin, double xmax);
  void ParseUnits();
};

#endif
