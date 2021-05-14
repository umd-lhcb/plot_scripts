#include "core/table_row.hpp"

TableRow::TableRow(const std::string &label,
                   std::size_t lines_before,
                   std::size_t lines_after):
  label_(label),
  cuts_({"1"}),
  cut_("1"),
  weight_("1"),
  lines_before_(lines_before),
  lines_after_(lines_after),
  is_data_row_(false){
  }

TableRow::TableRow(const std::string &label,
                   const NamedFunc &cut,
                   std::size_t lines_before,
                   std::size_t lines_after,
                   const NamedFunc &weight):
  label_(label),
  cuts_({cut}),
  cut_(cut),
  weight_(weight),
  lines_before_(lines_before),
  lines_after_(lines_after),
  is_data_row_(true){
  
  }

TableRow::TableRow(const std::string &label,
                   const std::vector<NamedFunc> &cuts,
                   std::size_t lines_before,
                   std::size_t lines_after,
                   const NamedFunc &weight):
  label_(label),
  cuts_(cuts),
  cut_(cuts[0]),
  weight_(weight),
  lines_before_(lines_before),
  lines_after_(lines_after),
  is_data_row_(true){
  }
