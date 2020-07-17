#ifndef H_GENERATE_BABY
#define H_GENERATE_BABY

#include <string>
#include <map>
#include <set>
#include <vector>

struct SimpleVariable{
  SimpleVariable(const std::string &type, const std::string &name);
  std::string type_;//!<Type of variable (e.g., int, float, etc.)
  std::string name_;//!<Name of variable (e.g., ht, met, etc.)
};

class Variable{
public:
  Variable(const std::string &name);

  const std::string & Name() const;
  std::string & Name();

  std::string Type() const;
  std::string Type(const std::string &baby_type) const;
  std::string DecoratedType() const;
  std::string DecoratedType(const std::string &baby_type) const;

  bool HasEntry(const std::string &baby_type) const;
  void SetEntry(const std::string &baby_type,
                const std::string &type);

  bool ImplementInBase() const;
  bool VirtualInBase() const;
  bool ImplementIn(const std::string &baby_type) const;
  std::string VarIndex(const std::string &baby_type) const;
  bool MultipleTypes() const;
  bool EverythingIn(const std::string &baby_type) const;

  bool operator<(const Variable& other) const;

private:
  std::set<std::string> GetTypeSet() const;

  std::string name_;//!<Name of variable (e.g., ht, met, etc.)
  std::map<std::string, std::string> type_map_;//!<Map from Baby type (basic, full, etc.) to variable type (int, float, etc.)
};

std::set<Variable> GetVariables(const std::vector<std::string> &files, std::vector<std::string> &tree_names);

bool IsComment(const std::string &line);

SimpleVariable GetVariable(std::string line);

void RemoveExtraSpaces(std::string &line);

void WriteBaseHeader(const std::set<Variable> &vars,
                     const std::vector<std::string> &types);

void WriteBaseSource(const std::set<Variable> &vars,
                     const std::vector<std::string> &types);

void WriteSpecializedHeader(const std::set<Variable> &vars,
                            const std::string &type,
                            const std::vector<std::string> &baby_types);

void WriteSpecializedSource(const std::set<Variable> &vars,
                            const std::string &type, const std::string &treename,
                            const std::vector<std::string> &baby_types);

void WriteMergedHeader(const std::set<Variable> &vars,
                       const std::vector<std::string> &types);

#endif
