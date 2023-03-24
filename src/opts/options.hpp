#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <map>
#include <string>

namespace compchem {
  
class OptionList {
protected:
  std::map<std::string, bool> bool_opts;
  std::map<std::string, int> int_opts;
  std::map<std::string, double> double_opts;
  std::map<std::string, std::string> string_opts;
  
  // To upper, for indexing the map.
  static std::string to_upper(const std::string &str);
  static std::string to_upper(const char *str);
  
public:

  OptionList();
  OptionList(const OptionList &copy);
  
  virtual ~OptionList();

  /*
   * 0 if not option, 1 for bool, 2 for int, 3 for float, and 4 for string.
   */
  int isoption(const std::string &str);
  int isoption(const char *str);
  
  bool isoptionbool(const std::string &str);
  bool isoptionbool(const char *str);
  
  bool isoptionint(const std::string &str);
  bool isoptionint(const char *str);

  bool isoptiondouble(const std::string &str);
  bool isoptiondouble(const char *str);
  
  bool isoptionstring(const std::string &str);
  bool isoptionstring(const char *str);
  
  bool getbooloption(const std::string &str);
  bool getbooloption(const char *str);
  
  int getintoption(const std::string &str);
  int getintoption(const char *str);
  
  double getdoubleoption(const std::string &str);
  double getdoubleoption(const char *str);
  
  const std::string &getstringoption(const std::string &str);
  const std::string &getstringoption(const char *str);
  
  void setbooloption(const std::string &str, bool value);
  void setbooloption(const char *str, bool value);
  
  void setintoption(const std::string &str, int value);
  void setintoption(const char *str, int value);
  
  void setdoubleoption(const std::string &str, double value);
  void setdoubleoption(const char *str, double value);
  
  void setstringoption(const std::string &str, const std::string &value);
  void setstringoption(const char *str, const std::string &value);
  void setstringoption(const std::string &str, const char *value);
  void setstringoption(const char *str, const char *value);
};

class GlobalOptions : public OptionList {
private :
  GlobalOptions() : OptionList() {
    ;
  }
  GlobalOptions(const GlobalOptions &copy) = delete;
  GlobalOptions(const GlobalOptions &&move) = delete;

  ~GlobalOptions() = default;

public:

  static GlobalOptions &getsingleton();

};

}


#endif
