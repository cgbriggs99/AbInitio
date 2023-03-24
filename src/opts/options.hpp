#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <map>
#include <

namespace compchem {
  
class OptionList {
protected:
  std::map<string, bool> bool_opts;
  std::map<string, int> int_opts;
  std::map<string, double> double_opts;
  std::map<string, string> string_opts;
  
  // To upper, for indexing the map.
  static string to_upper(const string &str);
  static string to_upper(const char *str);
  
public:

  OptionList();
  OptionList(const OptionList &copy);
  
  virtual ~OptionList();

  /*
   * 0 if not option, 1 for bool, 2 for int, 3 for float, and 4 for string.
   */
  int isoption(const string &str);
  int isoption(const char *str);
  
  bool isoptionbool(const string &str);
  bool isoptionbool(const char *str);
  
  bool isoptionint(const string &str);
  bool isoptionint(const char *str);

  bool isoptiondouble(const string &str);
  bool isoptiondouble(const char *str);
  
  bool isoptionstring(const string &str);
  bool isoptionstring(const char *str);
  
  bool getbooloption(const string &str);
  bool getbooloption(const char *str);
  
  int getintoption(const string &str);
  int getintoption(const char *str);
  
  double getdoubleoption(const string &str);
  double getdoubleoption(const char *str);
  
  const string &getstringoption(const string &str);
  const string &getstringoption(const char *str);
  
  void setbooloption(const string &str, bool value);
  void setbooloption(const char *str, bool value);
  
  void setintoption(const string &str, int value);
  void setintoption(const char *str, int value);
  
  void setdoubleoption(const string &str, double value);
  void setdoubleoption(const char *str, double value);
  
  void setstringoption(const string &str, const string &value);
  void setstringoption(const char *str, const string &value);
  void setstringoption(const string &str, const char *value);
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
