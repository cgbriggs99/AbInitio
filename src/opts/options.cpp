#include <map>
#include <utility>
#include <cctype>
#include "options.hpp"
#include "default_options.hpp"
#include <stdexcept>

using namespace compchem;
using namespace std;

OptionList::OptionList() {
  this->bool_opts = std::map<string, bool>();
  this->int_opts = std::map<string, int>();
  this->double_opts = std::map<string, double>();
  this->string_opts = std::map<string, string>();
}

OptionList::OptionList(const OptionList &copy) {
  this->bool_opts = std::map<string, bool>();
  this->int_opts = std::map<string, int>();
  this->double_opts = std::map<string, double>();
  this->string_opts = std::map<string, string>();

  for(pair<string, bool> kv : copy.bool_opts) {
    this->bool_opts[kv.first] = kv.second;
  }

  for(pair<string, int> kv : copy.int_opts) {
    this->int_opts[kv.first] = kv.second;
  }

  for(pair<string, double> kv : copy.double_opts) {
    this->double_opts[kv.first] = kv.second;
  }

  for(pair<string, string> kv : copy.string_opts) {
    this->string_opts[kv.first] = kv.second;
  }
}

OptionList::~OptionList() {
  this->bool_opts.clear();
  this->int_opts.clear();
  this->double_opts.clear();
  this->string_opts.clear();
}

string OptionList::to_upper(const string &str) {
  string out;

  for(char c : str) {
    if(islower(c)) {
      out.push_back(toupper(c));
    } else {
      out.push_back(c);
    }
  }

  return out;
}

string OptionList::to_upper(const char *str) {
  string out;

  int pos = 0;
  while(str[pos] != 0) {
    if(islower(str[pos])) {
      out.push_back(toupper(str[pos]));
    } else {
      out.push_back(str[pos]);
    }
    pos++;
  }

  return out;
}

bool OptionList::isoptionbool(const string &str) const {
  return this->bool_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptionbool(const char *str) const {
  return this->bool_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptionint(const string &str) const {
  return this->int_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptionint(const char *str) const {
  return this->int_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptiondouble(const string &str) const {
  return this->double_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptiondouble(const char *str) const {
  return this->double_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptionstring(const string &str) const {
  return this->string_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::isoptionstring(const char *str) const {
  return this->string_opts.count(OptionList::to_upper(str)) != 0;
}

bool OptionList::getbooloption(const string &str) const {
  try {
    return this->bool_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not a boolean option.");
    throw(*new std::out_of_range(message));
  }
}

bool OptionList::getbooloption(const char *str) const {
  try {
    return this->bool_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not a boolean option.");
    throw(*new std::out_of_range(message));
  }
}

int OptionList::getintoption(const string &str) const {
  try {
    return this->int_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not an integer option.");
    throw(*new std::out_of_range(message));
  }
}

int OptionList::getintoption(const char *str) const {
  try {
    return this->int_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not an integer option.");
    throw(*new std::out_of_range(message));
  }
}

double OptionList::getdoubleoption(const string &str) const {
  try {
    return this->double_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not a double option.");
    throw(*new std::out_of_range(message));
  }
}

double OptionList::getdoubleoption(const char *str) const {
  try {
    return this->double_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not a double option.");
    throw(*new std::out_of_range(message));
  }
}

const string &OptionList::getstringoption(const string &str) const {
  try {
    return this->string_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not a string option.");
    throw(*new std::out_of_range(message));
  }
}

const string &OptionList::getstringoption(const char *str) const {
  try {
    return this->string_opts.at(OptionList::to_upper(str));
  } catch(std::out_of_range &e) {
    string message = "Error: key \"";
    message.append(str);
    message.append("\" is not a string option.");
    throw(*new std::out_of_range(message));
  }
}

void OptionList::setbooloption(const string &str, bool value) {
  this->bool_opts[OptionList::to_upper(str)] = value;
}

void OptionList::setbooloption(const char *str, bool value) {
  this->bool_opts[OptionList::to_upper(str)] = value;
}

void OptionList::setintoption(const string &str, int value) {
  this->int_opts[OptionList::to_upper(str)] = value;
}

void OptionList::setintoption(const char *str, int value) {
  this->int_opts[OptionList::to_upper(str)] = value;
}

void OptionList::setdoubleoption(const string &str, double value) {
  this->double_opts[OptionList::to_upper(str)] = value;
}

void OptionList::setdoubleoption(const char *str, double value) {
  this->double_opts[OptionList::to_upper(str)] = value;
}

void OptionList::setstringoption(const string &str, const string &value) {
  this->string_opts[OptionList::to_upper(str)] = string(value);
}

void OptionList::setstringoption(const char *str, const string &value) {
  this->string_opts[OptionList::to_upper(str)] = string(value);
}

void OptionList::setstringoption(const string &str, const char *value) {
  this->string_opts[OptionList::to_upper(str)] = string(value);
}

void OptionList::setstringoption(const char *str, const char *value) {
  this->string_opts[OptionList::to_upper(str)] = string(value);
}

GlobalOptions &GlobalOptions::getsingleton() {
  static GlobalOptions opts;
  static bool first = true;
  if(first) {
    DefaultOptionsFactory::initializeoptions(opts);
    first = false;
  }
  return opts;
}
