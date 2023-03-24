#include "test.h" 
#include "../src/opts/options.hpp"
#include "../src/opts/default_options.hpp"

using namespace compchem;
using namespace std;

int test_options() {
  int warns = 0;
  
  OptionList *opts = new OptionList();

  GlobalOptions &glob = GlobalOptions::getsingleton();

  ASSERT_WARN(&glob == &GlobalOptions::getsingleton(), warns);

  opts->setbooloption("Bool Option", true);
  ASSERT_WARN(opts->isoptionbool("bool option"), warns);
  ASSERT_WARN(opts->isoptionbool("BOOL OPTION"), warns);
  ASSERT_WARN(!opts->isoptionint("Bool Option"), warns);
  ASSERT_WARN(opts->getbooloption("bOoL oPtIoN"), warns);

  opts->setintoption("Int Option 123", 5);
  ASSERT_WARN(opts->isoptionint("int option 123"), warns);
  ASSERT_WARN(opts->isoptionint("INT OPTION 123"), warns);
  ASSERT_WARN(opts->getintoption("iNt OpTiOn 123") == 5, warns);

  opts->setdoubleoption("Double Option\n", 1);
  ASSERT_WARN(opts->isoptiondouble("double option\n"), warns);
  ASSERT_WARN(opts->isoptiondouble("DOUBLE OPTION\n"), warns);
  ASSERT_WARN(opts->getdoubleoption("dOuBlE oPtIoN\n") == 1, warns);

  opts->setstringoption("String Option", "string");
  ASSERT_WARN(opts->isoptionstring("string option"), warns);
  ASSERT_WARN(opts->getstringoption("STRING OPTION") == "string", warns);
  ASSERT_WARN(opts->getstringoption("sTrInG oPtIoN") == "string", warns);

  delete opts;

  return warns;
}

int main(void) {
  int warns = 0, errs = 0;

  int ret = test_options();
  if(ret == -1) {
    errs++;
  } else {
    warns += ret;
  }

  return warns + errs;
}
