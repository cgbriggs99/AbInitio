#ifndef TEI_ARRAY_HPP
#define TEI_ARRAY_HPP

namespace compchem {

class TEIArray {
private :
  double *data;
  int dim;
public :
  TEIArray(double * __restrict__ data, int dim);
  TEIArray(int dim);
  TEIArray(const TEIArray &copy);

  ~TEIArray();

  double at(int mu, int nu, int lam, int sig) const;
  double &at(int mu, int nu, int lam, int sig);
  double at_direct(int index) const;
  double &at_direct(int index);

  double operator()(int mu, int nu, int lam, int sig) const;
  double &operator()(int mu, int nu, int lam, int sig);

  const double *getdata() const;
  int getdim() const;
  int getindex(int mu, int nu, int lam, int sig) const;
  int getsize() const;

  void indextoquad(int index, int *mu, int *nu, int *lam, int *sig) const;
};
  
int biggest_trinum_leq(int val);

int biggest_triind_leq(int val);

int triangular_num(int index);
  
}


#endif
