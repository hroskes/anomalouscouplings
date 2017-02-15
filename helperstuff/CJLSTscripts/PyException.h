#include "TPyException.h"
#include "TString.h"

class PyException : public PyROOT::TPyException {
  private:
    TString message_;
  public:
    PyException(const TString& message);
    virtual const char* what() const throw();
};
