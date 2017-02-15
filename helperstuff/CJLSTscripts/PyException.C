#include "PyException.h"
#include "TString.h"

PyException::PyException(const TString& message) : message_(message) {}
const char* PyException::what() const throw() {
  return message_;
}
