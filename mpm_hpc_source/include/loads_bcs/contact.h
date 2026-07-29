#ifndef MPM_CONTACT_H_
#define MPM_CONTACT_H_

//! Alias for JSON
#include "json.hpp"
using Json = nlohmann::json;

namespace mpm {

//! Contact class to store the load on a set
//! \brief Contact class to store the load on a set
//! \details Contact stores the load on a set using mathematical functions, so
//! the load can vary dynamically with time
class Contact {
 public:
  // Constructor
  //! \param[setid] setid  set id
  //! \param[normal] contact surface outward normal vector
  //! contact
  Contact(int setid, unsigned dir, double normal)
      : setid_{setid}, dir_{dir}, normal_{normal} {};

  // Set id
  int setid() const { return setid_; }

  // Direction
  unsigned dir() const {return dir_; }
  
  // Normal
  double normal() const { return normal_; }

 private:
  // ID
  int setid_;
  //
  unsigned dir_;
  // nomal component
  double normal_;
};
}  // namespace mpm
#endif  // MPM_CONTACT_H_
