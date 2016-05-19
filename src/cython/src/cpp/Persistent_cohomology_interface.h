/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016  INRIA Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PERSISTENT_COHOMOLOGY_INTERFACE_H
#define	PERSISTENT_COHOMOLOGY_INTERFACE_H

#include <gudhi/Persistent_cohomology.h>

namespace Gudhi {

template<class FilteredComplex>
class Persistent_cohomology_interface : public
persistent_cohomology::Persistent_cohomology<FilteredComplex,persistent_cohomology::Field_Zp> {
 public:
  Persistent_cohomology_interface(FilteredComplex* stptr)
  : persistent_cohomology::Persistent_cohomology<FilteredComplex,persistent_cohomology::Field_Zp>(*stptr) { std::cout << "ctor" << std::endl; }
  void get_persistence(int homology_coeff_field, double min_persistence) {
    std::cout << "1" << std::endl; 
    persistent_cohomology::Persistent_cohomology<FilteredComplex,persistent_cohomology::Field_Zp>::init_coefficients(homology_coeff_field);
    std::cout << "2" << std::endl; 
    persistent_cohomology::Persistent_cohomology<FilteredComplex,persistent_cohomology::Field_Zp>::compute_persistent_cohomology(min_persistence);
    std::cout << "3" << std::endl; 
    persistent_cohomology::Persistent_cohomology<FilteredComplex,persistent_cohomology::Field_Zp>::output_diagram();
    std::cout << "4" << std::endl; 
  }

};

} // namespace Gudhi

#endif  // PERSISTENT_COHOMOLOGY_INTERFACE_H

