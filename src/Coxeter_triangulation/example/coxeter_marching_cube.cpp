// #define LEMINISCATE_OF_GERONO
// #define CIRCLE
// #define WAVY_CIRCLE
// #define SMILEY
// #define SPHERE
// #define TORUS
// #define DOUBLE_TORUS
// # define DEBUG_TRACES


#include <iostream>
#include <vector>
#include <fstream>

#include <gudhi/Simple_coxeter_system_remastered.h>
#include <gudhi/Coxeter_complex/Trie.h>
#include <gudhi/Hasse_diagram_persistence.h>
#include "output_hasse_to_medit.h"


// #include <gudhi/Points_off_io.h>
// #include <gudhi/Coxeter_system.h>
// #include <gudhi/Coxeter_complex.h>
// #include <gudhi/Coxeter_complex/Off_point_range.h>
#include <gudhi/Clock.h>

// #include "memory_usage.h"
// #include "cxx-prettyprint/prettyprint.hpp"

using Cell_id = typename Simple_coxeter_system::Alcove_id;
using Point_d = Eigen::VectorXd;
using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
using Hasse_boundary = std::vector<std::pair<Hasse_cell*, int> >;

// The vertices are always smaller than the non-vertices
struct Hasse_cell_comparator {
  bool operator() (Hasse_cell* l_it, Hasse_cell* r_it) {
    if (l_it->get_dimension() == 0)
      if (r_it->get_dimension() == 0)
        return l_it < r_it;
      else
        return true;
    else if (r_it->get_dimension() == 0)
      return false;
    else
      return l_it->get_boundary() < r_it->get_boundary();
  }
};
using Hasse_diagram = std::set<Hasse_cell*, Hasse_cell_comparator>;

struct Cell_id_comparator {
  bool operator()(const Cell_id& lhs, const Cell_id& rhs) {
    if (lhs.size() < rhs.size())
      return true;
    else if (lhs.size() > rhs.size())
      return false;
    for (std::size_t k = 0; k < lhs.size(); ++k) {
      if (lhs[k] < rhs[k])
        return true;
      else if (lhs[k] > rhs[k])
        return false;
      else if (!lhs.is_fixed(k) && rhs.is_fixed(k))
        return true;
      else if (lhs.is_fixed(k) && !rhs.is_fixed(k))
        return false;
    }
    return false;
  }
};
using VC_map = std::map<Cell_id, Hasse_cell*, Cell_id_comparator>;
using VP_map = std::map<Hasse_cell*, Point_d>;
using Full_cell_trie = Trie<Cell_id>;

#ifdef CIRCLE
const unsigned amb_d = 2; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double r = 5;
Eigen::Vector2d point1(r, 0.0);
std::vector<Point_d> seed_points = {point1};
std::string name = "circle";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0), y = p(1);
  Eigen::VectorXd coords(cod_d);
  coords(0) = x*x + y*y - r*r;
  return coords;
}
#endif

#ifdef LEMINISCATE_OF_GERONO
const unsigned amb_d = 2; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double a = 1;
Eigen::Vector2d point1(std::sqrt(a)+0.01, 0.01);
std::vector<Point_d> seed_points = {point1};
std::string name = "leminiscate_of_gerono";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0)-0.01, y = p(1)-0.01;
  Eigen::VectorXd coords(cod_d);
  coords(0) = x*x*x*x - a*x*x + a*y*y;
  return coords;
}
#endif

#ifdef WAVY_CIRCLE
const unsigned amb_d = 2; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double r = 1;
double sr = 0.1;
int n = 25;
Eigen::Vector2d point1(r, 0.0);
std::vector<Point_d> seed_points = {point1};
std::string name = "wavy_circle";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0), y = p(1);
  double pi = 3.14159265358;
  Eigen::VectorXd coords(cod_d);
  double theta = 0;
  if (x > 0)
    theta = std::atan(y/x);
  else if (x < 0 && y >= 0)
    theta = std::atan(y/x) + pi;
  else if (x < 0 && y < 0)
    theta = std::atan(y/x) - pi;
  else if (x == 0 && y > 0)
    theta = pi/2;
  else
    theta = -pi/2;
  double real_pos = r + sr*std::sin(n*theta);
  coords(0) = x*x + y*y - real_pos*real_pos;
  return coords;
}
#endif

#ifdef SMILEY
const unsigned amb_d = 2; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double r = 4;
double sr = 0.5;
double ir = 3.5;
double ey = 1.5, ex = 1.5;
double offset = 0.01;
Eigen::Vector2d point1(r+offset, offset);
Eigen::Vector2d point2(offset, offset);
Eigen::Vector2d point3(-ex+sr+offset, ey+offset);
Eigen::Vector2d point4(ex+sr+offset, ey+offset);
std::vector<Point_d> seed_points = {point1, point2, point3, point4};
std::string name = "smiley";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0) - offset, y = p(1) - offset;
  Eigen::VectorXd coords(cod_d);
  if (x*x + y*y > r*r)
    coords(0) = 1;
  else if (x*x + y*y == r*r)
    coords(0) = 0;
  else if ((x*x + y*y < ir*ir) && (y < 0))
    coords(0) = 1;
  else if ((x*x + y*y == ir*ir) && (y < 0))
    coords(0) = 0;
  else if ((x*x + y*y <= ir*ir) && (y == 0))
    coords(0) = 0;
  else if ((x+ex)*(x+ex) + (y-ey)*(y-ey) < sr*sr)
    coords(0) = 1;
  else if ((x+ex)*(x+ex) + (y-ey)*(y-ey) == sr*sr)
    coords(0) = 0;
  else if ((x-ex)*(x-ex) + (y-ey)*(y-ey) < sr*sr)
    coords(0) = 1;
  else if ((x-ex)*(x-ex) + (y-ey)*(y-ey) == sr*sr)
    coords(0) = 0;
  else
    coords(0) = -1;
  return coords;
}
#endif

#ifdef SPHERE
const unsigned amb_d = 3; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double r = 5;
Eigen::Vector3d point1(r, 0.0, 0.0);
std::vector<Point_d> seed_points = {point1};
std::string name = "sphere";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0), y = p(1), z = p(2);
  Eigen::VectorXd coords(cod_d);
  coords(0) = x*x + y*y + z*z - r*r;
  return coords;
}
#endif

#ifdef TORUS
const unsigned amb_d = 3; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double r = 5;
double sr = 1;
Eigen::Vector3d point1(r+sr, 0.0, 0.0);
std::vector<Point_d> seed_points = {point1};
std::string name = "torus";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0), y = p(1), z = p(2);
  Eigen::VectorXd coords(cod_d);
  coords(0) = (z*z + (std::sqrt(x*x + y*y) - r)*(std::sqrt(x*x + y*y) - r) - sr*sr);
  return coords;
}
#endif


#ifdef DOUBLE_TORUS
const unsigned amb_d = 3; // Ambient (domain) dimension
const unsigned cod_d = 1; // Codomain dimension
double r = 5;
double sr = 1;
double ofs = 7;
Eigen::Vector3d point1(r+sr, 0.0, 0.0);
Eigen::Vector3d point2(r+sr+ofs, 0.0, 0.0);
std::vector<Point_d> seed_points = {point1, point2};
std::string name = "double_torus";

Eigen::VectorXd f(Eigen::VectorXd p) {
  double x = p(0), y = p(1), z = p(2);
  Eigen::VectorXd coords(cod_d);
  coords(0) = (z*z + (std::sqrt(x*x + y*y) - r)*(std::sqrt(x*x + y*y) - r) - sr*sr)
    * (y*y + (std::sqrt((x-ofs)*(x-ofs) + z*z) - r)*(std::sqrt((x-ofs)*(x-ofs) + z*z) - r) - sr*sr);
  return coords;
}
#endif


// using Point_vector = std::vector< Point_d >;

// using Coxeter_complex = Gudhi::Coxeter_complex<Point_vector, Coxeter_system>;

void mark(const Cell_id& c_id, Full_cell_trie& trie) {
  trie.add(c_id);
  // std::cout << "Added " << c_id << ". Trie is now: " << trie << "\n";
#ifdef DEBUG_TRACES
  std::cout << "Size of the trie: " << trie.size() << "\n";
#endif
}

bool is_marked(const Cell_id& c_id, Full_cell_trie& trie) {
  return trie.contains(c_id);
}

void add_hasse_vertex(const Cell_id& f_id, Eigen::VectorXd& cart_coords, Hasse_diagram& hd, VC_map& vc_map, VP_map& vp_map) {
  Hasse_cell* new_cell = new Hasse_cell(0);
  auto res_pair = vc_map.emplace(f_id, new_cell);
  if (!res_pair.second) {
    delete new_cell;
    return;
  }
  hd.emplace(new_cell);
  vp_map.emplace(new_cell, cart_coords);
}

Hasse_cell* insert_hasse_subdiagram(const Cell_id& c_id,
                                    const std::vector<Cell_id>& meet_faces,
                                    Hasse_diagram& hd,
                                    VC_map& vc_map,
                                    unsigned cod_d,
                                    const Simple_coxeter_system& cs) {
  // if (c_id[0] == -2 && c_id[1] == -4 && c_id[2] == -5 && c_id[3] == -1 && c_id[4] == -4 && c_id[5] == -6)
  //   std::cout << "Problem!\n";
#ifdef DEBUG_TRACES
  std::cout << "  Insert_hasse_subdiagram for " << c_id << ". Meet_faces = " << meet_faces << "\n";
#endif
  if (c_id.dimension() == cod_d) {
    if (std::find(meet_faces.begin(), meet_faces.end(), c_id) != meet_faces.end())
      return vc_map.find(c_id)->second;
    else
      return 0;
  }
  else {
    Hasse_cell* new_cell = new Hasse_cell(c_id.dimension() - cod_d);
    Hasse_boundary& boundary = new_cell->get_boundary();
    for (auto f_id: cs.face_range(c_id, c_id.dimension() - 1)) {
      Hasse_cell* facet_cell = insert_hasse_subdiagram(f_id, meet_faces, hd, vc_map, cod_d, cs);
      if (facet_cell != 0)
        boundary.push_back(std::make_pair(facet_cell, 1));
    }
    if (boundary.empty()) {
      delete new_cell;
      return 0;
    }
    else {
#ifdef DEBUG_TRACES
      std::cout << "Boundary of the new cell = " << boundary << "\n";
      if (new_cell->get_dimension() == 1 && boundary.size() != 2) {
        std::cout << "Problem!\n";
      }
#endif
      auto res_pair = hd.emplace(new_cell);
      if (!res_pair.second) {
        delete new_cell;
        return *res_pair.first;
      }
      else
        return new_cell;
    }
  }
}

bool intersects(const Cell_id& f_id,
                Eigen::MatrixXd& li_matrix,
                Eigen::VectorXd& last_column,
                Hasse_diagram& hd,
                VC_map& vc_map,
                VP_map& vp_map,
		unsigned cod_d,
                const Simple_coxeter_system& cs) {
#ifdef DEBUG_TRACES
  std::cout << " Intersects called for face " << f_id << "\n";
#endif
  std::size_t amb_d = li_matrix.cols();
#ifdef DEBUG_TRACES
  std::cout << "amb_d = " << amb_d << ", cod_d = " << cod_d << "\n";
#endif
  std::size_t curr_row = cod_d;
  for (std::size_t k: cs.normal_basis_range(f_id)) {
    std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
    std::size_t i = (j*j + j - 2)/2 - k;
    for (std::size_t col = 0; col < amb_d; ++col)
      li_matrix(curr_row, col) = 0;
    for (std::size_t l = i; l < j; ++l)
      for (std::size_t col = 0; col < amb_d; ++col)
        li_matrix(curr_row, col) += cs.simple_root_matrix()(l, col);
    last_column(curr_row) = f_id[k] / f_id.level();
    curr_row++;
  }
#ifdef DEBUG_TRACES
  std::cout << "  lin_interpolation_matrix (after completion):\n" << li_matrix << "\n";
  std::cout << "  last_column:\n" << last_column << "\n";
#endif
  Eigen::VectorXd intersection = li_matrix.colPivHouseholderQr().solve(last_column);
#ifdef DEBUG_TRACES
  std::cout << "  intersection:\n" << intersection << "\n";
  std::cout << "  rel_error = " << (li_matrix * intersection - last_column).norm() / last_column.norm() << "\n";
#endif
  if ((li_matrix * intersection - last_column).norm() / last_column.norm() > 1e-5 / f_id.level())
    return false;
  Cell_id i_id = cs.query_point_location(intersection, f_id.level());
#ifdef DEBUG_TRACES
  std::cout << "  cell containing the intersection: " << i_id << "\n";
#endif
  for (std::size_t k = 0; k < i_id.size(); ++k)
    if (!f_id.is_fixed(k) && i_id[k] != f_id[k])
      return false;
  add_hasse_vertex(f_id, intersection, hd, vc_map, vp_map);
  return true;
}

template <typename Function>
void compute_hasse_diagram(std::vector<Point_d>& seed_points,
                           double level,
                           unsigned amb_d,
                           unsigned cod_d,
                           Hasse_diagram& hd,
                           VP_map& vp_map,
                           const Function& f) {
  const Simple_coxeter_system cs('A', amb_d);
  VC_map vc_map;
  Full_cell_trie trie(level, amb_d);
  Full_cell_trie visit_stack(level, amb_d);
#ifdef DEBUG_TRACES
  std::cout << "root_t_:\n" << cs.simple_root_matrix() << "\n";
#endif
  for (const Point_d& p: seed_points) {
    visit_stack.add(cs.query_point_location(p, level));
  }
  
  while (!visit_stack.empty()) {
    Cell_id c_id = visit_stack.pop();
#ifdef DEBUG_TRACES
    std::cout << "Simplex: "  << c_id << "\n";
#endif
    mark(c_id, trie);
    std::vector<Cell_id> meet_faces;
    Eigen::MatrixXd
      point_matrix(amb_d + 1, amb_d + 1),
      value_matrix(cod_d, amb_d + 1);
    int i = 0;
    for (Cell_id v_id: cs.face_range(c_id, 0)) {
      Eigen::VectorXd cart_coords = cs.cartesian_coordinates_of_vertex(v_id);
#ifdef DEBUG_TRACES
      std::cout << " Vertex: "  << v_id << "\n" << cart_coords << "\n";
#endif
      for (std::size_t j = 0; j < amb_d; ++j)
        point_matrix(j, i) = cart_coords(j);
      point_matrix(amb_d, i) = 1;
      Eigen::VectorXd val_vector = f(cart_coords);
      for (std::size_t j = 0; j < cod_d; ++j)
        value_matrix(j, i) = val_vector(j);
      ++i;
    }
#ifdef DEBUG_TRACES
    std::cout << " point_matrix:\n" << point_matrix << "\n";
#endif
    Eigen::MatrixXd point_matrix_inv = point_matrix.colPivHouseholderQr().inverse();
#ifdef DEBUG_TRACES
    std::cout << " point_matrix_inv:\n" << point_matrix_inv << "\n";
    std::cout << " value_matrix:\n" << value_matrix << "\n";
#endif
    Eigen::MatrixXd lin_interpolation_matrix = value_matrix * point_matrix_inv;
    Eigen::VectorXd last_column = lin_interpolation_matrix.col(amb_d);
#ifdef DEBUG_TRACES
    std::cout << " lin_interpolation_matrix (before completion):\n" << lin_interpolation_matrix << "\n";
#endif
    lin_interpolation_matrix.conservativeResize(amb_d, amb_d);
    last_column.conservativeResize(amb_d);
    last_column = -last_column;
    for (Cell_id f_id: cs.face_range(c_id, cod_d))
      if (intersects(f_id, lin_interpolation_matrix, last_column, hd, vc_map, vp_map, cod_d, cs)) {
#ifdef DEBUG_TRACES
        std::cout << " Result = true\n";
#endif
        meet_faces.push_back(f_id);
      }
#ifdef DEBUG_TRACES
      else
        std::cout << " Result = false\n";
    std::cout << " Size of vc_map = " << vc_map.size() << "\n";
#endif
    insert_hasse_subdiagram(c_id, meet_faces, hd, vc_map, cod_d, cs);
    for (const Cell_id& f_id: meet_faces)
      // if (!is_marked(f_id))
      for (Cell_id cf_id: cs.coface_range(f_id, amb_d))
        if (!is_marked(cf_id, trie))
          visit_stack.add(cf_id);    
  }
}

/** TEST CIRCLE */
void test_circle(double level) {
  const unsigned amb_d = 2; // Ambient (domain) dimension
  const unsigned cod_d = 1; // Codomain dimension
  double r = 5;
  Eigen::Vector2d point1(r, 0.0);
  std::vector<Point_d> seed_points = {point1};
  std::string name = "circle";
  std::cout << "Test 1: " << name << "...\n";

  struct Function {
    Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
      double x = p(0), y = p(1);
      Eigen::VectorXd coords(cod_d_);
      coords(0) = x*x + y*y - r_*r_;
      return coords;
    }

    Function(unsigned cod_d, double r) : cod_d_(cod_d), r_(r) {}
    unsigned cod_d_;
    double r_;
  } f(cod_d, r);
  Hasse_diagram hd;
  VP_map vp_map;
  Gudhi::Clock t;
  compute_hasse_diagram(seed_points, level, amb_d, cod_d, hd, vp_map, f);
  t.end();
  std::vector<unsigned> dimensions(amb_d-cod_d+1, 0);
  for (auto cell: hd) {
    dimensions[cell->get_dimension()]++;
  }
  std::cout << "Simplices by dimension: " << dimensions << "\n";
  std::cout << "Reconstruction time: " <<  t.num_seconds() << "s\n";
  output_hasse_to_medit(hd, vp_map, "marching_cube_output_"+name);
  std::cout << "Wrote the reconstruction in marching_cube_output_" << name << ".mesh\n";
}

/** TEST WAVY CIRCLE */
void test_wavy_circle(double level) {
  const unsigned amb_d = 2; // Ambient (domain) dimension
  const unsigned cod_d = 1; // Codomain dimension
  unsigned n = 25;
  double r = 5;
  double sr = 0.5;
  Eigen::Vector2d point1(r, 0.0);
  std::vector<Point_d> seed_points = {point1};
  std::string name = "wavy_circle";
  std::cout << "Test 2: " << name << "...\n";

  struct Function {
    Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
      double x = p(0), y = p(1);
      double pi = 3.14159265358;
      Eigen::VectorXd coords(cod_d_);
      double theta = 0;
      if (x > 0)
	theta = std::atan(y/x);
      else if (x < 0 && y >= 0)
	theta = std::atan(y/x) + pi;
      else if (x < 0 && y < 0)
	theta = std::atan(y/x) - pi;
      else if (x == 0 && y > 0)
	theta = pi/2;
      else
	theta = -pi/2;
      double real_pos = r_ + sr_*std::sin(n_*theta);
      coords(0) = x*x + y*y - real_pos*real_pos;
      return coords;
    }

    Function(unsigned cod_d, unsigned n, double r, double sr)
      : cod_d_(cod_d), n_(n), r_(r), sr_(sr) {}
    unsigned cod_d_, n_;
    double r_, sr_;
    
  } f(cod_d, n, r, sr);
  Hasse_diagram hd;
  VP_map vp_map;
  Gudhi::Clock t;
  compute_hasse_diagram(seed_points, level, amb_d, cod_d, hd, vp_map, f);
  t.end();
  std::vector<unsigned> dimensions(amb_d-cod_d+1, 0);
  for (auto cell: hd) {
    dimensions[cell->get_dimension()]++;
  }
  std::cout << "Simplices by dimension: " << dimensions << "\n";
  std::cout << "Reconstruction time: " <<  t.num_seconds() << "s\n";
  output_hasse_to_medit(hd, vp_map, "marching_cube_output_"+name);
  std::cout << "Wrote the reconstruction in marching_cube_output_" << name << ".mesh\n";
}

/** TEST SMILEY */
void test_smiley(double level) {
  const unsigned amb_d = 2; // Ambient (domain) dimension
  const unsigned cod_d = 1; // Codomain dimension
  double r = 4;
  double sr = 0.5;
  double ir = 3.5;
  double ey = 1.5, ex = 1.5;
  double offset = 0.01;
  Eigen::Vector2d point1(r+offset, offset);
  Eigen::Vector2d point2(offset, offset);
  Eigen::Vector2d point3(-ex+sr+offset, ey+offset);
  Eigen::Vector2d point4(ex+sr+offset, ey+offset);
  std::vector<Point_d> seed_points = {point1, point2, point3, point4};
  std::string name = "smiley";
  std::cout << "Test 3: " << name << "...\n";

  struct Function {
    Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
      double x = p(0) - offset, y = p(1) - offset;
      Eigen::VectorXd coords(cod_d);
      if (x*x + y*y > r*r)
	coords(0) = 1;
      else if (x*x + y*y == r*r)
	coords(0) = 0;
      else if ((x*x + y*y < ir*ir) && (y < 0))
	coords(0) = 1;
      else if ((x*x + y*y == ir*ir) && (y < 0))
	coords(0) = 0;
      else if ((x*x + y*y <= ir*ir) && (y == 0))
	coords(0) = 0;
      else if ((x+ex)*(x+ex) + (y-ey)*(y-ey) < sr*sr)
	coords(0) = 1;
      else if ((x+ex)*(x+ex) + (y-ey)*(y-ey) == sr*sr)
	coords(0) = 0;
      else if ((x-ex)*(x-ex) + (y-ey)*(y-ey) < sr*sr)
	coords(0) = 1;
      else if ((x-ex)*(x-ex) + (y-ey)*(y-ey) == sr*sr)
	coords(0) = 0;
      else
	coords(0) = -1;
      return coords;
    }

    Function(unsigned cod_d_,
	     double r_,
	     double sr_,
	     double ir_,
	     double ey_,
	     double ex_,
	     double offset_)
      : cod_d(cod_d_), r(r_), sr(sr_), ir(ir_), ey(ey_), ex(ex_), offset(offset_) {}
    unsigned cod_d;
    double r = 4;
    double sr = 0.5;
    double ir = 3.5;
    double ey = 1.5, ex = 1.5;
    double offset = 0.01;    
  } f(cod_d, r, sr, ir, ey, ex, offset);
  Hasse_diagram hd;
  VP_map vp_map;
  Gudhi::Clock t;
  compute_hasse_diagram(seed_points, level, amb_d, cod_d, hd, vp_map, f);
  t.end();
  std::vector<unsigned> dimensions(amb_d-cod_d+1, 0);
  for (auto cell: hd) {
    dimensions[cell->get_dimension()]++;
  }
  std::cout << "Simplices by dimension: " << dimensions << "\n";
  std::cout << "Reconstruction time: " <<  t.num_seconds() << "s\n";
  output_hasse_to_medit(hd, vp_map, "marching_cube_output_"+name);
  std::cout << "Wrote the reconstruction in marching_cube_output_" << name << ".mesh\n";
}


int main(int argc, char * const argv[]) {
  test_circle(1.5);
  test_wavy_circle(2.5);
  test_smiley(10.1);
}
