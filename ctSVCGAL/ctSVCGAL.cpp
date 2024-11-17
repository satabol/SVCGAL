/*
# This file is part of project Sverchok. It's copyrighted by the contributors
# recorded in the version control history of the file, available from
# its original location https://github.com/nortikin/sverchok/commit/master
#
# SPDX-License-Identifier: GPL3
# License-Filename: LICENSE
*/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/compute_outer_frame_margin.h>

#include <CGAL/Boolean_set_operations_2/Gps_polygon_validation.h>

//#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
//#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/draw_straight_skeleton_2.h>
#include <CGAL/draw_surface_mesh.h>

//#include "CGAL/input_helpers.h" // polygon reading, random polygon with weights generation
#include <CGAL/extrude_skeleton.h>

#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>

//#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <vector>
#include <memory>
//#include <list>

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>

namespace SS = CGAL::CGAL_SS_i;
//namespace PMP = CGAL::Polygon_mesh_processing;
namespace Internal = CGAL::Straight_skeleton_extrusion::internal;

// Kernel choice:
// EPICK: Robust and fast
// EPECK_with_sqrt: Exact and slow
// EPECK: More robust, and less slow than EPECK_with_sqrt

//using ExceptionWithFailedContours = CGAL::Straight_skeleton_extrusion::internal::ExceptionWithFailedContours<Point_2>;

#include <unordered_map>

// This actually works? 
// https://stackoverflow.com/a/25155315/2391876
#ifdef _WIN32 // This is what Visual Studio defines 
#define DLLEXPORT __declspec(dllexport)
// This resolves the "M_PI (pi) is not defined properly" error. 
// (https://stackoverflow.com/questions/6563810/m-pi-works-with-math-h-but-not-with-cmath-in-visual-studio)
// There may be better approaches to this, but I'm not sure what. 
// I also really feel like this should not be necessary, see OpenSubDiv -> CMakeLists.txt -> /D_USE_MATH_DEFINES. 
#define _USE_MATH_DEFINES
#include <math.h>
//#include "ctypes_SVCGAL.h"
//using namespace CGAL::Straight_skeleton_extrusion::internal;
#elif __linux__
#define DLLEXPORT 
#elif __APPLE__
#define DLLEXPORT 
#endif

#define VAL2STR(val) #val

namespace CGAL {
  namespace Straight_skeleton_extrusion {
    namespace internal {

      enum Err {
        _0001,
        _0002,
        _0003,
        _0004,
        _0005,
        _0006,
        _0007,
        _0008,
        _0009,
        _0010,
        _0011,
        _0012,
        _0013,
        _0014,
        _0015,
        _0016,
        _0017,
        _0018,
        _0019,
        _0020,
        _0021,
        _0022,
        _0023,
        _0024,
        _0025,
        _0026,
        _0027,
        _0028,
        _0029,
        _0030,
        _0031,
        _0032,
      };

      using K = CGAL::Exact_predicates_inexact_constructions_kernel;
      using FT = K::FT;
      using Point_2 = K::Point_2;
      using Segment_2 = K::Segment_2;
      using Line_2 = K::Line_2;
      using Point_3 = K::Point_3;
      using Vector_3 = K::Vector_3;

      using Polygon_2 = CGAL::Polygon_2<K>;
      using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

      using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
      using Straight_skeleton_2_ptr = std::shared_ptr<Straight_skeleton_2>;

      using Mesh = CGAL::Surface_mesh<Point_3>;

      typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
      typedef boost::graph_traits<Mesh>::edge_descriptor   edge_descriptor;


      //---------------- Vertex container implementation. ----------------
      struct Vertex {
        // Minimal required interface ----------------------
        Vertex() {}

        Vertex(Vertex const& src) {
          _position[0] = src._position[0];
          _position[1] = src._position[1];
          _position[2] = src._position[2];
        }

        void Clear(void* = 0) {
          _position[0] = _position[1] = _position[2] = 0.0f;
        }

        void AddWithWeight(Vertex const& src, float weight) {
          _position[0] += weight * src._position[0];
          _position[1] += weight * src._position[1];
          _position[2] += weight * src._position[2];
        }

        void SetPosition(float x, float y, float z) {
          _position[0] = x;
          _position[1] = y;
          _position[2] = z;
        }

        const float* GetPosition() const {
          return _position;
        }

      private:
        float _position[3];
      };

      //---------------- SVCGAL ----------------

      struct ObjVertex {
        uint32_t p = (uint32_t)-1;
        uint32_t n = (uint32_t)-1;
        uint32_t uv = (uint32_t)-1;

        ObjVertex() {
        }

        ObjVertex(uint32_t pi) {
          p = pi;
        }

        bool operator==(const ObjVertex& v) const {
          return v.p == p && v.n == n && v.uv == uv;
        }
      };

      struct ObjVertexHash {
        std::size_t operator()(const ObjVertex& v) const {
          size_t hash = std::hash<uint32_t>()(v.p);
          hash = hash * 37 + std::hash<uint32_t>()(v.uv);
          hash = hash * 37 + std::hash<uint32_t>()(v.n);
          return hash;
        }
      };

      typedef std::unordered_map<ObjVertex, uint32_t, ObjVertexHash> VertexMap;

      ///// <summary>
      ///// Структура обмена данными с Python
      ///// </summary>
      //typedef struct _MD {

      //  bool has_error = false; // Если в структуре содердаться ошибки, а не корректный вариант расчёта
      //  char* str_error = NULL; // Если есть ошибка, то тут должно быть описание как она возникла

      //  int polygon_id = -1; // id обрабатываемого полигона. Мало ли при возврате клиенту пригодится, т.к. расчёт многопоточный. -1 - не обрабатывался. Должен быть 0 и больше.

      //  // if has_error = false
      //  int nn_verts = NULL; // Result mesh count of verts (по одному значению на объект)
      //  int nn_edges = NULL; // Result mesh count of edges (по одному значению на объект)
      //  int nn_faces = NULL; // Result mesh count of faces (по одному значению на объект)

      //  float* vertices = NULL; // Result mesh vertices
      //  int* edges = NULL; // Result mesh edges
      //  int* faces = NULL; // Result mesh faces

      //  // if has_error = true
      //  // Информация по сбойным контурам:
      //  int ftcs_count = 0; // failed triangulated contours - количество сбойных контуров
      //  int* ftcs_vertices_counter = NULL; // Список количеств вершин отдельно по контурам
      //  const char** ftcs_vertices_description = NULL; // Description для контура (предназначен для ошибок).
      //  float* ftcs_vertices_list = NULL; // Общий список вершин контуров
      //} MESH_DATA;

      struct MESH_DATA2 {

        bool has_error = false; // Если в структуре содердаться ошибки, а не корректный вариант расчёта
        char* str_error = NULL; // Если есть ошибка, то тут должно быть описание как она возникла

        int polygon_id = -1; // id обрабатываемого полигона. Мало ли при возврате клиенту пригодится, т.к. расчёт многопоточный. -1 - не обрабатывался. Должен быть 0 и больше.

        // if has_error = false
        int  nn_objects = 0; // Количество объектов результата (меняется в зависимости от join_mode)
        int* nn_objects_indexes = NULL; // Индексы объектов из тех, что были переданы (по одному индексу на объект)
        int* nn_offsets_counts = NULL; // Список количеств индексов offsets, которые были использованы для получения результатов по объектам
        int* nn_offsets_indexes = NULL; // Сами индексы offsets, которые были использованы для получения результатов по объектам

        int* nn_verts = NULL; // Result mesh count of verts (по одному значению на объект)
        int* nn_edges = NULL; // Result mesh count of edges (по одному значению на объект)
        int* nn_faces = NULL; // Result mesh count of faces (по одному значению на объект)

        float* vertices = NULL; // Result mesh vertices
        int* edges = NULL; // Result mesh edges
        int* faces = NULL; // Result mesh faces

        // nn_errors_count:
        //        nn_description1_per_error1
        //        nn_contours_per_error1
        //            nn_vertices_per_contour1
        //                vertices_of_errors
        // 
        int   nn_source_objects_count = 0;    // Количество ошибок возвращается в соответствии с исходными объектами, по которым возвращается набор ошибок при получении результата
        int* nn_source_objects_indexes = NULL; // Индексы исходных объектов
        int* nn_errors_count = NULL; // Количество ошибок на объект (по одному элементу на объект, а количество элементов соответствует nn_objects)
        char** nn_description1_per_error1 = NULL; // По одному описанию на ошибку (количество записей на объект соответствует)
        int* nn_contours_per_error1 = NULL; // Количество контуров на ошибку. Количество контуров в ошибках (1-один контур, >1 обычно обычно все контуры объекта. Для понимания интерпретации результата)
        int* nn_vertices_per_contour1 = NULL; // Количество вершин на ошибку. 
        float* vertices_of_errors = NULL; // Все вершины с ошибками (их количество должно соответствовать сумме nn_vertices_per_contour1 для всех объектов)


        //// if has_error = true
        //// Информация по сбойным контурам:
        //int ftcs_count = 0; // failed triangulated contours - количество сбойных контуров
        //int* ftcs_vertices_counter = NULL; // Список количеств вершин отдельно по контурам
        //const char** ftcs_vertices_description = NULL; // Description для контура (предназначен для ошибок).
        //float* ftcs_vertices_list = NULL; // Общий список вершин контуров
      };

      //static int check_if_canceled(float progress,
      //  void (*update_cb)(void*, float progress, int* cancel),
      //  void* update_cb_data) {
      //  int cancel = 0;
      //  update_cb(update_cb_data, progress, &cancel);
      //  return cancel;
      //}

      ///// <summary>
      ///// Скопировать контуры в список со сбойными контурами.
      ///// </summary>
      ///// <param name="mesh_data"></param>
      ///// <param name="contours"></param>
      //template <typename MD>
      //void copy_contours_to_failed(MD* mesh_data, std::map<int, std::vector<Point_2>>& contours) {
      //  mesh_data->ftcs_count = contours.size();
      //  mesh_data->has_error = true;

      //  if (mesh_data->ftcs_count > 0) {

      //    // Сохранить количества вершин сбойных контуров (пока ещё не вершины)
      //    mesh_data->ftcs_vertices_counter = (int*)malloc(sizeof(int) * mesh_data->ftcs_count);
      //    mesh_data->ftcs_vertices_description = (const char**)malloc(sizeof(char*) * mesh_data->ftcs_count);
      //    int count_of_vertices = 0;
      //    const char* description = "unknown error. This contour has status failed. Investigate shape.";
      //    for (int I = 0; I <= contours.size() - 1; I++) {
      //      if (contours.find(I) == contours.end()) {
      //        continue;
      //      }
      //      std::vector<Point_2> contour1 = contours[I];
      //      int ftcs_id = I;
      //      int ftcs_vertices_counter = contour1.size();
      //      mesh_data->ftcs_vertices_counter[I] = ftcs_vertices_counter;
      //      mesh_data->ftcs_vertices_description[I] = description;
      //      count_of_vertices += ftcs_vertices_counter;
      //    }

      //    // Запомнить все вершины в сбойных контурах:
      //    mesh_data->ftcs_vertices_list = (float*)malloc(sizeof(float[3]) * count_of_vertices);
      //    int J = 0;
      //    for (int I = 0; I <= contours.size() - 1; I++) {
      //      if (contours.find(I) == contours.end()) {
      //        continue;
      //      }
      //      std::vector<Point_2> contour1 = contours[I];
      //      for (auto v_it = contour1.begin(); v_it != contour1.end(); ++v_it) {
      //        Point_2 p2 = *v_it;
      //        mesh_data->ftcs_vertices_list[J * 3 + 0] = p2.x();
      //        mesh_data->ftcs_vertices_list[J * 3 + 1] = p2.y();
      //        mesh_data->ftcs_vertices_list[J * 3 + 2] = 0.0; // p.z(); - Point_2 no z
      //        J++;
      //      }
      //    }
      //  }
      //}



      //template <typename Point_2>
      //class ContoursItem {
      //public:
      //  std::vector<Point_2> contour;
      //  char* description = NULL;
      //  ContoursItem() {
      //  }
      //  ContoursItem(std::vector<Point_2> _contour, char* _description) {
      //    // Просто присвоить, чтобы потом освободить. Либо в общей функции freemem, либо в set_description
      //    description = _description;
      //    contour = _contour;
      //  }

      //  ContoursItem(std::vector<Point_2> _contour, const char* _description) {
      //    contour = _contour;
      //    set_description(_description);
      //  }

      //  ContoursItem(std::vector<Point_2> _contour, const char* _description, int _value) {
      //    contour = _contour;
      //    set_description(_description, _value);
      //  }

      //  void set_description(char* _description) {
      //    if (description != NULL) {
      //      free(description);
      //      description = NULL;
      //    }
      //    description = _description;
      //  }

      //  // Вывести только description:
      //  void set_description(const char* _description) {
      //    if (description != NULL) {
      //      free(description);
      //      description = NULL;
      //    }

      //    // скопировать значение и не запоминать исходный адрес.
      //    if (_description != NULL) {
      //      int sz = std::strlen(_description) + 1;
      //      if (sz > 1000) {
      //        sz = 1000;
      //      }
      //      char* buf = (char*)malloc(sizeof(char) * sz);
      //      snprintf(buf, sz, "%s", _description);
      //      description = buf;
      //    }
      //  }

      //  // Вывести description с некоторым значением int:
      //  void set_description(const char* _description, int value) {
      //    if (description != NULL) {
      //      free(description);
      //      description = NULL;
      //    }

      //    // скопировать значение и не запоминать исходный адрес.
      //    if (_description != NULL) {
      //      int sz = std::strlen(_description) + 20;
      //      if (sz > 1000) {
      //        sz = 1000;
      //      }
      //      char* buf = (char*)malloc(sizeof(char) * sz);
      //      snprintf(buf, sz, _description, value);
      //      description = buf;
      //    }
      //  }
      //};

      //template <typename Point_2>
      //class ExceptionWithFailedContours : public std::exception {
      //public:
      //  /// <summary>
      //  /// General description.
      //  /// </summary>
      //  char* message = NULL;
      //  std::vector<ContoursItem<Point_2>> items;
      //  explicit ExceptionWithFailedContours(char* _message) : std::exception() {
      //    message = _message;
      //  }
      //  explicit ExceptionWithFailedContours(const char* _message) : std::exception() {
      //    set_message(_message);
      //  }
      //  const char* what() const noexcept override {
      //    return message;
      //  }

      //  //Перезаписать общее сообщение
      //  void set_message(char* _message) {
      //    if (message != NULL) {
      //      free(message);
      //      message = NULL;
      //    }
      //    message = _message;
      //  }

      //  //Перезаписать общее сообщение
      //  void set_message(const char* _message) {
      //    if (message != NULL) {
      //      free(message);
      //      message = NULL;
      //    }

      //    // скопировать значение и не запоминать исходный адрес.
      //    if (_message != NULL) {
      //      int sz = std::strlen(_message) + 1;
      //      if (sz > 1000) {
      //        sz = 1000;
      //      }
      //      char* buf = (char*)malloc(sizeof(char) * sz);
      //      snprintf(buf, sz, "%s", _message);
      //      message = buf;
      //    }
      //  }
      //};

      typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
      typedef K::Point_2                    Point;
      typedef CGAL::Polygon_2<K>            Polygon_2;
      typedef CGAL::Polygon_with_holes_2<K> PolygonWithHoles;
      typedef std::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr;
      typedef std::vector<PolygonWithHolesPtr> PolygonWithHolesPtrVector;

      typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

      typedef CGAL::Polygon_2<Kernel>  Contour;
      typedef std::shared_ptr<Contour> ContourPtr;
      typedef std::vector<ContourPtr>  ContourSequence;

      typedef CGAL::Straight_skeleton_2<Kernel> Ss;
      typedef CGAL::Straight_skeleton_builder_traits_2<Kernel>       SsBuilderTraits;
      typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss> SsBuilder;

      typedef CGAL::Polygon_offset_builder_traits_2<Kernel>                    OffsetBuilderTraits;
      typedef CGAL::Polygon_offset_builder_2<Ss, OffsetBuilderTraits, Contour> OffsetBuilder;



      template <typename PolygonWithHoles,
        typename PointRange, typename FaceRange,
        typename NamedParameters
      >
      bool construct_horizontal_faces(const PolygonWithHoles& p,
        const FT altitude,
        PointRange& points,
        FaceRange& faces,
        const NamedParameters& np,
        const bool invert_faces = false,
        const bool verbose = false
      ) {

        using Default_kernel = typename Kernel_traits<Point>::type;
        using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
          NamedParameters,
          Default_kernel>::type;
        using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Geom_traits>;
        //typedef CGAL::Constrained_triangulation_face_base_2<K>            Fb;
        //typedef CGAL::Triangulation_data_structure_2<Vb, Fb>               TDS;
        //typedef CGAL::Exact_predicates_tag                                Itag;
        //typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
        //typedef CDT::Face_handle                                          Face_handle;
        //typedef CDT::Point                                                Point;
        //typedef std::pair<Face_handle, int> Edge;
        //using CDT_Vertex_handle = typename CDT::Vertex_handle;
        //using CDT_Face_handle = typename CDT::Face_handle;

        using Fb = CGAL::Constrained_triangulation_face_base_2<Geom_traits>;
        using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
        using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
        using CDT = CGAL::Constrained_Delaunay_triangulation_2<Geom_traits, TDS, Itag>;
        using CDT_Vertex_handle = typename CDT::Vertex_handle;
        using CDT_Face_handle = typename CDT::Face_handle;
        typedef std::pair<CDT_Face_handle, int> Edge;


        CDT cdt;

        //std::map<int, Polygon_2> invalid_contours;

//        try {
          bool was_polygon_exception = false;
          //if (verbose == true) {
          //  printf("\nTriangulating boundary: outer_boundary size: %u", (unsigned int)p.outer_boundary().size());
          //}
          //ExceptionWithFailedContours<Point_2> ex_failed_contours("");
          try {
            cdt.insert_constraint(p.outer_boundary().begin(), p.outer_boundary().end(), true /*close*/);
          } catch (const std::exception& ex) {
            was_polygon_exception = true;
            // general exception. Save contour with exception for future debug. No extrude skeleton after exception.
            //invalid_contours.insert({ -1, p.outer_boundary() });
            //ex_failed_contours.items[-1] = ContoursItem<Point_2>({ std::vector<Point_2>(p.outer_boundary().vertices_begin(), p.outer_boundary().vertices_end()), (char*)ex.what() });
            // Сохранить все остальные контуры тоже:
            int I = 0;
            for (auto h_it = p.holes_begin(); h_it != p.holes_end(); ++h_it, I++) {
              //invalid_contours.insert({ I, *h_it });
              //ex_failed_contours.items[I] = ContoursItem<Point_2>({ std::vector<Point_2>((*h_it).vertices_begin(), (*h_it).vertices_end()), (char*)ex.what() });
            }
            //throw ex_failed_contours; // no one contour can by added as hole if outer border is invalid
          }
          int I = 0;
          for (auto h_it = p.holes_begin(); h_it != p.holes_end(); ++h_it, I++) {
            Polygon_2 polygon = *h_it;
            try {
              cdt.insert_constraint(h_it->begin(), h_it->end(), true /*close*/);
            } catch (const std::exception& ex) {
              was_polygon_exception = true;
              //invalid_contours.insert({ I, *h_it });
              //ex_failed_contours.items[I] = ContoursItem<Point_2>({ std::vector<Point_2>((*h_it).vertices_begin(), (*h_it).vertices_end()), (char*)ex.what() });
            }
          }

          //if (ex_failed_contours.items.size() > 0) {
          //  throw ex_failed_contours;
          //} else {
          //}
        //} catch (ExceptionWithFailedContours<Point_2>& ex) {
        //  throw ex;
        //} catch (const std::exception& ex) {
        //  // По идее сюда не попасть. Стоит для страховки
        //  //printf("\n>>>>> You are not allowed be here. %s", ex.what());
        //  throw ex;
        //}

        //printf("\nContinue construct_horizontal_face...");

          if (was_polygon_exception == false) {

            std::size_t id = points.size(); // point ID offset (previous faces inserted their points)
            for (CDT_Vertex_handle vh : cdt.finite_vertex_handles()) {
              points.emplace_back(cdt.point(vh).x(), cdt.point(vh).y(), altitude);
              vh->info() = id++;
            }

#ifdef CGAL_SLS_DEBUG_DRAW
            // CGAL::draw(cdt);
#endif

            std::unordered_map<CDT_Face_handle, bool> in_domain_map;
            boost::associative_property_map< std::unordered_map<CDT_Face_handle, bool> > in_domain(in_domain_map);

            CGAL::mark_domain_in_triangulation(cdt, in_domain);

            for (CDT_Face_handle f : cdt.finite_face_handles()) {
              if (!get(in_domain, f))
                continue;

              int v0_idx = f->vertex(0)->info();
              int v1_idx = f->vertex(1)->info();
              int v2_idx = f->vertex(2)->info();
              Point_3 p3_0(points[v0_idx]);
              Point_3 p3_1(points[v1_idx]);
              Point_3 p3_2(points[v2_idx]);
              Point_2 p2_0(p3_0.x(), p3_0.y());
              Point_2 p2_1(p3_1.x(), p3_1.y());
              Point_2 p2_2(p3_2.x(), p3_2.y());
              FT tri_area = CGAL::area(p2_0, p2_1, p2_2);

              // invert faces for the z=0 plane (bottom face)
              if (invert_faces)
                faces.push_back({ f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info() });
              else
                faces.push_back({ f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info() });
            }
            //printf("\nExit construct_horizontal_face...");
          } else {
            //printf("\n construct_horizontal_faces: was_polygon_exception==true!!!!");
          }
          bool res = !was_polygon_exception;
          return res;
      }

      extern "C" {

        /// <summary>
        /// Применяется при топологической сортировке контуров
        /// </summary>
        struct SRT {
          Polygon_2 polygon;
          enum ContourType {
            OUTER, HOLE
          } ct;
          double area; // Площадь контура (независимо от типа внутренний или внешний контур площадь положительна)
          std::vector<std::vector<Point_2>> contours; // Вершины контуров. Первый контур - внешний. В него будут добавляться holes.
          SRT(Polygon_2 _polygon, ContourType _ct, double _area, std::vector<std::vector<Point_2>> _contours)
            :polygon(_polygon), ct(_ct), contours(_contours) {
            area = _area;
          }
        };

        /// <summary>
        /// Преобразовать PolygonWithHoles2 в Mesh (полигональный mesh)
        /// </summary>
        /// <param name="pwh1"></param>
        /// <param name="altitude"></param>
        /// <param name="sm1"></param>
        static bool ConvertPolygonWithHoles2IntoMesh(CGAL::Straight_skeleton_extrusion::internal::Polygon_with_holes_2& pwh1, double altitude, Mesh& sm1) {
          std::vector<Point_3> polygon_points;
          std::vector<std::vector<std::size_t> > polygon_faces;
          bool res = construct_horizontal_faces(pwh1, altitude, polygon_points, polygon_faces, CGAL::parameters::maximum_height((double)0.0), true);

          if (res == true) {
            namespace PMP = ::CGAL::Polygon_mesh_processing;
            PMP::merge_duplicate_points_in_polygon_soup(polygon_points, polygon_faces);
            if (!PMP::is_polygon_soup_a_polygon_mesh(polygon_faces)) {
              PMP::orient_polygon_soup(polygon_points, polygon_faces);
            }

            PMP::polygon_soup_to_polygon_mesh(polygon_points, polygon_faces, sm1);
          } else {

          }
          return res;
        }

        /// <summary>
        /// Преобразовать набор контуров в набор Mesh. Функция определяет вложенность контуров друг в друга и на основе этой вложенности строит полигоны с отверстиями. Потом из этих полигонов строит несколько Mesh.
        /// Boundaries и holes должны быть корректно направлены (boundaries - ccw, против часовой, hols - cw, по часовой.
        /// </summary>
        /// <param name="cs_need_sort"></param>
        /// <param name="vect_pwh"></param>
        static void ConvertContourSequenceIntoPolygonsWithHoles(
          CGAL::Straight_skeleton_extrusion::internal::ContourSequence& _cs_in,
          std::vector<CGAL::Straight_skeleton_extrusion::internal::Polygon_with_holes_2>& vect_pwh_out
        ) {
          // Отсортировать полигоны по площади 
          sort(_cs_in.begin(), _cs_in.end(), [](const ContourPtr& c1, const ContourPtr& c2) {return abs(c1->area()) > abs(c2->area()); });

          // Зассортировать контуры друг в друга:
          std::vector<ContourSequence> stack_cs; // vect_collect_for_polygons
          for (auto& c1 : _cs_in) {
            ContourSequence cs1;
            cs1.push_back(c1);
            stack_cs.push_back(cs1);
          }

          // Результаты сортировки складывать в этот вектор. Каждый остров в этом векторе будет иметь один внешний контур и
          // ноль или больше количество holes.
          std::vector<ContourSequence> vect_cs_contours_for_polygons;

          while (stack_cs.size() > 0) {
            for (int I = stack_cs.size() - 1; I >= 0; I--) {
              ContourSequence& vcfp_I = stack_cs[I];
              ContourPtr& vcfp_I0 = vcfp_I[0];
              if (vcfp_I0->is_counterclockwise_oriented() == true) {
                vect_cs_contours_for_polygons.push_back(std::move(vcfp_I));
                stack_cs.erase(stack_cs.begin() + I);
              } else {
                break;
              }
            }

            // Искать с конца в списке holes и найти в какой boundary они попадают:
            for (int I = stack_cs.size() - 1; I >= 0; I--) {
              ContourSequence& vcfp_I = stack_cs[I];
              ContourPtr& vcfp_I0 = stack_cs[I][0];
              // Если найдена hole то искать boundary
              if (vcfp_I0->is_clockwise_oriented() == true) {
                bool is_parent_hole_found = false;
                if (stack_cs.size() >= 2) {
                  Point_2 p0_I0 = vcfp_I0->vertices().at(0);
                  for (int IJ = stack_cs.size() - 2; IJ >= 0; IJ--) {
                    ContourSequence& vcfp_IJ = stack_cs[IJ];
                    ContourPtr& vcfp_IJ0 = stack_cs[IJ][0];
                    // Если boundary найдена, то проверить подходит ли эта hole в него? Если подходит, то вставить его в этот boundary.
                    if (vcfp_IJ0->is_counterclockwise_oriented() == true) {
                      // https://stackoverflow.com/questions/64046074/is-there-a-cgal-function-that-checks-if-a-point-is-inside-a-linear-polygon-with
                      // Если хотя бы одна точка находится внутри, то считаем. что и весь контур находится внутри (Предыдущие проверки должны били исключить пересечения до этого)
                      Contour p2(vcfp_IJ0.get()->begin(), vcfp_IJ0.get()->end()); // Не может быть clockwise(), т.к. этот контур и был выбран как counterclockwise()
                      if (CGAL::POSITIVE == CGAL::oriented_side(p0_I0, p2)) {
                        for (auto& c1 : vcfp_I) {
                          vcfp_IJ.push_back(std::move(c1));
                        }
                        stack_cs.erase(stack_cs.begin() + I);
                        is_parent_hole_found = true;
                        break;
                      }
                    }
                  }
                }
                if (is_parent_hole_found == false) {
                  // Такого быть не должно. Надо проверить.
                  vect_cs_contours_for_polygons.push_back(std::move(vcfp_I));
                  stack_cs.erase(stack_cs.begin() + I);
                }
              } else {
                // Не должно происходить, т.к. все не is_counterclockwise_oriented должны были быть удалены с конца стека на прошлом цикле, но не проверить нельзя.
                // Update - может происходить, когда это не первый цикл. Надо выйти, чтобы пойти в while на следующий круг и обработать такой контур.
                break;
              }
            }
          }

          for (auto& cs : vect_cs_contours_for_polygons) {
            Polygon_with_holes_2 pwh(Polygon_2(cs[0]->begin(), cs[0]->end()));
            for (auto contour_it = cs.begin() + 1; contour_it != cs.end(); ++contour_it) {
              auto& hole = *contour_it;
              pwh.add_hole(*hole.get());
            }
            vect_pwh_out.push_back(pwh);
          }
        }


        // Object Index, Offset, Altitude, Polygon. Сопоставленные objects, planes, offsets, altitudes и помнить порядок в котором заданы offsets:
        struct OIOA {
          int object_index;
          int offset_index;

          double offset;
          double altitude;
          // Этот тип имеет смысл только при описании полигона с отрицательным offset
          enum NEGATIVE_TYPE {
            WITH_EXTERNAL_FRAME,
            INVERTED_HOLE,
          } neg_type;

          // Если будет получен результат для заданного offset, то сюда будут записаны контуры по результатам расчёта.
          //std::vector<Contour> polygon_contours;
          ContourSequence polygon_contours;

          OIOA(int _object_index, int _offset_index, double _offset, double _altitude) {
            object_index = _object_index;
            offset_index = _offset_index;
            offset = _offset;
            altitude = _altitude;
            neg_type = NEGATIVE_TYPE::INVERTED_HOLE;
          }

          OIOA(int _object_index, int _offset_index, double _offset, NEGATIVE_TYPE _neg_type, double _altitude) {
            object_index = _object_index;
            offset_index = _offset_index;
            offset = _offset;
            neg_type = _neg_type;
            altitude = _altitude;
          }
        }; // OOAI

        // Список массив PWH и набор атрибутов OIOA (позже поменяю Altitude на Matrix)
        struct VECT_OIOA_PWHS {
          OIOA oioa1;

          std::vector<Polygon_with_holes_2> vect_pwh;
          VECT_OIOA_PWHS(OIOA& _oioa1, std::vector<Polygon_with_holes_2>& _vect_pwh)
            :oioa1(_oioa1), vect_pwh(_vect_pwh) {
          }
        };

        /// <summary>
        /// Конвертировать набор Polygon_with_holes_2 в наборы контуров в виде полигональной сетки vertices, edges, faces, но faces отсутствуют, но этот параметр нужен полигональной сетке.
        /// </summary>
        /// <param name="_vect_sm_in"></param>
        /// <param name="vect_object1_verts"></param>
        /// <param name="vect_object1_edges"></param>
        /// <param name="vect_object1_faces"></param>
        static void ConvertPWH2IntoVerticesEdgesNoFaces(
          std::vector<VECT_OIOA_PWHS> _in_vect_pwh_offset_index_order,
          std::vector<std::vector<float>>& _vect_object1_verts_out,
          std::vector<std::vector<unsigned int>>& _vect_object1_edges_out,
          std::vector<std::vector<unsigned int>>& _vect_object1_faces_out
        ) {
          //throw std::exception("not realized");

          int contour_cursor = 0;
          for (auto& pwh1_offset : _in_vect_pwh_offset_index_order) {
            //pwh1_offset.oioa1.altitude;
            std::vector<std::vector<Point_2>> cs;
            for (auto& pwh1 : pwh1_offset.vect_pwh) {
              cs.push_back(std::vector<Point_2>(pwh1.outer_boundary().begin(), pwh1.outer_boundary().end()));
              for (auto& hole1 : pwh1.holes()) {
                cs.push_back(std::vector<Point_2>(hole1.begin(), hole1.end()));
              }
            }

            unsigned int point_cursor = 0;
            for (auto& cs1 : cs) {
              std::vector<std::vector<float>> vect_points;
              std::vector<std::vector<unsigned int>> vect_edges;
              unsigned int p_index = 0;
              for (auto& p1 : cs1) {
                vect_points.push_back({ (float)p1.x(), (float)p1.y(), (float)pwh1_offset.oioa1.altitude });
                vect_edges.push_back({ contour_cursor + point_cursor + p_index, contour_cursor + point_cursor + p_index + 1 });
                p_index++;
              }
              vect_edges[vect_edges.size() - 1][1] = contour_cursor + point_cursor + 0;
              point_cursor += vect_edges.size();

              _vect_object1_verts_out.insert(_vect_object1_verts_out.end(), vect_points.begin(), vect_points.end());
              _vect_object1_edges_out.insert(_vect_object1_edges_out.end(), vect_edges.begin(), vect_edges.end());
            }
            contour_cursor += point_cursor;
          }
        }

        /// <summary>
        /// Конвертировать набор Meshes в наборы vertices, edges, faces
        /// </summary>
        /// <param name="_vect_sm_in"></param>
        /// <param name="vect_object1_verts"></param>
        /// <param name="vect_object1_edges"></param>
        /// <param name="vect_object1_faces"></param>
        static void ConvertMeshesIntoVerticesEdgesFaces(
          std::vector<Mesh> _vect_sm_in,
          std::vector<std::vector<float>>& vect_object1_verts,
          std::vector<std::vector<unsigned int>>& vect_object1_edges,
          std::vector<std::vector<unsigned int>>& vect_object1_faces
        ) {
          unsigned int start_shape_vert_index = 0;
          // Если равно 0, то результатов расчёта нет и тогда и не нужно считать mesh, т.к. вместе с ними считается и set_object1_offsets, через который позже считается mesh_data->nn_offsets_counts, который должен быть равен 0, а не количеству oioa
          if (_vect_sm_in.size() > 0) {
            for (auto& sm1 : _vect_sm_in) {
              // Загрузить вершины sm1
              for (auto& p1 : sm1.points()) {
                vect_object1_verts.push_back(std::vector<float>{ (float)p1.x(), (float)p1.y(), (float)p1.z() });
              }

              for (const auto& edge : sm1.edges()) {
                CGAL::SM_Halfedge_index e = edge.halfedge();
                CGAL::SM_Vertex_index vd1 = source(e, sm1);
                CGAL::SM_Vertex_index vd2 = target(e, sm1);
                vect_object1_edges.push_back({ start_shape_vert_index + vd1, start_shape_vert_index + vd2 });
              }

              {
                for (Mesh::Face_index face_index : sm1.faces()) {
                  // Определение индексов faces: https://stackoverflow.com/questions/46808246/cgal-get-face-data-from-surface-mesh
                  CGAL::Vertex_around_face_circulator<Mesh> vcirc(sm1.halfedge(face_index), sm1), done(vcirc);
                  std::vector<int> verts_indices;
                  do {
                    verts_indices.push_back(*vcirc++);
                  } while (vcirc != done);
                  // Reverse faces normals:
                  vect_object1_faces.push_back({ start_shape_vert_index + verts_indices[2], start_shape_vert_index + verts_indices[1], start_shape_vert_index + verts_indices[0] });
                }
              }
              start_shape_vert_index += sm1.vertices().size();
            }
          }
        }

        /// <summary>
        /// Содержание ошибки
        /// </summary>
        struct ObjectError {
          int object_index{ -1 };
          // Набор контуров, которые содержат ошибку. Это может быть как один контур, так и вся фигура из нескольких контуров.
          ContourSequence cs1;
          // Текстовое описание ошибки
          char* message{ NULL };

          void set_message(char* _message) {
            if (message != NULL) {
              free(message);
              message = NULL;
            }
            message = _message;
          }

          //Перезаписать общее сообщение
          void set_message(const char* _message) {
            if (message != NULL) {
              free(message);
              message = NULL;
            }

            // скопировать значение и не запоминать исходный адрес.
            if (_message != NULL) {
              int sz = std::strlen(_message) + 1;
              if (sz > 1000) {
                sz = 1000;
              }
              char* buf = (char*)malloc(sizeof(char) * sz);
              snprintf(buf, sz, "%s", _message);
              message = buf;
            }
          }

          ObjectError(int _object_index, std::vector<Point_2> _points, char* _message)
            :object_index{ _object_index } {
            cs1.push_back(std::make_shared<Contour>(_points.begin(), _points.end()));
            set_message(_message);
          }

          ObjectError(int _object_index, std::vector<Point_2> _points, const char* _message)
            :object_index{ _object_index } {
            cs1.push_back(std::make_shared<Contour>(_points.begin(), _points.end()));
            set_message(_message);
          }

          ObjectError(int _object_index, ContourPtr c1, char* _message)
            :object_index{ _object_index } {
            cs1.push_back(c1);
            set_message(_message);
          }

          ObjectError(int _object_index, ContourPtr c1, const char* _message)
            :object_index{ _object_index } {
            cs1.push_back(c1);
            set_message(_message);
          }

          ObjectError(int _object_index, Contour c1, char* _message)
            :object_index{ _object_index } {
            cs1.push_back(std::make_shared<Contour>(c1.begin(), c1.end()));
            set_message(_message);
          }

          ObjectError(int _object_index, Contour c1, const char* _message)
            :object_index{ _object_index } {
            cs1.push_back(std::make_shared<Contour>(c1.begin(), c1.end()));
            set_message(_message);
          }

          ObjectError(int _object_index, ContourSequence _cs1, char* _message)
            :object_index{ _object_index } {
            cs1.assign(_cs1.begin(), _cs1.end());
            set_message(_message);
          }

          ObjectError(int _object_index, ContourSequence _cs1, const char* _message)
            :object_index{ _object_index } {
            cs1.assign(_cs1.begin(), _cs1.end());
            set_message(_message);
          }

          ObjectError(int _object_index, const Polygon_with_holes_2& _pwh, char* _message)
            :object_index{ _object_index } {
            cs1.push_back(std::make_shared<Contour>(_pwh.outer_boundary().begin(), _pwh.outer_boundary().end()));
            for (auto& h1 : _pwh.holes()) {
              cs1.push_back(std::make_shared<Contour>(h1.begin(), h1.end()));
            }
            set_message(_message);
          }

          ObjectError(int _object_index, const Polygon_with_holes_2& _pwh, const char* _message)
            :object_index{ _object_index } {
            cs1.push_back(std::make_shared<Contour>(_pwh.outer_boundary().begin(), _pwh.outer_boundary().end()));
            for (auto& h1 : _pwh.holes()) {
              cs1.push_back(std::make_shared<Contour>(h1.begin(), h1.end()));
            }
            set_message(_message);
          }
        };

        /// <summary>
        /// CGAL Bidirection offsets of Polygon_2 with holes. Многопоточность сначала считает в параллель Straight Skeleton, затем offsets.
        /// Перед рассчётом производится подготовка исходных данных - все объекты деляться на две группы: с отрицательными offsets и положительными offsets.
        /// Отрицательные offsets предварительно обрабатываются - инвертируются holes, а внешние контура самых внешних объектов объединяются в несколько
        /// holes и производится расчёт внешнего frame, в который они впишут общий offset. Потом этот frame объединяется с holes и записывается обратно в
        /// рассчётные объекты, но уже с новыми параметрами. После этого производится расчёт в параллель всех Skeleton-ов, потом производится расщепление на
        /// цикла на рассчёт независимых offset для каждого Straight Skeleton. Результат возвращается пользователю в виде meshes (Edges или Faces), соответствующих исходным объектам.
        /// <image url="..\code_images\file_0031.png" scale="1.0"/>
        /// </summary>
        /// <param name="in_polygon_id"></param>
        /// <param name="in_offset"></param>
        /// <param name="in_altitude"></param>
        /// <param name="in_count_of_contours"></param>
        /// <param name="in_lens_of_contours"></param>
        /// <param name="in_vertices"></param>
        /// <param name="only_tests_for_valid"></param>
        /// <param name="res_type"></param>
        /// <param name="verbose"></param>
        /// <returns></returns>
        DLLEXPORT MESH_DATA2* straight_skeleton_2d_offset(
          // К чему относятся входные данные. <image url="..\code_images\file_0012.png" scale="1.0"/>
          // Контуры объектов считаются независимо друт от друга
          // Контуры планов могут влиять друг на друга при отрицательных offset-ах. Контуры с положительными offset будут рассчитываться независимо.
          // Также при отрицательных offset инвертируются holes и рассчитываются как положительные offset. Требуется избегать расчётов Straight Skeleton-ов
          // при вложенных контурах, т.к. это приводит к зависанию расчёта (есть такой глюк). Контуры в Straight Skeleton добавляются независимо,
          // поэтому можно случайно добавить контуры внутренних объектов (что с типами Polygon_with_holes_2 не происходит, т.к. у них есть деление на outer_boundary
          // и holes.
          int in_count_of_objects,        // Количество независимых исходных объектов
          // Shape mode: <image url="..\code_images\file_0020.png" scale=".1"/>
           int in_shapes_mode[],          // Тип обработки контуров объекта: 0-FULL_MODE (полная обработка), 1-EXCLUDE_HOLES (только внешние контуры, отверстия не учитываются), 2-INVERT_HOLES (отключить внешние контуры)
           int in_count_of_offsets[],     // Количества офсетов, которые надо посчитать (размер массива соответствует in_count_of_objects)
         float in_offsets[],              // Массив value offsets, которые надо посчитать
           int in_count_of_altitudes[],   // Количества высот, которые надо учитывать при расстановке результатов. (размер массива соответствует in_count_of_objects)
         float in_altitudes[],            // Массив value offsets, которые надо посчитать
           int in_count_of_planes[],      // Количество планов в объектах (размер массива соответствует in_count_of_objects)
           int in_contours_in_planes[],   // Количество контуров в планах. Первый контур плана всегда outer, остальные - holes. Вложенность планов одного объекта в упаковке не учитывается.
           int in_vertices_in_contours[], // Количество вершин в контуре.
         float in_vertices[][3],          // вершины контуров (являются одной цепочкой, которую надо поделить на части из параметра lens_of_contours
          bool only_tests_for_valid,      // Проверка всех контуров на валидность без расчёта. (true - только проверить на валидность, false - проверить на валидность и посчитать)
           int in_res_type,               // Тип результата. 0 - контуры, 1 - mesh
          bool force_z_zero,              // Обязательное преобразование Z=0 при любых обстоятельствах (иногда при редактировании модели можно двигать курсор в плоскости XY b получать Z "слегка" отличный от 0).
        // <image url="..\code_images\file_0030.png" scale=".3"/>
           int source_objects_join_mode,    // Предобработка исходных Meshes: # 0 - split - разделить на отдельные независимые-объекты, 1 - keep - оставить как есть, 2 - merge all meshes - объеденить в один большой объект.
           int results_join_mode,           // Что сделать с полученными mesh: # 0 - split - разделить на отдельные mesh-объекты, 1 - keep - оставить в своих объектах, 2 - merge all meshes - объеденить в один большой объект. В этом ноде дополнительные параметры к объектам не нужны
          bool verbose                     // Подробный вывод
        ) {
          if (verbose == true) {
#ifdef _DEBUG
            printf("\nstraight_skeleton_2d_offset in DEBUG MODE");
#endif
            printf("\nstraight_skeleton_2d_offset_internal");
            printf("\nmode verbose is on");
          }

          // TODO:: Переделать altitude на преобразование Matrix
          // https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Aff__transformation__3.html
          CGAL::Aff_transformation_3<Kernel> ac3(
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            //0.0, 0.0, 0.0, 1.0,
            1.0
          );

          MESH_DATA2* mesh_data = new MESH_DATA2();

          CGAL::Real_timer timer;
          timer.start();

          // Типы результатов, которые будут возвращены пользователю после расчётов - контуры или Faces
          // <image url="..\code_images\file_0027.png" scale=".2"/>
          enum ResType {
            EDGES, FACES
          }result_type;

          if (in_res_type == ResType::EDGES) {
            result_type = ResType::EDGES;
          } else { //if (in_res_type == ResType::FACES) {
            result_type = ResType::FACES;
          }
          //else {
          //  // TODO: Сообщить об ошибке
          //  throw std::exception("Unknows result_type. Allowed only 0 (for edges) or 1 (for faces)");
          //}

          // Параметры одного исходного полигонального объекта
          struct ClsObject {
            // Исходный индекс объекта
            int object_index;
            // Список offsets и altitudes
            std::vector<OIOA> vect_oioa;
            // Контуры исходных планов
            std::vector<ContourSequence> source_planes;
            // Полигоны для расчётов (позиции полигонов могут отличаться от позиций контуров в source_planes)
            std::vector<std::shared_ptr<Polygon_with_holes_2>> vect_polygons_with_holes;

            ClsObject(int _object_index, std::vector<OIOA> _vect_ooai, std::vector<ContourSequence> _source_planes)
              :source_planes(_source_planes), vect_oioa(_vect_ooai) {
              object_index = _object_index;
            }
          };

          std::vector<OIOA> res_contours; // Результаты расчёта всех контуров. Результаты будут сортироваться с учётом вложенности контуров по offset_index, в которой offset подавались на вход (по объектам) и с учётом параметра join_mode.
          std::vector<ObjectError> res_errors; // Собирать ошибки в объектах при загрузке исходных данных и когда попытка использовать контур привела к его исключению из расчёта. Относится только к исходным объектам.

          int v_pos = 0;
          int general_count_of_vertices = 0;
          std::map<int, std::vector<Point_2>> contours;
          std::vector<ClsObject> objects; // Информация по объектам будет нужна и после расчёта (нужно помнить начальные индексы входных объектов, потому что иногда объект полностью выпадает из расчёта, а нужно передать инфу, что в объекте ничего нет и его ошибки
          //try 
          {
            {
              double intersection_timer = 0;
              // Загрузка исходных данных по объектам с некоторыми проверками.
              int plane_cursor = 0;
              int contour_cursor = 0;
              int vert_cursor = 0;
              int offset_cursor = 0;

              for (int I = 0; I <= in_count_of_objects - 1; I++) {
                std::vector<ContourSequence> planes;
                int count_of_planes1 = in_count_of_planes[I];
                for (int IJ = 0; IJ <= count_of_planes1 - 1; IJ++, plane_cursor++) {
                  ContourSequence plane1;
                  for (int IJK = 0; IJK <= in_contours_in_planes[plane_cursor] - 1; IJK++, contour_cursor++) {
                    ContourPtr contour1 = std::make_shared<Contour>();
                    for (int IJKL = 0; IJKL <= in_vertices_in_contours[contour_cursor] - 1; IJKL++, vert_cursor++) {
                      bool is_z_not_zero = false;
                      if (force_z_zero == false) {
                        if (is_zero(in_vertices[vert_cursor][2]) == false) {
                          // Проверка if Z==0
                          is_z_not_zero = true;
                        }
                      }
                      Point_2 p2(in_vertices[vert_cursor][0], in_vertices[vert_cursor][1]);
                      contour1->push_back(p2);
                    }
                    if (contour1->size() < 3) {
                      // Контур не должен содержать в себе меньше трёх точек. Запомнить некорректный контур и исключиться его из расчёта
                      res_errors.push_back(ObjectError(I, contour1, VAL2STR(Err::_0009)". Polygon contains less 3 points."));
                      continue;
                    }
                    if (contour1->size() > 0) {
                      plane1.push_back(contour1);
                    }
                  }
                  if (plane1.size() > 0) {
                    // Проверить plane, который получился при загрузке нескольких контуров.
                    // Его внешний контур (boundary) и все holes должены быть is_simple. Если boundary не is_simple, то такой plane1 полностью исключается.
                    // Если holes не is_simple, то он исключается.
                    // Если что-то из внутренних holes пересекается друг с другом, то первый такой hole отменяется.
                    auto& boundary = plane1[0];
                    Polygon_2 p_boundary(boundary->begin(), boundary->end());
                    if (boundary->is_simple() == false) {
                      res_errors.push_back(ObjectError(I, plane1[0], VAL2STR(Err::_0010)". Boundary polygon is not simple."));
                      continue;
                    } else {
                      // Проверить что holes is_simple.
                      ContourSequence allowed_contours;

                      for (auto c1_it = plane1.begin(); c1_it != plane1.end(); ++c1_it) {
                        if (c1_it == plane1.begin()) {
                          continue;
                        }

                        Contour c1(c1_it->get()->begin(), c1_it->get()->end());
                        // Проверку на is_simple надо делать раньше, чем is_clockwise_oriented (в отладке тоже проверяется is_clockwise_oriented и, если нет, то падает с исключением ).
                        if (c1_it->get()->is_simple() == false) {
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Err::_0010)". Hole is not simple."));
                          continue;
                        }
                        // Hole должен быть ориентирован по часовой стрелке.
                        if (c1_it->get()->is_clockwise_oriented() == true) {
                          c1.reverse_orientation();
                        }
                        Polygon_2 p_hole1 = Polygon_2(c1.begin(), c1.end());

                        // Проверить, что все точки p_hole1 находятся внутри boundary. (есть исключение, которое рассматривается немного ниже, но точку снаружи всё ранво надо проверять. Также это может быть быстрее, чем проверять пересечения.)
                        bool is_hole_inside = true;
                        for (auto& hole1_point : p_hole1) {
                          // https://stackoverflow.com/questions/64046074/is-there-a-cgal-function-that-checks-if-a-point-is-inside-a-linear-polygon-with
                          if (CGAL::POSITIVE != CGAL::oriented_side(hole1_point, p_boundary)) {
                            is_hole_inside = false;
                            break;
                          }
                        }
                        if (is_hole_inside == false) {
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Err::_0023)". Hole is not inside boundary."));
                          continue;
                        }

                        // В продолжении проверки что точки holes находятся внутри объекта. Исключение, когда точки находятся внутри, но сами контуры пересекаются:
                        // https://stackoverflow.com/questions/65021923/is-there-a-cgal-function-for-finding-all-intersecting-points-between-a-2d-line
                        // https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
                        // <image url="..\code_images\file_0028.png" scale=".2"/>
                        // Пример отображения такого пересечения <image url="F:\Enternet\2024\24.11.10\SVCGAL\ctSVCGAL\code_images\file_0033.png" scale=".2"/>
                        CGAL::Real_timer timer1;
                        timer1.start();

                        std::vector<Point_2> c_intersect;
                        for (auto& boundary_edge_it = p_boundary.edges_begin(); boundary_edge_it != p_boundary.edges_end(); ++boundary_edge_it) {
                          for (auto& hole_edge_it = p_hole1.edges_begin(); hole_edge_it != p_hole1.edges_end(); ++hole_edge_it) {
                            const auto res = CGAL::intersection(*boundary_edge_it, *hole_edge_it);
                            if (res) {
                              if (const Segment_2* s = std::get_if<Segment_2>(&*res)) {
                                c_intersect.push_back(s->point(0));
                                c_intersect.push_back(s->point(1));
                              } else {
                                const Point_2* p = std::get_if<Point_2 >(&*res);
                                c_intersect.push_back(*p);
                              }
                            }
                          }
                        }
                        if (c_intersect.size() > 0) {
                          res_errors.push_back(ObjectError(I, c_intersect, VAL2STR(Err::_0030)". Hole intersects with boundary."));
                          continue;
                        }
                        timer1.stop();
                        intersection_timer += timer1.time();

                        // Проверить, что hole не пересекается с другими holes: <image url="..\code_images\file_0029.png" scale=".2"/>
                        bool is_objects_intersects = false;
                        try {
                          for (auto& ap1 : allowed_contours) {
                            Polygon_2 pap1(ap1.get()->begin(), ap1.get()->end());
                            pap1.reverse_orientation(); // Иногда такая проверка выдаёт исключение, но пока не знаю почему, поэтому она обёрнута в try/catch
                            is_objects_intersects = CGAL::do_intersect(p_hole1, pap1);
                            if (is_objects_intersects == true) {
                              break;
                            }
                          }
                        } catch (std::exception _ex) {
                          // Пока буду считать, что при любых неудачных попытках расчёта пересечения они пересекаются.
                          is_objects_intersects = true;
                        }
                        if (is_objects_intersects == true) {
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Err::_0024)". Hole intersect with hole."));
                          continue;
                        }

                        // Дополнительных проверок по пересечению внутренних контуров дополнительно пока не делается.
                        // Запомнить корректные контура.
                        allowed_contours.push_back(*c1_it);
                      }
                      // Добавить Boundary контур в начало списка
                      allowed_contours.insert(allowed_contours.begin(), boundary);
                      // Если имеются допустимые контуры, то загрузить их в план
                      if (allowed_contours.size() > 0) {
                        plane1.clear();
                        plane1.assign(allowed_contours.begin(), allowed_contours.end());
                      }
                    }
                    // Если в плане что-то есть, то загрузить его в массив допустимых планов для текущего object1
                    if (plane1.size() > 0) {
                      planes.push_back(plane1);
                    }
                  }
                }
                // Даже если plane нет, но нужно загрузить сопутствующие параметры offsets для объекта:
                {
                  int in_count_of_offsets1 = in_count_of_offsets[I];
                  std::vector<OIOA> vect_ooai;
                  for (int IJ = 0; IJ <= in_count_of_offsets1 - 1; IJ++) {
                    float in_offset1 = in_offsets[offset_cursor];
                    float in_altitude1 = in_altitudes[offset_cursor];
                    vect_ooai.push_back(OIOA(I, IJ, in_offset1, in_altitude1));
                    offset_cursor++;
                  }
                  objects.push_back(ClsObject(I, vect_ooai, planes));
                }
              }
              printf("\n" VAL2STR(Err::_0032) ". SS 2D Offset. Intersections boundary holes time: %.5g", intersection_timer);
            }

            {
              // Представить список объектов в зависимости от параметра source_objects_join_mode
              std::vector<ClsObject> vect_split_mode;
              ClsObject objects_merge_mode(0, std::vector<OIOA>(), std::vector<ContourSequence>());
              int object_index = 0;
              // Continue with good contours to load into Polygon_with_holes_2:
              for (auto& object1 : objects) {
                for (auto& plane1 : object1.source_planes) {
                  if (plane1.size() > 0) {

                    if (plane1[0]->is_clockwise_oriented()) {
                      plane1[0]->reverse_orientation();
                    }
                    std::shared_ptr<Polygon_with_holes_2> ppwh(new Polygon_with_holes_2(Polygon_2(plane1[0]->begin(), plane1[0]->end())));

                    for (int I = 1; I <= plane1.size() - 1; I++) {
                      if (plane1[I]->is_counterclockwise_oriented() == true) {
                        plane1[I]->reverse_orientation();
                      }
                      ppwh->add_hole(*plane1[I]);
                    }
                    //bool is_valid_pwh = true;
                    //try {
                    //  if (ppwh.get()->has_holes() == false) {
                    //    // https://github.com/CGAL/cgal/issues/5278
                    //    // По факту проверяется, что полигон без holes простой и замкнутый. Замкнутый тут при загрузке, простой проверяется раньше,
                    //    // поэтому полигон без holes можно пропустить.
                    //    //typedef Exact_predicates_exact_constructions_kernel kernel_type;
                    //    //typedef Polygon_2<kernel_type> polygon_type;
                    //    //bool valid = CGAL::is_valid_polygon(Polygon_2(ppwh.get()->outer_boundary().begin(), ppwh.get()->outer_boundary().end() ), typename Polygon_set_2<kernel_type>::Traits_2());
                    //  } else {
                    //    // Неофициальный метод проверки валидности полигона с holes: https://github.com/CGAL/cgal/issues/3667
                    //    // Может дать ложные срабатывания в E:\Enternet\2024\24.09.27\cgal\Arrangement_on_surface_2\include\CGAL\Arr_segment_traits_2.h:611 (2 и 5)
                    //    using PolygonSet = CGAL::Polygon_set_2<K>;
                    //    auto traits = std::make_unique<PolygonSet::Traits_2>();
                    //    // Иногда в режиме отладки эта функция выдаёт исключение.
                    //    is_valid_pwh = is_valid_polygon_with_holes(*ppwh.get(), *traits);
                    //  }
                    //} catch (std::exception _ex) {
                    //  is_valid_pwh = false;
                    //}
                    //if (is_valid_pwh == false) {
                    //  /* A valid polygon with holes is :
                    //   * 1 - Has empty or closed boundary and all the holes are closed
                    //   * 2 - The PWH is relatively simple polygon (holes are simple...)
                    //   * 3 - Has it's boundary oriented counterclockwise and the holes oriented
                    //   *     clockwise
                    //   * 4 - All the segments (boundary and holes) do not cross or intersect in their
                    //   *     relative interior
                    //   * 5 - The holes are on the interior of the boundary polygon if the boundary
                    //   *     is not empty
                    //   */
                    //  res_errors.push_back(ObjectError(object1.object_index, *ppwh.get(), VAL2STR(Err::_0022) ". Holes of the Polygon intersect amongst themselves or with outer boundary ()"));
                    //} else 
                    if (source_objects_join_mode == 0 /* 0 - split */) {
                      std::vector<ContourSequence> vect_cs;
                      vect_cs.push_back( plane1 );
                      std::vector<OIOA> vect_oioa1;
                      for (auto& oioa1 : object1.vect_oioa) {
                        OIOA _oioa1(object_index, oioa1.offset_index, oioa1.offset, oioa1.altitude);
                        vect_oioa1.push_back(_oioa1);
                      }
                      auto obj1 = ClsObject(object_index, vect_oioa1, vect_cs );
                      obj1.vect_polygons_with_holes.push_back(ppwh);
                      vect_split_mode.push_back(obj1);
                      object_index++;
                    } else if (source_objects_join_mode == 1 /* 1 - Keep source */) {
                      object1.vect_polygons_with_holes.push_back(ppwh);
                    } else /* 2 - merge */ {
                      if (objects_merge_mode.vect_oioa.size() == 0) {
                        // При merge использовать только параметры первого объекта. После загрузки данных oioa первого объекта остальные данные объектов oioa не используются
                        objects_merge_mode.vect_oioa.assign(object1.vect_oioa.begin(), object1.vect_oioa.end());
                      }
                      objects_merge_mode.source_planes.push_back(plane1);
                      objects_merge_mode.vect_polygons_with_holes.push_back(ppwh);
                    }
                  }
                }
              }
              if (source_objects_join_mode == 0 /* 0 - split */) {
                objects.clear();
                objects = vect_split_mode;
              } else if (source_objects_join_mode == 1 /* 1 - Keep source */) {
                // none. Use "objects"
              } else /* 2 - merge */ {
                objects.clear();
                objects.push_back(objects_merge_mode);
              }
              if (verbose == true) {
                printf("\n" VAL2STR(Err::_0028) ". SS 2D Offset. Source objects join mode: %i, count of objects: %zu. ", source_objects_join_mode, objects.size() );
              }
            }


            // Проверить, что теперь полигоны каждого объекта между собой не пересекаются (раньше были только holes отдельных полигонов):
            for (auto& object1 : objects) {
              {
                std::set<int> set_intersected_polygons_indexes;
                for (int I = object1.vect_polygons_with_holes.size() - 1; I >= 1; I--) {
                  const Polygon_with_holes_2& pwh1 = *object1.vect_polygons_with_holes[I].get();
                  for (int IJ = I - 1; IJ >= 0; IJ--) {
                    if (set_intersected_polygons_indexes.find(IJ) != set_intersected_polygons_indexes.end()) {
                      continue;
                    }
                    const Polygon_with_holes_2& pwh2 = *object1.vect_polygons_with_holes[IJ].get();
                    bool is_objects_intersects = false;
                    try {
                      is_objects_intersects = CGAL::do_intersect(pwh1, pwh2);
                    } catch (std::exception _ex) {
                      // Пока буду считать, что при любых неудачных попытках расчёта пересечения они пересекаются.
                      is_objects_intersects = true;
                    }
                    if (is_objects_intersects == true) {
                      // Удалить оба объекта как некорректные и запомнить их vertices для отображения как ошибки.
                      set_intersected_polygons_indexes.insert(I);
                      set_intersected_polygons_indexes.insert(IJ);

                      // Лучше бы как-то отметить, что полигон пересекается с другим полигоном, чем удалять полигон.
                      // TODO: подумать, как это лучше сделать (в дальнейшем этот параметр надо будет учитывать).
                      // Update: всё равно исходная позиция полигона вроде как не определяема и в разных системах может меняться, поэтому можно оставить пока так)
                      // <image url="..\code_images\file_0019.png" scale="1.0"/>
                      if (verbose == true) {
                        printf("\n" VAL2STR(Err::_0015) ". Polygons (%u, %u) INTERSECTS. Both marked as error", I, IJ);
                      }
                      break;
                    }
                  }
                }
                std::vector<int> vect_intersected_polygon_indexes;
                for (auto& idx : set_intersected_polygons_indexes) {
                  vect_intersected_polygon_indexes.push_back(idx);
                }
                sort(vect_intersected_polygon_indexes.begin(), vect_intersected_polygon_indexes.end(), [](const int v1, const int& v2) {return v1 > v2; });
                for (auto& I : vect_intersected_polygon_indexes) {
                  res_errors.push_back(ObjectError(object1.object_index, *object1.vect_polygons_with_holes[I].get(), VAL2STR(Err::_0011)". Polygons intersects"));
                  object1.vect_polygons_with_holes.erase(object1.vect_polygons_with_holes.begin() + I);
                }
              }
            }

            if (verbose) {
              printf("\nSS Offset objects to process: %zu", objects.size());
            }

            {
              int object1_index = 0;
              // Привести контуры объекта в соответствии с настройкой shape_mode.
              // Тип обработки контуров объекта: 
              // 0-ORIGINAL_BOUNDARIES (полная обработка), 
              // 1-EXCLUDE_HOLES (Keep only outer boundary), 
              // 2-INVERT_HOLES (Exclude outer boundary and fill holes).
              // <image url="..\code_images\file_0026.png" scale=".1"/>
              for (auto& object1 : objects) {
                int shapes_mode1 = in_shapes_mode[object1_index];
                if (shapes_mode1 == 0) {
                  // Оставить исходные объекты без изменений
                } else {
                  // Обработать их в соответствии с настройкой shapes_mode1

                  ContourSequence contours_sorted_by_area;
                  for (auto& pwh0 : object1.vect_polygons_with_holes) {
                    contours_sorted_by_area.push_back(std::make_shared<Contour>(pwh0->outer_boundary().begin(), pwh0->outer_boundary().end()));
                    for (auto& hole1 : pwh0->holes()) {
                      contours_sorted_by_area.push_back(std::make_shared<Contour>(hole1.begin(), hole1.end()));
                    }
                  }
                  // Сортировку по площади выполнить независимо от направления контура:
                  sort(contours_sorted_by_area.begin(), contours_sorted_by_area.end(), [](const ContourPtr& p1, const ContourPtr& p2) {return abs(p1->area()) > abs(p2->area()); });
                  std::vector<ContourSequence> stack_cs; // стек для выполнения сортировки вложения контуров.
                  for (auto& c1 : contours_sorted_by_area) {
                    ContourSequence cs1;
                    cs1.push_back(c1);
                    stack_cs.push_back(cs1);
                  }

                  std::vector<ContourSequence> external_contours;
                  std::vector<ContourSequence> vect_cs_inverted_holes;
                  while (stack_cs.size() > 0) {
                    // Перенести последние контура, являющиеся отверстиями, в список результирующих контуров (они могут содержать в себе бывшие внешние контура, кто-то уже добавился)
                    // <image url="..\code_images\file_0004.png" scale=".5"/>
                    for (int I = stack_cs.size() - 1; I >= 0; I--) {
                      ContourSequence& vcfp_I = stack_cs[I];
                      ContourPtr& vcfp_I0 = vcfp_I[0];
                      if (vcfp_I0->is_clockwise_oriented() == true) {
                        vect_cs_inverted_holes.push_back(std::move(vcfp_I));
                        stack_cs.erase(stack_cs.begin() + I);
                      } else {
                        break;
                      }
                    }

                    // Если последний контур является Outer Boundary, то поискать в кого он может быть вложен.
                    // <image url="..\code_images\file_0005.png" scale="1.0"/>
                    for (int I = stack_cs.size() - 1; I >= 0; I--) {
                      ContourSequence& vcfp_I = stack_cs[I];
                      ContourPtr& vcfp_I0 = stack_cs[I][0];
                      // В OUTER контуры ничего добавлять на нужно. Только они добавляются в HOLE-контуры.
                      if (vcfp_I0->is_counterclockwise_oriented() == true) {
                        bool is_parent_hole_found = false;
                        if (stack_cs.size() > 2) {
                          Point_2 p0_I0 = vcfp_I0->vertices().at(0);
                          for (int IJ = stack_cs.size() - 2; IJ >= 0; IJ--) {
                            ContourSequence& vcfp_IJ = stack_cs[IJ];
                            ContourPtr& vcfp_IJ0 = stack_cs[IJ][0];
                            if (vcfp_IJ0->is_clockwise_oriented() == true) {
                              // https://stackoverflow.com/questions/64046074/is-there-a-cgal-function-that-checks-if-a-point-is-inside-a-linear-polygon-with
                              // Если хотя бы одна точка находится внутри, то считаем. что и весь контур находится внутри (Предыдущие проверки должны били исключить пересечения до этого)
                              const Contour& c_IJ = *vcfp_IJ0.get();
                              const Contour& p2 = c_IJ.is_clockwise_oriented() ? Contour(c_IJ.vertices().rbegin(), c_IJ.vertices().rend()) : c_IJ;
                              if (CGAL::POSITIVE == CGAL::oriented_side(p0_I0, p2)) {
                                for (auto& c1 : vcfp_I) {
                                  vcfp_IJ.push_back(std::move(c1));
                                }
                                stack_cs.erase(stack_cs.begin() + I);
                                is_parent_hole_found = true;
                                break;
                              }
                            }
                          }
                        }
                        if (is_parent_hole_found == false) {
                          // <image url="..\code_images\file_0006.png" scale=".2"/>
                          // Если не найдена HOLE, в которой он мог бы находится, значит это общий внешний контур.
                          external_contours.push_back(std::move(vcfp_I));
                          stack_cs.erase(stack_cs.begin() + I);
                        }
                      } else {
                        // Не должно происходить, т.к. все не is_counterclockwise_oriented должны были быть удалены с конца стека на прошлом цикле, но не проверить нельзя.
                        // update: очень может быть, если это не первый элемент в цикле и перед этим мог отработать is_counterclockwise_oriented, он был соотнесён со своим holes
                        // и теперь настала очередь is_counterclockwise_oriented.
                        break;
                      }
                    }
                  }
                  // Удалить исходные полигоны объекта и в зависимости от shapes_mode1 выбрать либо inverted_holes или external_contours
                  object1.vect_polygons_with_holes.clear();

                  if (shapes_mode1 == 1 /*EXCLUDE_HOLES*/) {
                    for (auto& ec1 : external_contours) {
                      Polygon_2 plg(ec1[0]->begin(), ec1[0]->end());
                      if (ec1[0]->is_clockwise_oriented() == true) {
                        plg.reverse_orientation();
                      }
                      object1.vect_polygons_with_holes.push_back(std::make_shared<Polygon_with_holes_2>(plg));
                    }
                  } else if (shapes_mode1 == 2 /*INVERT_HOLES*/) {
                    for (auto& ec1 : vect_cs_inverted_holes) {
                      Polygon_2 plg(ec1[0]->begin(), ec1[0]->end());
                      if (ec1[0]->is_clockwise_oriented() == true) {
                        plg.reverse_orientation();
                      }
                      Polygon_with_holes_2 pwh = Polygon_with_holes_2(plg);
                      for (auto h_it = ec1.begin() + 1; h_it != ec1.end(); ++h_it) {
                        Contour hole1((*h_it)->begin(), (*h_it)->end());
                        if (hole1.is_counterclockwise_oriented() == true) {
                          hole1.reverse_orientation();
                        }
                        pwh.add_hole(hole1);
                      }
                      object1.vect_polygons_with_holes.push_back(std::make_shared<Polygon_with_holes_2>(pwh));
                    }
                  } else {
                    // TODO:: Исключение
                  }
                }
                object1_index++;
              }
            }

            // TODO: уже не актуально проверять параметр тут. Он должен быть поближе к расчётам.
            if (only_tests_for_valid == false) {

              {
                {
                  // Пройтись по объектам и запустить расчёты по отрицательным отступам.
                  // Получить список объектов с отрицательными offsets:
                  struct ObjectNegativeOffsetSS {
                    ClsObject object1;
                    //SsBuilder ssb; <- Непонятно по какой причине, но если ssb будет определён тут, то он работать не будет, а вызов функции construct_skeleton будет завершаться исключением access violation.
                    std::shared_ptr<Ss> ss = nullptr;
                  };

                  std::vector<ObjectNegativeOffsetSS> vect_objectsNegativeOffsetSS;
                  {
                    int I = 0;
                    std::vector<int> vectI_to_erase;
                    for (auto& object1 : objects) {
                      for (auto& oioa1 : object1.vect_oioa) {
                        if (oioa1.offset < 0) {
                          ClsObject negObject(object1.object_index, object1.vect_oioa, std::vector<ContourSequence>());
                          for (auto& pwh : object1.vect_polygons_with_holes) {
                            negObject.vect_polygons_with_holes.push_back(std::move(pwh));
                          }
                          for (auto& c1 : object1.source_planes) {
                            negObject.source_planes.push_back(std::move(c1));
                          }
                          vect_objectsNegativeOffsetSS.push_back({ negObject });
                          vectI_to_erase.push_back(I);
                          break;
                        }
                      }
                      I++;
                    }
                    for (auto idx_it = vectI_to_erase.rbegin(); idx_it != vectI_to_erase.rend(); ++idx_it) {
                      int& idx = *idx_it;
                      objects.erase(objects.begin() + idx); // Обязательно удалить объект с отрицательным offsets, чтобы CGAL не пытался его рассчитывать напрямую. Ниже будет отдельное ответвление для подготовки объектов с отрицательным offsets к расчёту в общем стеке.
                    }
                  }

                  if (vect_objectsNegativeOffsetSS.size() > 0)
                    // Дополнительно собрать информацию по инвертированным holes и передать их на обработку в positive offsets
                  {
                    const int threadNumbers = boost::thread::hardware_concurrency();
                    boost::asio::thread_pool pool1(threadNumbers);
                    boost::mutex mtx_;
                    for (auto& object1NegativeOffsetSS : vect_objectsNegativeOffsetSS) {
                      boost::asio::post(pool1, [&object1NegativeOffsetSS, &objects, &mtx_] {
                        // Отрицательные отступы считаются сразу на весь object1

                        {
                          // Краткое описание алгоритма расчёта внешнего контура
                        // <image url="..\code_images\file_0001.png" scale="1.0"/>
                        // 1. Определить общий максимальный frame, который охватит все контуры при максимальном offset (CGAL::compute_outer_frame_margin: https://doc.cgal.org/latest/Straight_skeleton_2/group__PkgStraightSkeleton2OffsetFunctions.html#ga80c0848e0145bbd531b1fc178fd07d33).
                        // 2. Инвертировать все Holes. При этом все Outer Boundaries станут границами Holes в пределах внешнего Frame
                        // 3. Соединить Frame и все старые Outer Boundaries как один Polygon_with_holes_2 и рассчитать его как внутренний offset (но уже давать положительные значения offset как исходные) (удалить после такого рассчёта все внешние offset, т.к. они не должны фигурировать в результате)
                        // 4. Рассчитать остальные инвертированные holes как острова (и тоже уже с положительными offset)
                        // 5. Объеденить результаты (3) и (4)
                        // Дополнительно - запоминать кто из контуров использовался в качестве определяющего внешний контур тут уже не важно (может в будущем и понадобиться, но прямо сейчас не нужно)

                        // Если при расчёте отрицательного offset буду обнаружены инвертированные holes, то для рассчёта таких инвертированных holes можно будет использовать хак:
                        // Контуры можно посчитать только положительные. Но раз инвертированные holes являются замкнутыми, то их теперь можно посчитать на равне с положительными,
                        // для этого нужно инвертировать негативное значение отступа, а altitude оставить таким же, как было. И всё это "упаковать" как новый объект.
                        // <image url="..\code_images\file_0013.png" scale=".5"/>

                        // Начало работы по алгоритму:
                        // Выполнить декомпозицию всех контуров и отсортировать их по площадям. (при этом их ориентация останется не тронутой)
                        // <image url="..\code_images\file_0003.png" scale=".5"/>
                          ContourSequence contours_sorted_by_area; // Список контуров, отсортированных по площади
                          for (auto& polygon1_with_holes : object1NegativeOffsetSS.object1.vect_polygons_with_holes) {
                            // <!-- ТАК ДЕЛАТЬ НЕЛЬЗЯ!!! В этом случае сохраняется связь с исходным объектом и если исходный объект удаляется,
                            // то его деструкторы УНИЧТОЖАТ объекты, на которые ссылается "умный" указатель! Поэтому нужно делать как указано ниже -
                            // СОЗДАВАТЬ отдельные объекты прямо под умные указатели, чтобы кроме умных указателей к ним никто не имел доступ!!! -->
                            //contours_sorted_by_area.push_back(std::shared_ptr<Contour>(&polygon1_with_holes->outer_boundary()));
                            //for (auto& hole1 : polygon1_with_holes->holes()) {
                            //  contours_sorted_by_area.push_back(std::shared_ptr<Contour>(&hole1));
                            //}
                            // <!-- -->

                            // НАДО ДЕЛАТЬ ТАК!!!
                            contours_sorted_by_area.push_back(std::make_shared<Contour>(polygon1_with_holes->outer_boundary().begin(), polygon1_with_holes->outer_boundary().end()));
                            for (auto& hole1 : polygon1_with_holes->holes()) {
                              contours_sorted_by_area.push_back(std::make_shared<Contour>(hole1.begin(), hole1.end()));
                            }
                          }
                          // Сортировку по площади выполнить независимо от направления контура:
                          sort(contours_sorted_by_area.begin(), contours_sorted_by_area.end(), [](const ContourPtr& p1, const ContourPtr& p2) {return abs(p1->area()) > abs(p2->area()); });
                          std::vector<ContourSequence> stack_cs; // vect_collect_for_polygons
                          for (auto& c1 : contours_sorted_by_area) {
                            ContourSequence cs1;
                            cs1.push_back(c1);
                            stack_cs.push_back(cs1);
                          }

                          // Отсортировать полигоны в топологической последовательности. 
                          // Перед расчётом важно убедиться, что если какой-то полигон вложен в другой полигон, то такие полигоны не могут использоваться вместе для рассчёта external offset
                          // (Во-первых внешний и так больше, а во вторых функция construct_skeleton входит в бесконечный цикл, если передать объект с такой топологией).
                          // Они должны быть преобразованы в инвертированные holes (т.е. стать островами)
                          // Если предположить, что все объекты отделены друг от друга условной "водой", т.е. между внутренней границей одного объекта и внутренней границей другого всегда есть расстояние,
                          // то достаточно определить, что первая точка внешнего контура одного объекта является вложенной в hole другого объекта.
                          // Сначала нужно получить сортированный список контуров по площадям. По этому списку будет работать алгоритм инвертирования holes. Все Holes, которые будут в итоге инвертированы, даже с
                          // островами, могут быть рассчитаны отдельно и независимо друг от друга. Когда в цикле сортировки останутся только старые outer boundaries, то именно они и пойдут в расчёт общего
                          // внешнего контура.

                          std::vector<ContourSequence> external_contours;
                          std::vector<ContourSequence> vect_cs_inverted_holes;
                          while (stack_cs.size() > 0) {

                            // Перенести последние контура, являющиеся отверстиями, в список результирующих контуров (они могут содержать в себе бывшие внешние контура, кто-то уже добавился)
                            // <image url="..\code_images\file_0004.png" scale=".5"/>
                            for (int I = stack_cs.size() - 1; I >= 0; I--) {
                              ContourSequence& vcfp_I = stack_cs[I];
                              ContourPtr& vcfp_I0 = vcfp_I[0];
                              if (vcfp_I0->is_clockwise_oriented() == true) {
                                vect_cs_inverted_holes.push_back(std::move(vcfp_I));
                                stack_cs.erase(stack_cs.begin() + I);
                              } else {
                                break;
                              }
                            }

                            // Если последний контур является Outer Boundary, то поискать в кого он может быть вложен.
                            // <image url="..\code_images\file_0005.png" scale=".2"/>
                            for (int I = stack_cs.size() - 1; I >= 0; I--) {
                              ContourSequence& vcfp_I = stack_cs[I];
                              ContourPtr& vcfp_I0 = stack_cs[I][0];
                              // В OUTER контуры ничего добавлять на нужно. Только они добавляются в HOLE-контуры.
                              if (vcfp_I0->is_counterclockwise_oriented() == true) {
                                bool is_parent_hole_found = false;
                                if (stack_cs.size() > 2) {
                                  Point_2 p0_I0 = vcfp_I0->vertices().at(0);
                                  for (int IJ = stack_cs.size() - 2; IJ >= 0; IJ--) {
                                    ContourSequence& vcfp_IJ = stack_cs[IJ];
                                    ContourPtr& vcfp_IJ0 = stack_cs[IJ][0];
                                    if (vcfp_IJ0->is_clockwise_oriented() == true) {
                                      // https://stackoverflow.com/questions/64046074/is-there-a-cgal-function-that-checks-if-a-point-is-inside-a-linear-polygon-with
                                      // Если хотя бы одна точка находится внутри, то считаем. что и весь контур находится внутри (Предыдущие проверки должны били исключить пересечения до этого)
                                      const Contour& c_IJ = *vcfp_IJ0.get();
                                      const Contour& p2 = c_IJ.is_clockwise_oriented() ? Contour(c_IJ.vertices().rbegin(), c_IJ.vertices().rend()) : c_IJ;
                                      if (CGAL::POSITIVE == CGAL::oriented_side(p0_I0, p2)) {
                                        for (auto& c1 : vcfp_I) {
                                          vcfp_IJ.push_back(std::move(c1));
                                        }
                                        stack_cs.erase(stack_cs.begin() + I);
                                        is_parent_hole_found = true;
                                        break;
                                      }
                                    }
                                  }
                                }
                                if (is_parent_hole_found == false) {
                                  // <image url="..\code_images\file_0006.png" scale=".2"/>
                                  // Если не найдена HOLE, в которой он мог бы находится, значит это общий внешний контур.
                                  external_contours.push_back(std::move(vcfp_I));
                                  stack_cs.erase(stack_cs.begin() + I);
                                }
                              } else {
                                // Не должно происходить, т.к. все не is_counterclockwise_oriented должны были быть удалены с конца стека на прошлом цикле, но не проверить нельзя.
                                // UPDATE: очень может быть, если это не первый элемент в цикле и перед этим мог отработать is_counterclockwise_oriented, он был соотнесён со своим holes
                                // и теперь настала очередь is_counterclockwise_oriented.
                                break;
                              }
                            }
                          }

                          // Рассчитать вспомогательный внешний Frame для расчёта наибольшего отрицательного offset, в который войдут все внешние объекты.
                          // <image url="..\code_images\file_0009.png" scale=".2"/>
                          {
                            double max_abs_negative_offset = abs(std::min_element(object1NegativeOffsetSS.object1.vect_oioa.begin(), object1NegativeOffsetSS.object1.vect_oioa.end(), [](const OIOA& o1, const OIOA& o2) { return o1.offset < o2.offset; })->offset);
                            std::vector<double> vxmin, vymin, vxmax, vymax;
                            ContourSequence allowed_polygons; // Для каких полигонов удалось вычислить margin
                            for (auto& cs : external_contours) {
                              CGAL::Bbox_2 bbox = cs[0]->bbox();
                              // <image url="..\code_images\file_0010.png" scale=".2"/>
                              std::optional<double> margin = CGAL::compute_outer_frame_margin(cs[0]->begin(), cs[0]->end(), max_abs_negative_offset);
                              if (margin) {
                                allowed_polygons.push_back(std::move(cs[0]));
                                vxmin.push_back(bbox.xmin() - *margin);
                                vymin.push_back(bbox.ymin() - *margin);
                                vxmax.push_back(bbox.xmax() + *margin);
                                vymax.push_back(bbox.ymax() + *margin);
                              } else {
                                // TODO: Запомнить ошибку, почему-то не удалось посчитать margin
                              }
                            }
                            if (allowed_polygons.size() > 0) {
                              // Определить суммарные максимальные/минимальные координаты (https://en.cppreference.com/w/cpp/algorithm/min_element):

                              double xmin = *std::min_element(vxmin.begin(), vxmin.end());
                              double ymin = *std::min_element(vymin.begin(), vymin.end());
                              double xmax = *std::max_element(vxmax.begin(), vxmax.end());
                              double ymax = *std::max_element(vymax.begin(), vymax.end());
                              ContourPtr frame_outer_boundary = std::make_shared<Contour>();
                              frame_outer_boundary->push_back(Point_2(xmin, ymin));
                              frame_outer_boundary->push_back(Point_2(xmax, ymin));
                              frame_outer_boundary->push_back(Point_2(xmax, ymax));
                              frame_outer_boundary->push_back(Point_2(xmin, ymax));
                              // <image url="..\code_images\file_0011.png" scale=".2"/>

                              ContourSequence cs_frame;
                              cs_frame.push_back(frame_outer_boundary);
                              std::shared_ptr<Polygon_with_holes_2> main_outer_negative_offset_polygon = std::make_shared<Polygon_with_holes_2>(Polygon_2(frame_outer_boundary->begin(), frame_outer_boundary->end()));

                              for (auto& hole1 : allowed_polygons) {
                                hole1->reverse_orientation();
                                main_outer_negative_offset_polygon->add_hole(Polygon_2(hole1->begin(), hole1->end()));
                                cs_frame.push_back(std::move(hole1));
                              }

                              ClsObject object1FrameForNegativeOffsets(object1NegativeOffsetSS.object1.object_index, std::vector<OIOA>(), std::vector<ContourSequence>());
                              object1FrameForNegativeOffsets.source_planes.push_back(cs_frame);
                              for (auto& oioa1 : object1NegativeOffsetSS.object1.vect_oioa) {
                                if (oioa1.offset < 0) {
                                  object1FrameForNegativeOffsets.vect_oioa.push_back(OIOA(oioa1.object_index, oioa1.offset_index, oioa1.offset, OIOA::NEGATIVE_TYPE::WITH_EXTERNAL_FRAME, oioa1.altitude));
                                }
                              }
                              object1FrameForNegativeOffsets.vect_polygons_with_holes.push_back(main_outer_negative_offset_polygon); // <- Не забывать подгружать получающиеся полигоны для рассчёта. Так-то можно было бы сделать их получение после рассчёта отрицательных offsets, но раньше пришлось проверять пересечения контуров, а это требовало наличия полигонов. Поэтому, раз получены новые контуры для расчёта отрицательных полигонов, то нужно нужно их преобразовать в полигон.
                              mtx_.lock();
                              objects.push_back(object1FrameForNegativeOffsets);
                              mtx_.unlock();
                              // Расчёт Straight Skeleton для отрицательного offset выполняется вместе с расчётами положительного offset позже.
                            }
                          }

                          if (vect_cs_inverted_holes.size() > 0) {
                            // Если имеются инвертированные holes, то нужно добавить их в общий стек расчёта:
                            ClsObject object1WithInvertedHoles(object1NegativeOffsetSS.object1.object_index, std::vector<OIOA>(), vect_cs_inverted_holes);
                            for (OIOA& oap1 : object1NegativeOffsetSS.object1.vect_oioa) {
                              if (oap1.offset < 0) {
                                object1WithInvertedHoles.vect_oioa.push_back(OIOA(object1NegativeOffsetSS.object1.object_index, oap1.offset_index, oap1.offset /*abs берётся позже*/, OIOA::NEGATIVE_TYPE::INVERTED_HOLE, oap1.altitude));
                              }
                            }
                            if (object1WithInvertedHoles.vect_oioa.size() > 0) {
                              // Добавить в объект с инвертированными holes сами инвертированные полигоны:
                              for (auto& plane1 : vect_cs_inverted_holes) {
                                if (plane1.size() > 0) {
                                  if (plane1[0]->is_clockwise_oriented()) {
                                    plane1[0]->reverse_orientation();
                                  }
                                  Polygon_with_holes_2 pwh(Polygon_2(plane1[0]->begin(), plane1[0]->end()));
                                  for (int I = 1; I <= plane1.size() - 1; I++) {
                                    if (plane1[I]->is_counterclockwise_oriented() == true) {
                                      plane1[I]->reverse_orientation();
                                    }
                                    pwh.add_hole(*plane1[I]);
                                  }
                                  object1WithInvertedHoles.vect_polygons_with_holes.push_back(std::make_shared<Polygon_with_holes_2>(pwh));
                                }
                              }
                              // Добавить объект с инвертированными holes в общий список объектов для рассчёта его как "обыкновенные" положительные offset:
                              mtx_.lock();
                              vect_cs_inverted_holes.clear();
                              objects.push_back(object1WithInvertedHoles);
                              mtx_.unlock();
                            }
                          }

                          {
                            // Если в объекте есть не только отрицательные offsets, но и положительные, то нужно и задать параметры для их расчёта тоже:
                            std::vector<OIOA> vect_oioa;
                            for (auto& oioa1 : object1NegativeOffsetSS.object1.vect_oioa) {
                              if (oioa1.offset >= 0) {
                                vect_oioa.push_back(OIOA(oioa1.object_index, oioa1.offset_index, oioa1.offset, oioa1.altitude));
                              }
                            }
                            if (vect_oioa.size() > 0) {
                              ClsObject object1WithPositiveOffsets(object1NegativeOffsetSS.object1.object_index, vect_oioa, object1NegativeOffsetSS.object1.source_planes);
                              for (auto& pwh1 : object1NegativeOffsetSS.object1.vect_polygons_with_holes) {
                                object1WithPositiveOffsets.vect_polygons_with_holes.push_back(pwh1);
                              }
                              mtx_.lock();
                              objects.push_back(object1WithPositiveOffsets);
                              mtx_.unlock();
                            }
                          }
                        }
                        }
                      ); // boost::asio
                    }
                    pool1.join();
                  }
                }

                // Расчёт offsets. Напоминаю, что тут находятся не только исходные параметры, но и полученные путём преобразования параметров при подготовке данных для расчёта отрицательных offsets.
                // Для начала надо получить один скелетон на plane.
                /* Positive Offset Multithread Struct. Keep single polygon with holes and a list of offsets and altitudes. */
                struct POMS {
                  int I = 0;
                  std::shared_ptr<Polygon_with_holes_2> polygon1_with_holes;
                  // List of
                  std::vector<OIOA> oioa;
                  // Keep Struct Skeleton (or zero if failed to calc)
                  std::shared_ptr<Ss> ss{ nullptr };
                  POMS(int _object_index, std::shared_ptr<Polygon_with_holes_2> _polygon1_with_holes, std::vector<OIOA> _oioa)
                    :polygon1_with_holes(_polygon1_with_holes), oioa(_oioa) {
                    I = _object_index;
                  }
                };

                std::vector<POMS> vect_polygon1_oioa;
                for (auto& object1 : objects) {
                  for (auto& polygon_with_holes1 : object1.vect_polygons_with_holes) {
                    // Отсортировать по offsets. Будет удобно считать, если при увеличении отступа пропадёт контур,то и считать следующие бОльшие контуры не надо
                    // - в многопоточности это уже не актуально. Может и можно в итоге прервать, но надо поискать другой метод. sort(positive_offsets_data1.vect_oap.begin(), positive_offsets_data1.vect_oap.end(), [](const OAP& oap1, const OAP& oap2) {return oap1.offset < oap2.offset; });
                    // Положительные отступы считаются по одному
                    POMS p1(object1.object_index, std::shared_ptr<Polygon_with_holes_2>(), object1.vect_oioa);
                    p1.polygon1_with_holes = std::move(polygon_with_holes1);
                    vect_polygon1_oioa.push_back(p1);
                  }
                }

                int threadNumbers = boost::thread::hardware_concurrency();
                boost::asio::thread_pool pool3(threadNumbers);
                boost::mutex mtx_;

                int threads_counts = 0;
                for (auto& polygon1_oioa : vect_polygon1_oioa) {
                  boost::asio::post(pool3, [&polygon1_oioa, &res_errors, &res_contours, &mtx_, threads_counts, verbose] {
                    SsBuilder ssb;
                    int verts_count = polygon1_oioa.polygon1_with_holes->outer_boundary().size();
                    {
                      // Загрузить outer boundary и holes в ssb.

                      // Загрузить внешний контур:
                      ssb.enter_contour(polygon1_oioa.polygon1_with_holes->outer_boundary().begin(), polygon1_oioa.polygon1_with_holes->outer_boundary().end());
                      // Загрузить holes
                      for (auto& hole1 : polygon1_oioa.polygon1_with_holes->holes()) {
                        verts_count += hole1.size();
                        ssb.enter_contour(hole1.begin(), hole1.end());
                      }
                    }
                    // Construct the skeleton
                    CGAL::Real_timer timer1;
                    timer1.start();
                    try {
                      polygon1_oioa.ss = ssb.construct_skeleton();
                    } catch (std::exception _ex) {
                    }
                    timer1.stop();

                    if (verbose == true) {
                      mtx_.lock();
                      printf("\n " VAL2STR(Err::_0026) ". SS 2D Offset. Thread: %u, SsBuilder verts: %d, build time: %.5g", threads_counts, verts_count, timer1.time());
                      mtx_.unlock();
                    }
                    // Proceed only if the skeleton was correctly constructed.
                    mtx_.lock();
                    if (polygon1_oioa.ss) {
                      // Всё норм, ничего не делать
                    } else {
                      try {
                        // Не получилось посчитать Straight Skeleton. Надо сообщить о проблеме пользователю:
                        res_errors.push_back(ObjectError(polygon1_oioa.oioa[0].object_index, *polygon1_oioa.polygon1_with_holes.get(), VAL2STR(Err::_0006)". Offset. Error Build Straight Skeleton. Outer Boundary"));
                      } catch (std::exception _ex) {
                        int error = 0;
                      }
                    }
                    mtx_.unlock();
                    }  // boost::asio::post
                  );
                  threads_counts++;
                }
                pool3.join();

                boost::asio::thread_pool pool4(threadNumbers);
                for (POMS& polygon1_oioa : vect_polygon1_oioa) {
                  if (polygon1_oioa.ss) {
                    for (auto& oioa1 : polygon1_oioa.oioa) {
#ifdef _DEBUG
                      // Пропускать offset == 0 при отладке. Срабатывает проверка на 0 в offset skeleton, хотя сам skeleton при 0 строится нормально
                      if (is_zero(oioa1.offset)) {
                        continue;
                      }
#endif
                      boost::asio::post(pool4, [&polygon1_oioa, &oioa1, &mtx_, &res_contours] {
                        ContourSequence offset_contours; // Instantiate the container of offset contours
                        OffsetBuilder ob(*polygon1_oioa.ss); // Instantiate the offset builder with the skeleton

                        ob.construct_offset_contours(abs(oioa1.offset)/*Здесь считаются и отрицательные offsets для внешних контуров, но данные подготовлены для инверсии*/, std::back_inserter(offset_contours)); // Obtain the offset contours

                        // Ещё вариант конфигурации контуров, которые требуется идентифицировать при сопоставлении контуров (т.е. и количество контуров меняется и результирующие фигуры тоже меняются):
                        // TODO: Нужно разобраться где разбираться с конфигурацией. Тут разобраться не получится, потому что расчёт с положительным offset (хоть он тут и берётся через abs)
                        // может являтся частью расчёта отрицательного offset, но тут их сопоставить нельзя из-за нехватки всех данных для сопоставления (внешний контур отрицательного offset
                        // считается отдельно от внутренних inverted holes). Собирать их надо в другом месте.

                        if (offset_contours.size() > 0) {
                          // Load results contours (may be several contours if holes)
                          // Сохранить рассчитанные по offsets горизонтальные контуры как результат (их может быть несколько, когда в объекте имеются holes):
                          // <image url="..\code_images\file_0014.png" scale=".3"/>
                          if (oioa1.offset >= 0) {
                            // Если offset положительный, то загрузить все контуры.
                            // TODO: В документации возвращаемый результат не определяет какую-то гарантированную последовательность 
                            // возвращаемых контуров. Известно, что исходный полигон всегда один, но может иметь Holes, поэтому
                            // при возврате результат может содержать несколько контуров: один внешний (обязательно) и несколько внутренних (опционально, если изначально
                            // в исходном контуре были holes). Для отображения результата в режиме EDGES последовательность контуров не важна,
                            // а для режиме FACES это важно и пока это не реализовано.
                            // Update: Есть некоторая зацепка - направления контуров. Есть функции c1->is_counterclockwise_oriented() и c1->is_clockwise_oriented(). Они определяют
                            // какой из контуров является Hole, а какой Outer Boundary.
                            // <image url="..\code_images\file_0017.png" scale=".5"/>
                            int count_of_contours = offset_contours.size();
                            // Проверил последовательность выдачи контуров на сложной фигуре. Вывод - нельзя предположить,
                            // что контуры выдаются последовательно по определённому алгоритму. Например, сначала идёт 
                            // внешний контур, потом holes, потом опять внешний контор, потом опять его holes.
                            // Но раз есть корректные значения CW и CCW, то можно будет проще отсортировать контуры по вложенности.
                            // Складывается впечатление, что контуры считаются в последовательности, в которой они добавлены в объект или как они
                            // распознаны системой подготовки исходных данных, но на это нельзя опираться.
                            // <image url="..\code_images\file_0018.png" scale=".3"/>

                            for (ContourSequence::iterator i = offset_contours.begin(); i != offset_contours.end(); ++i) {
                              const ContourPtr& c1 = *i;
                              bool cw_oriented = c1->is_clockwise_oriented();
                              bool ccw_oriented = c1->is_counterclockwise_oriented();
                              bool is_simple = c1->is_simple();

                              if (c1->size() > 0) {
                                oioa1.polygon_contours.push_back(c1);
                              }
                            }
                          } else {
                            if (oioa1.neg_type == OIOA::NEGATIVE_TYPE::WITH_EXTERNAL_FRAME) {
                              // Если offset отрицательный, и отмечен признаком WITH_EXTERNAL_FRAME, то у такого offset не загружать внешний контур, 
                              // т.к. он является вспомогательным для рассчётов и не должен отображаться.
                              // Locate the offset contour that corresponds to the frame
                              // That must be the outmost offset contour, which in turn must be the one
                              // with the largetst unsigned area.
                              ContourSequence::iterator f = offset_contours.end();
                              double lLargestArea = 0.0;
                              for (ContourSequence::iterator i = offset_contours.begin(); i != offset_contours.end(); ++i) {
                                double lArea = CGAL_NTS abs((*i)->area()); //Take abs() as  Polygon_2::area() is signed.
                                if (lArea > lLargestArea) {
                                  f = i;
                                  lLargestArea = lArea;
                                }
                              }
                              // Remove the offset contour that corresponds to the frame.
                              offset_contours.erase(f);
                            }

                            // Для расчёта отрицательного offset все контуры инвертируются, поэтому при возврате результатов
                            // требуется инвертировать результаты повторно (чтобы вернуть исходную, корректную форму вложенности контуров).
                            for (auto& c1 : offset_contours) {
                              if (c1->size() > 0) {
                                c1->reverse_orientation();
                                oioa1.polygon_contours.push_back(c1);
                              } else {
                                // TODO: проверить, бывают ли такие ситуации, когда после получения результата контур не содержит точек?
                              }
                            }
                          }
                          mtx_.lock();
                          // Загрузить получившиеся данные в массив с результатами:
                          res_contours.push_back(oioa1);
                          mtx_.unlock();

                        }
                        }); // boost::asio::post
                    }
                  }
                }
                pool4.join();

                // Восстановить последовательность результатов расчётов offset по исходным позициям offset (как их задал пользователь):
                sort(res_contours.begin(), res_contours.end(), [](const OIOA& oap1, const OIOA& oap2) {return oap1.offset_index < oap2.offset_index; });
                sort(res_contours.begin(), res_contours.end(), [](const OIOA& oap1, const OIOA& oap2) {return oap1.object_index < oap2.object_index; });
                sort(res_errors.begin(), res_errors.end(), [](const ObjectError& oap1, const ObjectError& oap2) {return oap1.object_index < oap2.object_index; });
              }
            } else {
              if (verbose == true) {
                unsigned int errors_size = res_errors.size();
                if (errors_size == 0) {
                  printf("\nStraight Skeleton Offsets tests finished. Errors not found.");
                } else {
                  printf("\nStraight Skeleton Offsets tests finished. Errors count %u", errors_size);
                }
              }
            }
            mesh_data->has_error = false;
            mesh_data->str_error = NULL;
            // Расчёт offset заканчивается тут.
          }

          if (timer.is_running() == true) { // it can be stopped in catch
            timer.stop();
          }
          if (verbose == true) {
            printf("\nOffset computation of object took %.5g sec.", timer.time());
          }

          // Отмапить результаты offset-ов на индексы исходных объектов.
          std::map<int /*индекс исходного объекта*/, std::vector<OIOA> /*данные для полигонов для offset_index текущего объекта*/> map_res_by_objectindex;
          {
            for (auto& object1 : objects) {
              if (map_res_by_objectindex.find(object1.object_index) == map_res_by_objectindex.end()) {
                map_res_by_objectindex[object1.object_index] = std::vector<OIOA>();
              }
            }
            // 1. Сгруппировать контуры результатов по объектам
            for (auto& shape1 : res_contours) {
              map_res_by_objectindex[shape1.object_index].push_back(shape1);
            }
            // 2. Отсортировать индексы по возрастанию в каждой группе
            for (auto& map_res1 : map_res_by_objectindex) {
              // Сортировать по offset_index
              sort(map_res1.second.begin(), map_res1.second.end(), [](const OIOA& oioa1, const OIOA& oioa2) {return oioa1.offset_index < oioa2.offset_index; });
            }
          }
          mesh_data->nn_source_objects_count = map_res_by_objectindex.size(); // Общее значение, используется при подготовке всех выходных данных
          {
            // Собрать все зарегистрированные ошибки из объекта res_errors:
            // 1. Сгруппировать ошибки по объектам:
            std::map<int, std::vector<ObjectError>> map_res_errors;
            // У каждого объекта должен быть список ошибок, даже если он нулевой длины, поэтому нужно сначала создать списки,
            // а потом их заполнить. Если res_errors окажется пустым, то не получится сообщить КАЖДОМУ объекту, что у него 0 ошибок...
            for (auto& map_res1 : map_res_by_objectindex) {
              map_res_errors[map_res1.first] = std::vector<ObjectError>();
            }
            // ... Если ошибок не будет, то map_res_errors будет пустой
            for (auto& error1 : res_errors) {
              map_res_errors[error1.object_index].push_back(error1);
            }
            std::vector<int>                vect_errors_count; // Количество ошибок на объект
            std::vector<char*>              vect_description1_per_error1; // описания всех ошибок
            std::vector<int>                vect_contours_per_error1; // Количество контуров в ошибках (1-один контур, >1 обычно обычно все контуры объекта. Для понимания интерпретации результата)
            std::vector<int>                vect_vertices_per_contour1; // Количество вершин на контур
            std::vector<std::vector<double>> vect_vertices_of_errors;
            int error_counts_summa = 0;
            for (auto& kv : map_res_errors) {
              int errors_count = kv.second.size();
              error_counts_summa += errors_count;
              vect_errors_count.push_back(errors_count);
              for (auto& error1 : kv.second) {
                vect_description1_per_error1.push_back(error1.message);
                vect_contours_per_error1.push_back(error1.cs1.size());
                for (auto& c1 : error1.cs1) {
                  vect_vertices_per_contour1.push_back(c1.get()->size());
                  for (auto& v : *c1.get()) {
                    vect_vertices_of_errors.push_back(std::vector<double>{v.x(), v.y()});
                  }
                }
              }
            }

            // Определить количество памяти для выгрузки информации об ошибках:
            mesh_data->nn_errors_count = (int*)malloc(sizeof(int) * vect_errors_count.size());
            for (int I = 0; I <= vect_errors_count.size() - 1; I++) {
              mesh_data->nn_errors_count[I] = vect_errors_count[I];
            }

            if (error_counts_summa == 0) {
              mesh_data->nn_description1_per_error1 = NULL;
              mesh_data->nn_contours_per_error1 = NULL;
              mesh_data->nn_vertices_per_contour1 = NULL;
              mesh_data->vertices_of_errors = NULL;
            } else {

              mesh_data->nn_description1_per_error1 = (char**)malloc(sizeof(char*) * vect_description1_per_error1.size());
              mesh_data->nn_contours_per_error1 = (int*)malloc(sizeof(int) * vect_contours_per_error1.size());
              mesh_data->nn_vertices_per_contour1 = (int*)malloc(sizeof(int) * vect_vertices_per_contour1.size());
              mesh_data->vertices_of_errors = (float*)malloc(sizeof(float[3]) * vect_vertices_of_errors.size());
              for (int I = 0; I <= vect_description1_per_error1.size() - 1; I++) {
                mesh_data->nn_description1_per_error1[I] = vect_description1_per_error1[I];
              }
              for (int I = 0; I <= vect_contours_per_error1.size() - 1; I++) {
                mesh_data->nn_contours_per_error1[I] = vect_contours_per_error1[I];
              }
              for (int I = 0; I <= vect_vertices_per_contour1.size() - 1; I++) {
                mesh_data->nn_vertices_per_contour1[I] = vect_vertices_per_contour1[I];
              }
              for (int I = 0; I <= vect_vertices_of_errors.size() - 1; I++) {
                mesh_data->vertices_of_errors[I * 3 + 0] = vect_vertices_of_errors[I][0];
                mesh_data->vertices_of_errors[I * 3 + 1] = vect_vertices_of_errors[I][1];
                mesh_data->vertices_of_errors[I * 3 + 2] = 0.0;
              }
            }
          }

          {
            std::set<int> set_source_objects_index; // Исходные индексы объектов, пришедших на обработку
            mesh_data->nn_source_objects_indexes = (int*)malloc(sizeof(int) * mesh_data->nn_source_objects_count);

            for (auto& object1 : objects) {
              set_source_objects_index.insert(object1.object_index);
            }
            std::vector<int> vect_source_objects_indexes;
            for (auto& object1_index : set_source_objects_index) {
              vect_source_objects_indexes.push_back(object1_index);
            }
            for (int I = 0; I <= vect_source_objects_indexes.size() - 1; I++) {
              mesh_data->nn_source_objects_indexes[I] = vect_source_objects_indexes[I];
            }
          }

          if (only_tests_for_valid == false) {
            // Generate Results
            {
              // Требуется вернуть контуры с негативными отступами в исходное состояние, Сейчас некоторые из них могут содержать holes
              // <image url="..\code_images\file_0015.png" scale=".2"/>
              // Ещё вариант конфигурации контуров, которые требуется идентифицировать при сопоставлении контуров (т.е. и количество контуров меняется и результирующие фигуры тоже меняются):
              // <image url="..\code_images\file_0017.png" scale=".5"/>
              // Данные для возврата пользователю.
              //std::vector<int> vect_objects_index; // Исходные индексы объектов, отправленных на обработку
              std::vector<int> vect_objects_count_of_offsets; // Количество offsets объекта
              std::vector<int> vect_objects_count_of_verts; // Количество vertices объекта
              std::vector<int> vect_objects_count_of_edges; // Количество edges объекта
              std::vector<int> vect_objects_count_of_faces; // Количество faces объекта
              std::vector<int> vect_objects_offsets; // Использованные offsets
              std::vector<std::vector<float>> vect_objects_verts; // Координаты vertices
              std::vector<std::vector<unsigned int>> vect_objects_edges; // indexes of edges (per object1)
              std::vector<std::vector<unsigned int>> vect_objects_faces; // indexes of faces (per object1)

              std::vector<VECT_OIOA_PWHS> vect_pwhs_offset_index_order;

              for (auto& [object_index, res_by_object1] : map_res_by_objectindex) {
                // Собрать список контуров на одном уровне offset и сложить их в нужной последовательности,
                // чтобы они выстроились в валидную фигуру с внешним контуром и отверстиями:
                std::map<int, std::vector<OIOA>> object1_group_contours_by_offset_index;
                for (auto& oioa1 : res_by_object1) {
                  if (object1_group_contours_by_offset_index.find(oioa1.offset_index) == object1_group_contours_by_offset_index.end()) {
                    object1_group_contours_by_offset_index[oioa1.offset_index] = std::vector<OIOA>();
                  }
                  object1_group_contours_by_offset_index[oioa1.offset_index].push_back(oioa1);
                }

                for (auto& kv /*Насколько я понимаю, map работает как упорядоченная коллекция и выдаёт результат в порядке offset_index, что является исходным заданным порядком*/
                  : object1_group_contours_by_offset_index
                  ) {
                  std::vector<OIOA> vect_oioa = kv.second;

                  ContourSequence cs_for_mesh;
                  for (auto& oioa1 : vect_oioa) {
                    if (oioa1.polygon_contours.size() > 0) {
                      for (auto& cs : oioa1.polygon_contours) {
                        cs_for_mesh.push_back(cs);
                      }
                    }
                  }

                  if (cs_for_mesh.size() == 0) {

                  } else {
                    std::vector<Polygon_with_holes_2> vect_pwh; // набор pwh соответствующих одному offset для текущего объекта
                    ConvertContourSequenceIntoPolygonsWithHoles(cs_for_mesh, vect_pwh);
                    vect_pwhs_offset_index_order.push_back(VECT_OIOA_PWHS(vect_oioa[0], vect_pwh));
                  }
                }
              }

              // Определить способ обработки объектов по join_mode и сформировать новый список данных для выгрузки результатов (новый список результирующих объектов)
              std::map<int /*индекс результирующего объекта*/, std::vector<VECT_OIOA_PWHS> /*данные для полигонов для offset_index результирующего объекта*/> map_join_mesh;
              {
                // Как работает join_mode <image url="..\code_images\file_0024.png" scale=".2"/>. Для режима EDGES всё тоже самое, только выводятся только контуры соответствующего face без триангуляции (с учётом holes)
                // и с учётом направленности - внешние контуры закручены против часовой, а внутренние holes закручены по часовой: <image url="..\code_images\file_0025.png" scale=".2"/>
                int split_object_index = 0; // для join_mode==SPLIT индекс нового объекта назначается группе контуров с одним index_offset (а не каждому контуру, чтобы количество объектов в режиме Edges и Faces совпадало)
                int merge_object_index = 0; // Для join_mode=Merge индекс результирующего объекта один - 0, т.к. объект будет единственным.
                int  keep_object_index = 0; // Для join_mode=KEEP индексы результирующих объектов не меняются
                for (auto& pwhs_offset_index_order1 : vect_pwhs_offset_index_order) {
                  if (results_join_mode == 0) {
                    // Сгруппировать результирующие объекты по отдельным offset_index
                    // Разделить текущий набор pwhs на отдельные pwh1 и сделать каждый из них отдельным вектором (для удобства дельнейшей обработки)
                    for (auto& pwh1 : pwhs_offset_index_order1.vect_pwh) {
                      std::vector<Polygon_with_holes_2> pwh1_as_vector;
                      pwh1_as_vector.push_back(pwh1);
                      std::vector<VECT_OIOA_PWHS> vect_PWHS_ALTITUDE_as_vector;
                      vect_PWHS_ALTITUDE_as_vector.push_back(VECT_OIOA_PWHS(pwhs_offset_index_order1.oioa1, pwh1_as_vector));
                      map_join_mesh[split_object_index] = vect_PWHS_ALTITUDE_as_vector;
                      split_object_index++;
                    }
                  } else if (results_join_mode == 1) {
                    keep_object_index = pwhs_offset_index_order1.oioa1.object_index;
                    // Оставить результирующие объекты как в исходных объектах
                    if (map_join_mesh.find(keep_object_index) == map_join_mesh.end()) {
                      map_join_mesh[keep_object_index] = std::vector<VECT_OIOA_PWHS>();
                    }
                    map_join_mesh[keep_object_index].push_back(pwhs_offset_index_order1);
                  } else { // 2 - для всех остальных случаев
                    // Объеденить все результаты в один объект
                    if (map_join_mesh.find(merge_object_index) == map_join_mesh.end()) {
                      map_join_mesh[merge_object_index] = std::vector<VECT_OIOA_PWHS>();
                    }
                    map_join_mesh[merge_object_index].push_back(pwhs_offset_index_order1);
                  }
                }
              }

              mesh_data->nn_objects = map_join_mesh.size();
              mesh_data->nn_objects_indexes = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
              mesh_data->nn_offsets_counts = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
              mesh_data->nn_verts = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
              mesh_data->nn_edges = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
              mesh_data->nn_faces = (int*)malloc(sizeof(int) * mesh_data->nn_objects);

              {
                // Предполагаю использовать эту структуру в будущем для распараллеливания расчёта триангуляции.
                struct OBJECT1_MESH {
                  std::vector<std::vector<       float>> vect_res_object1_verts; // Координаты vertices (per object1)
                  std::vector<std::vector<unsigned int>> vect_res_object1_edges; // indexes of edges (per object1)
                  std::vector<std::vector<unsigned int>> vect_res_object1_faces; // indexes of faces (per object1)
                  std::set<int> set_object1_offsets_indexes; // Индексы offset из которых получился mesh для этого объекта
                };
              }

              // Время на триангуляцию тратится незначительное. Можно пока не делать распараллеливание. <image url="..\code_images\file_0023.png" scale="1.0"/>
              for (auto& [object_index, vect_pwhs_with_offset_index] : map_join_mesh) {
                //int object_index = KV.first;
                //auto& vect_pwhs_with_offset_index = KV.second;

                std::vector<std::vector<       float>> vect_res_object1_verts; // Координаты vertices (per object1)
                std::vector<std::vector<unsigned int>> vect_res_object1_edges; // indexes of edges (per object1)
                std::vector<std::vector<unsigned int>> vect_res_object1_faces; // indexes of faces (per object1)

                std::set<int> set_object1_offsets_indexes; // Индексы offset из которых получился mesh для этого объекта
                for (auto& pwhs_offset_index_order1 : vect_pwhs_with_offset_index) {
                  set_object1_offsets_indexes.insert(pwhs_offset_index_order1.oioa1.offset_index);
                }

                //  Рассчитать результат уже независимо от join_mode:
                if (result_type == ResType::EDGES) {
                  ConvertPWH2IntoVerticesEdgesNoFaces(vect_pwhs_with_offset_index, vect_res_object1_verts, vect_res_object1_edges, vect_res_object1_faces);
                } else if (result_type == ResType::FACES) {
                  std::vector<Mesh> vect_sm_offset_index_order;
                  for (auto& pwh_altitude : vect_pwhs_with_offset_index) {
                    for (auto& pwh1 : pwh_altitude.vect_pwh) {
                      Mesh sm1;
                      ConvertPolygonWithHoles2IntoMesh(pwh1, pwh_altitude.oioa1.altitude, sm1);
                      vect_sm_offset_index_order.push_back(sm1);
                    }
                  }
                  // Из набора Mesh надо сформировать полигональную сетку
                  ConvertMeshesIntoVerticesEdgesFaces(vect_sm_offset_index_order, vect_res_object1_verts, vect_res_object1_edges, vect_res_object1_faces);
                }

                mesh_data->nn_objects_indexes[object_index] = object_index; // Индекс объекта можно взять из первого результата
                mesh_data->nn_offsets_counts[object_index] = set_object1_offsets_indexes.size();
                mesh_data->nn_verts[object_index] = vect_res_object1_verts.size();
                mesh_data->nn_edges[object_index] = vect_res_object1_edges.size();
                mesh_data->nn_faces[object_index] = vect_res_object1_faces.size();

                //vect_objects_index         .push_back(vect_oioa[0].object_index); // Индекс объекта можно взять из первого результата
                vect_objects_count_of_verts.push_back(vect_res_object1_verts.size());
                vect_objects_count_of_edges.push_back(vect_res_object1_edges.size());
                vect_objects_count_of_faces.push_back(vect_res_object1_faces.size());

                vect_objects_offsets.insert(vect_objects_offsets.end(), set_object1_offsets_indexes.begin(), set_object1_offsets_indexes.end());
                vect_objects_verts.insert(vect_objects_verts.end(), vect_res_object1_verts.begin(), vect_res_object1_verts.end());
                vect_objects_edges.insert(vect_objects_edges.end(), vect_res_object1_edges.begin(), vect_res_object1_edges.end());
                vect_objects_faces.insert(vect_objects_faces.end(), vect_res_object1_faces.begin(), vect_res_object1_faces.end());
              }

              // записать результаты vertices, edges и faces
              if (vect_objects_verts.size() == 0) {
                mesh_data->nn_offsets_indexes = NULL;
                mesh_data->vertices = NULL;
                mesh_data->edges = NULL;
                mesh_data->faces = NULL;
              } else {
                mesh_data->nn_offsets_indexes = NULL;
                mesh_data->vertices = NULL;
                mesh_data->edges = NULL;
                mesh_data->faces = NULL;

                if (vect_objects_offsets.size() > 0) {
                  mesh_data->nn_offsets_indexes = (int*)malloc(sizeof(int) * vect_objects_offsets.size());
                }

                if (vect_objects_verts.size() > 0) {
                  mesh_data->vertices = (float*)malloc(sizeof(float[3]) * vect_objects_verts.size());
                }

                if (vect_objects_edges.size() > 0) {
                  mesh_data->edges = (int*)malloc(sizeof(int[2]) * vect_objects_edges.size());
                }

                if (vect_objects_faces.size() > 0) {
                  mesh_data->faces = (int*)malloc(sizeof(int[3]) * vect_objects_faces.size());
                }

                // Загрузить данные объектов в выходной результат
                for (int I = 0; I <= (int)vect_objects_offsets.size() - 1; I++) {
                  mesh_data->nn_offsets_indexes[I] = vect_objects_offsets[I];
                }
                for (int I = 0; I <= (int)vect_objects_verts.size() - 1; I++) {
                  mesh_data->vertices[I * 3 + 0] = vect_objects_verts[I][0];
                  mesh_data->vertices[I * 3 + 1] = vect_objects_verts[I][1];
                  mesh_data->vertices[I * 3 + 2] = vect_objects_verts[I][2];
                }
                for (int I = 0; I <= (int)vect_objects_edges.size() - 1; I++) {
                  mesh_data->edges[I * 2 + 0] = vect_objects_edges[I][0];
                  mesh_data->edges[I * 2 + 1] = vect_objects_edges[I][1];
                }
                for (int I = 0; I <= (int)vect_objects_faces.size() - 1; I++) {
                  mesh_data->faces[I * 3 + 0] = vect_objects_faces[I][0];
                  mesh_data->faces[I * 3 + 1] = vect_objects_faces[I][1];
                  mesh_data->faces[I * 3 + 2] = vect_objects_faces[I][2];
                }
              }
            }
          } else {
            mesh_data = mesh_data;
          }
          return mesh_data;
        }


        /*
        Background for Sverchok Node Straight Skeleton 2D Extrude.
        //<image url="..\code_images\file_0032.png" scale="1.0"/>
        */
        DLLEXPORT MESH_DATA2* straight_skeleton_2d_extrude(
          // К чему относятся входные данные. <image url="..\code_images\file_0012.png" scale="1.0"/>
          // Контуры объектов считаются независимо друт от друга
          // Контуры планов могут влиять друг на друга при отрицательных offset-ах. Контуры с положительными offset будут рассчитываться независимо.
          // Также при отрицательных offset инвертируются holes и рассчитываются как положительные offset. Требуется избегать расчётов Straight Skeleton-ов
          // при вложенных контурах, т.к. это приводит к зависанию расчёта (есть такой глюк). Контуры в Straight Skeleton добавляются независимо,
          // поэтому можно случайно добавить контуры внутренних объектов (что с типами Polygon_with_holes_2 не происходит, т.к. у них есть деление на outer_boundary
          // и holes.
           int in_count_of_objects,        // Количество независимых исходных объектов
          // Shape mode: <image url="..\code_images\file_0020.png" scale=".1"/>
           int in_shapes_mode[],                // Тип обработки контуров объекта: 0-FULL_MODE (полная обработка), 1-EXCLUDE_HOLES (только внешние контуры, отверстия не учитываются), 2-INVERT_HOLES (отключить внешние контуры)
          bool in_restrict_height_of_extrude[], // Ограничивать ли объект по высоте. true - ограничивать, false - не ограничивать
        double in_height_of_extrude[],          // высоты объектов
        double in_angles[],                     // углы по вершинам контуров (являются одной цепочкой, которую надо поделить на части из параметра vertices_in_contours)
           int in_count_of_planes[],            // Количество планов в объектах (размер массива соответствует in_count_of_objects)
           int in_contours_in_planes[],         // Количество контуров в планах. Первый контур плана всегда outer, остальные - holes. Вложенность планов одного объекта в упаковке не учитывается.
           int in_vertices_in_contours[],       // Количество вершин в контуре.
         float in_vertices[][3],                // вершины контуров (являются одной цепочкой, которую надо поделить на части из параметра vertices_in_contours
          bool only_tests_for_valid,            // Проверка всех контуров на валидность без расчёта. (true - только проверить на валидность, false - проверить на валидность и посчитать)
          bool force_z_zero,                    // Обязательное преобразование Z=0 при любых обстоятельствах (иногда при редактировании модели можно двигать курсор в объёме в плоскости XY и получать Z "слегка" отличный от 0).
           int result_join_mode,                // Что сделать с полученными mesh: # 0 - split - разделить на отдельные mesh-объекты, 1 - keep - оставить в своих объектах, 2 - merge all meshes - объеденить в один большой объект. В этом ноде дополнительные параметры к объектам не нужны
          bool verbose                          // Подробный вывод
        ) {
          if (verbose == true) {
#ifdef _DEBUG
            printf("\nstraight_skeleton_2d_extrude in DEBUG MODE");
#endif
            printf("\nstraight_skeleton_2d_extrude");
            printf("\nmode verbose is on");
          }

          // TODO:: Переделать altitude на преобразование Matrix
          // https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Aff__transformation__3.html
          CGAL::Aff_transformation_3<Kernel> ac3(
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            //0.0, 0.0, 0.0, 1.0,
            1.0
          );

          MESH_DATA2* mesh_data = new MESH_DATA2();

          CGAL::Real_timer timer;
          timer.start();

          // Структура для хранения соответствия контура и его углов
          // Это позволит делать инверсию объекта с сохранением углов на всех контурах.
          // В некоторых случаях эффект может быть странным, но само сохранение важнее.
          struct C1A1 /* Contour1, Angles1 */ {
            ContourPtr c1;
            std::vector<FT> angles;
            C1A1(ContourPtr _c1, std::vector<FT> _angles)
              : c1(_c1), angles(_angles) {
            }
            void reverse() {
              c1->reverse_orientation();
              std::reverse(angles.begin(), angles.end());
            }
            C1A1 reversed() {
              ContourPtr _c1(std::make_shared<Contour>(c1.get()->begin(), c1.get()->end()));
              std::vector<FT> _angles(angles);
              C1A1 c1a1(_c1, _angles);
              c1a1.reverse();
              return c1a1;
            }
          };

          struct PWH_CA /*Polygon_with_holes_2, Contours, Angles*/ {
            Polygon_with_holes_2& pwh;
            std::vector<C1A1>& vect_c1a1;
            PWH_CA(Polygon_with_holes_2& _pwh, std::vector<C1A1>& _vect_c1a1)
              :pwh(_pwh), vect_c1a1(_vect_c1a1) {
            }
          };

          struct ClsObject {
            // Исходный индекс объекта
            int object_index;
            int shape_mode;
            bool restrict_height_of_extrude;
            float height_of_extrude;

            // Валиден ли результат работы функции Extrude
            bool ss_extrude_is_valid;

            // Исходные точки, переданные на обработку
            std::vector<Point_2> vect_points;

            // Полигоны для расчётов (позиции полигонов могут отличаться от позиций контуров в source_planes)
            std::vector<std::shared_ptr<Polygon_with_holes_2>> vect_polygons_with_holes;

            // Соответствие контуров и их углов
            std::vector<std::vector<C1A1>> planes;

            std::vector<Mesh> vect_sm; // результат расчёта объекта в виде список Mesh (контуров объекта может быть несколько по числу Polygon_with_holes_2 или может получаться несколько при инверсии объекта)

            ClsObject(int _object_index, int _shape_mode, bool _restrict_height, float _height, std::vector<Point_2> _vect_points, std::vector<std::vector<C1A1>> _source_planes)
              :vect_points(_vect_points), planes(_source_planes) {
              object_index = _object_index;
              shape_mode = _shape_mode;
              restrict_height_of_extrude = _restrict_height;
              height_of_extrude = _height;
              ss_extrude_is_valid = false;
            }
          };

          std::vector<ObjectError> res_errors; // Собирать ошибки в объектах по контурам, когда попытка использовать контур привела к его исключению из расчёта.

          int v_pos = 0;
          int general_count_of_vertices = 0;
          std::map<int, std::vector<Point_2>> contours;
          std::vector<ClsObject> objects; // Информация по объектам будет нужна и после расчёта (нужно помнить начальные индексы входных объектов, потому что иногда объект полностью выпадает из расчёта, а нужно передать инфу, что в объекте ничего нет и его ошибки
          {
            {
              int plane_cursor = 0;
              int contour_cursor = 0;
              int vert_cursor = 0;
              int angle_cursor = 0;

              for (int I = 0; I <= in_count_of_objects - 1; I++) {
                double* object_angles = &in_angles[vert_cursor];
                std::vector<std::vector<C1A1>> planes;
                std::vector<Point_2> vect_points; // Запомнить все точки, которые переданы в этот объект
                int count_of_planes1 = in_count_of_planes[I];
                for (int IJ = 0; IJ <= count_of_planes1 - 1; IJ++, plane_cursor++) {
                  std::vector<C1A1> plane1;
                  for (int IJK = 0; IJK <= in_contours_in_planes[plane_cursor] - 1; IJK++, contour_cursor++) {
                    ContourPtr contour1_points = std::make_shared<Contour>();
                    std::vector<FT> contour1_angles;
                    for (int IJKL = 0; IJKL <= in_vertices_in_contours[contour_cursor] - 1; IJKL++, vert_cursor++, angle_cursor++) {
                      bool is_z_not_zero = false;
                      if (force_z_zero == false) {
                        if (is_zero(in_vertices[vert_cursor][2]) == false) {
                          // Проверка if Z==0
                          is_z_not_zero = true;
                        }
                      }
                      Point_2 p2(in_vertices[vert_cursor][0], in_vertices[vert_cursor][1]);
                      contour1_points->push_back(p2);
                      vect_points.push_back(p2);
                      contour1_angles.push_back(in_angles[angle_cursor]);
                    }
                    if (contour1_points->size() < 3) {
                      // Collect wrong contour:
                      res_errors.push_back(ObjectError(I, contour1_points, VAL2STR(Err::_0009)". Polygon contains less 3 points."));
                      continue;
                    }
                    if (contour1_points->size() > 0) {
                      plane1.push_back(C1A1(contour1_points, contour1_angles));
                    }
                  }
                  if (plane1.size() > 0) {
                    // Проверить plane.
                    // Его внешний контур и все holes должены быть is_simple. Если он не is_simple, то такой план полностью исключается - это проверяется тут,
                    // Если holes не is_simple, то он исключается.
                    // Если что-то из внутренних пересекается с внешним, то внутренний контур отменяется. - это проверяется при расчёте Straight Skeleton позже.
                    auto& boundary = plane1[0];
                    Polygon_2 p_boundary(boundary.c1->begin(), boundary.c1->end());
                    if (boundary.c1->is_simple() == false) {
                      res_errors.push_back(ObjectError(I, plane1[0].c1, VAL2STR(Err::_0010)". Boundary polygon is not simple."));
                      continue;
                    } else {
                      // Проверить что holes is_simple.
                      std::vector<C1A1> allowed_contours;

                      for (auto c1_it = plane1.begin(); c1_it != plane1.end(); ++c1_it) {
                        if (c1_it == plane1.begin()) {
                          continue;
                        }

                        Contour c1(c1_it->c1->begin(), c1_it->c1->end());
                        // Проверку на is_simple надо делать раньше, чем is_clockwise_oriented (в отладке тоже проверяется is_clockwise_oriented и, если нет, то падает с исключением ).
                        if (c1.is_simple() == false) {
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Err::_0010)". Hole is not simple."));
                          continue;
                        }
                        if (c1.is_clockwise_oriented() == true) {
                          c1.reverse_orientation();
                        }
                        Polygon_2 p_hole1 = Polygon_2(c1.begin(), c1.end());
                        // Проверить, что все точки p_hole1 находятся внутри boundary. (есть исключение, которое рассматривается немного ниже, но точку снаружи всё ранво надо проверять. Также это может быть быстрее, чем проверять пересечения.)
                        // <image url="..\code_images\file_0028.png" scale=".2"/>
                        bool is_hole_inside = true;
                        for (auto& hole1_point : p_hole1) {
                          // https://stackoverflow.com/questions/64046074/is-there-a-cgal-function-that-checks-if-a-point-is-inside-a-linear-polygon-with
                          if (CGAL::POSITIVE != CGAL::oriented_side(hole1_point, p_boundary)) {
                            is_hole_inside = false;
                            break;
                          }
                        }
                        if (is_hole_inside == false) {
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Err::_0023)". Hole is not inside boundary."));
                          continue;
                        }

                        // В продолжении проверки что точки holes находятся внутри объекта. Исключение, когда точки находятся внутри, но сами контуры пересекаются:
                        // https://stackoverflow.com/questions/65021923/is-there-a-cgal-function-for-finding-all-intersecting-points-between-a-2d-line
                        // https://doc.cgal.org/latest/Kernel_23/group__intersection__linear__grp.html
                        // <image url="..\code_images\file_0028.png" scale=".2"/>
                        // Пример отображения такого пересечения <image url="F:\Enternet\2024\24.11.10\SVCGAL\ctSVCGAL\code_images\file_0033.png" scale=".2"/>
                        std::vector<Point_2> c_intersect;
                        for (auto& boundary_edge_it = p_boundary.edges_begin(); boundary_edge_it != p_boundary.edges_end(); ++boundary_edge_it) {
                          for (auto& hole_edge_it = p_hole1.edges_begin(); hole_edge_it != p_hole1.edges_end(); ++hole_edge_it) {
                            const auto res = CGAL::intersection(*boundary_edge_it, *hole_edge_it);
                            if (res) {
                              if (const Segment_2* s = std::get_if<Segment_2>(&*res)) {
                                c_intersect.push_back(s->point(0));
                                c_intersect.push_back(s->point(1));
                              } else {
                                const Point_2* p = std::get_if<Point_2 >(&*res);
                                c_intersect.push_back(*p);
                              }
                            }
                          }
                        }
                        if (c_intersect.size() > 0) {
                          res_errors.push_back(ObjectError(I, c_intersect, VAL2STR(Err::_0031)". Hole intersects with boundary."));
                          continue;
                        }

                        // Проверить, что hole не пересекается с другими holes: <image url="..\code_images\file_0029.png" scale=".2"/>
                        bool is_objects_intersects = false;
                        try {
                          for (auto& ap1 : allowed_contours) {
                            Polygon_2 pap1(ap1.c1->begin(), ap1.c1->end());
                            pap1.reverse_orientation(); // Иногда такая проверка выдаёт исключение, но пока не знаю почему, поэтому она обёрнута в try/catch
                            is_objects_intersects = CGAL::do_intersect(p_hole1, pap1);
                            if (is_objects_intersects == true) {
                              break;
                            }
                          }
                        } catch (std::exception _ex) {
                          // Пока буду считать, что при любых неудачных попытках расчёта пересечения они пересекаются.
                          is_objects_intersects = true;
                        }
                        if (is_objects_intersects == true) {
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Err::_0024)". Hole intersect with hole."));
                          continue;
                        }

                        // Дополнительных проверок по пересечению внутренних контуров можно тут не делать, т.к. ещё будут проверки при построении skeleton
                        // и там могут быть выявлены всякие дополнительные пересечения holes как с внешним контуром, так и между собой.
                        allowed_contours.push_back(*c1_it);
                      }
                      allowed_contours.insert(allowed_contours.begin(), boundary);
                      // Если имеются допустимые контуры, то загрузить ими план
                      if (allowed_contours.size() > 0) {
                        plane1.clear();
                        plane1.assign(allowed_contours.begin(), allowed_contours.end());
                      }
                    }
                    // Если в плане что-то есть, то загрузить его в массив допустимых планов для текущего object1
                    if (plane1.size() > 0) {
                      planes.push_back(plane1);
                    }
                  }
                }
                // Даже если planes нет, но нужно загрузить сопутствующие параметры offsets:
                {
                  int shape_mode = in_shapes_mode[I];
                  bool restrict_height_of_extrude = in_restrict_height_of_extrude[I];
                  float height_of_extrude = in_height_of_extrude[I];
                  objects.push_back(ClsObject(I, shape_mode, restrict_height_of_extrude, height_of_extrude, vect_points, planes));
                }
              }
            }

            // Continue with good contours to load into Polygon_with_holes_2:
            for (auto& object1 : objects) {
              // По идее все planes тут >0, т.к. раньше была на это проверка.
              for (auto& plane1 : object1.planes) {
                if (plane1.size() > 0) {

                  if (plane1[0].c1->is_clockwise_oriented()) {
                    plane1[0].reverse();
                  }
                  std::shared_ptr<Polygon_with_holes_2> ppwh(new Polygon_with_holes_2(Polygon_2(plane1[0].c1->begin(), plane1[0].c1->end())));

                  for (int I = 1; I <= plane1.size() - 1; I++) {
                    if (plane1[I].c1->is_counterclockwise_oriented() == true) {
                      plane1[I].reverse();

                    }
                    ppwh->add_hole(*(plane1[I].c1));
                  }
                  //bool is_valid_pwh = true;
                  //try {
                  //  if (ppwh.get()->has_holes() == false) {
                  //    // https://github.com/CGAL/cgal/issues/5278
                  //    // По факту проверяется, что полигон без holes простой и замкнутый. Замкнутый тут при загрузке, простой проверяется раньше,
                  //    // поэтому полигон без holes можно пропустить.
                  //    //typedef Exact_predicates_exact_constructions_kernel kernel_type;
                  //    //typedef Polygon_2<kernel_type> polygon_type;
                  //    //bool valid = CGAL::is_valid_polygon(Polygon_2(ppwh.get()->outer_boundary().begin(), ppwh.get()->outer_boundary().end() ), typename Polygon_set_2<kernel_type>::Traits_2());
                  //  } else {
                  //    // Неофициальный метод проверки валидности полигона с holes: https://github.com/CGAL/cgal/issues/3667
                  //    // Может дать ложные срабатывания в E:\Enternet\2024\24.09.27\cgal\Arrangement_on_surface_2\include\CGAL\Arr_segment_traits_2.h:611 (2 и 5)
                  //    using PolygonSet = CGAL::Polygon_set_2<K>;
                  //    auto traits = std::make_unique<PolygonSet::Traits_2>();
                  //    // Иногда в режиме отладки эта функция выдаёт исключение.
                  //    is_valid_pwh = is_valid_polygon_with_holes(*ppwh.get(), *traits);
                  //  }
                  //} catch (std::exception _ex) {
                  //  is_valid_pwh = false;
                  //}
                  //if (is_valid_pwh == false) {
                  //  /* A valid polygon with holes is :
                  //   * 1 - Has empty or closed boundary and all the holes are closed
                  //   * 2 - The PWH is relatively simple polygon (holes are simple...)
                  //   * 3 - Has it's boundary oriented counterclockwise and the holes oriented
                  //   *     clockwise
                  //   * 4 - All the segments (boundary and holes) do not cross or intersect in their
                  //   *     relative interior
                  //   * 5 - The holes are on the interior of the boundary polygon if the boundary
                  //   *     is not empty
                  //   */
                  //  res_errors.push_back(ObjectError(object1.object_index, *ppwh.get(), VAL2STR(Err::_0022) ". Holes of the Polygon intersect amongst themselves or with outer boundary ()"));
                  //} else 
                  object1.vect_polygons_with_holes.push_back(ppwh);  // -- Подумать, как сопоставить pwh и исходные планы
                }
              }
            }

            // Проверить, что объекты между собой не пересекаются:
            for (auto& object1 : objects) {
              {
                std::set<int> intersected_polygons;
                for (int I = object1.vect_polygons_with_holes.size() - 1; I >= 1; I--) {
                  const Polygon_with_holes_2& pwh1 = *object1.vect_polygons_with_holes[I].get();
                  for (int IJ = I - 1; IJ >= 0; IJ--) {
                    if (intersected_polygons.find(IJ) != intersected_polygons.end()) {
                      continue;
                    }
                    const Polygon_with_holes_2& pwh2 = *object1.vect_polygons_with_holes[IJ].get();
                    bool is_objects_intersects = false;
                    try {
                      is_objects_intersects = CGAL::do_intersect(pwh1, pwh2);
                    } catch (std::exception _ex) {
                      // Пока буду считать, что при любых неудачных попытках расчёта пересечения они пересекаются.
                      is_objects_intersects = true;
                    }
                    if (is_objects_intersects == true) {
                      intersected_polygons.insert(I);
                      intersected_polygons.insert(IJ);
                      // Лучше бы как-то отметить, что полигон пересекается с другим полигоном, чем удалять полигон.
                      // TODO: подумать, как это лучше сделать (в дальнейшем этот параметр надо будет учитывать).
                      // Update: всё равно исходная позиция полигона вроде как не определяема и в разных системах может меняться, поэтому можно оставить пока так)
                      // <image url="..\code_images\file_0019.png" scale="1.0"/>
                      if (verbose == true) {
                        printf("\n" VAL2STR(Err::_0015) ". Polygons (%u, %u) INTERSECTS. Both marked as error", I, IJ);
                      }
                      break;
                    }
                  }
                }
                std::vector<int> intersected_indexes;
                for (auto& idx : intersected_polygons) {
                  intersected_indexes.push_back(idx);
                }
                sort(intersected_indexes.begin(), intersected_indexes.end(), [](const int v1, const int& v2) {return v1 > v2; });
                for (auto& I : intersected_indexes) {
                  res_errors.push_back(ObjectError(object1.object_index, *object1.vect_polygons_with_holes[I].get(), VAL2STR(Err::_0011)". Polygons intersects"));
                  object1.vect_polygons_with_holes.erase(object1.vect_polygons_with_holes.begin() + I);
                  object1.planes.erase(object1.planes.begin() + I);
                }
              }
            }

            if (verbose) {
              printf("\nSS Extrude objects to process: %zu", objects.size());
            }

            {
              int object1_index = 0;
              // Привести объект в соответствии с настройкой shape_mode.
              // Тип обработки контуров объекта: 
              // 0-ORIGINAL_BOUNDARIES (полная обработка), 
              // 1-EXCLUDE_HOLES (Keep only outer boundary), 
              // 2-INVERT_HOLES (Exclude outer boundary and fill holes),
              // Обратить внимание на УГЛЫ - УГЛЫ ЗАГРУЖЕНЫ для контуров, поэтому смена режиме не должна сказывать на обработке углов по контурам
              for (auto& object1 : objects) {
                int shapes_mode1 = in_shapes_mode[object1_index];
                if (shapes_mode1 == 0) {
                  // Оставить исходные объекты без изменений
                } else {
                  // Обработать их в соответствии с настройкой shapes_mode1

                  //ContourSequence contours_sorted_by_area;
                  std::vector<C1A1> contours_sorted_by_area;
                  for (auto& plane1 : object1.planes) {
                    for (auto& c1 : plane1) {
                      contours_sorted_by_area.push_back(c1);
                    }
                  }
                  // Сортировку по площади выполнить независимо от направления контура:
                  sort(contours_sorted_by_area.begin(), contours_sorted_by_area.end(), [](const C1A1& p1, const C1A1& p2) {return abs(p1.c1->area()) > abs(p2.c1->area()); });
                  std::vector<std::vector<C1A1>> stack_cs; // vect_collect_for_polygons
                  for (auto& c1 : contours_sorted_by_area) {
                    std::vector<C1A1> cs1;
                    cs1.push_back(c1);
                    stack_cs.push_back(cs1);
                  }

                  std::vector<std::vector<C1A1>> external_contours;
                  std::vector<std::vector<C1A1>> vect_cs_inverted_holes;
                  while (stack_cs.size() > 0) {
                    // Перенести последние контура, являющиеся отверстиями, в список результирующих контуров (они могут содержать в себе бывшие внешние контура, кто-то уже добавился)
                    // <image url="..\code_images\file_0004.png" scale=".5"/>
                    for (int I = stack_cs.size() - 1; I >= 0; I--) {
                      std::vector<C1A1>& vcfp_I = stack_cs[I];
                      C1A1& vcfp_I0 = vcfp_I[0];
                      if (vcfp_I0.c1->is_clockwise_oriented() == true) {
                        vect_cs_inverted_holes.push_back(std::move(vcfp_I));
                        stack_cs.erase(stack_cs.begin() + I);
                      } else {
                        break;
                      }
                    }

                    // Если последний контур является Outer Boundary, то поискать в кого он может быть вложен.
                    // <image url="..\code_images\file_0005.png" scale="1.0"/>
                    for (int I = stack_cs.size() - 1; I >= 0; I--) {
                      std::vector<C1A1>& vcfp_I = stack_cs[I];
                      C1A1& vcfp_I0 = stack_cs[I][0];
                      // В OUTER контуры ничего добавлять на нужно. Только они добавляются в HOLE-контуры.
                      if (vcfp_I0.c1->is_counterclockwise_oriented() == true) {
                        bool is_parent_hole_found = false;
                        if (stack_cs.size() > 2) {
                          Point_2 p0_I0 = vcfp_I0.c1->vertices().at(0);
                          for (int IJ = stack_cs.size() - 2; IJ >= 0; IJ--) {
                            std::vector<C1A1>& vcfp_IJ = stack_cs[IJ];
                            C1A1& vcfp_IJ0 = stack_cs[IJ][0];
                            if (vcfp_IJ0.c1->is_clockwise_oriented() == true) {
                              // https://stackoverflow.com/questions/64046074/is-there-a-cgal-function-that-checks-if-a-point-is-inside-a-linear-polygon-with
                              // Если хотя бы одна точка находится внутри, то считаем. что и весь контур находится внутри (Предыдущие проверки должны били исключить пересечения до этого)
                              C1A1& c_IJ = vcfp_IJ0;
                              const C1A1& p2 = c_IJ.c1->is_clockwise_oriented() ? c_IJ.reversed() : c_IJ;
                              if (CGAL::POSITIVE == CGAL::oriented_side(p0_I0, *(p2.c1.get()))) {
                                for (auto& c1 : vcfp_I) {
                                  vcfp_IJ.push_back(std::move(c1));
                                }
                                stack_cs.erase(stack_cs.begin() + I);
                                is_parent_hole_found = true;
                                break;
                              }
                            }
                          }
                        }
                        if (is_parent_hole_found == false) {
                          // <image url="..\code_images\file_0006.png" scale="1.0"/>
                          // Если не найдена HOLE, в которой он мог бы находится, значит это общий внешний контур.
                          external_contours.push_back(std::move(vcfp_I));
                          stack_cs.erase(stack_cs.begin() + I);
                        }
                      } else {
                        // Не должно происходить, т.к. все не is_counterclockwise_oriented должны были быть удалены с конца стека на прошлом цикле, но не проверить нельзя.
                        // UPDATE: очень может быть, если это не первый элемент в цикле и перед этим мог отработать is_counterclockwise_oriented, он был соотнесён со своим holes
                        // и теперь настала очередь is_counterclockwise_oriented.
                        break;
                      }
                    }
                  }
                  // Удалить исходные полигоны объекта и в зависимости от shapes_mode1 выбрать либо inverted_holes или external_contours
                  object1.vect_polygons_with_holes.clear();
                  object1.planes.clear();

                  if (shapes_mode1 == 1 /*EXCLUDE_HOLES*/) {
                    for (auto& ec1 : external_contours) {
                      if (ec1[0].c1->is_clockwise_oriented() == true) {
                        ec1[0].reverse();
                      }
                      Polygon_2 pg2(ec1[0].c1->begin(), ec1[0].c1->end());
                      object1.vect_polygons_with_holes.push_back(std::make_shared<Polygon_with_holes_2>(pg2));
                      std::vector<C1A1> vect_c1a1(ec1.begin(), ec1.begin() + 1);
                      object1.planes.push_back(vect_c1a1);
                    }
                  } else if (shapes_mode1 == 2 /*INVERT_HOLES*/) {
                    for (auto& ec1 : vect_cs_inverted_holes) {
                      if (ec1[0].c1->is_clockwise_oriented() == true) {
                        ec1[0].reverse();
                      }
                      //Повнимательней тут с загрузкой контуров
                      Polygon_2 plg(ec1[0].c1->begin(), ec1[0].c1->end());
                      std::vector<C1A1> vect_c1a1(ec1.begin(), ec1.begin() + 1);
                      Polygon_with_holes_2 pwh = Polygon_with_holes_2(plg);

                      for (auto ec_it = ec1.begin() + 1; ec_it != ec1.end(); ++ec_it) {
                        if (ec_it->c1->is_counterclockwise_oriented() == true) {
                          ec_it->reverse();
                        }
                        Contour hole1(ec_it->c1->begin(), ec_it->c1->end());
                        pwh.add_hole(hole1);
                        vect_c1a1.push_back(*ec_it);
                      }
                      object1.planes.push_back(vect_c1a1);
                      object1.vect_polygons_with_holes.push_back(std::make_shared<Polygon_with_holes_2>(pwh));
                    }
                  } else {
                    // TODO:: Исключение
                  }
                }
                object1_index++;
              }
            }

            // TODO: уже не актуально проверять параметр тут. Он должен быть поближе к расчётам.
            if (only_tests_for_valid == false) {
              {
                // Рассчёт Straight Skeleton Extrude
                int threadNumbers = boost::thread::hardware_concurrency();
                boost::asio::thread_pool pool3(threadNumbers);
                boost::mutex mtx_;

                int threads_counts = 0;
                for (auto& object1 : objects) {
                  if (object1.vect_polygons_with_holes.size() == 0) {
                    res_errors.push_back(ObjectError(object1.object_index, object1.vect_points, "Non valid contour of object"));
                  } else {
                    for (auto& pwh : object1.vect_polygons_with_holes) {
                      object1.vect_sm.push_back(Mesh());
                    }
                    for (int pwh_index = 0; pwh_index <= (int)object1.vect_polygons_with_holes.size() - 1; pwh_index++) {
                      threads_counts++;
                      boost::asio::post(pool3, [&object1, pwh_index, &res_errors, verbose, threads_counts, &mtx_] {
                        auto& pwh1 = object1.vect_polygons_with_holes[pwh_index];
                        CGAL::Real_timer timer1;
                        timer1.start();

                        bool skeleton_is_valid = false;
                        std::vector<std::vector<FT>> vect_angles;
                        for (auto& c1a1 : object1.planes[pwh_index]) {
                          vect_angles.push_back(c1a1.angles);
                        }

                        CGAL::Real_timer timer3;
                        timer3.start();
                        int vertices_count = 0;
                        try {
                          Polygon_with_holes_2 p2 = *(pwh1.get());
                          if (verbose) {
                            vertices_count += p2.outer_boundary().size();
                            for (auto& hole1 : p2.holes()) {
                              vertices_count += hole1.size();
                            }
                          }
                          if (object1.restrict_height_of_extrude == true) {
                            object1.ss_extrude_is_valid = CGAL::extrude_skeleton(p2, object1.vect_sm[pwh_index], CGAL::parameters::angles(vect_angles).maximum_height(object1.height_of_extrude).verbose(verbose));
                          } else {
                            object1.ss_extrude_is_valid = CGAL::extrude_skeleton(p2, object1.vect_sm[pwh_index], CGAL::parameters::angles(vect_angles).verbose(verbose));
                          }
                        } catch (std::exception _ex) {
                          object1.ss_extrude_is_valid = false;
                        }
                        timer1.stop();
                        if (verbose) {
                          mtx_.lock();
                          printf("\n" VAL2STR(Err::_0027) ". SS Extrude: Thread: %i, Object: %i, ss_valid - %d, contours vertices: %u, extrude time - %.5g sec.", threads_counts, object1.object_index, object1.ss_extrude_is_valid, vertices_count, timer1.time());
                          mtx_.unlock();
                        }
                        if (object1.ss_extrude_is_valid == false) {
                          res_errors.push_back(ObjectError(object1.object_index, object1.vect_points, "Straight Skeleton Extrude failed"));
                        }

                        if (timer1.is_running() == true) {
                          timer1.stop();
                        }
                        }  // boost::asio::post
                      );
                    }
                  }
                }
                pool3.join();

                // Восстановить последовательность результатов расчётов offset по исходным позициям offset (как их задал пользователь):
                sort(res_errors.begin(), res_errors.end(), [](const ObjectError& oap1, const ObjectError& oap2) {return oap1.object_index < oap2.object_index; });
              }
            }

            {
              if (verbose == true) {
                unsigned int errors_size = res_errors.size();
                if (errors_size == 0) {
                  printf("\nStraight Skeleton Offsets tests finished. Errors not found.");
                } else {
                  printf("\nStraight Skeleton Offsets tests finished. Errors count %u", errors_size);
                }
              }
            }
            mesh_data->has_error = false;
            mesh_data->str_error = NULL;
            // Расчёт offset заканчивается тут.

          }

          if (timer.is_running() == true) { // it can be stopped in catch
            timer.stop();
          }
          if (verbose == true) {
            printf("\nAll Extrude computation of objects took %.5g sec.", timer.time());
          }

          // Преобразовать результат Extrude в полигональный Mesh и подготовить ошибки для возврата пользователю.

          mesh_data->nn_source_objects_count = objects.size(); // Значение, используется при подготовке данных по ошибкам исходных объектов
          {
            // Собрать все зарегистрированные ошибки из объекта res_errors:
            // 1. Сгруппировать ошибки по объектам:
            std::map<int, std::vector<ObjectError>> map_res_errors;
            // У каждого объекта должен быть список ошибок, даже если список нулевой длины, поэтому нужно сначала создать списки,
            // а потом их заполнить. Если res_errors будет пустым, то не получится сообщить КАЖДОМУ объекту, что у него 0 ошибок...
            for (auto& object1 : objects) {
              map_res_errors[object1.object_index] = std::vector<ObjectError>();
            }
            // ... Если ошибок не будет, то map_res_errors будет пустой
            for (auto& error1 : res_errors) {
              map_res_errors[error1.object_index].push_back(error1);
            }
            std::vector<int>                vect_errors_count; // Количество ошибок на объект
            std::vector<char*>              vect_description1_per_error1; // описания всех ошибок
            std::vector<int>                vect_contours_per_error1; // Количество контуров в ошибках (1-один контур, >1 обычно обычно все контуры объекта. Для понимания интерпретации результата)
            std::vector<int>                vect_vertices_per_contour1; // Количество вершин на контур
            std::vector<std::vector<double>> vect_vertices_of_errors;
            int error_counts_summa = 0;
            for (auto& kv : map_res_errors) {
              int errors_count = kv.second.size();
              error_counts_summa += errors_count;
              vect_errors_count.push_back(errors_count);
              for (auto& error1 : kv.second) {
                vect_description1_per_error1.push_back(error1.message);
                vect_contours_per_error1.push_back(error1.cs1.size());
                for (auto& c1 : error1.cs1) {
                  vect_vertices_per_contour1.push_back(c1.get()->size());
                  for (auto& v : *c1.get()) {
                    vect_vertices_of_errors.push_back(std::vector<double>{v.x(), v.y()});
                  }
                }
              }
            }

            // Определить количество памяти для выгрузки информации об ошибках:
            mesh_data->nn_errors_count = (int*)malloc(sizeof(int) * vect_errors_count.size());
            for (int I = 0; I <= vect_errors_count.size() - 1; I++) {
              mesh_data->nn_errors_count[I] = vect_errors_count[I];
            }

            if (error_counts_summa == 0) {
              mesh_data->nn_description1_per_error1 = NULL;
              mesh_data->nn_contours_per_error1 = NULL;
              mesh_data->nn_vertices_per_contour1 = NULL;
              mesh_data->vertices_of_errors = NULL;
            } else {

              mesh_data->nn_description1_per_error1 = (char**)malloc(sizeof(char*) * vect_description1_per_error1.size());
              mesh_data->nn_contours_per_error1 = (int*)malloc(sizeof(int) * vect_contours_per_error1.size());
              mesh_data->nn_vertices_per_contour1 = (int*)malloc(sizeof(int) * vect_vertices_per_contour1.size());
              mesh_data->vertices_of_errors = (float*)malloc(sizeof(float[3]) * vect_vertices_of_errors.size());
              for (int I = 0; I <= vect_description1_per_error1.size() - 1; I++) {
                mesh_data->nn_description1_per_error1[I] = vect_description1_per_error1[I];
              }
              for (int I = 0; I <= vect_contours_per_error1.size() - 1; I++) {
                mesh_data->nn_contours_per_error1[I] = vect_contours_per_error1[I];
              }
              for (int I = 0; I <= vect_vertices_per_contour1.size() - 1; I++) {
                mesh_data->nn_vertices_per_contour1[I] = vect_vertices_per_contour1[I];
              }
              for (int I = 0; I <= vect_vertices_of_errors.size() - 1; I++) {
                mesh_data->vertices_of_errors[I * 3 + 0] = vect_vertices_of_errors[I][0];
                mesh_data->vertices_of_errors[I * 3 + 1] = vect_vertices_of_errors[I][1];
                mesh_data->vertices_of_errors[I * 3 + 2] = 0.0;
              }
            }
          }

          std::vector<int> vect_source_objects_index; // Исходные индексы объектов, пришедших на обработку
          mesh_data->nn_offsets_counts = NULL; // Эти данные в функции ss_extrude не используются и клиентом не читаются
          mesh_data->nn_source_objects_indexes = (int*)malloc(sizeof(int) * mesh_data->nn_source_objects_count);

          for (auto& object1 : objects) {
            vect_source_objects_index.push_back(object1.object_index);
          }
          for (int I = 0; I <= vect_source_objects_index.size() - 1; I++) {
            if (mesh_data->nn_source_objects_indexes) {
              mesh_data->nn_source_objects_indexes[I] = vect_source_objects_index[I];
            }
          }

          if (only_tests_for_valid == false) {
            // Generate Results
            {
              // Данные для возврата пользователю.
#ifdef _DEBUG
              printf("\nSS Extrude, convert data: ");
#endif
              std::vector<int> vect_results_objects_index; // Исходные индексы объектов, пришедших на обработку
              std::vector<int> vect_objects_count_of_offsets; // Количество offsets объекта
              std::vector<int> vect_objects_count_of_verts; // Количество vertices объекта
              std::vector<int> vect_objects_count_of_edges; // Количество edges объекта
              std::vector<int> vect_objects_count_of_faces; // Количество faces объекта
              //std::vector<int> vect_objects_offsets; // Использованные offsets
              std::vector<std::vector<float>> vect_objects_verts; // Координаты vertices
              std::vector<std::vector<unsigned int>> vect_objects_edges; // indexes of edges (per object1)
              std::vector<std::vector<unsigned int>> vect_objects_faces; // indexes of faces (per object1)

              // Определить способ обработки объектов по join_mode
              std::map<int, std::vector<Mesh>> map_join_mesh;
              int split_index = 0;
              for (auto& object1 : objects) {
                if (result_join_mode == 0) {
                  for (auto& m1 : object1.vect_sm) {
                    map_join_mesh[split_index] = std::vector<Mesh>();
                    map_join_mesh[split_index].push_back(m1);
                    split_index++;
                  }
                } else if (result_join_mode == 1) {
                  map_join_mesh[object1.object_index] = object1.vect_sm;
                } else { // 2 - для всех остальных случаев
                  int index_of_single_object = 0;
                  if (map_join_mesh.find(index_of_single_object) == map_join_mesh.end()) {
                    map_join_mesh[index_of_single_object] = std::vector<Mesh>();
                  }
                  for (auto& m1 : object1.vect_sm) {
                    map_join_mesh[index_of_single_object].push_back(m1);
                  }
                }
              }

              mesh_data->nn_objects = map_join_mesh.size();
              mesh_data->nn_objects_indexes = (int*)malloc(sizeof(int) * mesh_data->nn_objects);

              mesh_data->nn_verts = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
              mesh_data->nn_edges = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
              mesh_data->nn_faces = (int*)malloc(sizeof(int) * mesh_data->nn_objects);

              int result_index = 0;
              for (auto& [object_index, vect_sm] : map_join_mesh) {
                mesh_data->nn_objects_indexes[result_index] = object_index;

                //std::set   <                     int>  set_object1_offsets; // offsets (per object1)
                std::vector<std::vector<       float>> vect_object1_verts; // Координаты vertices (per object1)
                std::vector<std::vector<unsigned int>> vect_object1_edges; // indexes of edges (per object1)
                std::vector<std::vector<unsigned int>> vect_object1_faces; // indexes of faces (per object1)

                CGAL::Real_timer timer2;
                timer2.start();

                ConvertMeshesIntoVerticesEdgesFaces(vect_sm, vect_object1_verts, vect_object1_edges, vect_object1_faces);
                timer2.stop();
                if (verbose) {
                  printf("\nObject %i: verts - %zu, edges - %zu, faces - %zu; triangulated time - %.5g sec.", object_index, vect_object1_verts.size(), vect_object1_edges.size(), vect_object1_faces.size(), timer2.time());
                }
                mesh_data->nn_verts[result_index] = vect_object1_verts.size();
                mesh_data->nn_edges[result_index] = vect_object1_edges.size();
                mesh_data->nn_faces[result_index] = vect_object1_faces.size();

                vect_objects_count_of_verts.push_back(vect_object1_verts.size());
                vect_objects_count_of_edges.push_back(vect_object1_edges.size());
                vect_objects_count_of_faces.push_back(vect_object1_faces.size());

                vect_objects_verts.insert(vect_objects_verts.end(), vect_object1_verts.begin(), vect_object1_verts.end());
                vect_objects_edges.insert(vect_objects_edges.end(), vect_object1_edges.begin(), vect_object1_edges.end());
                vect_objects_faces.insert(vect_objects_faces.end(), vect_object1_faces.begin(), vect_object1_faces.end());

                result_index++;
              }

              // записать результаты vertices, edges и faces
              if (vect_objects_verts.size() == 0) {
                mesh_data->nn_offsets_indexes = NULL;
                mesh_data->vertices = NULL;
                mesh_data->edges = NULL;
                mesh_data->faces = NULL;
              } else {
                mesh_data->vertices = (float*)malloc(sizeof(float[3]) * vect_objects_verts.size());
                mesh_data->edges = (int*)malloc(sizeof(int[2]) * vect_objects_edges.size());
                mesh_data->faces = (int*)malloc(sizeof(int[3]) * vect_objects_faces.size());

                // Загрузить данные объектов в выходной результат
                for (int I = 0; I <= vect_objects_verts.size() - 1; I++) {
                  mesh_data->vertices[I * 3 + 0] = vect_objects_verts[I][0];
                  mesh_data->vertices[I * 3 + 1] = vect_objects_verts[I][1];
                  mesh_data->vertices[I * 3 + 2] = vect_objects_verts[I][2];
                }
                for (int I = 0; I <= vect_objects_edges.size() - 1; I++) {
                  mesh_data->edges[I * 2 + 0] = vect_objects_edges[I][0];
                  mesh_data->edges[I * 2 + 1] = vect_objects_edges[I][1];
                }
                for (int I = 0; I <= vect_objects_faces.size() - 1; I++) {
                  // Normals orientation <image url="..\code_images\file_0022.png" scale=".3"/>
                  mesh_data->faces[I * 3 + 0] = vect_objects_faces[I][0];
                  mesh_data->faces[I * 3 + 1] = vect_objects_faces[I][2];
                  mesh_data->faces[I * 3 + 2] = vect_objects_faces[I][1]; // Не 0,1,2, 0,
                }
              }
            }
          } else {
            mesh_data = mesh_data;
          }
          return mesh_data;
        }

        DLLEXPORT void free_mem2(MESH_DATA2* md) {
#ifdef _DEBUG
          printf("\nfree_mem2 - start");
#endif
          if (md->str_error != NULL) {
            free((char*)md->str_error); // https://stackoverflow.com/questions/2819535/unable-to-free-const-pointers-in-c
            md->str_error = NULL;
          }

          if (md->nn_objects_indexes != NULL) {
            free(md->nn_objects_indexes);
            md->nn_objects_indexes = NULL;
          }
          if (md->nn_offsets_counts != NULL) {
            free(md->nn_offsets_counts);
            md->nn_offsets_counts = NULL;
          }
          if (md->nn_offsets_indexes != NULL) {
            free(md->nn_offsets_indexes);
            md->nn_offsets_indexes = NULL;
          }

          if (md->vertices != NULL) {
            free(md->vertices);
            md->vertices = NULL;
          }
          if (md->faces != NULL) {
            free(md->faces);
            md->faces = NULL;
          }
          if (md->edges != NULL) {
            free(md->edges);
            md->edges = NULL;
          }

          md->nn_verts = 0;
          md->nn_faces = 0;

          if (md->nn_source_objects_count > 0) {
            if (md->nn_source_objects_indexes != NULL) {
              free(md->nn_source_objects_indexes);
              md->nn_source_objects_indexes = NULL;
            }
            int nn_description1_per_error1_cursor = 0;
            for (int I = 0; I <= md->nn_source_objects_count - 1; I++) {
              for (int IJ = 0; IJ <= md->nn_errors_count[I] - 1; IJ++) {
                if (md->nn_description1_per_error1[nn_description1_per_error1_cursor] != NULL) {
                  free(md->nn_description1_per_error1[nn_description1_per_error1_cursor]);
                  nn_description1_per_error1_cursor++;
                }
              }
            }
            if (md->nn_description1_per_error1 != NULL) {
              free(md->nn_description1_per_error1);
              md->nn_description1_per_error1 = NULL;
            }

            if (md->nn_errors_count != NULL) {
              free(md->nn_errors_count);
              md->nn_errors_count = NULL;
            }
            if (md->nn_contours_per_error1 != NULL) {
              free(md->nn_contours_per_error1);
              md->nn_contours_per_error1 = NULL;
            }
            if (md->nn_vertices_per_contour1 != NULL) {
              free(md->nn_vertices_per_contour1);
              md->nn_vertices_per_contour1 = NULL;
            }
            if (md->vertices_of_errors != NULL) {
              free(md->vertices_of_errors);
              md->vertices_of_errors = NULL;
            }
          }
#ifdef _DEBUG
          printf("\nfree_mem2 - finish");
#endif
        }
      }
    }
  }
}