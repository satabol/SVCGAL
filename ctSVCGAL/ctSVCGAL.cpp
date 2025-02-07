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
#include <CGAL/Polygon_mesh_processing/connected_components.h>

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
#include "extrude_skeleton.h"

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
#include <algorithm>
//#include <list>

#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
//#include <boost/property_map/property_map.hpp> - исследовал возможность задания id для halfedges для Polygon_2, чтобы они сохранялить в SS - пока неудачно.

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
        _0033,
        _0034,
        _0035,
        _0036,
        _0037,
        _0038,
        _0039,
        _0040,
        _0041,
      };

#ifdef _DEBUG
      // Во время отладки требуется использовать тип не double/float, а rational. 
      // Тогда разные assertion, срабатывающие при _DEBUG не будут срабатывать и бросать исключения, 
      using K = CGAL::Exact_predicates_exact_constructions_kernel;
#else
      // При использовании этого типа в режиме _DEBUG возникает много "странных" исключений типа таких:
      // <image url="..\code_images\file_0038.png" scale="1.0"/>
      // Но в _RELEASE всё работает норм и можно использовать более простую математику float/double, также это должно ускорять рассчёты
      using K = CGAL::Exact_predicates_inexact_constructions_kernel;
#endif
      using FT = typename K::FT;
      using Point_2 = K::Point_2;
      using Segment_2 = K::Segment_2;
      using Line_2 = K::Line_2;
      using Point_3 = K::Point_3;
      using Vector_3 = K::Vector_3;

      using Polygon_2 = CGAL::Polygon_2<K>;
      using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;

      using Straight_skeleton_2 = CGAL::Straight_skeleton_2<K>;
      using Straight_skeleton_2_ptr = std::shared_ptr<Straight_skeleton_2>;

      using HDS = typename Straight_skeleton_2::Base;
      using HDS_Vertex_const_handle = typename HDS::Vertex_const_handle;
      using HDS_Halfedge_const_handle = typename HDS::Halfedge_const_handle;
      using HDS_Face_handle = typename HDS::Face_handle;

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
        int* nn_faces_indexes_counters = NULL; // Счётчик индексов вершин для каждой face. Количество вершин в faces является переменным (но не меньше 3-х)

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

      //typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
      typedef K::Point_2                    Point;
      typedef CGAL::Polygon_2<K>            Polygon_2;
      typedef CGAL::Polygon_with_holes_2<K> PolygonWithHoles;
      typedef std::shared_ptr<PolygonWithHoles> PolygonWithHolesPtr;
      typedef std::vector<PolygonWithHolesPtr> PolygonWithHolesPtrVector;

      //typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

      typedef CGAL::Polygon_2<K>  Contour;
      typedef std::shared_ptr<Contour> ContourPtr;
      typedef std::vector<ContourPtr>  ContourSequence;

      typedef CGAL::Straight_skeleton_2<K> Ss;
      typedef CGAL::Straight_skeleton_builder_traits_2<K>       SsBuilderTraits;
      typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits, Ss> SsBuilder;

      typedef CGAL::Polygon_offset_builder_traits_2<K>                    OffsetBuilderTraits;
      typedef CGAL::Polygon_offset_builder_2<Ss, OffsetBuilderTraits, Contour> OffsetBuilder;

//// *********** Данные для edgeId
//      // Тип для хранения id
//      typedef int IdType;
//
//      // Определяем структуру для хранения id
//      template<typename T>
//      struct EdgeIdProperty {
//        typedef IdType value_type;
//        typedef value_type& reference;
//        typedef T key_type;
//        typedef boost::writable_property_map_tag category;
//      };
//
//      // Функция для присвоения уникальных id рёбрам
//      void assign_edges_ids(const Polygon_2& polygon, std::map<typename Polygon_2::Edge_const_iterator, IdType>& edge_id_map, IdType& unique_is_edges) {
//        //IdType id = 0;
//        for (auto it = polygon.edges_begin(); it != polygon.edges_end(); ++it) {
//          edge_id_map[it] = unique_is_edges++;
//        }
//      }
//// *********** Данные для edgeId

//************  Данные для face в Mesh
      /// <summary>
      /// Пользовательские данные face исходных граней SS у результирующего Mesh
      /// </summary>
      struct MeshFacePropertyStruct {
        /// <summary>
        /// Какому ss_id принадлежит этот face (задаётся)
        /// </summary>
        int ss_id;

        /// <summary>
        /// face_id исходного из SS с ss_id. (задаётся)
        /// </summary>
        int ss_face_id;

        /// <summary>
        /// Знак SS - определяет к какой части области по отношению к исходной фигуре находится SS 
        /// <image url="..\code_images\file_0078.png" scale="1.0"/>
        /// </summary>
        enum SS_SIGN {
          /// <summary>
          /// SS построен снаружи исходного контура на отрицательных offset
          /// </summary>
          NEGATIVE = -1,

          /// <summary>
          /// Не определено, обычно в процессе.
          /// </summary>
          UNDEFINED=0,

          /// <summary>
          /// SS построен внутри исходного контура на положительных offset. 0 входит в POSITIVE.
          /// </summary>
          POSITIVE = +1,
        };

        /// <summary>
        /// Тип SS (отрицательный или положительный SS - если он получен с использованием не отрицательных offset, то SS положительный. Если использовались offset <0, то SS отрицательный. Примечание - в OIOA есть ещё типы отрицательных offset)
        /// </summary>
        SS_SIGN ss_sign;

        /// <summary>
        /// уникальный face_id (рассчитывается)
        /// </summary>
        int mesh_face_id;
        MeshFacePropertyStruct(int _ss_id, SS_SIGN _ss_sign, int _ss_face_id, int _mesh_face_id)
          :ss_id(_ss_id), ss_sign(_ss_sign), ss_face_id(_ss_face_id), mesh_face_id(_mesh_face_id) {
          
        }
        MeshFacePropertyStruct()
          :ss_id(-1 /*Не задано*/ ), ss_sign(MeshFacePropertyStruct::SS_SIGN::UNDEFINED), ss_face_id(-1 /*Не задано*/), mesh_face_id(-1 /*Не задано*/) {
        }
      };
      
      typedef CGAL::SM_Face_index face_descriptor;

      // Свойство для граней
      template<class T>
      using MeshFacePropertyMap = Mesh::Property_map<face_descriptor, T>;
//************ /Данные для face в Mesh

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
        using Fb = CGAL::Constrained_triangulation_face_base_2<Geom_traits>;
        using TDS = CGAL::Triangulation_data_structure_2<Vb, Fb>;
        using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
        using CDT = CGAL::Constrained_Delaunay_triangulation_2<Geom_traits, TDS, Itag>;
        using CDT_Vertex_handle = typename CDT::Vertex_handle;
        using CDT_Face_handle = typename CDT::Face_handle;
        typedef std::pair<CDT_Face_handle, int> Edge;


        CDT cdt;

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

      /// <summary>
      /// Сравнение умных указателей. Шаблон взят из: https://stackoverflow.com/questions/33650944/compare-shared-ptr-object-equality
      /// </summary>
      /// <typeparam name="T"></typeparam>
      /// <typeparam name="U"></typeparam>
      /// <param name="a"></param>
      /// <param name="b"></param>
      /// <returns></returns>
      //template<class T>
      //bool compare_shared_ptr(const std::shared_ptr<T>& a, const std::shared_ptr<T>& b) {
      //  bool res = false;
      //  if (a == b) {
      //    res = true;
      //  }
      //  //if (a && b) {
      //  //  res = (*a == *b);
      //  //}
      //  return res;
      //}


      template <typename SLSFacePoints, typename PointRange, typename FaceRange, typename NamedParameters>
      void triangulate_skeleton_face(
        SLSFacePoints& face_points,
        const bool invert_faces,
        PointRange& points,
        FaceRange& faces,
        const NamedParameters& np
      ) {

        using Default_kernel = typename Kernel_traits<Point>::type;
        using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
          NamedParameters,
          Default_kernel>::type;
        using PK = CGAL::Projection_traits_3<Geom_traits>;
        using PVbb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, PK>;
        using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
        using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
        using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
        using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
        using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
        using PCDT_Vertex_handle = typename PCDT::Vertex_handle;
        using PCDT_Face_handle = typename PCDT::Face_handle;

        CGAL_precondition(face_points.size() >= 3);

        // shift once to ensure that face_points[0] and face_points[1] are at z=0 and thus the normal is correct
        std::rotate(face_points.rbegin(), face_points.rbegin() + 1, face_points.rend());
        CGAL_assertion(face_points[0][2] == 0 && face_points[1][2] == 0);

        const Vector_3 n = CGAL::cross_product(face_points[1] - face_points[0], face_points[2] - face_points[0]);
        PK traits(n);
        PCDT pcdt(traits);

        try {
          pcdt.insert_constraint(face_points.begin(), face_points.end(), true /*close*/);
        } catch (const typename PCDT::Intersection_of_constraints_exception&) {
          std::cerr << "Warning: Failed to triangulate skeleton face" << std::endl;
          return;
        }

        std::size_t id = points.size(); // point ID offset (previous faces inserted their points);
        for (PCDT_Vertex_handle vh : pcdt.finite_vertex_handles()) {
          points.push_back(pcdt.point(vh));
          vh->info() = id++;
        }

#ifdef CGAL_SLS_DEBUG_DRAW
        // CGAL::draw(pcdt);
#endif

        std::unordered_map<PCDT_Face_handle, bool> in_domain_map;
        boost::associative_property_map< std::unordered_map<PCDT_Face_handle, bool> > in_domain(in_domain_map);

        CGAL::mark_domain_in_triangulation(pcdt, in_domain);

        for (PCDT_Face_handle f : pcdt.finite_face_handles()) {
          if (!get(in_domain, f))
            continue;

          // invert faces for exterior skeletons
          if (invert_faces)
            faces.push_back({ f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info() });
          else
            faces.push_back({ f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info() });
        }
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
        /// <param name="pwh1">[in] Входной полигон PWH (Polygon with Holes)</param>
        /// <param name="altitude">[in] Высота, на которой надо разместить mesh</param>
        /// <param name="sm1">[out] Результирующая Mesh</param>
        static bool ConvertPolygonWithHoles2IntoMesh(CGAL::Straight_skeleton_extrusion::internal::Polygon_with_holes_2& pwh1, FT altitude, Mesh& sm1) {
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

        /// <summary>
        /// Тип SS (NEGATIVE - построен на негативном отступе, 0 - не определено, POSITIVE - SS построен на положительном отступе)
        /// </summary>
        enum SS_TYPE {
          UNDEFINED = 0,
          
          /// <summary>
          /// Элемент получен из исходного PWH.
          /// </summary>
          POSITIVE = 1,
          
          /// <summary>
          /// Элемент получен добавлением виртуального внешнего 4-х гранного face.
          /// </summary>
          NEGATIVE_WITH_EXTERNAL_FRAME = -2,

          /// <summary>
          /// Используется при обработке отрицательного offset, когда элемент исходного контура является отверстием в исходном mesh. См. цикл рассчёта External Frame при отрицательном offset.
          /// </summary>
          NEGATIVE_INVERTED_HOLE = -1,
        };

        /// <summary>
        /// Object Index, Offset index, offset, Altitude, Polygon, ss, ss_id. Сопоставленные objects, planes, offsets, altitudes и помнить порядок в котором заданы offsets:
        /// </summary>
        struct OIOA {

          /// <summary>
          /// Исходный индекс объекта для которого строится offset
          /// </summary>
          int object_index;

          /// <summary>
          /// Индекс offset в списке offset (чтобы помнить порядок задания offset). Иногда может принимать значение "-1", что означает, что параметры offset недействительны,
          /// но по какой-то причине объект должен иметь этот instance, например, когда смешиваются точки пересечения offset-ов и точек объекта, то точки объекта тоже
          /// хранят в себе значение этого класса, но это только для совместимости в пределах списка)
          /// </summary>
          int offset_index;

          /// <summary>
          /// Величина рассматриваемого offset
          /// </summary>
          FT offset;

          /// <summary>
          /// Высота altitude для offset
          /// </summary>
          FT altitude;

          /// <summary>
          /// Straight Skeleton, соответствующий рассматриваемому offset-у (offset) (выставляется только после рассчёта SS).
          /// </summary>
          std::shared_ptr<Ss> ss{ nullptr };

          /// <summary>
          /// SS id, который соответствует текущему offset (не всегда одному offset соответствует один SS, бывает, что одному уровню соответствуют SS от разных подобъектов при положительных offset).
          /// -1 - значение SS ещё не было присвоено. Это нужно, чтобы при выполнении extrude между двумя offset с разными знаками можно было сопоставить соединяемые между собой рёбра. Такое возможно,
          /// потому что у соединяемых между собой нескольких SS-ов с разными знаками есть общие полурёбра, которые принадлежат разным SS. Именно через такие полурёбра и производится соединения SS с
          /// offset с разными знаками. Если offset-ы одного знака, то проблем с определением грани в общем-то нет, т.к. грань между двумя offsets общая, от одного SS.
          /// </summary>
          size_t ss_id = -1;

          /// <summary>
          /// Этот тип имеет смысл только при описании полигона с отрицательным offset. <image url="..\code_images\file_0044.png" scale="1.0"/>
          /// </summary>
          enum NEGATIVE_TYPE {
            /// <summary>
            /// Элемент получен добавлением виртуального внешнего 4-х гранного face.
            /// </summary>
            WITH_EXTERNAL_FRAME,
            /// <summary>
            /// Используется при обработке отрицательного offset, когда элемент исходного контура является отверстием в исходном mesh. См. цикл рассчёта External Frame при отрицательном offset.
            /// </summary>
            INVERTED_HOLE,
          } neg_type;

          /// <summary>
          /// Если будет получен результат для заданного offset, то сюда будут записаны контуры по результатам расчёта. (Это не PWH, а набор контуров)
          /// </summary>
          ContourSequence polygon_contours;

          /// <summary>
          /// Отображение HalfEdge контуров Straight Polygon (хранится снаружи этой структуры) для объекта с исходным индексом object_index на точку пересечения с offset
          /// </summary>
          std::unordered_map<HDS_Halfedge_const_handle, Point_2> map_halfedge_to_offset_points;

          OIOA(int _object_index, int _offset_index, FT _offset, FT _altitude)
            : object_index(_object_index), offset_index(_offset_index), offset(_offset), neg_type(NEGATIVE_TYPE::INVERTED_HOLE), altitude(_altitude){
            //object_index = _object_index;
            //offset_index = _offset_index;
            //offset = _offset;
            //altitude = _altitude;
            //neg_type = NEGATIVE_TYPE::INVERTED_HOLE;
          }

          OIOA(int _object_index, int _offset_index, FT _offset, NEGATIVE_TYPE _neg_type, FT _altitude)
            :object_index(_object_index), offset_index(_offset_index), offset(_offset), neg_type(_neg_type), altitude(_altitude) {
            //object_index = _object_index;
            //offset_index = _offset_index;
            //offset = _offset;
            //neg_type = _neg_type;
            //altitude = _altitude;
          }
        }; // OOAI

        /// <summary>
        /// Сопоставление SS, его полурёбер, точек пересечения offset с полурёбрами и способ получения SS (neg_type) при отрицательных offset. Используется для расчёта Extrude.
        /// </summary>
        struct SS__HalfEdge__Point_2 {

          /// <summary>
          /// Straight Skeleton, соответствующий параметрам offset-а oioa (выставляется только после рассчёта SS).
          /// </summary>
          std::shared_ptr<Ss> ss{ nullptr };
          /// <summary>
          /// id Straight Skeleton из ss. -1 - не установлено.
          /// </summary>
          size_t ss_id = -1;

          /// <summary>
          /// Этот тип имеет смысл только при описании полигона с отрицательным offset. Используется для удаления внешнего контура при построении extrude или контура Straight Skeleton.
          /// </summary>
          OIOA::NEGATIVE_TYPE neg_type;

          /// <summary>
          /// Сопоставление halfedge и точки пересечения halfedge с offset из oioa. Примечание: сумма halfedges не обязательно представляет собой какой-то один PWH. Это может быть и несколько PWH.
          /// </summary>
          std::unordered_map<HDS_Halfedge_const_handle, Point_2> map__halfedge__offset_points;

          SS__HalfEdge__Point_2(std::shared_ptr<Ss> _ss, size_t _ss_id, OIOA::NEGATIVE_TYPE _neg_type , std::unordered_map<HDS_Halfedge_const_handle, Point_2> _map__halfedge__offset_points)
            :ss{ _ss }, ss_id{ _ss_id }, neg_type {_neg_type}, map__halfedge__offset_points{_map__halfedge__offset_points}
          {

          }
        };

        /// <summary>
        /// Параметры offset. Из чего образуются контуры текущего offset.
        /// oioa1 - индексные атрибуты offset
        /// vect_pwh - вектор PWH из которых состоит offset
        /// vect__SS__HalfEdge__Point_2 - вектор SS, испрользуемых при построении контуров offset из oioa1
        /// </summary>
        struct OIOA_OFFSET_SS_PARAMS {

          /// <summary>
          /// Параметры отступа offset: Object Index, Offset index, offset, Altitude, Polygon. Сопоставленные objects, planes, offsets, altitudes и помнить порядок в котором заданы offsets:
          /// </summary>
          OIOA oioa1;

          /// <summary>
          /// Состав в полигонах текущего offset (с его параметрами). Имеет смысл, когда mode является edges или faces.
          /// </summary>
          std::vector<Polygon_with_holes_2> vect_pwh;

          /// <summary>
          /// Вектор точек для расчёта полигонов в mode Extrude и Straight_Skeleton.
          /// </summary>
          std::vector<std::vector<Point_3>> vect_vect_ss_face_points;

          /// <summary>
          /// результирующая mesh от вектора точек vect_points с дополнительными атрибутами для расчётов bevel и Extrude.
          /// </summary>
          Mesh mesh_ss_01_source; 

          /// <summary>
          /// Получается из mesh при merge всех vertices.
          /// </summary>
          Mesh mesh_ss_02_merged;
          
          /// <summary>
          /// Сопоставление вектора SS и его параметров текущему значению offset из oioa. Сделан вектор, потому что текущему значению offset могут соответствовать несколько SS
          /// с его точками пересечения (нескольким островам на одном offset могу соответствовать несколько Straight Skeleton-ов).
          /// </summary>
          std::vector<SS__HalfEdge__Point_2> vect__SS__HalfEdge__Point_2;

          ///// <summary>
          ///// Точки пересечения HalfEdge с offset всех граней Straight Skeleton (SS, задаётся снаружи).
          ///// </summary>
          //std::unordered_map<HDS_Halfedge_const_handle, Point_2> map_halfedge_to_offset_points;

          OIOA_OFFSET_SS_PARAMS(OIOA& _oioa1, std::vector<Polygon_with_holes_2>& _vect_pwh, std::vector<SS__HalfEdge__Point_2> _vect__SS__HalfEdge__Point_2)
            :oioa1(_oioa1), vect_pwh(_vect_pwh), vect__SS__HalfEdge__Point_2(_vect__SS__HalfEdge__Point_2) {
          }
        };


        /// <summary>
        /// Расчёт offsets. Напоминаю, что тут находятся не только исходные параметры, но и полученные путём преобразования параметров при подготовке данных для расчёта отрицательных offsets.
        /// Для начала надо получить один скелетон на plane.
        /// POMS - Positive Offset Multithread Struct. Keep single polygon with holes and a list of offsets and altitudes. Одиночный полигон с несколькими параметрами offset.
        /// </summary>
        struct POMS {
          /// <summary>
          /// Индекс исходного объекта. Иногда объект состоит из нескольких PWH, поэтому в векторе POMS может быть несколько экземпляром POMS с одинаковым индексом.
          /// </summary>
          int object_index = 0;

          /// <summary>
          /// PWH. Исходный полигон для которого рассчитывается Straight Skeleton (SS). Это только один из полигонов исходного объекта (хотя иногда и может и представлять целый объект)
          /// </summary>
          std::shared_ptr<Polygon_with_holes_2> polygon1_with_holes;


          /// <summary>
          /// Вектор параметров для нескольких offset для рассматриваемого PWH (даже если offset один, то его надо обернуть в вектор)
          /// </summary>
          std::vector<OIOA> vect_oioa;

          SS_TYPE ss_type;

          size_t ss_id = -1;

          /// <summary>
          /// Keep Struct Skeleton (or zero if failed to calc)
          /// </summary>
          std::shared_ptr<Ss> ss{ nullptr };

          /// <summary>
          /// результирующая mesh от вектора точек vect_points с дополнительными атрибутами для расчётов bevel и Extrude.
          /// </summary>
          Mesh mesh_ss_01_source;

          /// <summary>
          /// Получается из mesh при merge всех vertices.
          /// </summary>
          Mesh mesh_ss_02_merged;

          /// <summary>
          /// Сопоставление вектора SS и его параметров текущему значению offset из oioa. Сделан вектор, потому что текущему значению offset могут соответствовать несколько SS
          /// с его точками пересечения (нескольким островам на одном offset могу соответствовать несколько Straight Skeleton-ов).
          /// </summary>
          std::vector<SS__HalfEdge__Point_2> vect__SS__HalfEdge__Point_2;


          POMS(int _object_index, std::shared_ptr<Polygon_with_holes_2> _polygon1_with_holes, std::vector<OIOA> _vect_oioa, SS_TYPE _ss_type)
            :polygon1_with_holes(_polygon1_with_holes), vect_oioa(_vect_oioa), ss_type(_ss_type) {
            object_index = _object_index;
          }
        };

        /// <summary>
        /// Конвертировать набор Polygon_with_holes_2 в наборы контуров в виде полигональной сетки vertices, edges, faces, но faces отсутствуют, но этот параметр нужен полигональной сетке.
        /// </summary>
        /// <param name="_vect_sm_in">[in] Входной набор полигонов</param>
        /// <param name="vect_object1_verts">[out] Результат verts</param>
        /// <param name="vect_object1_edges">[out] Результат edges</param>
        /// <param name="vect_object1_faces">[out] Результат faces</param>
        static void ConvertPWH2IntoVerticesEdgesNoFaces(
          std::vector<OIOA_OFFSET_SS_PARAMS> _in_vect_pwh_offset_index_order,
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
                vect_points.push_back({ (float)CGAL::to_double(p1.x()), (float)CGAL::to_double(p1.y()), (float)CGAL::to_double(pwh1_offset.oioa1.altitude) });
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
                vect_object1_verts.push_back(std::vector<float>{ (float)CGAL::to_double(p1.x()), (float)CGAL::to_double(p1.y()), (float)CGAL::to_double(p1.z()) });
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
        /// Неочевидный корректный результат для mode faces:
        /// <image url="..\code_images\file_0033.png" scale="1.0"/>
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
          CGAL::Aff_transformation_3<K> ac3(
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
            EDGES,
            FACES,
            BEVEL,  // Последовательная заливка боковых поверхностей при переходе от одного offset к другому с учётом перепада высот (altitudes). По сути - почти Bevel, только с учётом алгоритма Straight Skeleton.
                    // см. функцию construct_lateral_faces. Тонкость в работе - возможны отрицательные altitudes.
            STRAIGHT_SKELETON, // Вывести схему Straight Skeleton (SS Faces). Если будут указаны offsets, то поделить всё на offset-ы (altitude не учитываются/игнорируются)
          }result_type;

          if (in_res_type == ResType::EDGES) {
            result_type = ResType::EDGES;
          } else if (in_res_type == ResType::BEVEL) {
            result_type = ResType::BEVEL;
          } else if (in_res_type == ResType::STRAIGHT_SKELETON) {
            result_type = ResType::STRAIGHT_SKELETON;
          } else {
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
            
            SS_TYPE ss_type;

            ClsObject(int _object_index, std::vector<OIOA> _vect_ooai, std::vector<ContourSequence> _source_planes, SS_TYPE _ss_type)
              :source_planes(_source_planes), vect_oioa(_vect_ooai), ss_type(_ss_type) {
              object_index = _object_index;
            }
          };

          std::vector<OIOA> res_contours; // Результаты расчёта всех контуров. Результаты будут сортироваться с учётом вложенности контуров по offset_index, в которой offset подавались на вход (по объектам) и с учётом параметра join_mode.
          std::vector<ObjectError> res_errors; // Собирать ошибки в объектах при загрузке исходных данных и когда попытка использовать контур привела к его исключению из расчёта. Относится только к исходным объектам.

          int v_pos = 0;
          int general_count_of_vertices = 0;
          std::map<int, std::vector<Point_2>> contours;
          std::vector<ClsObject> objects; // Информация по объектам будет нужна и после расчёта (нужно помнить начальные индексы входных объектов, потому что иногда объект полностью выпадает из расчёта, а нужно передать инфу, что в объекте ничего нет и его ошибки
          
          // Вектор параметров для рассчёта SS.
          std::vector<POMS> vect_polygon1_oioa;

          //// Карта уникальных идентификаторов исходных edges, используемых для построения SS (чтобы потом найти смежные faces для примыкающих друг к другу SS.
          //IdType unique_is_edges = 1;
          //std::map<typename Polygon_2::Edge_const_iterator, IdType> edge_id_map;

          //try 
          {
            {
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
                        // Проверить, что все точки p_hole1 находятся внутри boundary. Есть исключение, когда точки находятся внутри, но сами контуры пересекаются, но пока так, подумать над решением:
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
                      for (auto& cs : plane1) {
                        //assign_edges_ids( *cs.get(), edge_id_map, unique_is_edges);
                        auto edges_size = cs.get()->edges().size();
                        int II = 0;
                      }
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
                  objects.push_back(ClsObject(I, vect_ooai, planes, SS_TYPE::UNDEFINED) );
                }
              }
            }

            {
              // Представить список объектов в зависимости от параметра source_objects_join_mode
              std::vector<ClsObject> vect_split_mode;
              ClsObject objects_merge_mode(0, std::vector<OIOA>(), std::vector<ContourSequence>(), SS_TYPE::UNDEFINED);
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
                      auto obj1 = ClsObject(object_index, vect_oioa1, vect_cs, SS_TYPE::UNDEFINED );
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
                          ClsObject negObject(object1.object_index, object1.vect_oioa, std::vector<ContourSequence>(), SS_TYPE::UNDEFINED /*Пока как неопределённый тип. Он вычислиться позже*/);
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

                    // Остальные объекты отметить как POSITIVE:
                    for (auto& object1 : objects) {
                      object1.ss_type = SS_TYPE::POSITIVE;
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
                        // Контуры можно посчитать только положительные. Но раз инвертированные holes являются замкнутыми, то их теперь можно посчитать наравне с положительными,
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

                          // Отсортировать полигоны в топологической последовательности по вложенности. 
                          // Перед расчётом важно убедиться, что если какой-то полигон вложен в другой полигон, то такие вложенные полигоны не должны использоваться вместе для рассчёта external offset
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
                                bool is_parent_hole_found = false; // Изначально предполагается, что рассматриваемый контур является самым внешним, который никуда не вложен.
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
                          } // while (stack_cs.size() > 0) 

                          // Рассчитать вспомогательный внешний Frame для расчёта наибольшего отрицательного offset, в который войдут все внешние объекты.
                          // <image url="..\code_images\file_0009.png" scale=".2"/>
                          {
                            // std::min_element - потому что ищем максимум отрицательного отступа. Соответственно максимальное значение по модулю будет у минимального отрицательного элемента (-20, -5, -1, -10) - максимальное отрицательное числе =-20.
                            auto max_abs_negative_offset = abs(std::min_element(object1NegativeOffsetSS.object1.vect_oioa.begin(), object1NegativeOffsetSS.object1.vect_oioa.end(), [](const OIOA& o1, const OIOA& o2) { return o1.offset < o2.offset; })->offset);
                            std::vector<FT> vxmin, vymin, vxmax, vymax;
                            ContourSequence allowed_polygons; // Для каких полигонов удалось вычислить margin
                            for (auto& cs : external_contours) {
                              CGAL::Bbox_2 bbox = cs[0]->bbox();
                              // <image url="..\code_images\file_0010.png" scale=".2"/>
                              std::optional<FT> margin = CGAL::compute_outer_frame_margin(cs[0]->begin(), cs[0]->end(), max_abs_negative_offset);
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

                              FT xmin = *std::min_element(vxmin.begin(), vxmin.end());
                              FT ymin = *std::min_element(vymin.begin(), vymin.end());
                              FT xmax = *std::max_element(vxmax.begin(), vxmax.end());
                              FT ymax = *std::max_element(vymax.begin(), vymax.end());
                              // Внешний контур для охвата всех вложенных контуров (виртуальный контур, который позже будет использован для построения отрицательного offset)
                              ContourPtr frame_outer_boundary = std::make_shared<Contour>();
                              frame_outer_boundary->push_back(Point_2(xmin, ymin));
                              frame_outer_boundary->push_back(Point_2(xmax, ymin));
                              frame_outer_boundary->push_back(Point_2(xmax, ymax));
                              frame_outer_boundary->push_back(Point_2(xmin, ymax));
                              // <image url="..\code_images\file_0011.png" scale=".2"/>

                              // Собрать контуры вместе, чтобы потом построить из них SS с отрицательным offset.
                              ContourSequence cs_frame;
                              cs_frame.push_back(frame_outer_boundary);
                              std::shared_ptr<Polygon_with_holes_2> main_outer_negative_offset_polygon = std::make_shared<Polygon_with_holes_2>(Polygon_2(frame_outer_boundary->begin(), frame_outer_boundary->end()));

                              for (auto& hole1 : allowed_polygons) {
                                hole1->reverse_orientation();
                                main_outer_negative_offset_polygon->add_hole(Polygon_2(hole1->begin(), hole1->end()));
                                cs_frame.push_back(std::move(hole1));
                              }

                              ClsObject object1FrameForNegativeOffsets(object1NegativeOffsetSS.object1.object_index, std::vector<OIOA>(), std::vector<ContourSequence>(), SS_TYPE::NEGATIVE_WITH_EXTERNAL_FRAME);
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
                            ClsObject object1WithInvertedHoles(object1NegativeOffsetSS.object1.object_index, std::vector<OIOA>(), vect_cs_inverted_holes, SS_TYPE::NEGATIVE_INVERTED_HOLE);
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
                            // Если в объекте есть не только отрицательные offsets, но и положительные, то нужно задать параметры и для их расчёта тоже:
                            std::vector<OIOA> vect_oioa;
                            for (auto& oioa1 : object1NegativeOffsetSS.object1.vect_oioa) {
                              if (oioa1.offset >= 0) {
                                vect_oioa.push_back(OIOA(oioa1.object_index, oioa1.offset_index, oioa1.offset, oioa1.altitude));
                              }
                            }
                            //if (vect_oioa.size() > 0) - Это условие отменено, т.е. в будущем для рассчёта faces SS требуются все SS и положительные и отрицательные. vect_oioa.size==0, когда все offset сконцентрированы в отрицательной зоне (даже без 0)
                            {
                              ClsObject object1WithPositiveOffsets(object1NegativeOffsetSS.object1.object_index, vect_oioa, object1NegativeOffsetSS.object1.source_planes, SS_TYPE::POSITIVE);
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

                //std::vector<POMS> vect_polygon1_oioa;
                for (auto& object1 : objects) {
                  for (auto& polygon_with_holes1 : object1.vect_polygons_with_holes) {
                    // Отсортировать по offsets. Будет удобно считать, если при увеличении отступа пропадёт контур,то и считать следующие бОльшие контуры не надо
                    // - в многопоточности это уже не актуально. Может и можно в итоге прервать, но надо поискать другой метод. sort(positive_offsets_data1.vect_oap.begin(), positive_offsets_data1.vect_oap.end(), [](const OAP& oap1, const OAP& oap2) {return oap1.offset < oap2.offset; });
                    // Положительные отступы считаются по одному
                    POMS p1(object1.object_index, std::shared_ptr<Polygon_with_holes_2>(), object1.vect_oioa, object1.ss_type /*Изначально все создаваемые полигоны неопределены, т.к. отрицательные отступы рассчитываются разделением PWH-полигонов на вспомогательные объекты*/);
                    p1.polygon1_with_holes = std::move(polygon_with_holes1);
                    vect_polygon1_oioa.push_back(p1);
                  }
                }

                int threadNumbers = boost::thread::hardware_concurrency();
                boost::asio::thread_pool pool3(threadNumbers);
                boost::mutex mtx_;

                int threads_counts = 0;
                size_t ss_id_counter = 0;
                for (auto& polygon1_oioa : vect_polygon1_oioa) {
                  ss_id_counter++;
                  boost::asio::post(pool3, [ss_id_counter , &polygon1_oioa, &res_errors, &res_contours, &mtx_, threads_counts, verbose] {
                    SsBuilder ssb;
                    int verts_count = polygon1_oioa.polygon1_with_holes->outer_boundary().size();
                    {
                      // TODO: Временный параметр для равномерного offset. надо будет добавить weight во входные сокеты нода и заменить на получение weight из интерфейса нода
                      std::vector<std::vector<FT> > uniform_weights;
                      uniform_weights.reserve(polygon1_oioa.polygon1_with_holes->number_of_holes() + 1);

                      uniform_weights.push_back(std::vector<FT>(polygon1_oioa.polygon1_with_holes->outer_boundary().size(), FT(1)));
                      for (const auto& hole : polygon1_oioa.polygon1_with_holes->holes())
                        uniform_weights.push_back(std::vector<FT>(hole.size(), FT(1)));

                      // Загрузить outer boundary и holes в ssb.
                      {
                        // 1. Загрузить внешний контур:
                        ssb.enter_contour(polygon1_oioa.polygon1_with_holes->outer_boundary().begin(), polygon1_oioa.polygon1_with_holes->outer_boundary().end());
                        // 2. Загрузить holes
                        for (auto& hole1 : polygon1_oioa.polygon1_with_holes->holes()) {
                          verts_count += hole1.size();
                          ssb.enter_contour(hole1.begin(), hole1.end());
                        }
                        //ssb.enter_contour_weights(uniform_weights.begin(), uniform_weights.end());
                      }
                    }
                    // Рассчитать SS
                    CGAL::Real_timer timer1;
                    timer1.start();
                    try {
                      polygon1_oioa.ss = ssb.construct_skeleton();
                    } catch (std::exception _ex) {
                      // При ошибке исходный SS останется null
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
                      polygon1_oioa.ss_id = ss_id_counter;
                      // Запомнить ss_id в соответствующих ему oioa (ss_id не зависит от индескса объекта и является сквозным)
                      for (auto& oioa1 : polygon1_oioa.vect_oioa) {
                        oioa1.ss_id = ss_id_counter;
                      }
                    } else {
                      try {
                        // Не получилось посчитать Straight Skeleton. Надо сообщить о проблеме пользователю:
                        res_errors.push_back(ObjectError(polygon1_oioa.vect_oioa[0].object_index, *polygon1_oioa.polygon1_with_holes.get(), VAL2STR(Err::_0006)". Offset. Error Build Straight Skeleton. Outer Boundary"));
                      } catch (std::exception _ex) {
                        int error = 0;
                      }
                    }
                    mtx_.unlock();
                  }  // boost::asio::post
                  );
                  threads_counts++;
                } // for
                pool3.join();

                // Рассчитать offset-ы в многопоточном режиме для всех SS
                boost::asio::thread_pool pool4(threadNumbers);
                for (POMS& polygon1_oioa : vect_polygon1_oioa) {
                  if (polygon1_oioa.ss) {
                    for (auto& oioa1 : polygon1_oioa.vect_oioa) {
                      oioa1.ss = polygon1_oioa.ss; // Сохранить привязку offset и SS.
#ifdef _DEBUG
                      // Пропускать offset == 0 в Solution Configuration = _DEBUG, т.к. в _DEBUG срабатывает проверка на 0 в offset skeleton, хотя сам skeleton при 0 в Release строится нормально
                      if (is_zero(oioa1.offset)) {
                        continue;
                      }
#endif
                      boost::asio::post(pool4, [&polygon1_oioa, &oioa1, &mtx_, &verbose, &res_contours] {
                        ContourSequence offset_contours; // Instantiate the container of offset contours - Куда складывать контуры offset после вычисления offset. На один offset может быть несколько контуров.
                        //std::unordered_map<HDS_Halfedge_const_handle, Point_2> map_halfedge_to_offset_points; // Сопоставление точек пересечения offset с SS по halfedge.
                        //OffsetBuilder offsetBuilder(*polygon1_oioa.ss); // Instantiate the offset builder with the skeleton
                        //offsetBuilder.construct_offset_contours(abs(oioa1.offset)/*Здесь считаются и отрицательные offsets для внешних контуров, но данные подготовлены для инверсии, это учтётся позже*/, std::back_inserter(offset_contours)); // Obtain the offset contours
                        
                        CGAL::offset_params(oioa1.ss, oioa1.offset, CGAL::parameters::verbose(verbose), oioa1.map_halfedge_to_offset_points, std::back_inserter(offset_contours));

                        // Ещё вариант конфигурации контуров, которые требуется идентифицировать при сопоставлении контуров (т.е. и количество контуров меняется и результирующие фигуры тоже меняются):
                        // TODO: Нужно разобраться где разбираться с конфигурацией. Тут разобраться не получится, потому что расчёт с положительным offset (хоть он тут и берётся через abs)
                        // может являтся частью расчёта отрицательного offset, но тут их сопоставить нельзя из-за нехватки всех данных для сопоставления (внешний контур отрицательного offset
                        // считается отдельно от внутренних inverted holes). Собирать их надо в другом месте.

                        //if (offset_contours.size() > 0) - отмена проверки количество контуров для текущего offset. Загружать результат рассчёта offset, даже если при расчёте нет результирующих контуров для указанного offset.
                        // Требуется для последующих рассчётов mode Bevel и Straight Skeleton.
                        {
                          // Load results contours (may be several contours if holes)
                          // Сохранить рассчитанные по offsets горизонтальные контуры как результат (их может быть несколько, когда в объекте имеются holes):
                          // <image url="..\code_images\file_0014.png" scale=".3"/>
                          if (oioa1.offset >= 0) {
                            // Если offset положительный, то загрузить все контуры.
                            // TODO: В документации возвращаемый результат не определяет какую-то гарантированную последовательность 
                            // возвращаемых контуров. Известно, что исходный полигон всегда один, но может иметь Holes, поэтому
                            // при возврате результат может содержать некоторое количество контуров с отверстиями (опционально, если изначально
                            // в исходном контуре были holes). Если отступ указан большим, то контуров может и не быть.
                            // Для отображения результата в режиме EDGES последовательность контуров не важна, а для режима FACES это важно и пока это не реализовано.
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
                              // т.к. он является вспомогательным для рассчётов и не должен отображаться/выводится в результат.
                              // Locate the offset contour that corresponds to the frame
                              // That must be the outmost offset contour, which in turn must be the one
                              // with the largest unsigned area.
                              ContourSequence::iterator f = offset_contours.end();
                              FT lLargestArea = 0.0;
                              for (ContourSequence::iterator i = offset_contours.begin(); i != offset_contours.end(); ++i) {
                                FT lArea = CGAL_NTS abs((*i)->area()); //Take abs() as  Polygon_2::area() is signed.
                                if (lArea > lLargestArea) {
                                  f = i;
                                  lLargestArea = lArea;
                                }
                              }
                              // Remove the offset contour that corresponds to the frame max area.
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
          std::map<
            int /*индекс исходного объекта*/,
            std::vector<OIOA> /*данные для полигонов для offset_index текущего объекта*/
          > map_res_by_objectindex;
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
            // 2. Отсортировать offset по соответствующим им индексам по возрастанию в каждой группе (в результате расчёта в многопоточном режиме offset не обязательно будут посчитаны
            // в порядке их поступления, требуется расположить их в том порядке, в котором они поступили в расчёт)
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
            std::vector<std::vector<FT>> vect_vertices_of_errors;
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
                    vect_vertices_of_errors.push_back(std::vector<FT>{v.x(), v.y()});
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
                mesh_data->vertices_of_errors[I * 3 + 0] = (float)CGAL::to_double(vect_vertices_of_errors[I][0]);
                mesh_data->vertices_of_errors[I * 3 + 1] = (float)CGAL::to_double(vect_vertices_of_errors[I][1]);
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
              std::vector<unsigned int>              vect_objects_faces_indexes_counters; // Счётчики количества индексов каждой face
              std::vector<std::vector<unsigned int>> vect_objects_faces; // indexes of faces (per object1)

              // Вектор срезов. Каждый объект состоит из уровня одного offset, содержащего вектор PWH, вектором SS (из которых получен вектор PWH), и map HalfEdges,
              // соответствующих SS с соответствующим набором точек пересечения по offset.
              std::vector<OIOA_OFFSET_SS_PARAMS> vect_pwhs_offset_index_order;

              // Цикл обработки параметра vect_pwhs_offset_index_order (тут идёт распределение по object_index и offset_index):
              // После обработки в цикле в этом векторе будет содержаться следующая информация <image url="..\code_images\file_0042.png" scale=".1"/>
              for (auto& [object_index, vect_object1] : map_res_by_objectindex) {
                // Собрать список контуров на одном уровне offset и сложить их в нужной последовательности,
                // чтобы они выстроились в валидную фигуру с внешним контуром и отверстиями:
                std::map<int, std::vector<OIOA>> object1_group_contours_by_offset_index;
                for (auto& oioa1 : vect_object1) {
                  if (object1_group_contours_by_offset_index.find(oioa1.offset_index) == object1_group_contours_by_offset_index.end()) {
                      object1_group_contours_by_offset_index[oioa1.offset_index] = std::vector<OIOA>();
                  }
                  object1_group_contours_by_offset_index[oioa1.offset_index].push_back(oioa1);
                }

                for (auto& [idx, vect_oioa] /*Насколько я понимаю, map работает как упорядоченная коллекция и выдаёт результат в порядке offset_index, что является исходным заданным порядком*/
                  : object1_group_contours_by_offset_index
                  ) {
                  
                  ContourSequence cs_for_mesh;
                  std::vector<SS__HalfEdge__Point_2> vect__SS__HalfEdge__Point_2;
                  for (auto& oioa1 : vect_oioa) {
                    if (oioa1.polygon_contours.size() > 0) {
                      for (auto& cs : oioa1.polygon_contours) {
                        cs_for_mesh.push_back(cs);
                      }
                    }
                    // Запоминать SS__HalfEdge__Point_2 даже когда oioa1.polygon_contours.size()==0 (SS__HalfEdge__Point_2 нужны для рассчётов в mode Bevel и Straight Skeleton)
                    auto inst_SS_HalfEdge_Point_2 = SS__HalfEdge__Point_2(oioa1.ss, oioa1.ss_id, oioa1.neg_type, oioa1.map_halfedge_to_offset_points);
                    vect__SS__HalfEdge__Point_2.push_back(inst_SS_HalfEdge_Point_2);
                  }
                  // Упорядочить вектор по ss_id в сторону увеличения. Для удобства, вроде как не должно ни на что влиять, но при отладке удобнее смотреть на одинаковый порядок.
                  // Состав vect__SS__HalfEdge__Point_2: <image url="..\code_images\file_0039.png" scale=".025"/>
                  sort(vect__SS__HalfEdge__Point_2.begin(), vect__SS__HalfEdge__Point_2.end(), [](const auto& hhsp1, const auto& hhsp2) {return hhsp1.ss_id < hhsp2.ss_id; });

                  {
                    std::vector<Polygon_with_holes_2> vect_pwh; // набор pwh соответствующих одному offset для текущего объекта. Может не содержать контуров (требуется для рассчётов в mode Bevel и Straight Skeleton)
                    if (cs_for_mesh.size() > 0) {
                      ConvertContourSequenceIntoPolygonsWithHoles(cs_for_mesh, vect_pwh);
                    }
                    vect_pwhs_offset_index_order.push_back(OIOA_OFFSET_SS_PARAMS(
                      // <image url="..\code_images\file_0040.png" scale=".2"/>
                      vect_oioa[0] /* Все элементы в vect_oioa имеют одинаковые object_index и offset_index oioa, поэтому и берём 0-й oioa объект */,
                      // <image url="..\code_images\file_0041.png" scale=".1"/> - список PWH из которого состоит "срез" текущего offset по всему объекту (пример при разных offset, но при вызове этот параметр только при одном offset)
                      vect_pwh,
                      vect__SS__HalfEdge__Point_2 /*Запомнить, какие SS связаны с текущим уровнем offset*/
                    ));
                  }
                }
              }

              // Определить способ обработки объектов по join_mode и подготовить набор данных для результатов (новый список результирующих объектов)
              std::map<
                int /* индекс результирующего объекта */,
                std::vector<OIOA_OFFSET_SS_PARAMS> /* данные для полигонов для offset_index результирующего объекта */
              > map_join_mesh;
              {
                if (result_type == ResType::EDGES || result_type == ResType::FACES) {
                  // Как работает join_mode <image url="..\code_images\file_0024.png" scale=".2"/>. Для режима EDGES всё тоже самое, только выводятся только контуры соответствующего face без триангуляции (с учётом holes)
                  // и с учётом направленности - внешние контуры закручены против часовой, а внутренние holes закручены по часовой: <image url="..\code_images\file_0025.png" scale=".2"/>

                  int split_object_index = 0; // для join_mode==SPLIT индекс нового объекта назначается группе контуров с одним index_offset (а не каждому контуру, чтобы количество объектов в режиме Edges и Faces совпадало)
                  int merge_object_index = 0; // Для join_mode==MERGE индекс результирующего объекта один - 0, т.к. объект будет единственным.
                  int  keep_object_index = 0; // Для join_mode==KEEP  индексы результирующих объектов не меняются

                  // Сгруппировать или разгруппировать объекты по параметру results_join_mode
                  for (auto& pwhs_offset_index_order1 : vect_pwhs_offset_index_order) {
                    if (results_join_mode == 0) {
                      // Сгруппировать результирующие объекты по отдельным offset_index
                      // Разделить текущий набор pwhs на отдельные pwh1 и сделать каждый из них отдельным вектором (для удобства дельнейшей обработки)
                      for (auto& pwh1 : pwhs_offset_index_order1.vect_pwh) {
                        std::vector<Polygon_with_holes_2> pwh1_as_vector; // Один элемент в целом списке векторов (для результата количество элементов всё равно должно быть запомнено в векторах)
                        pwh1_as_vector.push_back(pwh1);
                        std::vector<OIOA_OFFSET_SS_PARAMS> vect_PWHS_ALTITUDE_as_vector;
                        std::vector<SS__HalfEdge__Point_2> vect__SS__HalfEdge__Point_2;
                        vect__SS__HalfEdge__Point_2.push_back(SS__HalfEdge__Point_2(pwhs_offset_index_order1.oioa1.ss, pwhs_offset_index_order1.oioa1.ss_id, pwhs_offset_index_order1.oioa1.neg_type, pwhs_offset_index_order1.oioa1.map_halfedge_to_offset_points));
                        vect_PWHS_ALTITUDE_as_vector.push_back(OIOA_OFFSET_SS_PARAMS(pwhs_offset_index_order1.oioa1, pwh1_as_vector, vect__SS__HalfEdge__Point_2 /*Из какого SS получен текущий контур*/));
                        map_join_mesh[split_object_index] = vect_PWHS_ALTITUDE_as_vector;
                        split_object_index++;
                      }
                    } else if (results_join_mode == 1) {
                      keep_object_index = pwhs_offset_index_order1.oioa1.object_index;
                      // Оставить результирующие объекты как в исходных объектах
                      if (map_join_mesh.find(keep_object_index) == map_join_mesh.end()) {
                        map_join_mesh[keep_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                      }
                      map_join_mesh[keep_object_index].push_back(pwhs_offset_index_order1);
                    } else { // 2 - для всех остальных случаев
                      // Объеденить все результаты в один объект
                      if (map_join_mesh.find(merge_object_index) == map_join_mesh.end()) {
                        map_join_mesh[merge_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                      }
                      map_join_mesh[merge_object_index].push_back(pwhs_offset_index_order1);
                    }
                  }
                } else if (result_type == ResType::BEVEL) {

                  // ПОКА НЕ ИСПОЛЬЗУЕТСЯ!!!

                  //// Тут посчитать Extrude и определить и количественный состав данных и структурный состав данных.
                  //// Расчёт делается на основе данных в vect_pwhs_offset_index_order
                  //vect_pwhs_offset_index_order;

                  //for (size_t I = 0; I <= (int)vect_pwhs_offset_index_order.size() - 2 /*bevel может быть только между двух offset*/; I++) {
                  //  OIOA_OFFSET_SS_PARAMS& oioa_offset_ss_params_I0 = vect_pwhs_offset_index_order[I + 0];
                  //  OIOA_OFFSET_SS_PARAMS& oioa_offset_ss_params_I1 = vect_pwhs_offset_index_order[I + 1];
                  //  
                  //  // Возможные варианты соотношения offset-ов разных уровней:
                  //  // 0. Оба значения одинаковые. Это ошибка. В принципе могут попадаться одинаковые значения, но не вместе.
                  //  // 1. оба offset с одним знаком или один offset знаковый, а второй 0. Все параметры можно брать из одного SS
                  //  // 2. один offset положительный, другой offset отрицательный.
                  //  FT val_offset0 = oioa_offset_ss_params_I0.oioa1.offset;
                  //  FT val_offset1 = oioa_offset_ss_params_I1.oioa1.offset;
                  //  auto sign_offset0 = sign(val_offset0);
                  //  auto sign_offset1 = sign(val_offset1);


                  //  // Если offset-ы разных знаков, то найти какими faces (из существующих ss) сопоставляются два ребра между собой?
                  //  if (sign_offset0!=sign_offset1 && sign_offset0!= CGAL::Sign::ZERO && sign_offset1!= CGAL::Sign::ZERO) {
                  //    
                  //    // Идти по рёбрам 0-го offset.
                  //    for (auto SS_HalfEdge_Point_2_iter_0 = oioa_offset_ss_params_I0.vect__SS__HalfEdge__Point_2.begin(); SS_HalfEdge_Point_2_iter_0 != oioa_offset_ss_params_I0.vect__SS__HalfEdge__Point_2.end(); ++SS_HalfEdge_Point_2_iter_0) {
                  //      
                  //      auto& SS_HalfEdge_Point_2_inst_0 = *oioa_offset_ss_params_I0.vect__SS__HalfEdge__Point_2.begin();
                  //      const HDS& hds_0 = static_cast<const HDS&>(*SS_HalfEdge_Point_2_iter_0->ss);
                  //      SS_HalfEdge_Point_2_iter_0->ss->size_of_border_halfedges();
                  //      auto size_of_border_edges_0 = SS_HalfEdge_Point_2_iter_0->ss->size_of_border_edges();
                  //      auto size_of_border_halfedges_0 = SS_HalfEdge_Point_2_iter_0->ss->size_of_border_halfedges();
                  //      auto size_of_halfedges_0 = SS_HalfEdge_Point_2_iter_0->ss->size_of_halfedges();
                  //      
                  //      // Цикл по faces полученного SS (SS состоит из faces)
                  //      bool ignore_frame_faces_0 = false;
                  //      size_t fc_0 = 0;
                  //      size_t fc_0_filtered = 0;
                  //      for (const HDS_Face_handle hds_f_0 : CGAL::faces(hds_0)) {
                  //        ignore_frame_faces_0 = (val_offset0 < 0) && (SS_HalfEdge_Point_2_iter_0->neg_type== OIOA::NEGATIVE_TYPE::WITH_EXTERNAL_FRAME);
                  //        if (ignore_frame_faces_0 && fc_0++ < 4) {
                  //          continue;
                  //        }
                  //        fc_0_filtered++;

                  //        HDS_Halfedge_const_handle hds_h_0 = hds_f_0->halfedge(), done = hds_h_0;
                  //        int id_hds_half_edge_0 = hds_h_0->id();
                  //        HDS_Vertex_const_handle hds_vert_current_0  = hds_h_0->vertex();
                  //        
                  //        auto p0 = hds_vert_current_0->point();
                  //        auto p0_x = CGAL::to_double(p0.x());
                  //        auto p0_y = CGAL::to_double(p0.y());
                  //        HDS_Vertex_const_handle hds_vert_opposite_0 = hds_h_0->opposite()->vertex();

                  //        // Найти SS, у которого будет face, смежная с текущим face от 0-го итератора рёбер offset
                  //        size_t fc_1 = 0;
                  //        size_t fc_1_filtered = 0;
                  //        bool adjacent_founded = false;
                  //        size_t ss_id_0 = -1;
                  //        for (auto SS_HalfEdge_Point_2_iter_1 = oioa_offset_ss_params_I1.vect__SS__HalfEdge__Point_2.begin(); SS_HalfEdge_Point_2_iter_1 != oioa_offset_ss_params_I1.vect__SS__HalfEdge__Point_2.end(); ++SS_HalfEdge_Point_2_iter_1) {
                  //          auto& SS_HalfEdge_Point_2_inst_1   = *oioa_offset_ss_params_I1.vect__SS__HalfEdge__Point_2.begin();
                  //          const HDS& hds_1 = static_cast<const HDS&>(*SS_HalfEdge_Point_2_iter_1->ss);

                  //          bool ignore_frame_faces_1 = false;
                  //          size_t fc_1 = 0;
                  //          size_t fc_1_filtered = 0;
                  //          ss_id_0 = SS_HalfEdge_Point_2_iter_0->ss_id;
                  //          for (const HDS_Face_handle hds_f_1 : CGAL::faces(hds_1)) {
                  //            ignore_frame_faces_1 = (val_offset1 < 0) && (SS_HalfEdge_Point_2_inst_1.neg_type == OIOA::NEGATIVE_TYPE::WITH_EXTERNAL_FRAME);
                  //            if (ignore_frame_faces_1 && fc_1++ < 4) {
                  //              continue;
                  //            }
                  //            fc_1_filtered++;

                  //            HDS_Halfedge_const_handle hds_h_1           = hds_f_1->halfedge(), done = hds_h_1;
                  //            int id_hds_half_edge_1 = hds_h_1->id();
                  //            HDS_Vertex_const_handle hds_vert_current_1  = hds_h_1->vertex();
                  //            auto p1 = hds_vert_current_1->point();
                  //            auto p1_x = CGAL::to_double(p1.x());
                  //            auto p1_y = CGAL::to_double(p1.y());
                  //            HDS_Vertex_const_handle hds_vert_opposite_1 = hds_h_1->opposite()->vertex();

                  //            if (p0_x == p1_x && p0_y == p1_y) {
                  //              if (verbose == true) {
                  //                size_t ss_id_1 = SS_HalfEdge_Point_2_iter_1->ss_id;
                  //                printf("\n " VAL2STR(Err::_0029) ". SS 2D Offset.     founded ss_id_0: %zu, edge_0 idx: %zu, edge_0 id: %u, ss_id_1: %zu, edge_1 idx: %zu, edge_1 id: %u ", ss_id_0, fc_0_filtered, id_hds_half_edge_0, ss_id_1, fc_1_filtered, id_hds_half_edge_1 );
                  //              }
                  //              adjacent_founded = true;
                  //              break;
                  //            }
                  //          }
                  //          if (adjacent_founded == true) {
                  //            break;
                  //          }
                  //        }
                  //        if (adjacent_founded == false) {
                  //          printf("\n " VAL2STR(Err::_0029) ". SS 2D Offset. not founded ss_id_0: %zu, edge_0 idx: %zu, edge_0 id: %u, offset: %.5g ", ss_id_0, fc_0_filtered, id_hds_half_edge_0, CGAL::to_double(val_offset0) );
                  //        }
                  //      }
                  //    }
                  //  }
                  //}

                } else if (result_type == ResType::STRAIGHT_SKELETON) {
                  // Рассчитать Straight Skeleton с учётом offsets в горизонтальной плоскости (offsets сортируются только по индексам в порядке, заданном пользователем, не по значениям)
                  // <image url="..\code_images\file_0043.gif" scale=".3"/>

                  CGAL::Real_timer timer1; // Тратится около 3% времени на merge и расчёт ajacent faces.
                  timer1.start();
                  // Вспомогательные переменные:
                  bool ignore_frame_faces = false;
                  bool invert_faces = false;

                  std::map<int /*object_index*/, std::map<int /*ss_id*/, std::vector<OIOA_OFFSET_SS_PARAMS>>> map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS;
                  for (auto& oioa_offset_ss_params : vect_pwhs_offset_index_order) {
                    if (map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS.find(oioa_offset_ss_params.oioa1.object_index) == map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS.end()) {
                      map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[oioa_offset_ss_params.oioa1.object_index] = std::map<int, std::vector<OIOA_OFFSET_SS_PARAMS>>();
                    }
                    for (auto& ss_halfedge_point_2 : oioa_offset_ss_params.vect__SS__HalfEdge__Point_2) {
                      if (map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[oioa_offset_ss_params.oioa1.object_index].find(ss_halfedge_point_2.ss_id) == map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[oioa_offset_ss_params.oioa1.object_index].end()) {
                        map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[oioa_offset_ss_params.oioa1.object_index][ss_halfedge_point_2.ss_id] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                      }
                      map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[oioa_offset_ss_params.oioa1.object_index][ss_halfedge_point_2.ss_id].push_back(
                        OIOA_OFFSET_SS_PARAMS(oioa_offset_ss_params.oioa1, oioa_offset_ss_params.vect_pwh, std::vector<SS__HalfEdge__Point_2>(1, ss_halfedge_point_2))
                      );
                    }
                  }

                  // 
                  std::map<int /*object_index*/, std::map<int /*ss_id*/, std::map<int /*ss_face_id*/, HDS_Face_handle>>> map__object_index__ss_id__ss_face_id__face_handle;

                  // Собрать всех ss от одного объекта в одном векторе:
                  std::map<int /*object_index*/, std::vector<POMS>> map__object_index__poms;
                  {
                    for (auto& poms : vect_polygon1_oioa) {
                      if (map__object_index__poms.find(poms.object_index) == map__object_index__poms.end()) {
                        map__object_index__poms[poms.object_index] = std::vector<POMS>();
                      }
                      map__object_index__poms[poms.object_index].push_back(poms);
                    }
                    for (auto& [object_index, vect_poms] : map__object_index__poms) {
                      if (map__object_index__ss_id__ss_face_id__face_handle.find(object_index) == map__object_index__ss_id__ss_face_id__face_handle.end()) {
                        map__object_index__ss_id__ss_face_id__face_handle[object_index] = std::map<int /*ss_id*/, std::map<int /*ss_face_id*/, HDS_Face_handle>>();
                      }

                      // Сквозная нумерация всех faces во всех faces в Straight Skeleton-ах в пределах одного объекта object_index.
                      size_t face_continuous_numbering = 0;

                      for (auto& poms : vect_poms) {
                        // Получить все faces каждого SS, пронумеровать их в сквозном порядке и в каждой face сохранить информацию о привязке к ss_id, исходной face_id от ss и присвоить новый сквозной номер face.
                        // Примечание: Во время получения faces пропускать главный отрицательный SS (если он есть) и не учитывать его внешний вспомогательный контур.
                        {
                          if (map__object_index__ss_id__ss_face_id__face_handle[object_index].find(poms.ss_id) == map__object_index__ss_id__ss_face_id__face_handle[poms.object_index].end()) {
                            map__object_index__ss_id__ss_face_id__face_handle[object_index][poms.ss_id] = std::map<int /*ss_face_id*/, HDS_Face_handle>();
                          }
                          {
                            struct MESH_SS_STRUCT {
                              /// <summary>
                              /// Точки Mesh SS с их исходными point id.
                              /// </summary>
                              std::map<int /*point_id*/, Point_2> map_points;

                              /// <summary>
                              /// faces Mesh SS с их исходными face id SS.
                              /// </summary>
                              std::map<int /*face_id*/, std::vector<int/*point_id*/>> map_vector_face_points_id;

                              /// <summary>
                              /// Набор игнорируемых face_id (для внешнего frame).
                              /// </summary>
                              std::set<int> set_ignore_faces_id; // Вектор игнорирования faces

                              /// <summary>
                              /// Набор игнорируемых point_id (для внешнего frame).
                              /// </summary>
                              std::set<int> set_ignore_vertices_id;
                            } mesh_data;

                            const HDS& hds = static_cast<const HDS&>(*(poms.ss));
                            std::size_t fc = 0;

                            // Цикл по faces полученного SS (SS состоит из faces)
                            for (const HDS_Face_handle ss_face_handle : CGAL::faces(hds)) {
                              int ss_face_id = ss_face_handle->id(); // текущий face_id
                              // Если offset отрицательный и ss получен при помощи вспомогательного внешнего frame,
                              // то удалить external frame только у ss, полученного с помощью вспомогательного frame и
                              // ничего не удалять у ss, являющихся inverted_holes (инвертированными отверстиями)
                              // <image url="..\code_images\file_0045.png" scale="1.0"/>
                              ignore_frame_faces = poms.ss_type==SS_TYPE::NEGATIVE_WITH_EXTERNAL_FRAME;
                              bool ignore_frame_vertices = false; // Игнорировать вершины из внешнего frame. Если frame будет не внешний, то удалить такую вершины их вектора игнорирования
                              if (ignore_frame_faces && fc++ < 4) {
                                //continue;
                                ignore_frame_vertices = true;
                                mesh_data.set_ignore_faces_id.insert(ss_face_id); // Признак, что данный face в будущем не нужно использовать
                              }
                              if (map__object_index__ss_id__ss_face_id__face_handle[object_index][poms.ss_id].find(ss_face_id) == map__object_index__ss_id__ss_face_id__face_handle[object_index][poms.ss_id].end()) {
                                map__object_index__ss_id__ss_face_id__face_handle[object_index][poms.ss_id][ss_face_id] = ss_face_handle;
                              }

                              // Точки одного полигона рассматриваемого SS <image url="..\code_images\file_0050.png" scale=".2"/>
                              std::vector<Point_3> vect_ss_face_points;
                              std::vector<int> vertex_face_indexes;

                              // Если я правильно понял, то hds_h - полуребро Straight Skeleton, которое рассматривается в настоящий момент.
                              HDS_Halfedge_const_handle hds_h = ss_face_handle->halfedge(), done = hds_h;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                              HDS_Halfedge_const_handle contour_h = hds_h->defining_contour_edge();
                              CGAL_assertion(hds_h == contour_h);
                              //const bool is_vertical = (contour_h->weight() == vertical_weight);
#endif
                              int verts_counter = 0;

                              do {
                                // <image url="..\code_images\file_0035.png" scale="1.0"/>
                                HDS_Vertex_const_handle hds_vert_opposite = hds_h->opposite()->vertex();
                                HDS_Vertex_const_handle hds_vert_current = hds_h->vertex();
                                verts_counter++;

                                const Point_2& off_p = hds_vert_current->point();
                                vect_ss_face_points.emplace_back(off_p.x(), off_p.y(), 0.0);

                                int vertex_id = hds_vert_current->id();
                                vertex_face_indexes.push_back(vertex_id);
                                mesh_data.map_points[vertex_id] = off_p;

                                // Немного хитрое условие. Если попался Vertex для внешнего frame, то пометить его на удаление.
                                // Если этот vertex попался позже в face, который не является внешним, то отложить его удаление
                                if (ignore_frame_vertices == true) {
                                  mesh_data.set_ignore_vertices_id.insert(vertex_id);
                                } else {
                                  mesh_data.set_ignore_vertices_id.erase(vertex_id);
                                }
                                hds_h = hds_h->next();
                              } while (hds_h != done);

                              if (vect_ss_face_points.size() < 3) {
                                std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
                                // При этом не запоминаются mesh_data.map_vector_face_points_id[face_id]. Возможно останутся висящие точки.
                                // TODO: потом надо будет проверить mesh на одиночные точки
                                continue;
                              }

                              mesh_data.map_vector_face_points_id[ss_face_id] = vertex_face_indexes; // Сюда уже не пропускается face, если точек меньше 3, хотя сами точки запоминаются.

                              // Сформировать из новых индексов новый face.
                              //auto& face_index = mesh_ss.add_face(vect_vertex_indexes);
                              // Добавить в новый face атрибут с параметрами ss_id и номер сквозной нумерации <image url="..\code_images\file_0051.png" scale=".2"/>
                              //mesh_ss_property_map[face_index] = MeshFacePropertyStruct(ss_halfedge_point_2.ss_id, face_continuous_numbering++);

                              //map_object_ss_geometry_points[object_index].push_back(face_points);
                              //triangulate_skeleton_face(face_points, invert_faces, points, faces, CGAL::parameters::verbose(verbose));

                            } // Цикл по faces

                            int s1 = mesh_data.map_points.size();
                            int s2 = mesh_data.map_vector_face_points_id.size();
                            int s3 = mesh_data.set_ignore_vertices_id.size();
                            int x = 0;

                            // Удалить неиспользуемые faces:
                            for (auto& face_id : mesh_data.set_ignore_faces_id) {
                              mesh_data.map_vector_face_points_id.erase(face_id);
                              map__object_index__ss_id__ss_face_id__face_handle[object_index][poms.ss_id].erase(face_id);
                            }
                            // Удалить неиспользуемые vertices
                            for (auto& point_id : mesh_data.set_ignore_vertices_id) {
                              mesh_data.map_points.erase(point_id);
                            }

                            // Добавить точки для новой mesh и запомнить их новые индексы:
                            Mesh& mesh_ss = poms.mesh_ss_01_source; // vect__oioa__offset_ss_params[0].mesh_ss_01_source;
                            // Добавить в новый face атрибут с параметрами ss_id и номер сквозной нумерации <image url="..\code_images\file_0051.png" scale=".2"/>
                            auto [mesh_ss_property_map, created] = mesh_ss.add_property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct", MeshFacePropertyStruct(-1, MeshFacePropertyStruct::SS_SIGN::UNDEFINED, -1, -1));

                            std::map<int, Mesh::vertex_index> map_renamed_point_id;
                            int new_point_counter = 0;
                            for (auto& [point_id, point2] : mesh_data.map_points) {
                              Point_3 p3(point2.x(), point2.y(), 0.0);
                              auto& new_point_id = mesh_ss.add_vertex(p3);
                              map_renamed_point_id[point_id] = new_point_id;
                            }
                            // Создать новые faces на основе новых индексов:
                            for (auto& [ss_face_id /*исходный face_id от SS*/, vect_points_id] : mesh_data.map_vector_face_points_id) {
                              std::vector<Mesh::vertex_index> face_new_point_id;
                              for (auto& prev_point_id : vect_points_id) {
                                face_new_point_id.push_back(map_renamed_point_id[prev_point_id]);
                              }
                              Mesh::Face_index face_index = mesh_ss.add_face(face_new_point_id);
                              mesh_ss_property_map[face_index] = MeshFacePropertyStruct(poms.ss_id, poms.ss_type==SS_TYPE::POSITIVE ? MeshFacePropertyStruct::SS_SIGN::POSITIVE : poms.ss_type==SS_TYPE::NEGATIVE_WITH_EXTERNAL_FRAME ||poms.ss_type==SS_TYPE::NEGATIVE_INVERTED_HOLE ? MeshFacePropertyStruct::SS_SIGN::NEGATIVE : MeshFacePropertyStruct::SS_SIGN::UNDEFINED, ss_face_id, face_continuous_numbering);
                              face_continuous_numbering++;
                            }
                            //vect__oioa__offset_ss_params[0].vect_vect_ss_face_points.push_back(vect_ss_face_points);
                            // Тут получена Mesh SS, у которой не надо делать merge points. Вроде как теперь дальнейший merge нескольких mesh_ss должно пройти быстрее?


                            // 
                            // Примечание: ввиду того, что offset пока не учитывается, то геометрия ss хранится только в первом элементе map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_id][ss_id][0]
                            // Примечание 2: подумать, как вообще вынести отсюда сохранение mesh_ss в map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS.
                          }
                        }
                      }
                    }
                  }

                  
                  // Связать Object_index и Mesh SS, чтобы в дальнейшем определить смежные faces SS для каждого объекта object_index целиком. Пример результата таблицы смежных faces <image url="..\code_images\file_0049.png" scale=".1"/>
                  
                  std::map<int, Mesh> map_object_index_mesh_ss_joined;
                  for (auto& [object_index, vect_poms] : map__object_index__poms){
                    // Просуммировать все SS Mesh, которые получены в текущем объекте с сохранением свойств faces (ss_id и face_index).
                    // Для этого нужно выполнить загрузку всех faces в один mesh и смержить точки.
                    // Создаем новую поверхность для объединения
                    Mesh result_ss_Mesh;
                    // Чтобы дополнительные свойства добавлялись в result_mesh_ss нужно и туда тоже добавить это свойство (без этого это свойство в исходных mesh копироваться не будет)
                    result_ss_Mesh.add_property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct", MeshFacePropertyStruct(-1, MeshFacePropertyStruct::SS_SIGN::UNDEFINED, -1, -1));
                    for (auto& poms : vect_poms) {
                      result_ss_Mesh.join(poms.mesh_ss_01_source); // - в исходниках вроде как делает merge vertices и faces -  update - нет, не делает
                    }
                    map_object_index_mesh_ss_joined[object_index] = result_ss_Mesh;
                  }

                  std::map<int, Mesh> map_object_index_mesh_ss_merged;
                  for (auto& [object_index, mesh_ss_joined] : map_object_index_mesh_ss_joined) {
                    // Смержить вершины результирующего mesh (faces оставить как есть).
                    // В результате одинаковые vertices объединяться, а faces останутся и при этом у них останется прежнее количество вершин. Однако в результате
                    // нормализации faces может измениться последовательность halfedges, но это не имеет никакого значения, т.к. соседство faces друг с другом от
                    // этого не поменяется.
                    // https://github.com/CGAL/cgal/issues/5039
                    std::vector<Point_3>                  result_ss_mesh_source_points; // Для сравнения с результатом merge
                    std::vector<std::vector<std::size_t>> result_ss_mesh_source_polygons; // Для сравнения с результатом merge
                    // Параметры для получения результатов merge:
                    std::vector<Point_3>                  result_ss_mesh_merged_points; 
                    std::vector<std::vector<std::size_t>> result_ss_mesh_merged_polygons;
                    // Пример сравнения result_ss_mesh_source_points и result_ss_mesh_merged_polygons <image url="..\code_images\file_0048.png" scale=".1"/>
                    
                    Mesh mesh_ss_merged; // Сюда положить результат merge
                    auto [mesh_ss_merged_property_map, created] = mesh_ss_merged.add_property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct", MeshFacePropertyStruct(-1, MeshFacePropertyStruct::SS_SIGN::UNDEFINED, -1, -1)); // Добавить структуру свойств. По умолчанию в новый объект свойства из join-объектов не копируются

                    {
                      // merge вершин. Полезная информация о Polygon Mesh: https://doc.cgal.org/latest/Polygon_mesh_processing/index.html, о Surface Mesh: https://doc.cgal.org/latest/Surface_mesh/index.html
                      CGAL::Real_timer timer_merge;
                      timer_merge.start();
                      //CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(result_ss_mesh_source, result_ss_mesh_source_points, result_ss_mesh_source_polygons);
                      CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(mesh_ss_joined, result_ss_mesh_merged_points, result_ss_mesh_merged_polygons);
                      CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons);
                      // Примечание: После merge иногда остаются точки, которые "кажутся" не слившимися, но это скорее всего погрешности исходных данных. Иногда
                      // задаётся многоугольник, у которого точки должны сойтись в центре, но из-за погрешности размещания углов появляются дополнительные точки
                      // которые отличаются друг от друга в 4-м знаке, 7-м, 9-м знаке и т.д. Вроде как не надо с этим бороться, хотя это и выглядит странно:
                      // <image url="..\code_images\file_0052.png" scale=".2"/> (на картинке есть ещё не помеченные близкие точки, где индекс ребра находится прямо на точке)

                      if (!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(result_ss_mesh_merged_polygons)) {
                        CGAL::Polygon_mesh_processing::orient_polygon_soup(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons);
                      }
                      //CGAL::Polygon_mesh_processing::merge_duplicate_polygons_in_polygon_soup(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons); // - merge polygons НЕ ДЕЛАТЬ!!!
                      CGAL::Polygon_mesh_processing::repair_polygon_soup(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons, CGAL::parameters::erase_all_duplicates(true).require_same_orientation(true));
                      //CGAL::Polygon_mesh_processing::orient_polygon_soup(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons);
                      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons, mesh_ss_merged);
                      CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycles(mesh_ss_merged);
                      //CGAL::Polygon_mesh_processing::stitch_borders(mesh_ss_merged); // - Сшивание. Вроде как это нужно только для объёмных фигур. https://doc.cgal.org/latest/Polygon_mesh_processing/index.html
                      timer_merge.stop();
                      if (verbose) {
                        printf("\n " VAL2STR(Err::_0032) ". SS 2D Offset. object_index %u, ss merge ajacent timer %.5g", object_index, timer_merge.time());
                      }
                    }

                    map_object_index_mesh_ss_merged[object_index] = mesh_ss_merged;

                    // Загрузка атрибутов faces из старого результирующего Mesh:
                    auto& _map_mesh_ss_property_map = mesh_ss_joined.property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct");
                    if (_map_mesh_ss_property_map.has_value() == true) {
                      auto& mesh_ss_source_property_map = _map_mesh_ss_property_map.value();
                      //int counter = 0;
                      //// Вывод property_map старого результирующего mesh:
                      //for (auto& face_index : result_ss_mesh.faces()) {
                      //  auto& val = mesh_ss_property_map[face_index];
                      //  if (verbose == true) {
                      //    printf("\n " VAL2STR(Err::_0030) ". SS 2D Offset. SS merged: %u. ss_id: %u, face index: %d", counter++, val.ss_id, val.face_index);
                      //  }
                      //}

                      // для тестов:
                      auto number_of_edges_source     = mesh_ss_joined.number_of_edges();
                      auto number_of_halfedges_source = mesh_ss_joined.number_of_halfedges();
                      auto number_of_faces_source     = mesh_ss_joined.number_of_faces();

                      auto number_of_edges_merged     = mesh_ss_merged.number_of_edges(); // Так странно - количество edges можно узнать, а пройтись по списке edges и получить их вершины на их рёбрах - не получилось пока найти.
                      auto number_of_halfedges_merged = mesh_ss_merged.number_of_halfedges();
                      auto number_of_faces_merged     = mesh_ss_merged.number_of_faces();

                      if (mesh_ss_merged.is_valid(true) == false) {
                        // TODO: Не знаю пока надо или нет это проверять (проверка mesh на валидность) - непонятно как реагировать.
                      }

                      // Расчёт смежных faces из разных Straight Skeleton текущего объекта object_index.
                      // Есть исключение при получении результата. Если все offset сконцентрированы только в положительных SS или только в отрицательных SS,
                      // то map_adjacent_mesh_faces остаётся пустым. Например:
                      // <image url="..\code_images\file_0064.png" scale="1.0"/>
                      std::map<int /*mesh_face_id*/,          int /*mesh_face_id*/                      >  map__adjacent_mesh_faces;
                      std::map<int /*mesh_face_id*/,                              MeshFacePropertyStruct>  map__mesh_face_id__face_info;
                      std::map<int /*       ss_id*/, std::map<int /*ss_face_id*/, MeshFacePropertyStruct>> map__ss_id__ss_face_id__info;

                      // Скопировать информацию о faces в смерженный новый mesh (ss_id и сквозной индекс face_index)
                      for (auto& f1 : mesh_ss_merged.faces()) {
                        // Скопировать faces property из старого результирующего mesh в новый результируюший mesh (проверил, что faces идут по очереди как и в исходном несмерженном результирующем mesh).
                        mesh_ss_merged_property_map[f1] = mesh_ss_source_property_map[f1];
                      }

                      // Тут непосредственно рассчитываются смежные faces на новом Mesh merged
                      for (auto& hd : mesh_ss_merged.halfedges()) {
                        // Пропускаем внешние границы Mesh
                        if (mesh_ss_merged.is_border(hd) == true) {
                          continue;
                        }
                        //// Иногда проверка halfedge на border не срабатывает, а при этом сама edge является border. - обновление - надо и opposite halfedge тоже надо проверять на is_border.
                        //auto& edge_index = result_ss_mesh_merged.edge(hd);
                        //if (result_ss_mesh_merged.is_border(edge_index) == true) {
                        //  continue;
                        //}

                        // Получаем первое полуребро
                        auto& h1 = hd;
                        // Работа verbose==true при ошибке (вывод в к консоль): <image url="..\code_images\file_0047.png" scale=".1"/>
                        if (mesh_ss_merged.is_valid(h1, verbose) == false || mesh_ss_merged.is_removed(h1) == true) {
                          continue;
                        }

                        // Получаем второе полуребро (противоположное)
                        auto& h2 = mesh_ss_merged.opposite(h1);
                        // И его тоже надо проверить, что оно не крайнее:
                        if (mesh_ss_merged.is_border(h2) == true) {
                          continue;
                        }
                        if (mesh_ss_merged.is_valid(h2, verbose) == false || mesh_ss_merged.is_removed(h2)==true) {
                          continue;
                        }

                        // Получаем первую грань
                        auto f1 = mesh_ss_merged.face(h1);

                        // Получаем вторую грань
                        auto f2 = mesh_ss_merged.face(h2);

                        // Проверить, что они валидные (по хорошему надо что-то вывести?)
                        if ( mesh_ss_merged.is_valid(f1, verbose) == false || mesh_ss_merged.is_removed(f1) == true) {
                          continue;
                        }
                        if (mesh_ss_merged.is_valid(f2, verbose) == false || mesh_ss_merged.is_removed(f2) == true ) {
                          continue;
                        }

                        auto& face1_info = mesh_ss_merged_property_map[f1];
                        auto& face2_info = mesh_ss_merged_property_map[f2];
                        // На этой схеме видно, что смежные faces только если ss_id у них разный.
                        // При одинаковом ss_id - это просто соседние faces.
                        // <image url="..\code_images\file_0051.png" scale=".2"/>
                        if (face1_info.ss_id != face2_info.ss_id) {
                          // Тонкий момент - нужно больше осмысления - точно ли не бывает такого, что пары всегда уникальны?
                          // Надо выяснить, может быть бывают случаи, когда один face имеет больше чем одну смежную пару?
                          // Например, один face длинный и граничит с несколькими faces от другого SS?
                          // Рассуждения: смежными бывают только faces из разнознаковых SS, а отрицательные SS
                          // получаются только из положительных SS. Положительные SS никогда не соприкасаются друг
                          // с другом. Следовательно смежные пары - это всегда производные от исходных контуров, типа inset,
                          // только наружу. Поэтому у одного face всегда только один смежный face.
                          // <image url="..\code_images\file_0046.png" scale=".2"/>
                          if (map__adjacent_mesh_faces.find(face1_info.mesh_face_id) == map__adjacent_mesh_faces.end()) {
                            if (map__adjacent_mesh_faces.find(face2_info.mesh_face_id) == map__adjacent_mesh_faces.end()) {
                              // Регистрация смежной пары faces между двумя разными SS в пределах одного объекта:
                              // Примечание: смежная пара всегда единственная. Если пара зарегистрирована, то ни одна из faces
                              // этой пары больше не будет смежной ни с какой другой faces. Напоминание - смежные пары существуют
                              // только у SS с разными знаками, а SS с одним знаком никогда на соединяются как смежные.
                              map__adjacent_mesh_faces[face1_info.mesh_face_id] = face2_info.mesh_face_id;
                              map__mesh_face_id__face_info[face1_info.mesh_face_id] = face1_info;
                              map__mesh_face_id__face_info[face2_info.mesh_face_id] = face2_info;
                            }
                          }
                        }
                      }

                      for (auto& fd : mesh_ss_merged.faces()) {
                        auto& face1_info = mesh_ss_merged_property_map[fd];
                        if (map__ss_id__ss_face_id__info.find(face1_info.ss_id) == map__ss_id__ss_face_id__info.end()) {
                          map__ss_id__ss_face_id__info[face1_info.ss_id] = std::map<int, MeshFacePropertyStruct>();
                        }
                        map__ss_id__ss_face_id__info[face1_info.ss_id][face1_info.ss_face_id] = face1_info;
                      }

                      if (verbose == true) {
                        // Вывести полученные пары: <image url="..\code_images\file_0049.png" scale=".2"/>
                        printf("\n faces pairs of object_index: %u :=========================", object_index);
                        for (auto& [f1_index, f2_index] : map__adjacent_mesh_faces) {
                          printf("\n " VAL2STR(Err::_0031) ". SS 2D Offset. object_index: % 2u, SS faces map: % 2u/% 2u. ss_id: % 2u/% 2u  >>>> % 2u/% 2u", object_index, map__mesh_face_id__face_info[f1_index].ss_id, f1_index, map__mesh_face_id__face_info[f2_index].ss_id, f2_index, f1_index, f2_index);
                        }
                      }
                      // Тут уже можно начать обрабатывать offset-ы текущего объекта object_id:
                      {
                        CGAL::Real_timer timer_bevel; // Тратится около 3% времени на merge и расчёт ajacent faces.
                        double t1 = 0;
                        //timer_bevel.start();
                        // Примечание: Все SS являются плоскими объектами, имеющими один или больше контуров: первый контур всегда внешний (если есть один, то это только внешний контур), остальные
                        // всегда внутренние. Рекурсии у SS нет! Если есть внутренние контуры, то топологически SS является двусторонним - внешняя и внутренняя стороны (не верх или низ как у плоскости,
                        // а именно контуры).

                        /// <summary>
                        /// Сопоставить определённым faces от SS набор точек пересечений с полурёбрами этого SS заданному offset.
                        /// Повтор другими словами: Есть набор offset с точками пересечений в SS, но в одном offset указывается вектор из нескольких SS,
                        /// с которыми этот offset пересекается. Требуется преобразовать этот набор в зависимость от SS, чтобы знать с какими offset пересекается этот ss_id.
                        /// </summary>
                        struct SS_Offset_HalfEdge_Intersections_Point_2{
                          /// <summary>
                          /// Информация об одном offset, относящемся к имеющемуся ss.
                          /// </summary>
                          //OIOA oioa;

                          FT oioa_offset;
                          int oioa_offset_index;
                          FT oioa_altitude;

                          /// <summary>
                          /// Какие точки пересечения с offset у текущего SS соответствуют halfedge текущего offset
                          /// </summary>
                          std::shared_ptr<Ss> ss;

                          /// <summary>
                          /// В каких точках пересекается ss с offset из oioa
                          /// </summary>
                          std::unordered_map<HDS_Halfedge_const_handle, Point_2> map__halfedge__offset_points;

                          //SS_Offset_HalfEdge_Intersections_Point_2(OIOA _oioa, std::shared_ptr<Ss> _ss, std::unordered_map<HDS_Halfedge_const_handle, Point_2> _map__halfedge__offset_points)
                          //  : oioa(_oioa), ss(_ss), map__halfedge__offset_points(_map__halfedge__offset_points){
                          //}
                          SS_Offset_HalfEdge_Intersections_Point_2(FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude, std::shared_ptr<Ss> _ss, std::unordered_map<HDS_Halfedge_const_handle, Point_2> _map__halfedge__offset_points)
                            : oioa_offset(_oioa_offset), oioa_offset_index(_oioa_offset_index), oioa_altitude(_oioa_altitude), ss(_ss), map__halfedge__offset_points(_map__halfedge__offset_points){
                          }
                        };

                        // Список всех offset у ss_id с информацией какие halfedge от исходного SS пересекаются с этим offset с сортировкой по offset_index (от меньшего к большему)
                        std::map<int /*ss_id*/, std::vector<SS_Offset_HalfEdge_Intersections_Point_2>> map__ss_id__he_intersections_point2;
                        std::map<int /*ss_id*/, std::shared_ptr<Ss> /*ss*/> map__ss_id__ss;

                        // Обработку контуров надо вести либо по map_adjacent_mesh_faces (если есть отрицательные SS), либо по всем faces от рассчитанных SS (если SS только положительные).
                        {
                          for (auto& [ss_id, oioa_offset_ss_params] : map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_index]) {
                            for (auto& offset_params : oioa_offset_ss_params) {
                              for (auto& he1 : offset_params.vect__SS__HalfEdge__Point_2 /*Какие точки пересечений offset и SS соответствуют halfedge*/) {
                                if (map__ss_id__he_intersections_point2.find(he1.ss_id) == map__ss_id__he_intersections_point2.end()) {
                                  map__ss_id__he_intersections_point2[he1.ss_id] = std::vector<SS_Offset_HalfEdge_Intersections_Point_2>();
                                }
                                map__ss_id__he_intersections_point2[he1.ss_id].push_back(SS_Offset_HalfEdge_Intersections_Point_2(offset_params.oioa1.offset, offset_params.oioa1.offset_index, offset_params.oioa1.altitude, he1.ss, he1.map__halfedge__offset_points));
                                map__ss_id__ss[he1.ss_id] = he1.ss;
                              }
                            }
                          }
                          // Отсортировать точки пересечения по offset_index (потому что соединение отступов происходит по offset_index, а не по величине отступа offset),
                          // чтобы позже при обходе контура выстроить точки пересечения в каждом halfedge по порядку offset_index (с учётом наклона/slope halfedge):
                          for (auto& [ss_id, vect_intersection_points] : map__ss_id__he_intersections_point2) {
                            std::sort(vect_intersection_points.begin(), vect_intersection_points.end(), [](const SS_Offset_HalfEdge_Intersections_Point_2& elem1, const SS_Offset_HalfEdge_Intersections_Point_2& elem2) {
                              bool res = elem1.oioa_offset_index < elem2.oioa_offset_index;
                              return res;
                            });
                          }

                          // offset-ы строятся на сегментах контуров SS. Если нет отрицательных SS, то нужно выбрать один face (любой из faces SS находится на контуре), если есть отрицательные SS, то нужно выбрать два смежных face.
                          int map__ajacent_faces_size = map__adjacent_mesh_faces.size();
                          int log_counter1 = 0;
                          // Вектор полигонов, которые должны получится после обхода:
                          struct OFFSET_POINT {
                            Point_3 point;
                            bool altitude_calculated;
                            OFFSET_POINT(Point_3 _point, bool _altitude_calculated)
                              : point(_point), altitude_calculated(_altitude_calculated) {

                            }
                          };

                          {
                            if (verbose == true) {
                              // Вывести список пересечений offset со всеми halfedge (для отладки). Потом удалить этот блок
                              printf("\n=====================\nIntersection points");
                              for (auto& [ss_id, vect_ss_offset_intersection_points /*относится только к одному ss_id*/] : map__ss_id__he_intersections_point2) {
                                printf("\nss_id=%u", ss_id);
                                for (auto& offset_intersections_point : vect_ss_offset_intersection_points) {
                                  printf("\n  offset=% 2.8lf: [", CGAL::to_double(offset_intersections_point.oioa_offset));
                                  for (auto& [he_handle, point] : offset_intersections_point.map__halfedge__offset_points) {
                                    printf("% 2u:(% 2.8lf, % 2.8lf) ", he_handle->slope(), CGAL::to_double(point.x()), CGAL::to_double(point.y()));
                                  }
                                  printf(" ]");
                                }
                              }
                            }
                          }
                          
                          /// <summary>
                          /// Точка, используемая для точек плана Beveled SS. Будет использоваться для рассчёта Beveled SS. Список этих точек представляет собой "план" по которому уже в 3D будет строится полная версия Beveled SS.
                          /// Beveled SS
                          /// </summary>
                          struct SHAPE_POINT {
                            /// <summary>
                            /// Индекс точки в общем массиве точек.
                            /// </summary>
                            int index = -1;

                            /// <summary>
                            /// Вектор точек типа POINT. Точка типа POINT не может напрямую использоваться при построении faces, ограниченной двумя offset, т.к. пары offset или цепочки offset могут пересекать эту точку несколько раз.
                            /// При этом сами offset определены однозначно и у них дубликатов быть не может, даже если они совпадают. Но использоваться этот вектор будет только когда
                            /// начнётся разделение SS сегмента на полигоны. Также возможны ситуации, что точка никогда не понадобится, если эта точка не будет находится между двумя offsets (хотя
                            /// в будущем точка может потребоваться для caps-ов).
                            /// <image url="..\code_images\file_0073.png" scale=".2"/>
                            /// </summary>
                            //std::vector<SHAPE_POINT> vect__POINTS;

                            /// <summary>
                            /// Сама точка. Только высота по Z у точки рассчитывается не всегда. Высота заранее известна у точек, лежащих на offset.
                            /// Если точка находится на контуре, то её высоту надо будет рассчитать позже. Это будет сделано только после обработки всех
                            /// точек.
                            /// </summary>
                            Point_3 point;

                            /// <summary>
                            /// Рассчитана ли высота у точки? (true - рассчитана, false - не рассчитана)
                            /// </summary>
                            bool is_altitude_calc;

                            enum POINT_TYPE {
                              /// <summary>
                              /// Точка из контура (не пересекающая edge). Одна и таже точка типа POINT может встречаться несколько раз при обработке переходов offset через эту точку при реверсных направлениях пар offset.
                              /// </summary>
                              POINT,

                              /// <summary>
                              /// Точка из offset. Точка пересекает offset вверх (в контуре. В faces это направление должно учитывать относительные положения edges, из которых строится face между 2-мя offset)
                              /// </summary>
                              UP,

                              /// <summary>
                              /// Точка из offset. Точка пересекает edge вниз
                              /// </summary>
                              DOWN
                            };

                            //OIOA oioa;
                            FT oioa_offset;
                            int oioa_offset_index;
                            FT oioa_altitude;

                            /// <summary>
                            /// Чтобы найти параметры offset.
                            /// </summary>
                            //int offset_index;

                            SHAPE_POINT(int _index, Point_3& _point, bool _is_altitude_calc, FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude)
                              : index(_index), point(_point), is_altitude_calc(_is_altitude_calc), oioa_offset(_oioa_offset), oioa_offset_index(_oioa_offset_index), oioa_altitude(_oioa_altitude) {

                            }

                            //SHAPE_POINT& operator=(const SHAPE_POINT& elem) {
                            //  this->index = elem.index;
                            //  this->point = elem.point;
                            //  this->is_altitude_calc = elem.is_altitude_calc;
                            //  this->oioa = elem.oioa;

                            //  return *this;
                            //}

                            SHAPE_POINT()
                              :index(-1), point(0.0, 0.0, 0.0), is_altitude_calc(false), oioa_offset(0.0), oioa_offset_index(-1), oioa_altitude(0.0)
                            {

                            }
                          };

                          //std::map<int, CALC_POINT> map__calc_points;

                          /// <summary>
                          /// Вектор результирующих точек SS-Beveled-проекции
                          /// </summary>
                          struct SHAPE_POINTS {
                            /// <summary>
                            /// Список точек с индексами.
                            /// </summary>
                            std::map<int /*point_index*/, SHAPE_POINT> map__point_index__calc_point;
                            
                            /// <summary>
                            /// Внутренний счётчик точек при создании новый точек
                            /// </summary>
                            int ss_points_counter=0;

                            /// <summary>
                            /// Сопоставление индексов результирующих точек пересечениям с полурёбрами he во время обхода segment-ов.
                            /// Параметр map__ss_id__he_handle и map__ss_id__vertex_id работают в паре при определении результирующих точек.
                            /// </summary>
                            std::map<std::tuple<int /*ss_id*/, HDS_Halfedge_const_handle, int /*offset_index*/>, int /*уникальный идентификатор точки*/> map__ss_id__he_handle;

                            /// <summary>
                            /// Сопоставление индексов результирующих точек с точками самого SS во время обхода segment-ов.
                            /// Параметр map__ss_id__he_handle и map__ss_id__vertex_id работают в паре при определении результирующих точек.
                            /// </summary>

                            std::map<std::tuple<int /*ss_id*/, int        /*vertex_id*/                       >, int /*уникальный идентификатор точки*/> map__ss_id__vertex_id;

                            int get_index_or_append_he__offset_index(int ss_id, HDS_Halfedge_const_handle& he_handle, int offset_index, Point_3& _point, bool _is_altitude_calc, FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude) {
                              auto& ss_id__he_handle__offset_index = std::tuple(ss_id, he_handle, offset_index);
                              auto& it = map__ss_id__he_handle.find(ss_id__he_handle__offset_index);
                              int res_counter = -1;
                              if (it == map__ss_id__he_handle.end()) {
                                res_counter = ss_points_counter++;
                                map__ss_id__he_handle[ss_id__he_handle__offset_index] = res_counter;
                                map__point_index__calc_point[res_counter] = SHAPE_POINT(res_counter, _point, _is_altitude_calc, _oioa_offset, _oioa_offset_index, _oioa_altitude);
                              } else {
                                auto& [tpl, _res_counter] = (*it);
                                res_counter = _res_counter;
                              }
                              return res_counter;
                              };

                            int get_index_or_append_vertex_id (int ss_id, int vertex_id, Point_3& _point, bool _is_altitude_calc, FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude) {
                              auto& ss_id__vertex_id = std::tuple(ss_id, vertex_id);
                              auto& it = map__ss_id__vertex_id.find(ss_id__vertex_id);
                              int res_counter = -1;
                              if (it == map__ss_id__vertex_id.end()) {
                                res_counter = ss_points_counter++;
                                map__ss_id__vertex_id[ss_id__vertex_id] = res_counter;
                                map__point_index__calc_point[res_counter] = SHAPE_POINT(res_counter, _point, _is_altitude_calc, _oioa_offset, _oioa_offset_index, _oioa_altitude);
                              } else {
                                auto& [tpl, _res_counter] = (*it);
                                res_counter = _res_counter;
                              }
                              return res_counter;
                            };
                          };

                          // Список точек, образующихся в результате рассчётов сегментов.
                          SHAPE_POINTS calc_points;

                          {
                            {
                              /// <summary>
                              /// Точки, используемые при обработке плана Project Beveled SS. Сами точки SHAPE_POINT не могут использоваться в плане, потому что
                              /// кроме индекса самой точки надо ещё помнить в каком направлении (slope) алгоритм через неё проходит при обходе faces против часовой стрелки.
                              /// Если SHAPE_POINT на месте пересечения находится на границе между двух рёбер, то в контур она попадает в двух экземплярах с учётом slope (типа UP или DOWN).
                              /// Если SHAPE_POINT типа POINT, то она попадает в контур в одном экземпляре и её slope не учитывается.
                              /// Пример 1: <image url="..\code_images\file_0075.png" scale="1.0"/>
                              /// Пример 2: <image url="..\code_images\file_0079.png" scale=".2"/>
                              /// </summary>
                              struct CONTOUR_SHAPE_POINT {
                                int point_index;
                                SHAPE_POINT::POINT_TYPE type;

                                CONTOUR_SHAPE_POINT(int _point_index, SHAPE_POINT::POINT_TYPE _type)
                                  :point_index(_point_index), type(_type){

                                }
                              };

                              // Подобрать данные по offset на каждый сегмент всех SS, но учесть, что если нет отрицательных SS, то все сегменты будут принадлежать соответствующим SS,
                              // но при наличии отрицательных SS требуется построить все сегменты только на отрицательных SS и привязать результаты подбора сегментов к ss_face_id с отрицательными SS. 
                              // Примечание: Если все offset сосредоточены только в положительных или только в отрицательных offset, то map_ajacent_faces будет пустым.
                              // Пример привязки offset к ss_face_id из отрицательного SS: <image url="..\code_images\file_0063.png" scale="1.0"/>

                              // Информация о сегментах для offset с привязкой к ss_id:
                              // <image url="..\code_images\file_0065.png" scale="1.0"/>
                              std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector<MeshFacePropertyStruct>>> map__ss_id__mesh_face_id__segment_faces; // Рассчитать segment_faces
                              std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector<CONTOUR_SHAPE_POINT>>>    map__ss_id__mesh_face_id__segment_contour; // Преобразовать segment_faces в segment_contour

                              // Обработать сегменты каждого SS (сегмент - это отрезок у "подножия" SS, граничащий либо с пустотой вокруг фигуры, если нет отрицательных SS либо с отрицательными SS):
                              for (auto& [ss_id, vect__he_intersections_point2 /*относится только к одному ss_id*/] : map__ss_id__he_intersections_point2) {

                                // Если есть смежные faces, то для расчёта нужно обрабатывать смежные faces по этой таблице.
                                // Примечание: если в таблице смежных faces есть записи, то все faces во всех SS имеют смежные faces.
                                // Если таблица смежных faces пустая, то можно обработать текущий ss_id в одиночку.
                                // Примечание: если в таблице смежных faces нет записей, то вообще не существует смежных faces (это возможно когда все offset сосредоточены в положительных SS
                                // или когда все offset сосредоточены в отрицательных offset)
                                if (map__ajacent_faces_size > 0) {
                                  // Проверить, есть ли в таблице смежных faces выбрать первый face выбранного ss_id (достаточно проверить только первый face, потому что все mesh_face_id в ключе будут от одного ss_id)
                                  auto& mesh_face_id = map__ss_id__ss_face_id__info[ss_id].begin()->second.mesh_face_id; /* выбирать первый элемент через begin() потому что если попадётся внешний отрицательный контур SS, полученный через виртуальный frame, то нумерация начнётся не с 0, а с 4, поэтому выбираем итератором */
                                  if (map__adjacent_mesh_faces.find(mesh_face_id) == map__adjacent_mesh_faces.end()) {
                                    // Не обрабатывать этот SS, т.к. faces этого SS будут обработаны с помощью смежных faces (либо уже были обработаны).
                                    continue;
                                  }
                                }

                                // Рассчитать какие faces надо использовать для сегментов:
                                for (auto& [ss_face_id, mesh_face_id0_info] : map__ss_id__ss_face_id__info[ss_id]) {

                                  // Рассчитать, какие faces нужны для построения набора offset на внешнем ребре SS (по одному face или парами смежных face?).
                                  std::vector<MeshFacePropertyStruct> segment_faces;
                                  segment_faces.push_back(mesh_face_id0_info); // Примечание - таблица map__ajacent_faces не гарантирует, что ключём будут только отрицательные SS, в value - +SS. Может быть и наоборот. Поэтому позже будет произведён реверс списка, чтобы на первом месте была ss_face_id из отрицательной SS.
                                  // Запомнить базовый mesh_face к которому будет привязка offset-ов
                                  MeshFacePropertyStruct& target_ss_face_info = mesh_face_id0_info;
                                  // Определить нужно ли обрабатывать только этот ss_face_id как один face или, если есть смежный face, то уже нужно обработать уже двоих вместе со смежным face?
                                  if (map__ajacent_faces_size > 0) {
                                    // Получить информацию о смежном face:
                                    auto& mesh_face_id1_info /*adjacent face*/ = map__mesh_face_id__face_info[map__adjacent_mesh_faces[mesh_face_id0_info.mesh_face_id]];

                                    segment_faces.push_back(mesh_face_id1_info);
                                    auto& vect_ss1_offset_intersection_points = map__ss_id__he_intersections_point2[mesh_face_id1_info.ss_id];

                                    if (mesh_face_id1_info.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE) {
                                      target_ss_face_info = mesh_face_id1_info;
                                      std::reverse(segment_faces.begin(), segment_faces.end()); // Чтобы отрицательный SS оказался на первом месте.
                                    }
                                  }

                                  if (map__ss_id__mesh_face_id__segment_faces.find(ss_id) == map__ss_id__mesh_face_id__segment_faces.end()) {
                                    map__ss_id__mesh_face_id__segment_faces[ss_id] = std::map<int /*mesh_face_id*/, std::vector<MeshFacePropertyStruct>>();
                                  }
                                  map__ss_id__mesh_face_id__segment_faces[ss_id][target_ss_face_info.mesh_face_id] = segment_faces;
                                }
                              }

                              {
                                // Преобразовать полученные данные в смесь точек исходного многоугольника и точек пересечения:
                                for (auto& [ss_id, mesh_face_id__segment_faces] : map__ss_id__mesh_face_id__segment_faces) {
                                  for (auto& [mesh_face_id, segment_faces] : mesh_face_id__segment_faces) {
                                    // Вектор halfedges для текущего сегмента (в текущем сегменте одна или две faces)
                                    std::vector<std::vector<HDS_Halfedge_const_handle>> vect_segment__halfedges_handle;
                                    // Все точки на контуре вокруг сегмента. <image url="..\code_images\file_0079.png" scale=".2"/>
                                    std::vector<CONTOUR_SHAPE_POINT> vect_contour_shape_point;
                                    for (int I = 0; I <= (int)segment_faces.size() - 1; I++) {
                                    //for (auto& segment_faces_I : segment_faces) {
                                      auto& segment_faces_I = segment_faces[I];
                                      std::vector<CGAL::Sign> face_slopes;
                                      auto& ss_face_handle = map__object_index__ss_id__ss_face_id__face_handle[object_index][segment_faces_I.ss_id][segment_faces_I.ss_face_id];
                                      vect_segment__halfedges_handle.push_back(std::vector<HDS_Halfedge_const_handle>());
                                      {
                                        HDS_Halfedge_const_handle hds_h = ss_face_handle->halfedge(), hds_h_start=hds_h, hds_h_finish = hds_h;
                                        // Обработать вместе точки контура и точки пересечения. Добавить их в контур в порядке, в котором они
                                        // встречаются при обходе контура Face.
                                        do {
                                          HDS_Halfedge_const_handle hds_offset_h = hds_h;
                                          int hds_h_slope = hds_h->slope();
                                          if (hds_h_slope == CGAL::NEGATIVE) { // ensure same geometric point on both sides
                                            hds_offset_h = hds_offset_h->opposite();
                                          }

                                          // Проверить все offset на предмет пересечения со всем offset. 
                                          struct POINT_INFO {
                                            /// <summary>
                                            /// Параметры offset, с которым у рассматриваемой halfedge возникло пересечение
                                            /// </summary>
                                            //OIOA oioa;
                                            FT  oioa_offset;
                                            int oioa_offset_index;
                                            FT  oioa_altitude;


                                            /// <summary>
                                            /// Точка пересечения с halfedge
                                            /// </summary>
                                            Point_2 point;

                                            //int intersection_index;

                                            //POINT_INFO(OIOA& _oioa, Point_2& _point )
                                            //  :oioa(_oioa), point(_point) {
                                            //}
                                            POINT_INFO(FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude, Point_2& _point )
                                              :oioa_offset(_oioa_offset), oioa_offset_index(_oioa_offset_index), oioa_altitude(_oioa_altitude), point(_point) {
                                            }
                                          };
                                          // Список offset с которыми произошли пересечения
                                          std::vector<POINT_INFO> vect__intersection_point_info;
                                          // В зависимости от уклона ребра SS и знака SS необходимо правильно определить последовательность добавления точек пересечения с ребром.
                                          // При движении "вверх" при положительном +SS индексы добавляются последовательно как они и объявлены клиентом,
                                          // а при движении вниз - в обратном порядке. В отрицательной области -SS надо делать наоборот - добавлять индексы
                                          // в обратной последовательности при движении (-)"наружу"/"вниз" и добавлять в прямой последовательности при
                                          // движении "вверх"/"внутрь" по направлении линии 0-0: <image url="..\code_images\file_0071.png" scale="1.0"/>
                                          // Добавлять точки пересечения с he в направлении обхода контура
                                          if ( hds_h_slope == CGAL::NEGATIVE && segment_faces_I.ss_sign==MeshFacePropertyStruct::SS_SIGN::POSITIVE /* +SS */
                                            ||
                                            hds_h_slope == CGAL::POSITIVE && segment_faces_I.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE /* -SS */
                                            ) {
                                            // <image url="..\code_images\file_0076.png" scale="1.0"/>
                                            //for (auto& intersection_points : map_ss_id_he_intersections_point2[segment_faces_I.ss_id] /*ss_id брать не из цикла, а из рассматриваемой segment_faces_I*/) 
                                            for (auto& elem_it = map__ss_id__he_intersections_point2[segment_faces_I.ss_id /* ss_id брать не из цикла, а из рассматриваемой segment_faces_I, чтобы искать пересечения offset-ов со своим he */ ].rbegin(); elem_it != map__ss_id__he_intersections_point2[segment_faces_I.ss_id].rend(); ++elem_it){
                                              auto& intersection_points = *elem_it;
                                              auto& offset0_he_intersect_iter = intersection_points.map__halfedge__offset_points.find(hds_offset_h);
                                              if (offset0_he_intersect_iter != intersection_points.map__halfedge__offset_points.end()) {
                                                auto& [HDS_Halfedge_const_handle, point] = *offset0_he_intersect_iter;
                                                vect__intersection_point_info.push_back(POINT_INFO(intersection_points.oioa_offset, intersection_points.oioa_offset_index, intersection_points.oioa_altitude, point));
                                              }
                                            }
                                          } else {
                                            // <image url="..\code_images\file_0077.png" scale="1.0"/>
                                            for (auto& elem_it = map__ss_id__he_intersections_point2[segment_faces_I.ss_id].begin(); elem_it != map__ss_id__he_intersections_point2[segment_faces_I.ss_id].end(); ++elem_it)
                                            {
                                              auto& intersection_points = *elem_it;
                                              auto& offset0_he_intersect_iter = intersection_points.map__halfedge__offset_points.find(hds_offset_h);
                                              if (offset0_he_intersect_iter != intersection_points.map__halfedge__offset_points.end()) {
                                                auto& [HDS_Halfedge_const_handle, point] = *offset0_he_intersect_iter;
                                                vect__intersection_point_info.push_back(POINT_INFO(intersection_points.oioa_offset, intersection_points.oioa_offset_index, intersection_points.oioa_altitude, point));
                                              }
                                            }
                                          }
                                          
                                          // Добавить точки в контур
                                          {
                                            {
                                              if (hds_h_slope == CGAL::ZERO) {
                                                // Это ортогональная линия SS. Предполагаем, что не бывает подряд две ортогональных линии в одном направлении SS в одном face (Update - всё-таки есть исключение, но пока не знаю на что это влияет).
                                                // update: Проверил, что CGAL не пересекает между собой горизонтальные линии и горизонтальные halfedges. (Вообще это очень
                                                // неприятная ситуация, потому что точка пересечения по сути превращается в линию пересечения)
                                                // TODO: Временно пропускаю добавление каки-либо пересечений с ортогональной линией, посмотрим, что будет.
                                                // <image url="..\code_images\file_0062.png" scale="1.0"/> 
                                                // Нашлось исключение: <image url="..\code_images\file_0072.png" scale="1.0"/> (Надо будет проверить может ли оно влиять на рассчёт?)
                                              } else {
                                                SHAPE_POINT::POINT_TYPE pt = SHAPE_POINT::POINT_TYPE::UP;
                                                if (
                                                  hds_h_slope == CGAL::NEGATIVE
                                                  ) {
                                                  pt = SHAPE_POINT::POINT_TYPE::DOWN;
                                                }
                                                // Если face от отрицательного -SS, то инвертировать направление точки пересечения
                                                // <image url="..\code_images\file_0069.png" scale=".2"/>
                                                if (segment_faces_I.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE) {
                                                  switch (pt) {
                                                    case SHAPE_POINT::POINT_TYPE::UP:
                                                      pt = SHAPE_POINT::POINT_TYPE::DOWN;
                                                      break;
                                                    case SHAPE_POINT::POINT_TYPE::DOWN:
                                                      pt = SHAPE_POINT::POINT_TYPE::UP;
                                                      break;
                                                  }
                                                }
                                                // TODO - так-то их можно сразу обрабатывать на месте обнаружения. Подумать, как сделать это прямо в точке обнаружения выше.
                                                for (auto& he_offset : vect__intersection_point_info) {
                                                  int point_index = calc_points.get_index_or_append_he__offset_index(ss_id, hds_offset_h, he_offset.oioa_offset_index, Point_3(he_offset.point.x(), he_offset.point.y(), he_offset.oioa_altitude), true, he_offset.oioa_offset, he_offset.oioa_offset_index, he_offset.oioa_altitude );
                                                  vect_contour_shape_point.push_back(CONTOUR_SHAPE_POINT(point_index, pt));
                                                }
                                              }
                                            }

                                            // При обработке смежных faces и текущем SS, если он отрицательный (а при segment_faces.size() > 1 один из SS точно отрицательный)
                                            // пропускать добавление первой и последней точек (но обработка выше должна быть).
                                            if (segment_faces.size() > 1 && segment_faces_I.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE && (hds_h == hds_h_start || hds_h->next() == hds_h_finish)) {
                                              // Сначала была идея, что при обработке двух смежных faces не использовать их первые точки, но
                                              // столкнулся с проблемой задвоения индексов, когда встречались конечная и начальная точки от 
                                              // двух смежных SS, которые визуально совпадают, а топологически - это две разные точки.
                                              // Поэтому переделал этот способ и исключил добавление первой и последней точки у -SS при смежных faces.
                                              // <image url="..\code_images\file_0074.png" scale="1.0"/>
                                            } else {
                                              // Добавить точку контура в список (если она не первая в списке при двух смежных faces. Если faces>1, то тут первая точка уже пропущена уже пропущена)
                                              int he_vertex_id = hds_h->vertex()->id();
                                              int point_index = calc_points.get_index_or_append_vertex_id(segment_faces_I.ss_id, he_vertex_id, Point_3(hds_h->vertex()->point().x(), hds_h->vertex()->point().y(), 0 /*Z не выставляется (это altitude), будет считаться позже*/), false /*altitude не рассчитана, будет считаться позже*/, 0.0, -1 /*Это точка контура, поэтому offset-а у неё нет*/, 0.0 );
                                              vect_contour_shape_point.push_back(CONTOUR_SHAPE_POINT(point_index, SHAPE_POINT::POINT_TYPE::POINT));
                                            }
                                          }
                                          vect_segment__halfedges_handle.back().push_back(hds_offset_h);
                                          face_slopes.push_back(hds_h->slope()); // Предположительно это указывает направление движения к или от центра SS. update - нет, если представить себе, что SS является объёмной фигурой типа гор с вершинами, то edges являются ребрами на пиках, а slope указывает как смещается это полуребро по вертикали -1 - вниз, 0 - горизонтальное, +1 - вверх. Соответственно противоположное полуребро является инверсией по направлению.
                                          hds_h = hds_h->next();
                                        } while (hds_h != hds_h_finish);

                                        // Проверка, что наклонные slope() на противоположных (opposite()) рёбрах ожидаются противоположных знаков (кроме 0 - инверсия от 0 - всегда 0).
                                        // <image url="..\code_images\file_0061.png" scale=".1"/>
                                        if (verbose) {
                                          printf("\n " VAL2STR(Err::_0034) " ss_id=% 2u, ss_face_id % 2u, mesh_face_id % 2u, slopes: ", segment_faces_I.ss_id, segment_faces_I.ss_face_id, segment_faces_I.mesh_face_id);
                                          for (auto& slope : face_slopes) {
                                            printf("% 2d ", slope);
                                          }
                                        }
                                      }
                                    }

                                    if (map__ss_id__mesh_face_id__segment_contour.find(ss_id) == map__ss_id__mesh_face_id__segment_contour.end()) {
                                        map__ss_id__mesh_face_id__segment_contour[ss_id] = std::map<int /*ss_face_id*/, std::vector<CONTOUR_SHAPE_POINT>>();
                                      }
                                    // Пример распределения offsets одного сегмента при смежных faces и одновременного нахождения offsets и в отрицательном -SS и в положительном +SS:
                                    // <image url="..\code_images\file_0070.png" scale="1.0"/>
                                    map__ss_id__mesh_face_id__segment_contour[ss_id][mesh_face_id] = vect_contour_shape_point;
                                  }
                                }
                              }

                              // Параметры курсора для обхода точек
                              struct COLLECT_CURSOR {
                                /// <summary>
                                /// (true) Собирать или (false) пропускать точки тип POINT
                                /// </summary>
                                bool do_collect_points;

                                /// <summary>
                                /// Создавать ли новые faces? (если offset-ы находятся по обеим сторонам 0-0, то возможен только один face и его
                                /// надо создать заранее и новые не добавлять. Тогда все корректные точки (которые выше и ниже offset) надо добавить в уже созданный face.
                                /// </summary>
                                bool do_create_new_faces;

                                /// <summary>
                                /// После окончания обработки пары offset инвертировать faces
                                /// </summary>
                                bool do_reverse;

                                /// <summary>
                                /// Кто начинает новый face?
                                /// Индекс для первого пересечения, который начинает новые faces (+ добавить эту точку)
                                /// </summary>
                                int offset_index0;
                                /// <summary>
                                /// Тип первого пересечения, который начинает новые faces, после которого надо добавить точки и добавить саму эту точку.
                                /// </summary>
                                SHAPE_POINT::POINT_TYPE offset_type0;

                                /// <summary>
                                /// Индекс для второго пересечения, после которого надо добавить новую точку POINT (+ добавить эту точку)
                                /// </summary>
                                int offset_index1;

                                /// <summary>
                                /// Тип второго пересечения, после которого нужно собирать точки и добавлять их в последний face
                                /// </summary>
                                SHAPE_POINT::POINT_TYPE offset_type1;

                                COLLECT_CURSOR(bool _do_collect_points, bool _do_create_new_faces, bool _do_reverse, int _index_offset0, SHAPE_POINT::POINT_TYPE _offset_type0, int _index_offset1, SHAPE_POINT::POINT_TYPE _offset_type1)
                                  :do_collect_points(_do_collect_points), do_create_new_faces(_do_create_new_faces), do_reverse(_do_reverse), offset_index0(_index_offset0), offset_type0(_offset_type0), offset_index1(_index_offset1), offset_type1(_offset_type1){

                                }
                              };

                              /// <summary>
                              /// Состав и исходные данные одной Face в наборе faces одного сегмента.
                              /// </summary>
                              struct FACE_INFO {
                                /// <summary>
                                /// Вектор индексов вершин для текущего face
                                /// </summary>
                                std::vector/*point indexes*/<int> face_verts_indexes;

                                /// <summary>
                                /// Первый отступ
                                /// </summary>
                                //OIOA oioa0;
                                int oioa0_offset_index;
                                FT oioa0_offset;

                                /// <summary>
                                /// Второй отступ
                                /// </summary>
                                //OIOA oioa1;
                                int oioa1_offset_index;
                                FT oioa1_offset;

                                /// <summary>
                                /// Выполнять ли reverse. (По хорошему этот параметр должен быть где-то снаружи, т.к. это однократная операция)
                                /// </summary>
                                bool do_reverse;

                                //FACE_INFO(OIOA& _oioa0, OIOA& _oioa1, bool& _do_reverse)
                                //: oioa0(_oioa0), oioa1(_oioa1), do_reverse(_do_reverse){
                                //  
                                //}

                                FACE_INFO(int _oioa0_offset_index, FT _oioa0_offset, int _oioa1_offset_index, FT _oioa1_offset, bool& _do_reverse)
                                  : oioa0_offset_index(_oioa0_offset_index), oioa0_offset(_oioa0_offset), oioa1_offset_index(_oioa1_offset_index), oioa1_offset(_oioa1_offset), do_reverse(_do_reverse) {
                                  
                                }
                              };

                              // На этот момент получены списки offset-ов на всех сегментах. Все offset распределены по индексам относительно линии 0-0 и знаку SS. Надо создать faces.
                              // faces создаются на основе offset_index
                              //std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector/*faces*/<std::vector/*point indexes*/<int>>>> map__ss_id__mesh_face_id__faces__indexes;
                              res_contours;
                              std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>> map__ss_id__mesh_face_id__faces__indexes; // Состав одного сегмента (сразу привязывается к mesh_face_id).
                              for (auto& [ss_id, map_mesh_face_id_segment_points] : map__ss_id__mesh_face_id__segment_contour) {
                                // Список offset-ов;
                                std::set<int> set__offset_index;
                                std::vector<OIOA> vect_oioa;
                                for (auto& [ss_id, oioa_offset_ss_params] : map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_index]) {
                                  for (auto& offset_params : oioa_offset_ss_params) {
                                    if (set__offset_index.find(offset_params.oioa1.offset_index) == set__offset_index.end()) {
                                      set__offset_index.insert(offset_params.oioa1.offset_index);
                                      vect_oioa.push_back(offset_params.oioa1);
                                    }
                                  }
                                }
                                std::sort(vect_oioa.begin(), vect_oioa.end(), [](const OIOA& o1, const OIOA& o2) {
                                  bool res = o1.offset_index < o2.offset_index;
                                  return res;
                                });

                                map__ss_id__mesh_face_id__faces__indexes[ss_id] = std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>();
                                for (auto& [mesh_face_id, vect_segment_points] : map_mesh_face_id_segment_points) {
                                  //map__ss_id__mesh_face_id__faces__indexes[ss_id][mesh_face_id] = std::vector/*faces*/<FACE_INFO>();
                                  
                                  // Список edges, которые надо обработать по offset. В будущем они будут заменены входными сокетами ноды.
                                  // Пока просто линейная последовательность (0,1),(1,2),(2,3)...,(n-2, n-1). Также эта последовательность должна определять
                                  // направления нормалей faces

                                  struct EDGE_INFO {
                                    //OIOA oioa0;
                                    int oioa0_offset_index;
                                    FT  oioa0_offset;

                                    //OIOA oioa1;
                                    int oioa1_offset_index;
                                    FT  oioa1_offset;

                                    //EDGE_INFO(OIOA& _oioa0, OIOA& _oioa1)
                                    //  : oioa0(_oioa0), oioa1(_oioa1) {

                                    //}
                                    EDGE_INFO(int _oioa0_offset_index, FT _oioa0_offset, int _oioa1_offset_index, FT _oioa1_offset)
                                      : oioa0_offset_index(_oioa0_offset_index), oioa0_offset(_oioa0_offset), oioa1_offset_index(_oioa1_offset_index), oioa1_offset(_oioa1_offset) {

                                    }
                                  };
                                  std::vector<EDGE_INFO> vect_edges;
                                  for (int I = 0; I <= (int)vect_oioa.size() - 2; I++) {
                                    vect_edges.push_back(EDGE_INFO(vect_oioa[I+0].offset_index,vect_oioa[I+0].offset, vect_oioa[I+1].offset_index, vect_oioa[I+1].offset_index));
                                  }

                                  std::vector<FACE_INFO>& vect_faces /*список faces на текущем сегменте*/ = std::vector/*faces*/<FACE_INFO>();
                                  for (auto& edge_info : vect_edges) {
                                    auto& _oioa0_offset_index = edge_info.oioa0_offset_index;
                                    auto& _oioa0_offset = edge_info.oioa0_offset;
                                    auto& _oioa1_offset_index = edge_info.oioa1_offset_index;
                                    auto& _oioa1_offset = edge_info.oioa1_offset;
                                    
                                    // Сначала определить какой сегмент по абсолютному значению ближе к линии 0-0?
                                    auto& offset0     = vect_oioa[_oioa0_offset_index].offset;
                                    auto& offset1     = vect_oioa[_oioa1_offset_index].offset;

                                    // Подсчёт faces всегда выполняется против часовой стрелки от первой точки после старта.
                                    COLLECT_CURSOR& cursor = COLLECT_CURSOR(false, false, false,-1, SHAPE_POINT::POINT_TYPE::POINT, -1, SHAPE_POINT::POINT_TYPE::POINT);
                                    // Рассчёт зависит от того в каком отношении к 0-0 находятся рассматриваемые offset. Мысленно представляем куда движется стерка обхода и
                                    // в какой последовательности она должна пересечь offset-ы.
                                    if (offset0 < 0 && offset1 < 0) {
                                      // оба offset ниже 0-0
                                      if (offset0 > offset1) {
                                        // 0
                                        // 1
                                        cursor = COLLECT_CURSOR( false, true,  true, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                      }else if (offset1 > offset0) {
                                        // 1
                                        // 0
                                        cursor = COLLECT_CURSOR( false, true,  false, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                      } else
                                        /* остальные, когда offset0 == offset1 */
                                      if (_oioa0_offset_index < _oioa1_offset_index) {
                                        // 1
                                        // 0
                                        cursor = COLLECT_CURSOR( false, true, false, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                      } else if (_oioa0_offset_index > _oioa1_offset_index) {
                                        // 0
                                        // 1
                                        cursor = COLLECT_CURSOR( false, true,  true, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                      } 
                                    } else if (offset0 > 0 && offset1 > 0) {
                                      // оба offset выше 0-0
                                      if (offset0 > offset1) {
                                        // 0
                                        // 1
                                        cursor = COLLECT_CURSOR( false, true,  true, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                      }else if (offset1 > offset0) {
                                        // 1
                                        // 0
                                        cursor = COLLECT_CURSOR( false, true, false, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                      } else
                                        /* остальные, когда offset0 == offset1 */
                                        if (_oioa0_offset_index > _oioa1_offset_index) {
                                          cursor = COLLECT_CURSOR(  false, true,  true, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                        } else if (_oioa1_offset_index > _oioa0_offset_index) {
                                          cursor = COLLECT_CURSOR(  false, true, false, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                        } 
                                    } else if (offset0 < 0 && offset1 > 0) {
                                      // Добавить новый face и разрешить добавлять в него точки типа POINT сразу
                                      cursor = COLLECT_CURSOR(true, false,false, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                      vect_faces.push_back(FACE_INFO(_oioa0_offset_index, _oioa0_offset, _oioa1_offset_index, _oioa1_offset, cursor.do_reverse));
                                    } else if (offset0 > 0 && offset1 < 0) {
                                      // Добавить новый face и разрешить добавлять в него точки типа POINT сразу
                                      cursor = COLLECT_CURSOR(true, false,  true, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                      vect_faces.push_back(FACE_INFO(_oioa0_offset_index, _oioa0_offset, _oioa1_offset_index, _oioa1_offset, cursor.do_reverse));
                                    } else {
                                      // Осталось, когда одна из точек в нуле
                                    }

                                    for (auto& contour_point : vect_segment_points) {
                                      auto& calc_point = calc_points.map__point_index__calc_point[contour_point.point_index];
                                      //calc_point.type;  .oioa.offset_index
                                      if (contour_point.type == SHAPE_POINT::POINT_TYPE::POINT) {
                                        if (cursor.do_collect_points == true) {
                                          vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                        }
                                        continue;
                                      }

                                      if (contour_point.type == cursor.offset_type0 && calc_point.oioa_offset_index== cursor.offset_index0) {
                                        if (cursor.do_create_new_faces == true) {
                                          // Создать новый face и запомнить эту точку:
                                          vect_faces.push_back(FACE_INFO(_oioa0_offset_index, _oioa0_offset, _oioa1_offset_index, _oioa1_offset, cursor.do_reverse));
                                        }
                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                        cursor.do_collect_points = true;
                                        continue;
                                      }
                                      if (contour_point.type != cursor.offset_type0 && calc_point.oioa_offset_index == cursor.offset_index0) {
                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                        cursor.do_collect_points = false;
                                        continue;
                                      }

                                      if (contour_point.type == cursor.offset_type1 && calc_point.oioa_offset_index == cursor.offset_index1) {
                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                        cursor.do_collect_points = true;
                                        continue;
                                      }

                                      if (contour_point.type != cursor.offset_type1 && calc_point.oioa_offset_index == cursor.offset_index1) {
                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                        cursor.do_collect_points = false;
                                        continue;
                                      }
                                    }
                                  }
                                  // Проверить vect_faces. Если там были faces с признаком reverse, то инвертировать порядок индексов:
                                  for (auto& face_info : vect_faces) {
                                    if (face_info.do_reverse == true) {
                                      std::reverse(face_info.face_verts_indexes.begin(), face_info.face_verts_indexes.end());
                                    }
                                  }
                                  map__ss_id__mesh_face_id__faces__indexes[ss_id][mesh_face_id] = vect_faces;
                                }
                              }

                              if (verbose) {
                                printf("\n ====================== LIST OF CONTOURS of object_index %u: =========================\n", object_index);
                                for (auto& [ss_id, map_mesh_face_id_segment_points] : map__ss_id__mesh_face_id__segment_contour) {
                                  for (auto& [mesh_face_id, vect_segment_points] : map_mesh_face_id_segment_points) {
                                    printf("\n " VAL2STR(Err::_0035) " contour: ss_id %u, mesh_face_id %u, length %u:", ss_id, mesh_face_id, (int)vect_segment_points.size());
                                    int counter = 0;
                                    for (auto& contour_point : vect_segment_points) {
                                      auto& calc_point = calc_points.map__point_index__calc_point[contour_point.point_index];
                                      // Временная сквозная нумерация точек пересечений и исходных точек в пределах одного segment. Позже надо будет разобраться как пронумеровать ВСЕ точки объекта
                                      printf("\n  % 2.8lf, % 2.8lf, % 2.8lf, type: %s, alt: %s:% 2.8lf, offset_index: % 3d, % 2.8lf; % 4u / % 3u", CGAL::to_double(calc_point.point.x()), CGAL::to_double(calc_point.point.y()), CGAL::to_double(calc_point.point.z()), contour_point.type == 0 ? "POINT" : contour_point.type == 1 ? "   UP" : " DOWN", calc_point.is_altitude_calc == true ? "    set" : "not set", CGAL::to_double(calc_point.oioa_altitude), calc_point.oioa_offset_index, CGAL::to_double(calc_point.oioa_offset), contour_point.point_index, counter++);
                                    }
                                    // Вывод только тех точек, которые образуют faces:
                                    int face_id = 0;
                                    for (auto& face_info : map__ss_id__mesh_face_id__faces__indexes[ss_id][mesh_face_id]) {
                                      printf("\n      face_id % 3u offsets: (% 2d, % 2d : % 2.8lf, % 2.8lf, reversed: %s); len: (% 3d): ", face_id, face_info.oioa0_offset_index, face_info.oioa1_offset_index, CGAL::to_double(face_info.oioa0_offset), CGAL::to_double(face_info.oioa1_offset), face_info.do_reverse == true ? "yes" : " no", (int)face_info.face_verts_indexes.size());
                                      for (auto& vert_index : face_info.face_verts_indexes) {
                                        printf(" % 3u", vert_index);
                                      }
                                      face_id++;
                                    }
                                  }
                                }
                                printf("\n ====================== /LIST OF CONTOURS of object_index %u: =========================\n", object_index);
                              }
                            }
                          }


                          // РАССЧЁТ геометрии SS для вывода (без offset-ов)

                          // Если есть смежные faces, то все данные для offset привязывать только к ss_id отрицательных SS (самый внешний + вложенные отрицательные SS)
                          // Примечание: Если все offset сосредоточены только в положительных или только в отрицательных offset, то map_ajacent_faces будет пустым.
                          // Пример привязки offset к ss_face_id из отрицательного SS: <image url="..\code_images\file_0063.png" scale="1.0"/>
                          std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector<SHAPE_POINT>>> map_ss_id_mesh_face_id_segment_points;

                          // Все результирующие точки по offset-ам в рассматриваемом SS:
                          std::vector<std::vector<std::vector<OFFSET_POINT>>> vect_ss_offset_faces;
                          // Обработать сегменты каждого SS (сегмент - это отрезок у "подножия" SS, граничащий либо с пустотой вокруг фигуры, если нет отрицательных SS либо с отрицательными SS):
                          for (auto& [ss_id, vect_ss_offset_intersection_points /*относится только к одному ss_id*/] : map__ss_id__he_intersections_point2) {

                            // Если есть смежные faces, то для расчёта нужно обрабатывать смежные faces по этой таблице.
                            // Примечание: если в таблице смежных faces есть записи, то все faces во всех SS имеют смежные faces.
                            // Если таблица смежных faces пустая, то можно обработать текущий ss_id в одиночку. Примечание: если в таблице смежных faces нет записей, то вообще не существует смежных faces (это только при положительных SS)
                            if (map__ajacent_faces_size > 0) {
                              // Проверить, есть ли в таблице смежных faces выбрать первый face выбранного ss_id (достаточно проверить только первый face, потому что все mesh_face_id в ключе будут от одного ss_id)
                              auto& mesh_face_id = map__ss_id__ss_face_id__info[ss_id].begin()->second.mesh_face_id; /* выбирать первый элемент через begin() потому что если попадётся внешний отрицательный контур SS, полученный через виртуальный frame, то нумерация начнётся не с 0, а с 4, поэтому выбираем итератором */
                              if (map__adjacent_mesh_faces.find(mesh_face_id) == map__adjacent_mesh_faces.end()) {
                                // Не обрабатывать этот SS, т.к. faces этого SS будут обработаны с помощью смежных faces (либо уже были обработаны).
                                continue;
                              }
                            }
                            
                            // Все точки в рассматриваемом SS:
                            // Пройтись и обработать faces у рассматриваемого ss_id
                            double contour_time = 0;
                            for (auto& [ss_face_id, face_info] : map__ss_id__ss_face_id__info[ss_id]) {

                              CGAL::Real_timer timer_contour_points_001;
                              timer_contour_points_001.start();

                              // Рассчитать, какие faces нужны для построения набора offset на внешнем ребре SS (по одному face или парами смежных face?).
                              std::vector<MeshFacePropertyStruct> segment_faces;

                              auto& mesh_face_id0_info /*source face*/ = map__ss_id__ss_face_id__info[ss_id][ss_face_id];
                              segment_faces.push_back(mesh_face_id0_info);
                              // Запомнить базовый mesh_face к которому будет привязка offset-ов
                              MeshFacePropertyStruct& target_ss_face_info = mesh_face_id0_info;
                              // Определить нужно ли обрабатывать только этот ss_face_id как один face или, если есть смежный face, то уже нужно обработать уже двоих вместе со смежным face?
                              if (map__ajacent_faces_size > 0) {
                                // Получить информацию о смежном face:
                                auto& mesh_face_id1_info /*adjacent face*/ = map__mesh_face_id__face_info[ map__adjacent_mesh_faces[mesh_face_id0_info.mesh_face_id] ];
                                // перед тем как добавлять смежный face надо проверить, а вдруг offset-ы находятся в одном ss_id?
                                if (mesh_face_id0_info.ss_id == mesh_face_id1_info.ss_id) {
                                  // Значит второй face не понадобиться.
                                } else {
                                  segment_faces.push_back(mesh_face_id1_info);
                                  auto& vect_ss1_offset_intersection_points = map__ss_id__he_intersections_point2[mesh_face_id1_info.ss_id];
                                  // Тут нужно просуммировать контуры смежных faces от положительного и отрицательного SS так, чтобы идти по halfedge в одном направлении,
                                  // вокруг этого смежного контура
                                  
                                }
                                if (mesh_face_id1_info.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE) {
                                  target_ss_face_info = mesh_face_id1_info;
                                }
                              }
//#ifdef DO_NOT_USE

                              /// <summary>
                              /// <summary>
                              struct CONTOUR_POINTS {
                                std::vector<SHAPE_POINT> points;
                              };

                              struct FACE_POINTS {
                                // Виртуально в этом классе считается, что точка входа в face через нижнюю границы - это UP, а выход через нижнюю границу  - это down,
                                // вход через верхнюю границу - DOWN, выход через UP. Именно поэтому offset_bottom должен быть указан внизу (с меньшим offset), а offset_top - должен быть указан наверху (с большим offset)
                                // Если оба offset равны (может быть при построении вертикальных стенок),  то нужно учитывать уже altitude. Если offset и alitude совпали (может так случится при анимации), то учитывается
                                // очерёдность этих offset в исходных данных.


                                /// <summary>
                                /// больший offset
                                /// </summary>
                                OIOA& oioa_top;

                                /// <summary>
                                /// нижний offset
                                /// </summary>
                                OIOA& oioa_bottom;

                                enum DIRECTION {
                                  /// <summary>
                                  /// Если изначально нижний и верхний offset также и являются нижним и верхним рассчётным offset
                                  /// </summary>
                                  COUNTER_CLOCKWISE,

                                  /// <summary>
                                  /// Если в исходных данных нижний и верхний offset указаны наоборот
                                  /// </summary>
                                  CLOCKWISE
                                };
                                /// <summary>
                                /// Направление обхода точек в зависимости от исходного отношения offset
                                /// </summary>
                                DIRECTION direction;

                                std::vector<SHAPE_POINT> points;
                              };
                              std::vector<SHAPE_POINT> segment_contour_points;
//#endif
                              // Получить face_handle для рассматриваемых faces (либо одна face, либо 2 смежные):
                              std::vector<HDS_Face_handle> vect_segment_faces_handle;
                              // Загрузить halfedges от заданных faces в этот вектор:
                              std::vector<std::vector<HDS_Halfedge_const_handle>> vect_segment_faces_edges_handle;
                              for (auto& segment_face : segment_faces) {
                                std::vector<CGAL::Sign> face_slopes;
                                auto& ss_face_handle = map__object_index__ss_id__ss_face_id__face_handle[object_index][segment_face.ss_id][segment_face.ss_face_id];
                                vect_segment_faces_handle.push_back(ss_face_handle);
                                vect_segment_faces_edges_handle.push_back(std::vector<HDS_Halfedge_const_handle>());
                                {
                                  HDS_Halfedge_const_handle hds_h = ss_face_handle->halfedge(), done = hds_h;
                                  do {
                                    HDS_Halfedge_const_handle hds_offset_h = hds_h;
                                    int hds_h_slope = hds_h->slope();
                                    if (hds_h_slope == CGAL::NEGATIVE) { // ensure same geometric point on both sides
                                      hds_offset_h = hds_offset_h->opposite();
                                    }

                                    // Проверить все offset на предмет пересечения со всем offset. 
                                    struct POINT_INFO {
                                      /// <summary>
                                      /// Параметры offset, с которым у рассматриваемой halfedge возникло пересечение
                                      /// </summary>
                                      //OIOA oioa;
                                      int oioa_offset_index;
                                      FT  oioa_offset;
                                      int oioa_altitude;

                                      /// <summary>
                                      /// Точка пересечения с halfedge
                                      /// </summary>
                                      Point_2 point;

                                      //POINT_INFO(OIOA& _oioa, Point_2& _point)
                                      //  :oioa(_oioa), point(_point) {
                                      //}
                                      POINT_INFO(FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude, Point_2& _point)
                                        :oioa_offset(_oioa_offset), oioa_offset_index(_oioa_offset_index), oioa_altitude(_oioa_altitude), point(_point) {
                                      }
                                    };
                                    // Список offset с которыми произошли пересечения
                                    std::vector<POINT_INFO> vect_he_offset;
                                    for (auto& intersection_points : vect_ss_offset_intersection_points) {
                                      auto& offset0_he_intersect_iter = intersection_points.map__halfedge__offset_points.find(hds_offset_h);
                                      if (offset0_he_intersect_iter != intersection_points.map__halfedge__offset_points.end()) {
                                        auto& [HDS_Halfedge_const_handle, point] = *offset0_he_intersect_iter;
                                        vect_he_offset.push_back(POINT_INFO(intersection_points.oioa_offset, intersection_points.oioa_offset_index, intersection_points.oioa_altitude, point));
                                      }
                                    }
                                    // Добавить точки в контур
                                    {
                                      std::sort(vect_he_offset.begin(), vect_he_offset.end(), [](const POINT_INFO& p1, const POINT_INFO& p2) {
                                        bool res = false;
                                        if (p1.oioa_offset > p2.oioa_offset) {
                                          res = true;
                                        }
                                        //else if (p1.oioa.offset == p2.oioa.offset) {
                                        //  if (p1.oioa.offset_index > p2.oioa.offset_index) {
                                        //    res = true;
                                        //  }
                                        //}
                                        return res;
                                      });
                                      {
                                        if (hds_h_slope == CGAL::ZERO) {
                                          // Это ортогональная линия SS. Предполагаем, что не бывает подряд две ортогональных линии в одном направлении SS в одном face.
                                          // update: Проверил, что CGAL не пересекает между собой горизонтальные линии и горизонтальные halfedges. (Вообще это очень
                                          // неприятная ситуация, потому что точка пересечения по сути превращается в линию пересечения)
                                          // TODO: Временно пропускаю добавление каки-либо пересечений с ортогональной линией, посмотрим, что будет.
                                          // <image url="..\code_images\file_0062.png" scale="1.0"/>
                                        }else{
                                          SHAPE_POINT::POINT_TYPE pt = SHAPE_POINT::POINT_TYPE::UP;
                                          if (hds_h_slope == CGAL::NEGATIVE) {
                                            pt = SHAPE_POINT::POINT_TYPE::DOWN;
                                          }
                                          for (auto& he_offset : vect_he_offset) {
                                            segment_contour_points.push_back(SHAPE_POINT(0, Point_3(he_offset.point.x(), he_offset.point.y(), he_offset.oioa_altitude), true, he_offset.oioa_offset, he_offset.oioa_offset_index, he_offset.oioa_altitude ));
                                          }
                                        }
                                      }
                                      // Добавить точку контура в список
                                      segment_contour_points.push_back(SHAPE_POINT(0, Point_3(hds_h->vertex()->point().x(), hds_h->vertex()->point().y(), 0 /*Z не выставляется (это altitude), будет считаться позже*/), false /*altitude не рассчитана, будет считаться позже*/, 0, /*Это точка контура, поэтому offset-а у неё нет*/ -1, 0.0 ));
                                    }
                                    vect_segment_faces_edges_handle.back().push_back(hds_offset_h);
                                    face_slopes.push_back(hds_h->slope()); // Предположительно это указывает направление движения к или от центра SS.
                                    hds_h = hds_h->next();
                                  } while (hds_h != done);

                                  // Проверка, что наклонные slope() на противоположных (opposite()) рёбрах ожидаются противоположных знаков (кроме 0 - инверсия от 0 - всегда 0).
                                  // <image url="..\code_images\file_0061.png" scale=".1"/>
                                  if (verbose) {
                                    printf("\n " VAL2STR(Err::_0037) " ss_id=% 2u, ss_face_id % 2u, mesh_face_id % 2u, slopes: ", segment_face.ss_id, segment_face.ss_face_id, segment_face.mesh_face_id);
                                    for (auto& slope : face_slopes) {
                                      printf("%d ", slope);
                                    }
                                  }
                                }
                              }

                              timer_contour_points_001.stop();
                              contour_time+=timer_contour_points_001.time();
                              
                              if (map_ss_id_mesh_face_id_segment_points.find(target_ss_face_info.ss_id) == map_ss_id_mesh_face_id_segment_points.end()) {
                                map_ss_id_mesh_face_id_segment_points[target_ss_face_info.ss_id] = std::map<int /*mesh_face_id*/, std::vector<SHAPE_POINT>>();
                              }
                              map_ss_id_mesh_face_id_segment_points[target_ss_face_info.ss_id][target_ss_face_info.mesh_face_id] = segment_contour_points;

                              // Преобразовать загруженные halfedges в один контур (из первого контура удаляется первый halfedge; из второго контура удаляется первый edge,
                              // затем контур переворачивается 
                              std::vector<HDS_Halfedge_const_handle> segment_combined_he;
                              segment_combined_he.insert(segment_combined_he.end(), vect_segment_faces_edges_handle[0].begin() + 1, vect_segment_faces_edges_handle[0].end());
                              if (vect_segment_faces_edges_handle.size() > 0) {
                                // <image url="..\code_images\file_0060.png" scale=".3"/>
                                segment_combined_he.insert(segment_combined_he.end(), vect_segment_faces_edges_handle[0].rbegin(), vect_segment_faces_edges_handle[0].rend() - 1);
                              }

                              // Обработать в одном сегменте все offset-ы:
                              std::vector<std::vector<OFFSET_POINT>> vect_segment_offsets_faces;
                              if (verbose == true) {
                                printf("\n" VAL2STR(Err::_0038) " object_index %u, Offsets: [", object_index);
                                for (int I = 0; I <= (int)vect_ss_offset_intersection_points.size() - 2; I++) {
                                  auto& vect_ss_offset_intersection_points_I0 = vect_ss_offset_intersection_points[I + 0];
                                  auto& offset = vect_ss_offset_intersection_points_I0.oioa_offset;
                                  auto& offset_index = vect_ss_offset_intersection_points_I0.oioa_offset_index;
                                  printf("%u: %.5g; ", offset_index, CGAL::to_double(offset));
                                }
                                printf("];");
                              }
                              timer_bevel.start();
                              for (int I = 0; I <= (int)vect_ss_offset_intersection_points.size() - 2; I++) {
                                double t1_0 = timer_bevel.time();
                                // Примечание: 
                                // 1. между двумя offset может образоваться несколько полигонов на одном сегменте
                                // пример, когда между offset может образоваться больше одного полигона: <image url="..\code_images\file_0053.png" scale=".2"/>
                                // 2. требуется сразу рассчитывать индексы точек при построении полигонов на одном сегменте, чтобы потом не делать merge точек
                                // соседних полигонов.

                                auto& vect_ss_offset_intersection_points_I0 = vect_ss_offset_intersection_points[I + 0];
                                auto& vect_ss_offset_intersection_points_I1 = vect_ss_offset_intersection_points[I + 1];
                                
                                auto *offset0_info = &vect_ss_offset_intersection_points_I0;
                                auto *offset1_info = &vect_ss_offset_intersection_points_I1;

                                t1_0 = timer_bevel.time() - t1_0;
                                t1 += t1_0;

                                // Координаты точек offset0:
                                //printf("\n");
                                //{
                                //  std::vector<std::unordered_map<HDS_Halfedge_const_handle, Point_2>> vect_offsets_points;
                                //  vect_offsets_points.push_back(offset0_info.map_halfedge_to_offset_points);
                                //  vect_offsets_points.push_back(offset1_info.map_halfedge_to_offset_points);
                                //  int IJ = 0;
                                //  for (auto& offset_points : vect_offsets_points) {
                                //    printf("\noffset intersections %u: [", I+IJ);
                                //    for (auto& [handle, point] : offset_points) {
                                //      printf("[%.5g, %.5g, 0.0],", CGAL::to_double(point.x()), CGAL::to_double(point.y()));
                                //    }
                                //    printf("], ");
                                //    IJ++;
                                //  }
                                //}

                                // Для рассчёта offset важно, чтобы offset0 был меньше offset1, поэтому, если это не так, то временно поменять их местами.
                                // После получения полигонов между этими offset нужно вернуть их ориентацию в исходное состояние (перевернуть нормали, может быть
                                // просто инвертировать индексы вершин?)
                                bool reverse_result_offset = false; // если 
                                if (offset0_info->oioa_offset > offset1_info->oioa_offset) {
                                  std::swap(offset0_info, offset1_info);
                                  reverse_result_offset = true;
                                }

                                if (vect_segment_faces_handle.size() == 1) {
                                  // Оба offset в одном ss_id (либо в верхнем, либо в нижнем).
                                  // Примечание: тут нужно помнить, что нижний, отрицательный ss face получен из положительного SS, поэтому направление обхода у него должно совпадать, а вот последовательность
                                  // пересечения с offset-ами будет обратная! (А вот это надо проверить!!!)
                                  // <image url="..\code_images\file_0054.png" scale=".2"/>

                                  std::vector<OFFSET_POINT> vect_offset01_polygon_points;
                                  /// <summary>
                                  /// Как повлиять на добавление точек исходного face при обработке offset0 и offset1
                                  /// </summary>
                                  enum AFFECT_MODE {
                                    SKIP=0, APPEND,
                                  };

                                  // если ss - положительный, то рассчёт идёт по контуру от offset0 к offset1 (достаточно проверить первый offset)
                                  if (offset0_info->oioa_offset >= 0) {
                                  } else {
                                    // если ss - отрицательный, то рассчёт идёт по контуру от offset1 к offset0
                                    std::swap(offset0_info, offset1_info);
                                  }

                                  for (auto& face_handle : vect_segment_faces_handle) {
                                    HDS_Halfedge_const_handle hds_h = face_handle->halfedge(), done = hds_h;
                                    AFFECT_MODE affect_mode = SKIP;
                                    //printf("\n edges: %u :=========================", object_index);
                                    do {
                                      HDS_Halfedge_const_handle hds_off_h = hds_h;
                                      // если проверяемый halfedge не пересекается с offset
                                      if (hds_h->slope() == CGAL::NEGATIVE) // ensure same geometric point on both sides
                                        hds_off_h = hds_off_h->opposite();

                                      auto& offset0_he_intersect_iter = offset0_info->map__halfedge__offset_points.find(hds_off_h);
                                      auto& offset1_he_intersect_iter = offset1_info->map__halfedge__offset_points.find(hds_off_h);
                                      HDS_Vertex_const_handle hds_vert_current = hds_h->vertex();
                                      const Point_2& p2_he = hds_vert_current->point();
                                      const Point_2& he_point = Point_2(p2_he.x(), p2_he.y());
                                      //printf("[%.5g, %.5g],", CGAL::to_double(he_point.x()), CGAL::to_double(he_point.y()) );

                                      if (offset0_he_intersect_iter == offset0_info->map__halfedge__offset_points.end() && offset1_he_intersect_iter == offset1_info->map__halfedge__offset_points.end()) {
                                        if (affect_mode == SKIP) {
                                          // пропустить добавление точки halfedge (либо полигон ещё не открыт и перебор идёт ниже offset0, либо перебор halfedge идёт выше offset1
                                          // <image url="..\code_images\file_0059.png" scale=".4"/>
                                        }else
                                          if (affect_mode == APPEND) {
                                            // Запомнить точку halfedge (перебор halfedge идёт между offset0 и offset1 или между offset1 и offset0):
                                            // <image url="..\code_images\file_0058.png" scale=".4"/> (нет пересечений ни с offset0, но с offset1)
                                            vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(he_point.x(), he_point.y(), 0.0), false /*высота z точки ещё не определена*/));
                                          }
                                      } else 
                                        if (offset0_he_intersect_iter != offset0_info->map__halfedge__offset_points.end() && offset1_he_intersect_iter == offset1_info->map__halfedge__offset_points.end()) {
                                          auto& offset0_he_intersect_point = offset0_he_intersect_iter->second;
                                          if (affect_mode == SKIP) {
                                            // В случае, если полигон открылся только одним offset0 (пересечения в offset1 нет),
                                            // то добавить и точку пересечения offset0-halfedge и точку halfedge: <image url="..\code_images\file_0057.png" scale=".3"/>
                                            //vect_offset01_polygon_points - должен быть сейчас пустым, потому что тут начинается новый полигон
                                            vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset0_he_intersect_point.x(), offset0_he_intersect_point.y(), offset0_info->oioa_altitude), true));
                                            vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(he_point.x(), he_point.y(), 0.0), false /*высота z точки ещё не определена*/));
                                            affect_mode = APPEND;
                                          } else if (affect_mode == APPEND) {
                                            vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset0_he_intersect_point.x(), offset0_he_intersect_point.y(), offset0_info->oioa_altitude), true));
                                            // Расчётный полигон копируется в список полигонов и очищается для дальнейшего наполнения:
                                            vect_segment_offsets_faces.push_back(std::vector<OFFSET_POINT>(vect_offset01_polygon_points.begin(), vect_offset01_polygon_points.end())); // Скопировать полученный полигон в список полигонов
                                            vect_offset01_polygon_points.clear(); // Начать новый вектор точек полигона
                                            affect_mode = SKIP;
                                          }
                                        } else 
                                          if (offset0_he_intersect_iter == offset0_info->map__halfedge__offset_points.end() && offset1_he_intersect_iter != offset1_info->map__halfedge__offset_points.end()) {
                                            auto& offset1_he_intersect_point = offset1_he_intersect_iter->second;
                                            if (affect_mode == SKIP) {
                                              // Если при пересечении offset1 произошёл вход в полигон, без выхода в offset0, то нужно запомнить
                                              // и точку пересечения offset1 и точку halfedge. При этом разрешить добавление следующих точек he, если они встретятся до точку пересечения offset1
                                              // <image url="..\code_images\file_0055.png" scale=".4"/>
                                              vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset1_he_intersect_point.x(), offset1_he_intersect_point.y(), offset1_info->oioa_altitude), true));
                                              vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(he_point.x(), he_point.y(), 0.0), false /*высота z точки ещё не определена*/));
                                              affect_mode = APPEND;
                                            } else if (affect_mode == APPEND) {
                                              // Если при переходе через offset1 произошёл выход из полигона, то запомнить точку пересечения offset1 и halfedge:
                                              // то нужно запомнить только точку пересечения halfedge <image url="..\code_images\file_0056.png" scale=".4"/>
                                              vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset1_he_intersect_point.x(), offset1_he_intersect_point.y(), offset1_info->oioa_altitude), true));
                                              affect_mode = SKIP;
                                            }
                                          } else 
                                            if (offset0_he_intersect_iter != offset0_info->map__halfedge__offset_points.end() && offset1_he_intersect_iter != offset1_info->map__halfedge__offset_points.end()) {
                                              auto& offset0_he_intersect_point = offset0_he_intersect_iter->second;
                                              auto& offset1_he_intersect_point = offset1_he_intersect_iter->second;
                                              if (affect_mode == SKIP) {
                                                if (vect_offset01_polygon_points.size() == 0) {
                                                  // Если полигон между offset-ами ещё пустой, то добавить в него в начало сразу обе точки в последовательности offset0,offset1 и оставить режим affect_mode=SKIP,
                                                  // потому что промежуточные точки halfedge не нужно добавлять при перечислении точек halfedge выше offset1
                                                  vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset0_he_intersect_point.x(), offset0_he_intersect_point.y(), offset0_info->oioa_altitude), true));
                                                  vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset1_he_intersect_point.x(), offset1_he_intersect_point.y(), offset1_info->oioa_altitude), true));
                                                } else {
                                                  // Если halfedge пересекается с двумя точками offset1 и offset0, значит это пересечение закрывает этот полигон. Добавить точки пересечения offset1 и offset0,
                                                  // скопировать точки полигона в список полигонов и очистить буфер точек
                                                  // affect_mode=SKIP, чтобы если начнётся новый полигон, то начать новый вектор для его точек
                                                  vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset1_he_intersect_point.x(), offset1_he_intersect_point.y(), offset1_info->oioa_altitude), true));
                                                  vect_offset01_polygon_points.push_back(OFFSET_POINT(Point_3(offset0_he_intersect_point.x(), offset0_he_intersect_point.y(), offset0_info->oioa_altitude), true));
                                                  vect_segment_offsets_faces.push_back(std::vector<OFFSET_POINT>(vect_offset01_polygon_points.begin(), vect_offset01_polygon_points.end())); // Скопировать полученный полигон в список полигонов
                                                  vect_offset01_polygon_points.clear(); // Начать новый вектор точек полигона
                                                }
                                                affect_mode = SKIP; // affect_mode не меняется
                                              } else {
                                                // Кроме как в mode SKIP две точки offset  не могут встретиться на одном halfedge.
                                              }
                                            }
                                      hds_h = hds_h->next();
                                    } while (hds_h != done);
                                  }

                                } else /*size==2*/ {
                                  // offset-ы в разных id
                                }
                              }
                              timer_bevel.stop();
                              //if (verbose == true) {
                              //  // Вывести точки faces:
                              //  printf("\n " VAL2STR(Err::_0039) " %u) faces: %u :", log_counter1++, object_index);
                              //  for (auto& face : vect_segment_offsets_faces) {
                              //    printf("\n[");
                              //    for (auto& vert : face) {
                              //      printf("[%.3g, %.3g, %.3g],", CGAL::to_double(vert.point.x()), CGAL::to_double(vert.point.y()), CGAL::to_double(vert.point.z()));
                              //    }
                              //    printf("],");
                              //  }
                              //}
                              vect_ss_offset_faces.push_back(vect_segment_offsets_faces);
                            }
                            printf("\n " VAL2STR(Err::_0041) ". SS 2D Offset. ss_id=% 2d, timer_contour_points_001 %.5g", ss_id, contour_time);
                          }

                          if (verbose == true) {
                            printf("\n " VAL2STR(Err::_0040) " ====================== GENERAL LIST OF POINTS OF object_index %u: =========================\n", object_index);
                            std::string str1;
                            
                            for (auto& vect_offset_faces : vect_ss_offset_faces) {
                              for (auto& face : vect_offset_faces) {
                                printf("[");
                                for (auto& vert : face) {
                                  printf("[%.3g, %.3g, %.3g],", CGAL::to_double(vert.point.x()), CGAL::to_double(vert.point.y()), CGAL::to_double(vert.point.z()));
                                }
                                printf("],");
                              }
                            }
                            printf("\n ====================== /GENERAL LIST OF POINTS object_index %u: =========================\n", object_index);
                          }
                        }
                        //timer_bevel.stop();
                        if (verbose) {
                          printf("\n " VAL2STR(Err::_0033) ". SS 2D Offset. object_index %u, ss timer_bevel %.5g, %.5g", object_index, timer_bevel.time(), t1);
                        }
                      }
                    }
                  }
                  timer1.stop();
                  if (verbose) {
                    printf("\n " VAL2STR(Err::_0031) ". SS 2D Offset. ss faces ajacent timer %.5g", timer1.time());
                  }

                  int split_object_index = 0; // для join_mode==SPLIT индекс нового объекта назначается группе контуров с одним index_offset (а не каждому контуру, чтобы количество объектов в режиме Edges и Faces совпадало)
                  int  keep_object_index = 0; // Для join_mode==KEEP  индексы результирующих объектов не меняются
                  int merge_object_index = 0; // Для join_mode==MERGE индекс результирующего объекта один - 0, т.к. объект будет единственным.
                  Mesh mesh_merged_all_objects; // При join_mode==MERGE выполнить join всех mesh

                  // Сгруппировать или разгруппировать объекты по параметру results_join_mode
                  // Примечание: все offset тут теряются
                  for (auto& [object_index, ss_id_OIOA_OFFSET_SS_PARAMS] : map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS) {
                    if (results_join_mode == 0) {
                      std::vector<Mesh> components;
                      // Разделить текущий набор точек по ss_id. Напоминаю, что ss_id сквозные.
                      for (auto& [ss_id, oioa_offset_ss_params] : ss_id_OIOA_OFFSET_SS_PARAMS) {
                        if (map_join_mesh.find(split_object_index) == map_join_mesh.end()) {
                          map_join_mesh[split_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                        }
                        //vect__oioa__offset_ss_params[0].mesh_ss_02_merged = vect__oioa__offset_ss_params[0].mesh_ss_01_source;

                        OIOA_OFFSET_SS_PARAMS offset_info(OIOA(merge_object_index, 0, 0.0, 0.0), std::vector<Polygon_with_holes_2>(), std::vector<SS__HalfEdge__Point_2>());
                        offset_info.mesh_ss_02_merged = oioa_offset_ss_params[0].mesh_ss_01_source;
                        map_join_mesh[split_object_index].push_back(offset_info);
                        split_object_index++;
                      }
                    } else if (results_join_mode == 1) {
                      for (auto& [ss_id, oioa_offset_ss_params] : ss_id_OIOA_OFFSET_SS_PARAMS) {
                        keep_object_index = object_index;
                        if (map_join_mesh.find(keep_object_index) == map_join_mesh.end()) {
                          map_join_mesh[keep_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                          // Запомнить только общий merged mesh в самом первом параметре oioa (выход с помощью break).
                          //vect__oioa__offset_ss_params[0].mesh_ss_02_merged = map_object_index_result_ss_Mesh[object_index];
                          //map_join_mesh[keep_object_index].push_back(vect__oioa__offset_ss_params[0]);
                          OIOA_OFFSET_SS_PARAMS offset_info(OIOA(merge_object_index, 0, 0.0, 0.0), std::vector<Polygon_with_holes_2>(), std::vector<SS__HalfEdge__Point_2>());
                          offset_info.mesh_ss_02_merged = map_object_index_mesh_ss_merged[object_index];
                          map_join_mesh[keep_object_index].push_back(offset_info);
                          break;
                        }
                        // Напоминание: offset пока не учитывается, поэтому только в первом элементе хранятся полигоны SS.
                      }
                    } else if (results_join_mode == 2) { // 2 - для всех остальных случаев
                      // Этот метод работает немного иначе, чем в режиме Straight Skeleton: объединение Mesh происходит сразу здесь,
                      // чтобы потом не пересчитывать индексы.
                      // Объеденить все результаты в один объект
                      for (auto& [ss_id, oioa_offset_ss_params] : ss_id_OIOA_OFFSET_SS_PARAMS) {
                        if (map_join_mesh.find(merge_object_index) == map_join_mesh.end()) {
                          map_join_mesh[merge_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                        }
                        // TODO: Нужно сделать Merge по вершинам. Сейчас его нет.
                        mesh_merged_all_objects.join(oioa_offset_ss_params[0].mesh_ss_01_source);
                        // Напоминание: offset пока не учитывается.
                        //map_join_mesh[merge_object_index].push_back(vect__oioa__offset_ss_params[0]);
                      }
                    }
                  }
                  if (results_join_mode == 2) {
                    OIOA_OFFSET_SS_PARAMS offset_info(OIOA(merge_object_index, 0, 0.0, 0.0), std::vector<Polygon_with_holes_2>(), std::vector<SS__HalfEdge__Point_2>());
                    offset_info.mesh_ss_02_merged = mesh_merged_all_objects;
                    map_join_mesh[merge_object_index].push_back(offset_info);
                  }

                } // else if (result_type == ResType::STRAIGHT_SKELETON)
              } // Конец рассчёта map_join_mesh. В настоящий момент известен количественный состав результата.

              mesh_data->nn_objects = map_join_mesh.size();
              if (map_join_mesh.size() == 0) {
              }else{
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

                // Выполнить обработку контуров. Время на триангуляцию тратится незначительное. Можно пока не делать распараллеливание.
                // <image url="..\code_images\file_0023.png" scale="1.0"/>
                for (auto& KV : map_join_mesh) {
                  int object_index = KV.first;
                  auto& vect_pwhs_with_offset_index = KV.second;

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
                    // Преобразовать вектор PWH полигонов в соответствующий вектор mesh
                    for (auto& pwh_altitude : vect_pwhs_with_offset_index) {
                      for (auto& pwh1 : pwh_altitude.vect_pwh) {
                        Mesh sm1;
                        ConvertPolygonWithHoles2IntoMesh(pwh1, pwh_altitude.oioa1.altitude, sm1);
                        vect_sm_offset_index_order.push_back(sm1);
                      }
                    }
                    // Из вектора Mesh-ей надо сформировать полигональную сетку
                    ConvertMeshesIntoVerticesEdgesFaces(vect_sm_offset_index_order, vect_res_object1_verts, vect_res_object1_edges, vect_res_object1_faces);
                  } else if (result_type == ResType::BEVEL) {
                  } else if (result_type == ResType::STRAIGHT_SKELETON) {
                    // TODO: На этом уровне данные от Extrude должны быть подготовлены. Тут их нужно получить только "извлечь" для отправки.
                    std::vector<Mesh> vect_sm_offset_index_order;
                    bool invert_faces = false;
                    // Триангулировать полигоны всех faces:
                    Mesh sm1; // About Mesh: https://doc.cgal.org/latest/Surface_mesh/index.html
                    
                    unsigned int start_shape_vert_index = 0;
                    for (auto& pwh_altitude : vect_pwhs_with_offset_index) {
                      
                      // Get vertices coords
                      for (auto& vertex : pwh_altitude.mesh_ss_02_merged.points()) {
                        vect_res_object1_verts.push_back(std::vector<float>{ (float)CGAL::to_double(vertex.x()), (float)CGAL::to_double(vertex.y()), (float)CGAL::to_double(vertex.z()) });
                      }
                      
                      // get edges vertex indexes
                      std::set<std::tuple<unsigned int, unsigned int>> set_tuple_edges;
                      for (auto& ei : pwh_altitude.mesh_ss_02_merged.edges()) {
                        auto& vi0 = pwh_altitude.mesh_ss_02_merged.vertex(ei, 0);
                        auto& vi1 = pwh_altitude.mesh_ss_02_merged.vertex(ei, 1);
                        //set_tuple_edges.insert({ vi0, vi1 });
                        vect_res_object1_edges.push_back({ vi0, vi1 });
                      }
                      //for (auto& [first_index, second_index] : set_tuple_edges) {
                      //  vect_res_object1_edges.push_back({ first_index, second_index });
                      //}

                      // get faces vertex indexes:
                      for (auto& fi : pwh_altitude.mesh_ss_02_merged.faces() ) {
                        std::vector<unsigned int> vect_face_vertex_id;
                        auto& hf = pwh_altitude.mesh_ss_02_merged.halfedge(fi);
                        for (auto& hi : halfedges_around_face(hf, pwh_altitude.mesh_ss_02_merged)) {
                          auto& vi = target(hi, pwh_altitude.mesh_ss_02_merged);
                          vect_face_vertex_id.push_back( vi.idx() );
                        }
                        vect_res_object1_faces.push_back(vect_face_vertex_id);
                      }

                      //for (auto& points : pwh_altitude.vect_vect_ss_face_points) {
                      //  if (points.size() >= 3) {
                      //    for (auto& p1 : points) {
                      //    }

                      //    vect_res_object1_verts.push_back(std::vector<float>{ (float)CGAL::to_double(points[0].x()), (float)CGAL::to_double(points[0].y()), (float)CGAL::to_double(points[0].z()) });
                      //    auto& face_indexes = std::vector<unsigned int>();
                      //    face_indexes.push_back(start_shape_vert_index+0);
                      //    for (int I = 1; I <= points.size() - 1; I++) {
                      //      auto& points_I = points[I];
                      //      vect_res_object1_verts.push_back(std::vector<float>{ (float)CGAL::to_double(points_I.x()), (float)CGAL::to_double(points_I.y()), (float)CGAL::to_double(points_I.z()) });
                      //      vect_res_object1_edges.push_back({ start_shape_vert_index + (I-1), start_shape_vert_index + (I-0) });
                      //      face_indexes.push_back(start_shape_vert_index + I);
                      //    }
                      //    vect_res_object1_edges.push_back({ (unsigned int)(start_shape_vert_index + points.size()-1), start_shape_vert_index + 0});
                      //    vect_res_object1_faces.push_back(face_indexes);
                      //    start_shape_vert_index += face_indexes.size();
                      //    //sm1.add_face(vect_idx);
                      //    // Триангулировать полигоны, полученные при расчёте соединений между offset-ами:
                      //    //triangulate_skeleton_face(points, invert_faces, faces_points, faces_indexes, CGAL::parameters::verbose(verbose));
                      //  }
                      //}

                      //if (sm1.num_vertices() > 0) {
                      //  namespace PMP = ::CGAL::Polygon_mesh_processing;
                      //  //PMP::merge_duplicate_points_in_polygon_soup(faces_points, faces_indexes);
                      //  //if (!PMP::is_polygon_soup_a_polygon_mesh(faces_indexes))
                      //  //  PMP::orient_polygon_soup(faces_points, faces_indexes);
                      //  //CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces_indexes));
                      //  //Mesh sm1;
                      //  //PMP::polygon_soup_to_polygon_mesh(faces_points, faces_indexes, sm1);
                      //  vect_sm_offset_index_order.push_back(sm1);
                      //}
                    }
                    // Из вектора Mesh-ей надо сформировать полигональную сетку
                    //if (vect_sm_offset_index_order.size() > 0) {
                    //  ConvertMeshesIntoVerticesEdgesFaces(vect_sm_offset_index_order, vect_res_object1_verts, vect_res_object1_edges, vect_res_object1_faces);
                    //}
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
                  // Список количества индексов в каждой face добавить последовательно в конец
                  for (int I = 0; I <= (int)vect_res_object1_faces.size() - 1; I++) {
                    vect_objects_faces_indexes_counters.push_back(vect_res_object1_faces[I].size());
                  }
                  vect_objects_faces.insert(vect_objects_faces.end(), vect_res_object1_faces.begin(), vect_res_object1_faces.end());
                }

                // записать результаты vertices, edges и faces
                if (vect_objects_verts.size() == 0) {
                  mesh_data->nn_offsets_indexes = NULL;
                  mesh_data->vertices = NULL;
                  mesh_data->edges = NULL;
                  mesh_data->faces = NULL;
                  mesh_data->nn_faces_indexes_counters = NULL;
                } else {
                  mesh_data->nn_offsets_indexes = NULL;
                  mesh_data->vertices = NULL;
                  mesh_data->edges = NULL;
                  mesh_data->faces = NULL;
                  mesh_data->nn_faces_indexes_counters = NULL;

                  if (vect_objects_offsets.size() > 0) {
                    mesh_data->nn_offsets_indexes = (int*)malloc(sizeof(int) * vect_objects_offsets.size());
                  }

                  // Загрузка данных от вершинах
                  if (vect_objects_verts.size() > 0) {
                    mesh_data->vertices = (float*)malloc(sizeof(float[3]) * vect_objects_verts.size());
                  }

                  // Загрузка данных о рёбрах
                  if (vect_objects_edges.size() > 0) {
                    mesh_data->edges = (int*)malloc(sizeof(int[2]) * vect_objects_edges.size());
                  }

                  // Инициализировать счётчики индексов для faces (ещё не сами индексы для faces):
                  if (vect_objects_faces.size() > 0) {
                    // Количество счётчиков индексов faces для каждой faces равно количеству faces:
                    mesh_data->nn_faces_indexes_counters = (int*)malloc(sizeof(int) * vect_objects_faces.size());
                    // Посчитать количество индексов для faces (faces имеют разную длину)
                    size_t faces_size = 0;
                    for (int I = 0; I <= (int)vect_objects_faces.size() - 1; I++) {
                      faces_size += vect_objects_faces[I].size();
                    }
                    mesh_data->faces = (int*)malloc(sizeof(int) * faces_size );
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
                  int face_idx = 0;
                  for (int I = 0; I <= (int)vect_objects_faces.size() - 1; I++) {
                    auto& vect_objects_faces_I = vect_objects_faces[I];
                    mesh_data->nn_faces_indexes_counters[I] = vect_objects_faces_I.size();
                    for (int IJ = 0; IJ <= vect_objects_faces_I.size() - 1; IJ++) {
                      mesh_data->faces[face_idx] = vect_objects_faces_I[IJ];
                      face_idx++;
                      //mesh_data->faces[I * 3 + 0] = vect_objects_faces[I][0];
                      //mesh_data->faces[I * 3 + 1] = vect_objects_faces[I][1];
                      //mesh_data->faces[I * 3 + 2] = vect_objects_faces[I][2];
                    }
                  }
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
          CGAL::Aff_transformation_3<K> ac3(
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
                        // Проверить, что все точки p_hole1 находятся внутри boundary. Есть исключение, когда точки находятся внутри, но сами контуры пересекаются, но пока так, подумать над решением:
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

                        // Проверить, что он не пересекается с другими holes
                        bool is_objects_intersects = false;
                        try {
                          for (auto& ap1 : allowed_contours) {
                            Polygon_2 pap1(ap1.c1->begin(), ap1.c1->end());
                            pap1.reverse_orientation();
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
            std::vector<std::vector<FT>> vect_vertices_of_errors;
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
                    vect_vertices_of_errors.push_back(std::vector<FT>{v.x(), v.y()});
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
                mesh_data->vertices_of_errors[I * 3 + 0] = (float)CGAL::to_double(vect_vertices_of_errors[I][0]);
                mesh_data->vertices_of_errors[I * 3 + 1] = (float)CGAL::to_double(vect_vertices_of_errors[I][1]);
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
              std::vector<unsigned int>              vect_objects_faces_indexes_counters; // Счётчики количества индексов каждой face
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
              for (auto& KV : map_join_mesh) {
                mesh_data->nn_objects_indexes[result_index] = KV.first;

                //std::set   <                     int>  set_object1_offsets; // offsets (per object1)
                std::vector<std::vector<       float>> vect_object1_verts; // Координаты vertices (per object1)
                std::vector<std::vector<unsigned int>> vect_object1_edges; // indexes of edges (per object1)
                std::vector<std::vector<unsigned int>> vect_object1_faces; // indexes of faces (per object1)

                CGAL::Real_timer timer2;
                timer2.start();

                ConvertMeshesIntoVerticesEdgesFaces(KV.second, vect_object1_verts, vect_object1_edges, vect_object1_faces);
                timer2.stop();
                if (verbose) {
                  printf("\nObject %i: verts - %zu, edges - %zu, faces - %zu; triangulated time - %.5g sec.", KV.first, vect_object1_verts.size(), vect_object1_edges.size(), vect_object1_faces.size(), timer2.time());
                }
                mesh_data->nn_verts[result_index] = vect_object1_verts.size();
                mesh_data->nn_edges[result_index] = vect_object1_edges.size();
                mesh_data->nn_faces[result_index] = vect_object1_faces.size();

                vect_objects_count_of_verts.push_back(vect_object1_verts.size());
                vect_objects_count_of_edges.push_back(vect_object1_edges.size());
                vect_objects_count_of_faces.push_back(vect_object1_faces.size());

                vect_objects_verts.insert(vect_objects_verts.end(), vect_object1_verts.begin(), vect_object1_verts.end());
                vect_objects_edges.insert(vect_objects_edges.end(), vect_object1_edges.begin(), vect_object1_edges.end());
                // Список количества индексов в каждой face добавить последовательно в конец
                for (int I = 0; I <= (int)vect_object1_faces.size() - 1; I++) {
                  vect_objects_faces_indexes_counters.push_back(vect_object1_faces[I].size());
                }
                vect_objects_faces.insert(vect_objects_faces.end(), vect_object1_faces.begin(), vect_object1_faces.end());

                result_index++;
              }

              // записать результаты vertices, edges и faces
              if (vect_objects_verts.size() == 0) {
                mesh_data->nn_offsets_indexes = NULL;
                mesh_data->vertices = NULL;
                mesh_data->edges = NULL;
                mesh_data->faces = NULL;
                mesh_data->nn_faces_indexes_counters = NULL;
              } else {
                mesh_data->nn_offsets_indexes = NULL;
                mesh_data->vertices = NULL;
                mesh_data->edges = NULL;
                mesh_data->faces = NULL;
                mesh_data->nn_faces_indexes_counters = NULL;

                mesh_data->vertices = (float*)malloc(sizeof(float[3]) * vect_objects_verts.size());
                mesh_data->edges = (int*)malloc(sizeof(int[2]) * vect_objects_edges.size());
                if (vect_objects_faces.size() > 0) {
                  // Количество счётчиков индексов faces для каждой faces равно количеству faces:
                  mesh_data->nn_faces_indexes_counters = (int*)malloc(sizeof(int) * vect_objects_faces.size());
                  // Посчитать количество индексов для faces (faces имеют разную длину, но в extrude они должны быть все триангулированы)
                  size_t faces_size = 0;
                  for (int I = 0; I <= (int)vect_objects_faces.size() - 1; I++) {
                    faces_size += vect_objects_faces[I].size();
                  }
                  mesh_data->faces = (int*)malloc(sizeof(int) * faces_size);
                  //mesh_data->faces = (int*)malloc(sizeof(int[3]) * vect_objects_faces.size());
                }

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
                //for (int I = 0; I <= vect_objects_faces.size() - 1; I++) {
                //  // Normals orientation <image url="..\code_images\file_0022.png" scale=".3"/>
                //  mesh_data->faces[I * 3 + 0] = vect_objects_faces[I][0];
                //  mesh_data->faces[I * 3 + 1] = vect_objects_faces[I][2];
                //  mesh_data->faces[I * 3 + 2] = vect_objects_faces[I][1]; // Не 0,1,2, 0,
                //}
                
                // Загрузка данных по индексам faces
                int face_idx = 0;
                for (int I = 0; I <= (int)vect_objects_faces.size() - 1; I++) {
                  auto& vect_objects_faces_I = vect_objects_faces[I];
                  mesh_data->nn_faces_indexes_counters[I] = vect_objects_faces_I.size();
                  for (int IJ = 0; IJ <= vect_objects_faces_I.size() - 1; IJ++) {
                    mesh_data->faces[face_idx] = vect_objects_faces_I[IJ];
                    face_idx++;
                  }
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
          if (md->edges != NULL) {
            free(md->edges);
            md->edges = NULL;
          }
          if (md->nn_faces_indexes_counters != NULL) {
            free(md->nn_faces_indexes_counters);
            md->nn_faces_indexes_counters = NULL;
          }
          if (md->faces != NULL) {
            free(md->faces);
            md->faces = NULL;
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