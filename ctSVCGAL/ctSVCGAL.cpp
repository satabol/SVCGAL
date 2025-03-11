/*
# This file is part of project Sverchok. It's copyrighted by the contributors
# recorded in the version control history of the file, available from
# its original location https://github.com/nortikin/sverchok/commit/master
#
# SPDX-License-Identifier: GPL3
# License-Filename: LICENSE
# Author of this file is satabol@yandex.ru. You can find me on the Telegram messenger channel https://t.me/sverchok_3d as 'Satabol'
# Documentation: https://nortikin.github.io/sverchok/docs/nodes/CAD/straight_skeleton_2d_offset.html
# This code has many images inline. To see them in the Visual Studio 2022 IDE use Visual Studio Addon "ImageComments2022": https://marketplace.visualstudio.com/items?itemName=ImageComments2022.ImageComments2022
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

//#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Plane_3.h>

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
#include <boost/crc.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <ctime>
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

      enum Msg {
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
        _0042,
        _0043,
        _0044,
        _0045,
        _0046,
        _0047,
        _0048,
        _0049,
        _0050,
        _0051,
        _0052,
        _0053,
        _0054,
        _0055,
        _0056,
        _0057,
        _0058,
        _0059,
        _0060,
        _0061,
        _0062,
        _0063,
        _0064,
        _0065,
        _0066,
        _0067,
        _0068,
        _0069,
        _0070,
        _0071,
        _0072,
        _0073,
        _0074,
        _0075,
        _0076,
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
        int* nn_objects_original_indexes = NULL; // предыдущие индексы объектов из тех (индексы исходных объектов)
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
          int object_original_index;

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

          OIOA(int _object_index, int _object_original_index, int _offset_index, FT _offset, FT _altitude)
            : object_index(_object_index), object_original_index(_object_original_index), offset_index(_offset_index), offset(_offset), neg_type(NEGATIVE_TYPE::INVERTED_HOLE), altitude(_altitude){
          }

          OIOA(int _object_index, int _object_original_index, int _offset_index, FT _offset, NEGATIVE_TYPE _neg_type, FT _altitude)
            :object_index(_object_index), object_original_index(_object_original_index), offset_index(_offset_index), offset(_offset), neg_type(_neg_type), altitude(_altitude) {
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

          OIOA_OFFSET_SS_PARAMS(OIOA _oioa1, std::vector<Polygon_with_holes_2> _vect_pwh, std::vector<SS__HalfEdge__Point_2> _vect__SS__HalfEdge__Point_2)
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
          int object_original_index = -1;

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


          POMS(int _object_index, int _object_original_index, std::shared_ptr<Polygon_with_holes_2> _polygon1_with_holes, std::vector<OIOA> _vect_oioa, SS_TYPE _ss_type)
            :polygon1_with_holes(_polygon1_with_holes), object_original_index(_object_original_index), vect_oioa(_vect_oioa), ss_type(_ss_type) {
            object_index = _object_index;
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
          std::vector/*point indexes*/<int> beveled_verts_indexes;

          /// <summary>
          /// Первый отступ
          /// </summary>
          //OIOA oioa0;
          int oioa0_offset_index;
          FT  oioa0_offset;
          FT  oioa0_altitude;

          /// <summary>
          /// Второй отступ
          /// </summary>
          //OIOA oioa1;
          int oioa1_offset_index;
          FT  oioa1_offset;
          FT  oioa1_altitude;

          /// <summary>
          /// Выполнять ли reverse. (По хорошему этот параметр должен быть где-то снаружи, т.к. это однократная операция)
          /// </summary>
          bool do_reverse;

          //FACE_INFO(OIOA& _oioa0, OIOA& _oioa1, bool& _do_reverse)
          //: oioa0(_oioa0), oioa1(_oioa1), do_reverse(_do_reverse){
          //  
          //}

          FACE_INFO(int _oioa0_offset_index, FT _oioa0_offset, FT _oioa0_altitude, int _oioa1_offset_index, FT _oioa1_offset, FT _oioa1_altitude, bool& _do_reverse)
            : oioa0_offset_index(_oioa0_offset_index), oioa0_offset(_oioa0_offset), oioa0_altitude(_oioa0_altitude), oioa1_offset_index(_oioa1_offset_index), oioa1_offset(_oioa1_offset), oioa1_altitude(_oioa1_altitude), do_reverse(_do_reverse) {

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

        static void FT_to_bytes(const FT& value, unsigned char* bytes) {
          std::memcpy(bytes, &value, sizeof(FT));
        }

        DLLEXPORT int fill_vertices(float src[][3], float* target, int size) {
          for (int I = 0; I <= size - 1; I++) {
            float f0 = src[I][0];
            float f1 = src[I][1];
            float f2 = src[I][2];
            target[I*3+0] = f0;
            target[I*3+1] = f1;
            target[I*3+2] = f2;
          }
          return 0;
        }
        DLLEXPORT int fill_edges(int* src[2], int* target[2], int size) {
          for (int I = 0; I <= size - 1; I++) {
            target[I][0] = src[I][0];
            target[I][1] = src[I][1];
          }
          return 0;
        }
        DLLEXPORT int fill_faces(int* src, int* target, int size) {
          for (int I = 0; I <= size - 1; I++) {
            target[I] = src[I];
          }
          return 0;
        }

        /// <summary>
        /// Кэширование SS на основе crc координат, по которым строится SS. <image url="..\code_images\file_0081.png" scale=".3"/>
        /// </summary>
        std::map<unsigned int /*crc - контрольная сумма координат исходных данных для построения skeleton*/, std::shared_ptr<Ss> /*ссылка на skeleton по crc*/> map__crc__ss;

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
          // К чему относятся входные данные. <image url="..\code_images\file_0012.png" scale=".3"/>
          // Контуры объектов считаются независимо друт от друга
          // Контуры планов могут влиять друг на друга при отрицательных offset-ах. Контуры с положительными offset будут рассчитываться независимо.
          // Также при отрицательных offset инвертируются holes и рассчитываются как положительные offset. Требуется избегать расчётов Straight Skeleton-ов
          // при вложенных контурах, т.к. это приводит к зависанию расчёта (есть такой глюк). Контуры в Straight Skeleton добавляются независимо,
          // поэтому можно случайно добавить контуры внутренних объектов (что с типами Polygon_with_holes_2 не происходит, т.к. у них есть деление на outer_boundary
          // и holes.
          int in_count_of_objects,        // Количество независимых исходных объектов. Shape mode: <image url="..\code_images\file_0020.png" scale=".1"/>
           int in_shapes_mode[],          // Тип обработки контуров объекта: 0-FULL_MODE (полная обработка), 1-EXCLUDE_HOLES (только внешние контуры, отверстия не учитываются), 2-INVERT_HOLES (отключить внешние контуры)
           int in_count_of_offsets[],     // Количества офсетов, которые надо посчитать (размер массива соответствует in_count_of_objects)
         float in_offsets[],              // Массив value offsets, которые надо посчитать
           int in_count_of_altitudes[],   // Количества высот, которые надо учитывать при расстановке результатов. (размер массива соответствует in_count_of_objects)
         float in_altitudes[],            // Массив value offsets, которые надо посчитать
           int in__profile_faces__counters[],           // Массив количества faces в профиле в каждом объекте
           int in__profile_faces__indexes_counters[],   // Количество индексов на каждый face в профиле
           int in__profile_faces__indexes[],            // Индексы для faces
           int in__profile_faces__close_modes_array[],  // Соответствие каждому профильному face свойства открыт он или закрыт 0 - открыт, 1 - закрыт, 2 - по парам как указано в списке (например: список 1,2,2,3,3,4,6,7,7,8 => преобразуется в: [1,2],[2,3],[3,4],[6,7],[7,8]), . Размер массива соответствует размеру in__profile_faces__indexes_counters (не сумме значений, а только длине) <image url="..\code_images\file_0089.png" scale=".5"/>
          //int in_count_of_offset_edges[],// 
          //int in_offset_edges[][2],      // 
           int in_count_of_planes[],      // Количество планов в объектах (размер массива соответствует in_count_of_objects)
           int in_contours_in_planes[],   // Количество контуров в планах. Первый контур плана всегда outer, остальные - holes. Вложенность планов одного объекта в упаковке не учитывается.
           int in_vertices_in_contours[], // Количество вершин в контуре.
         float in_vertices[][3],          // вершины контуров (являются одной цепочкой, которую надо поделить на части из параметра lens_of_contours
          bool only_tests_for_valid,      // Проверка всех контуров на валидность без расчёта. (true - только проверить на валидность, false - проверить на валидность и посчитать)
           int in_res_type,               // Тип результата. 0 - контуры, 1 - mesh
          bool force_z_zero,              // Обязательное преобразование Z=0 при любых обстоятельствах (иногда при редактировании модели можно двигать курсор в плоскости XY b получать Z "слегка" отличный от 0).
           int source_objects_join_mode,    // Предобработка исходных Meshes: # 0 - split - разделить на отдельные независимые-объекты, 1 - keep - оставить как есть, 2 - merge all meshes - объеденить в один большой объект.<image url="..\code_images\file_0030.png" scale=".3"/>
           int results_join_mode,           // Что сделать с полученными mesh: # 0 - split - разделить на отдельные mesh-объекты, 1 - keep - оставить в своих объектах, 2 - merge all meshes - объеденить в один большой объект. В этом ноде дополнительные параметры к объектам не нужны
          bool verbose,                     // Подробный вывод
          bool use_cache_of_straight_skeleton, // Применять кэш при рассчёте Offsets of Straight Skeleton
          bool bevel_more_split // Вдруг захочется нестандартного разбиения результата на более мелкие части в режиме beveled (разбиение по ss_id и profile_face_index)
          // <image url="..\code_images\file_0086.png" scale=".1"/>
        ) {

          if (verbose == true) {
#ifdef _DEBUG
            printf("\n " VAL2STR(Msg::_0073) ". straight_skeleton_2d_offset in DEBUG MODE");
#endif
            printf("\n " VAL2STR(Msg::_0057) ". straight_skeleton_2d_offset_internal. verbose_mode is on");
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

          struct OFFSET_EDGE {
            int start_index;
            int end_index;
            OFFSET_EDGE(int _start_index, int _end_index) 
            : start_index(_start_index), end_index(_end_index) {

            }
          };

          // Параметры одного исходного полигонального объекта
          struct ClsObject {
            // Исходный индекс объекта
            int object_index;
            // Исходный индекс объекта, на случай, если была применена операция split
            int object_original_index;
            // Список offsets и altitudes
            std::vector<OIOA> vect_oioa;
            std::vector<OFFSET_EDGE> vect_offset_edges;
            // Контуры исходных планов
            std::vector<ContourSequence> source_planes;
            // Полигоны для расчётов (позиции полигонов могут отличаться от позиций контуров в source_planes)
            std::vector<std::shared_ptr<Polygon_with_holes_2>> vect_polygons_with_holes;
            
            SS_TYPE ss_type;

            ClsObject(int _object_index, int _object_original_index, std::vector<OIOA> _vect_ooai, std::vector<ContourSequence> _source_planes, SS_TYPE _ss_type)
              :object_index(_object_index), object_original_index(_object_original_index), source_planes(_source_planes), vect_oioa(_vect_ooai), ss_type(_ss_type) {
            }
            //ClsObject(int _object_index, int _object_original_index, std::vector<OIOA> _vect_ooai, std::vector<ContourSequence> _source_planes, SS_TYPE _ss_type)
            //  :object_index(_object_index), object_original_index(_object_original_index), source_planes(_source_planes), vect_oioa(_vect_ooai), ss_type(_ss_type) {
            //}
          };

          std::vector<OIOA> res_contours; // Результаты расчёта всех контуров. Результаты будут сортироваться с учётом вложенности контуров по offset_index, в которой offset подавались на вход (по объектам) и с учётом параметра join_mode.
          std::vector<ObjectError> res_errors; // Собирать ошибки в объектах при загрузке исходных данных и когда попытка использовать контур привела к его исключению из расчёта. Относится только к исходным объектам.

          int v_pos = 0;
          int general_count_of_vertices = 0;
          std::map<int, std::vector<Point_2>> contours;
          std::vector<ClsObject> objects; // Информация по объектам будет нужна и после расчёта (нужно помнить начальные индексы входных объектов, потому что иногда объект полностью выпадает из расчёта, а нужно передать инфу, что в объекте ничего нет и его ошибки
          //std::map<int, std::vector<OFFSET_EDGE>> map__object_index__offset_edges;
          std::map<int, std::vector/*faces*/<std::vector<int /*индексы offset для соединения*/>>> map__object_index__profile_faces;
          std::map<int, std::vector/*faces*/<int>> map__object_index__profile_faces__open_modes;
          
          // Вектор параметров для рассчёта SS.
          std::vector<POMS> vect_polygon1_oioa;

          //// Карта уникальных идентификаторов исходных edges, используемых для построения SS (чтобы потом найти смежные faces для примыкающих друг к другу SS.
          //IdType unique_is_edges = 1;
          //std::map<typename Polygon_2::Edge_const_iterator, IdType> edge_id_map;

          //try 
          {
            {
              CGAL::Real_timer timer1;
              timer1.start();
              // Загрузка исходных данных по объектам с некоторыми проверками.
              int plane_cursor = 0;
              int contour_cursor = 0;
              int vert_cursor = 0;
              int offset_cursor = 0;
              int offset_faces_cursor = 0;
              int offset_faces_indexes_cursor = 0;
              int profile__faces_open_mode_cursor = 0;

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
                      res_errors.push_back(ObjectError(I, contour1, VAL2STR(Msg::_0009)". Polygon contains less 3 points."));
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
                      res_errors.push_back(ObjectError(I, plane1[0], VAL2STR(Msg::_0010)". Boundary polygon is not simple."));
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
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Msg::_0010)". Hole is not simple."));
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
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Msg::_0023)". Hole is not inside boundary."));
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
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Msg::_0024)". Hole intersect with hole."));
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
                    vect_ooai.push_back(OIOA(I, I, IJ, in_offset1, in_altitude1));
                    offset_cursor++;
                  }
                  
                  int in_count_of_profile_faces1 = in__profile_faces__counters[I]; // Количество faces в одном текущем объекте
                  std::vector<std::vector<int>> vect__profile_faces__face_indexes;
                  std::vector<int> vect__profile_faces__open_mode;
                  for (int IJ = 0; IJ <= in_count_of_profile_faces1 - 1; IJ++) {
                    int face_indexes_counter = in__profile_faces__indexes_counters[offset_faces_cursor++];
                    std::vector<int> vect__face_indexes;
                    for (int IJK = 0; IJK <= face_indexes_counter - 1; IJK++) {
                      int index = in__profile_faces__indexes[offset_faces_indexes_cursor++];
                      vect__face_indexes.push_back(index);
                    }
                    vect__profile_faces__face_indexes.push_back(vect__face_indexes);
                    int profile_face__open_mode = in__profile_faces__close_modes_array[profile__faces_open_mode_cursor++];
                    vect__profile_faces__open_mode.push_back(profile_face__open_mode);

                  }
                  map__object_index__profile_faces            [I] = vect__profile_faces__face_indexes;
                  map__object_index__profile_faces__open_modes[I] = vect__profile_faces__open_mode;

                  objects.push_back(ClsObject(I, I, vect_ooai, planes, SS_TYPE::UNDEFINED) );
                }
              }
              timer1.stop();
              if (verbose == true) {
                printf("\n " VAL2STR(Msg::_0057) ". SS 2D Offset. Load source data. Count of objects: %3zu, time of loading: % 10.5f", objects.size(), timer1.time());
                for (auto& object : objects) {
                  printf("\n " VAL2STR(Msg::_0067) ". SS 2D Offset. Object index: % 3d", object.object_index);
                  int plane_index = 0;
                  for (auto& plane : object.source_planes) {
                    printf("\n " VAL2STR(Msg::_0068) ".                       mesh: % 3d. contours : % 4d", plane_index++, (int)plane.size());
                  }
                }
              }
            }

            {
              if (verbose == true){
                printf("\n " VAL2STR(Msg::_0069) ". SS 2D Offset. Convert source data with source_objects_join_mode='%s'", source_objects_join_mode == 0 ? "SPLIT" : source_objects_join_mode == 1 ? "KEEP" : "MERGE");
              }
              CGAL::Real_timer timer1;
              timer1.start();
              // Представить список объектов в зависимости от параметра source_objects_join_mode
              std::vector<ClsObject> vect_split_mode;
              ClsObject objects_merge_mode(0, 0, std::vector<OIOA>(), std::vector<ContourSequence>(), SS_TYPE::UNDEFINED);
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
                    //  res_errors.push_back(ObjectError(object1.object_index, *ppwh.get(), VAL2STR(Msg::_0022) ". Holes of the Polygon intersect amongst themselves or with outer boundary ()"));
                    //} else 
                    if (source_objects_join_mode == 0 /* 0 - split */) {
                      std::vector<ContourSequence> vect_cs;
                      vect_cs.push_back( plane1 );
                      std::vector<OIOA> vect_oioa1;
                      for (auto& oioa1 : object1.vect_oioa) {
                        OIOA _oioa1(object_index, oioa1.object_original_index, oioa1.offset_index, oioa1.offset, oioa1.altitude);
                        vect_oioa1.push_back(_oioa1);
                      }
                      ClsObject obj1(object_index, object1.object_original_index, vect_oioa1, vect_cs, SS_TYPE::UNDEFINED );
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
              timer1.stop();
              if (verbose == true) {
                printf("\n " VAL2STR(Msg::_0028) ". SS 2D Offset. Time of conversion with join mode: %s : % 10.5f", source_objects_join_mode == 0 ? "SPLIT" : source_objects_join_mode == 1 ? "KEEP" : "MERGE", timer1.time() );
                printf("\n " VAL2STR(Msg::_0070) ". SS 2D Offset. After conversion. Count of objects: %3zu", objects.size());
                for (auto& object : objects) {
                  printf("\n " VAL2STR(Msg::_0071) ". SS 2D Offset. Object index: % 3d", object.object_index);
                  int plane_index = 0;
                  for (auto& pwh : object.vect_polygons_with_holes) {
                    printf("\n " VAL2STR(Msg::_0072) ".         polygon with holes: % 3d. contours : % 3d [external: % 6d verts, verts of holes: [", plane_index, 1/*external contour*/ + (int)pwh->number_of_holes(), (int)pwh->outer_boundary().size() );
                    for (auto& hole : pwh->holes()) {
                      printf(" % 6d", (int)hole.size());
                    }
                    printf(" ]");
                  }
                }
              }

            }

            {
              CGAL::Real_timer timer1;
              timer1.start();
              // Проверить, что теперь полигоны каждого объекта между собой не пересекаются (раньше были только holes отдельных полигонов):
              for (auto& object1 : objects) {
                {
                  std::set<int> set_intersected_polygons_indexes;
                  for (int I = (int)object1.vect_polygons_with_holes.size() - 1; I >= 1; I--) {
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
                          printf("\n ERROR. " VAL2STR(Msg::_0015) ". Polygons (%u, %u) INTERSECTS. Both marked as error", I, IJ);
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
                    res_errors.push_back(ObjectError(object1.object_index, *object1.vect_polygons_with_holes[I].get(), VAL2STR(Msg::_0011)". Polygons intersects"));
                    object1.vect_polygons_with_holes.erase(object1.vect_polygons_with_holes.begin() + I);
                  }
                }
              }

              timer1.stop();
              if (verbose == true) {
                printf("\n " VAL2STR(Msg::_0061) ". SS 2D Offset. Test of inner intersection. Time: % 10.5f", timer1.time());
                printf("\n " VAL2STR(Msg::_0058) ". SS Offset objects to process: %zu", objects.size());
              }
            }

            {
              CGAL::Real_timer timer1;
              timer1.start();
              int object1_index = 0;
              // Привести контуры объекта в соответствии с настройкой shape_mode.
              // Тип обработки контуров объекта: 
              // 0-ORIGINAL_BOUNDARIES (полная обработка), 
              // 1-EXCLUDE_HOLES (Keep only outer boundary), 
              // 2-INVERT_HOLES (Exclude outer boundary and fill holes).
              // <image url="..\code_images\file_0026.png" scale=".1"/>
              for (auto& object1 : objects) {
                // Видно, что не все буквы обработаны, потому что тут после split исходного объекта при увеличении количества объектов в параметре in_shapes_mode
                // "новые" индексы выходят за допустимый диапазон и тут или access violation или неопределённое поведение!
                // Если использовать object1_index напрямую, то возможна такая ошибка: <image url="..\code_images\file_0094.png" scale=".2"/>.

                //int shape_mode_index = object1_index <= in_count_of_objects ? object1_index : in_count_of_objects - 1;
                //int shapes_mode1 = in_shapes_mode[shape_mode_index];
                int shapes_mode1 = in_shapes_mode[object1.object_original_index];
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

              timer1.stop();
              if (verbose == true) {
                printf("\n " VAL2STR(Msg::_0060) ". SS 2D Offset. Sort contours. Time: % 10.5f", timer1.time());
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
                          ClsObject negObject(object1.object_index, object1.object_original_index, object1.vect_oioa, std::vector<ContourSequence>(), SS_TYPE::UNDEFINED /*Пока как неопределённый тип. Он вычислиться позже*/);
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
                        // 3. Соеденить Frame и все старые Outer Boundaries как один Polygon_with_holes_2 и рассчитать его как внутренний offset (но уже давать положительные значения offset как исходные) (удалить после такого рассчёта все внешние offset, т.к. они не должны фигурировать в результате)
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
                            auto max_abs_negative_offset_d = CGAL::to_double(max_abs_negative_offset);
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

                              ClsObject object1FrameForNegativeOffsets(object1NegativeOffsetSS.object1.object_index, object1NegativeOffsetSS.object1.object_original_index, std::vector<OIOA>(), std::vector<ContourSequence>(), SS_TYPE::NEGATIVE_WITH_EXTERNAL_FRAME);
                              object1FrameForNegativeOffsets.source_planes.push_back(cs_frame);
                              for (auto& oioa1 : object1NegativeOffsetSS.object1.vect_oioa) {
                                if (oioa1.offset < 0) {
                                  object1FrameForNegativeOffsets.vect_oioa.push_back(OIOA(oioa1.object_index, oioa1.object_original_index, oioa1.offset_index, oioa1.offset, OIOA::NEGATIVE_TYPE::WITH_EXTERNAL_FRAME, oioa1.altitude));
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
                            ClsObject object1WithInvertedHoles(object1NegativeOffsetSS.object1.object_index, object1NegativeOffsetSS.object1.object_original_index, std::vector<OIOA>(), vect_cs_inverted_holes, SS_TYPE::NEGATIVE_INVERTED_HOLE);
                            for (OIOA& oap1 : object1NegativeOffsetSS.object1.vect_oioa) {
                              if (oap1.offset < 0) {
                                object1WithInvertedHoles.vect_oioa.push_back(OIOA(object1NegativeOffsetSS.object1.object_index, object1NegativeOffsetSS.object1.object_original_index, oap1.offset_index, oap1.offset /*abs берётся позже*/, OIOA::NEGATIVE_TYPE::INVERTED_HOLE, oap1.altitude));
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
                                vect_oioa.push_back(OIOA(oioa1.object_index, oioa1.object_original_index, oioa1.offset_index, oioa1.offset, oioa1.altitude));
                              }
                            }
                            //if (vect_oioa.size() > 0) - Это условие отменено, т.е. в будущем для рассчёта faces SS требуются все SS и положительные и отрицательные. vect_oioa.size==0, когда все offset сконцентрированы в отрицательной зоне (даже без 0)
                            {
                              ClsObject object1WithPositiveOffsets(object1NegativeOffsetSS.object1.object_index, object1NegativeOffsetSS.object1.object_original_index, vect_oioa, object1NegativeOffsetSS.object1.source_planes, SS_TYPE::POSITIVE);
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

                for (auto& object1 : objects) {
                  if (object1.vect_polygons_with_holes.size() > 0) {
                    for (auto& polygon_with_holes1 : object1.vect_polygons_with_holes) {
                      // Отсортировать по offsets. Будет удобно считать, если при увеличении отступа пропадёт контур,то и считать следующие бОльшие контуры не надо
                      // - в многопоточности это уже не актуально. Может и можно в итоге прервать, но надо поискать другой метод. sort(positive_offsets_data1.vect_oap.begin(), positive_offsets_data1.vect_oap.end(), [](const OAP& oap1, const OAP& oap2) {return oap1.offset < oap2.offset; });
                      // Положительные отступы считаются по одному
                      POMS p1(
                        object1.object_index,
                        object1.object_original_index,
                        std::shared_ptr<Polygon_with_holes_2>(),
                        object1.vect_oioa,
                        object1.ss_type /*Изначально все создаваемые полигоны неопределены, т.к. отрицательные отступы рассчитываются разделением PWH-полигонов на вспомогательные объекты*/
                      );
                      p1.polygon1_with_holes = std::move(polygon_with_holes1);
                      vect_polygon1_oioa.push_back(p1);
                    }
                  } else {
                    // Если полигонов нет, то оставить этот вектор polygon1_with_holes пустым. Вернуть информацию по объекту всё равно надо.
                    POMS p1(
                      object1.object_index,
                      object1.object_original_index,
                      std::shared_ptr<Polygon_with_holes_2>(),
                      object1.vect_oioa,
                      object1.ss_type /*Изначально все создаваемые полигоны неопределены, т.к. отрицательные отступы рассчитываются разделением PWH-полигонов на вспомогательные объекты*/
                    );
                    vect_polygon1_oioa.push_back(p1);
                  }
                }

                int threadNumbers = boost::thread::hardware_concurrency();
                boost::mutex mtx_;

                {
                  CGAL::Real_timer timer1; // Общее время рассчёта всех SS (вместе с загрузкой данных в SS)
                  timer1.start();
                  boost::asio::thread_pool pool3(threadNumbers);
                  int threads_counts = 0;
                  size_t ss_id_counter = 0;
                  auto vect_polygon1_oioa_size = vect_polygon1_oioa.size();
                  auto vect_polygon1_oioa_rest = vect_polygon1_oioa_size;
                  for (auto& polygon1_oioa : vect_polygon1_oioa) {
                    ss_id_counter++;
                    boost::asio::post(pool3, [ss_id_counter, &polygon1_oioa, &res_errors, &res_contours, threads_counts, verbose, vect_polygon1_oioa_size, &vect_polygon1_oioa_rest, use_cache_of_straight_skeleton, &mtx_] {
                      unsigned int crc_val = 0;
                      bool is_crc_calculated = false;
                      if (polygon1_oioa.polygon1_with_holes == nullptr) {

                      } else {
                        SsBuilder ssb;
                        int verts_count = polygon1_oioa.polygon1_with_holes->outer_boundary().size();
                        {
                          // TODO: Временный параметр для равномерного offset. надо будет добавить weight во входные сокеты нода и заменить на получение weight из интерфейса нода
                          //std::vector<std::vector<FT> > uniform_weights;
                          //uniform_weights.reserve(polygon1_oioa.polygon1_with_holes->number_of_holes() + 1);

                          //uniform_weights.push_back(std::vector<FT>(polygon1_oioa.polygon1_with_holes->outer_boundary().size(), FT(1)));
                          //for (const auto& hole : polygon1_oioa.polygon1_with_holes->holes()) {
                          //  uniform_weights.push_back(std::vector<FT>(hole.size(), FT(1)));
                          //}

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
                          {
                            // Посчитать контрольные суммы:

                            // Загрузить все координаты контуров в вектор, по которому будет считаться CRC:
                            std::vector<FT> vect__points;
                            for (auto& point : polygon1_oioa.polygon1_with_holes->outer_boundary()) {
                              FT X = point.x();
                              FT Y = point.y();
                              vect__points.push_back(X);
                              vect__points.push_back(Y);
                            }
                            // 2. Загрузить holes
                            for (auto& hole1 : polygon1_oioa.polygon1_with_holes->holes()) {
                              verts_count += hole1.size();
                              for (auto& point : hole1) {
                                FT X = point.x();
                                FT Y = point.y();
                                vect__points.push_back(X);
                                vect__points.push_back(Y);
                              }
                            }

                            std::stringstream str_stream;
                            boost::archive::binary_oarchive oa(str_stream);
                            for (auto& p : vect__points) {
#ifdef _DEBUG
                              auto& num = p.exact();
#else
                              auto& num = p;
#endif // DEBUG

                              oa << num;
                            }
                            std::string serialized_data = str_stream.str();
                            const char* byte_buffer = serialized_data.c_str();
                            size_t length = serialized_data.length();

                            // Создаем объект для вычисления CRC-32
                            boost::crc_32_type crc;
                            crc.process_bytes(&byte_buffer[0], length);
                            crc_val = crc.checksum();
                            is_crc_calculated = true;
                            //if (verbose == true) {
                            //  mtx_.lock();
                            //  printf("\n " VAL2STR(Msg::_0055) ". ss_id=%2zu; crc=%2u", ss_id_counter, crc_val);
                            //  mtx_.unlock();
                            //}
                          }
                        }
                        // Рассчитать SS
                        CGAL::Real_timer timer2;
                        timer2.start();
                        bool is_ss_cached = false;
                        bool is_ss_used_as_chached = false;
                        try {
                          // Если рассчёт не требует кэша или в кэше ничего нет, то выполнить рассчёт SS в любом случае.
                          if (use_cache_of_straight_skeleton == false ||
                            map__crc__ss.find(crc_val) == map__crc__ss.end()
                            ) {
                            polygon1_oioa.ss = ssb.construct_skeleton();
                            // ... и не запоминать (вдруг утечка памяти какая-то?). Примечание: сюда можно попасть и при use_cache_of_straight_skeleton == true, если в кэше ничего нет.
                            // Update: решил пока перезаписывать кэш в любом случае.
                            //if (use_cache_of_straight_skeleton == true) {
                            map__crc__ss[crc_val] = polygon1_oioa.ss;
                            //}
                          } else {
                            polygon1_oioa.ss = map__crc__ss[crc_val];
                            is_ss_cached = true;
                            is_ss_used_as_chached = true;
                            //if (verbose == true) {
                            //  mtx_.lock();
                            //  printf("\n " VAL2STR(Msg::_0056) ". ss_id=%2zu; used as cached with crc=%2u", ss_id_counter, crc_val);
                            //  mtx_.unlock();
                            //}
                          }
                        } catch (std::exception _ex) {
                          // При ошибке исходный SS останется null
                        }

                        timer2.stop();
                        if (verbose == true) {
                          mtx_.lock();
                          vect_polygon1_oioa_rest--;
                          std::time_t t = std::time(nullptr);
                          std::tm* now = std::localtime(&t);

                          printf("\n " VAL2STR(Msg::_0026) ". SS 2D Offset. Calc SS. || ss_id: %4zu Thread: %4u/ %4zu/ %4zu, SsBuilder verts: % 6d, build time: % 12.5f, cached: %s, used from cache: %s, crc: %10u, time: %02d:%02d:%2d", ss_id_counter, threads_counts, vect_polygon1_oioa_rest /*осталось*/, vect_polygon1_oioa_size/*Количество объектов во время многопоточного рассчёта*/, verts_count, timer2.time(), is_ss_cached == true ? "YES" : " NO", is_ss_used_as_chached == true ? "YES" : " NO", crc_val, now->tm_hour, now->tm_min, now->tm_sec  /*Вывод текущего времени, чтобы видеть сколько времени прошло с момента задержки, когда большие скелетоны*/);
                          mtx_.unlock();
                        }
                      }
                      mtx_.lock();
                      // Proceed only if the skeleton was correctly constructed.
                      if (polygon1_oioa.ss) {
                        polygon1_oioa.ss_id = ss_id_counter;
                        // Запомнить ss_id в соответствующих ему oioa (ss_id не зависит от индескса объекта и является сквозным)
                        for (auto& oioa1 : polygon1_oioa.vect_oioa) {
                          oioa1.ss_id = ss_id_counter;
                        }
                      } else {
                        try {
                          // Не получилось посчитать Straight Skeleton. Надо сообщить о проблеме пользователю:
                          const Polygon_with_holes_2& _pwh = polygon1_oioa.polygon1_with_holes==nullptr ? Polygon_with_holes_2() : *polygon1_oioa.polygon1_with_holes.get();
                          res_errors.push_back(ObjectError(polygon1_oioa.vect_oioa[0].object_index, _pwh, VAL2STR(Msg::_0006)". Offset. Error Build Straight Skeleton. Outer Boundary"));
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
                  timer1.stop();
                  if (verbose == true) {
                    printf("\n " VAL2STR(Msg::_0045) ". SS 2D Offset. Count of SS: %d,                                                          build time: % 12.5f", threads_counts, timer1.time() );
                  }
                }

                {
                  CGAL::Real_timer timer1; // Общее время рассчёта всех offset
                  timer1.start();
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

                          CGAL::Real_timer timer2; // Общее время рассчёта всех SS (вместе с загрузкой данных)
                          timer2.start();
                          CGAL::offset_params(oioa1.ss, oioa1.offset, CGAL::parameters::verbose(verbose), oioa1.map_halfedge_to_offset_points, std::back_inserter(offset_contours));
                          timer2.stop();
                          if (verbose == true) {
                            mtx_.lock();
                            printf("\n " VAL2STR(Msg::_0046) ". SS 2D Offset. building offset. || ss_id: % 4zd, % 3d, % 10.5f, build time: % 10.5f, cntrs: %3zu ", oioa1.ss_id, oioa1.offset_index, CGAL::to_double(oioa1.offset), timer2.time(), offset_contours.size());
                            if (offset_contours.size() > 0) {
                              printf("[");
                              for (auto& c1 : offset_contours) {
                                printf(" %6zu", c1->vertices().size() );
                              }
                              printf(" ]");
                            }
                            mtx_.unlock();
                          }

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
                  timer1.stop();
                  if (verbose == true) {
                    printf("\n " VAL2STR(Msg::_0047) ". SS 2D Offset.                                  calc of offsets general time: % 10.5f", timer1.time());
                  }
                }

                // Восстановить последовательность результатов расчётов offset по исходным позициям offset (как их задал пользователь):
                sort(res_contours.begin(), res_contours.end(), [](const OIOA& oap1, const OIOA& oap2) {return oap1.offset_index < oap2.offset_index; });
                sort(res_contours.begin(), res_contours.end(), [](const OIOA& oap1, const OIOA& oap2) {return oap1.object_index < oap2.object_index; });
                sort(res_errors.begin(), res_errors.end(), [](const ObjectError& oap1, const ObjectError& oap2) {return oap1.object_index < oap2.object_index; });
              }
            } else {
              if (verbose == true) {
                unsigned int errors_size = res_errors.size();
                if (errors_size == 0) {
                  printf("\n " VAL2STR(Msg::_0074) ". Straight Skeleton Offsets tests finished. Errors not found.");
                } else {
                  printf("\n " VAL2STR(Msg::_0075) ". Straight Skeleton Offsets tests finished. Errors count %u", errors_size);
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
            printf("\n " VAL2STR(Msg::_0048) ". Offset computation of object took % 10.5f sec.", timer.time());
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

              // Сопоставление нового и старого индексов результирующих mesh
              // Если режим SPLIT, то в результате объекты могут разделиться на части, но части будут помнить кем они были раньше
              // Если режим  KEEP, то объекты не будут делиться, но тоже будут помнить, кем они были.
              // Если режим MERGE, то объекты объединяться и будут помнить тольк индекс первого объекта.
              std::map<int, int> map__new_object_index__original_object_index;

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
                        // Запомнить привязку индекса нового объекта к индексу оригинального объекта
                        map__new_object_index__original_object_index[split_object_index] = pwhs_offset_index_order1.oioa1.object_original_index;
                        split_object_index++;
                      }
                    } else if (results_join_mode == 1) {
                      keep_object_index = pwhs_offset_index_order1.oioa1.object_index;
                      // Оставить результирующие объекты как в исходных объектах
                      if (map_join_mesh.find(keep_object_index) == map_join_mesh.end()) {
                        map_join_mesh[keep_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                      }
                      map_join_mesh[keep_object_index].push_back(pwhs_offset_index_order1);
                      // Запомнить привязку индекса нового объекта к индексу оригинального объекта
                      map__new_object_index__original_object_index[keep_object_index] = pwhs_offset_index_order1.oioa1.object_original_index;
                    } else { // 2 - для всех остальных случаев
                      // Объеденить все результаты в один объект
                      if (map_join_mesh.find(merge_object_index) == map_join_mesh.end()) {
                        map_join_mesh[merge_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                      }
                      map_join_mesh[merge_object_index].push_back(pwhs_offset_index_order1);
                      if (map__new_object_index__original_object_index.find(merge_object_index) == map__new_object_index__original_object_index.end()) {
                        // Запомнить привязку индекса нового объекта к индексу оригинального объекта (для MERGE это надо делать один раз)
                        map__new_object_index__original_object_index[merge_object_index] = pwhs_offset_index_order1.oioa1.object_original_index;
                      }
                    }
                  }
                } else if (result_type == ResType::BEVEL || result_type == ResType::STRAIGHT_SKELETON) {
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
                        if (poms.ss == nullptr) {
                          // Иногда ss не рассчитан (может быть по причине ошибок в исходных данных), поэтому и mesh посчитать тоже не получится.
                          continue;
                        }
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
                              const auto& new_point_id = mesh_ss.add_vertex(p3);
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
                            // Тут получена Mesh SS, у которой не надо делать merge points. Вроде как теперь дальнейший merge нескольких mesh_ss должно пройти быстрее?
                            // 
                            // Примечание: ввиду того, что offset пока не учитывается, то геометрия ss хранится только в первом элементе map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_id][ss_id][0]
                            // Примечание 2: подумать, как вообще вынести отсюда сохранение mesh_ss в map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS.
                          }
                        }
                      }
                    }
                  }

                  
                  std::map<int, std::vector<Mesh>> map__object_index__mesh_ss_merged;

                  if (result_type == ResType::STRAIGHT_SKELETON) {
                    for (auto& [object_index, vect_poms] : map__object_index__poms) {
                      // Собрать все meshes от полученных сеток SS в один вектор
                      std::vector<Mesh> vect__object__meshes;
                      for (auto& poms : vect_poms) {
                        vect__object__meshes.push_back(poms.mesh_ss_01_source);
                      }
                      map__object_index__mesh_ss_merged[object_index] = vect__object__meshes;
                    }
                  } else if (result_type == ResType::BEVEL) {

                    // Связать Object_index и Mesh SS, чтобы в дальнейшем определить смежные faces SS для каждого объекта object_index целиком. Пример результата таблицы смежных faces <image url="..\code_images\file_0049.png" scale=".1"/>

                    std::map<int, Mesh> map__object_index__mesh_ss_joined;
                    for (auto& [object_index, vect_poms] : map__object_index__poms) {
                      // Просуммировать все SS Mesh, которые получены в текущем объекте с сохранением свойств faces (ss_id и face_index).
                      // Для этого нужно выполнить загрузку всех faces в один mesh и смержить точки.
                      // Создаем новую поверхность для объединения
                      Mesh result_ss_Mesh;
                      // Чтобы дополнительные свойства добавлялись в result_mesh_ss нужно и туда тоже добавить это свойство (без этого это свойство в исходных mesh копироваться не будет)
                      result_ss_Mesh.add_property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct", MeshFacePropertyStruct(-1, MeshFacePropertyStruct::SS_SIGN::UNDEFINED, -1, -1));
                      for (auto& poms : vect_poms) {
                        result_ss_Mesh.join(poms.mesh_ss_01_source); // - в исходниках вроде как делает merge vertices и faces -  update - нет, не делает
                      }
                      map__object_index__mesh_ss_joined[object_index] = result_ss_Mesh;
                    }

                    /*Примечание. доступ к элементам objects не нужно делать через этот object_index. Objects назначаются по одному object на ss_id, если есть отрицательные offsets.*/
                    for (auto& [object_index, mesh_ss_joined] : map__object_index__mesh_ss_joined) {
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

                      Mesh mesh_ss_merged; // Сюда положить результат merge для получения Mesh рассматриваемого SS.
                      auto [mesh_ss_merged__property_map, created] = mesh_ss_merged.add_property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct", MeshFacePropertyStruct(-1, MeshFacePropertyStruct::SS_SIGN::UNDEFINED, -1, -1)); // Добавить структуру свойств. По умолчанию в новый объект свойства из join-объектов не копируются

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
                          printf("\n " VAL2STR(Msg::_0032) ". SS 2D Offset. object_index %u, ss merge ajacent timer % 10.5f", object_index, timer_merge.time());
                        }
                      }



                      // Загрузка атрибутов faces из старого результирующего Mesh:
                      const auto& _map_mesh_ss_property_map = mesh_ss_joined.property_map<Mesh::Face_index, MeshFacePropertyStruct>("f:MeshFacePropertyStruct");
                      //if (_map_mesh_ss_property_map.has_value() == true) 
                      {
                        auto& mesh_ss_source_property_map = _map_mesh_ss_property_map.value();
                        //int counter = 0;
                        //// Вывод property_map старого результирующего mesh:
                        //for (auto& face_index : result_ss_mesh.faces()) {
                        //  auto& val = mesh_ss_property_map[face_index];
                        //  if (verbose == true) {
                        //    printf("\n " VAL2STR(Msg::_0030) ". SS 2D Offset. SS merged: %u. ss_id: %u, face index: %d", counter++, val.ss_id, val.face_index);
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
                        // <image url="..\code_images\file_0064.png" scale=".3"/>
                        std::map<int /*mesh_face_id*/, int /*mesh_face_id*/                                > map__adjacent_mesh_faces;
                        std::map<int /*mesh_face_id*/, MeshFacePropertyStruct                              > map__mesh_face_id__face_info;
                        std::map<int /*       ss_id*/, std::map<int /*ss_face_id*/, MeshFacePropertyStruct>> map__ss_id__ss_face_id__face_info;

                        // Скопировать информацию о faces в смерженный новый mesh (ss_id и сквозной индекс face_index)
                        for (const auto& f1 : mesh_ss_merged.faces()) {
                          // Скопировать faces property из старого результирующего mesh в новый результируюший mesh (проверил, что faces идут по очереди как и в исходном несмерженном результирующем mesh).
                          mesh_ss_merged__property_map[f1] = mesh_ss_source_property_map[f1];
                        }

                        // Тут непосредственно рассчитываются смежные faces на новом Mesh merged перебором всех halfedges:
                        for (const auto& hd : mesh_ss_merged.halfedges()) {
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
                          auto& he1 = hd;
                          // Работа verbose==true, при некорректной mesh, так выглядят ошибки (вывод в к консоль): <image url="..\code_images\file_0047.png" scale=".1"/>. 
                          if (mesh_ss_merged.is_valid(he1, verbose) == false || mesh_ss_merged.is_removed(he1) == true) {
                            continue;
                          }

                          // Получаем второе полуребро (противоположное)
                          const auto& he2 = mesh_ss_merged.opposite(he1);
                          // И его тоже надо проверить, что оно не крайнее:
                          if (mesh_ss_merged.is_border(he2) == true) {
                            continue;
                          }
                          if (mesh_ss_merged.is_valid(he2, verbose) == false || mesh_ss_merged.is_removed(he2) == true) {
                            continue;
                          }

                          // Получаем первую грань
                          auto face1 = mesh_ss_merged.face(he1);

                          // Получаем вторую грань
                          auto face2 = mesh_ss_merged.face(he2);

                          // Проверить, что они валидные (по хорошему надо что-то вывести?)
                          if (mesh_ss_merged.is_valid(face1, verbose) == false || mesh_ss_merged.is_removed(face1) == true) {
                            continue;
                          }
                          if (mesh_ss_merged.is_valid(face2, verbose) == false || mesh_ss_merged.is_removed(face2) == true) {
                            continue;
                          }

                          auto& face1_info = mesh_ss_merged__property_map[face1];
                          auto& face2_info = mesh_ss_merged__property_map[face2];
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

                        for (const auto& fd : mesh_ss_merged.faces()) {
                          auto& face1_info = mesh_ss_merged__property_map[fd];
                          map__ss_id__ss_face_id__face_info[face1_info.ss_id] = std::map<int, MeshFacePropertyStruct>();
                        }
                        for (const auto& fd : mesh_ss_merged.faces()) {
                          auto& face1_info = mesh_ss_merged__property_map[fd];
                          map__ss_id__ss_face_id__face_info[face1_info.ss_id][face1_info.ss_face_id] = face1_info;
                        }

#ifdef _DEBUG
                        if (verbose == true) {
                          // Вывести полученные пары: <image url="..\code_images\file_0049.png" scale=".2"/>
                          printf("\n faces pairs of object_index: %u :=========================", object_index);
                          for (auto& [f1_index, f2_index] : map__adjacent_mesh_faces) {
                            printf("\n " VAL2STR(Msg::_0031) ". SS 2D Offset. object_index: % 2u, SS faces map: % 2u/% 2u. ss_id: % 2u/% 2u  >>>> % 2u/% 2u", object_index, map__mesh_face_id__face_info[f1_index].ss_id, f1_index, map__mesh_face_id__face_info[f2_index].ss_id, f2_index, f1_index, f2_index);
                          }
                        }
#endif
                        // Тут уже можно начать обрабатывать offset-ы текущего объекта object_index:
                        {
                          CGAL::Real_timer timer_bevel; // Тратится около 3% времени на merge и расчёт ajacent faces.
                          double t1 = 0;
                          // Примечание: Все SS являются плоскими объектами, имеющими один или больше контуров: первый контур всегда внешний (если есть один, то это только внешний контур), остальные
                          // всегда внутренние. Рекурсии у SS нет! Если есть внутренние контуры, то топологически SS является двусторонним - внешняя и внутренняя стороны (не верх или низ как у плоскости,
                          // а именно контуры).

                          /// <summary>
                          /// Сопоставить определённым faces от SS набор точек пересечений с полурёбрами этого SS заданному offset.
                          /// Повтор другими словами: Есть набор offset с точками пересечений в SS, но в одном offset указывается вектор из нескольких SS,
                          /// с которыми этот offset пересекается. Требуется преобразовать этот набор в зависимость от SS, чтобы знать с какими offset пересекается этот ss_id.
                          /// </summary>
                          struct SS_Offset_HalfEdge_Intersections_Point_2 {
                            /// <summary>
                            /// Информация об одном offset, относящемся к имеющемуся ss.
                            /// </summary>
                            FT  oioa_offset;
                            int oioa_offset_index;
                            FT  oioa_altitude;

                            /// <summary>
                            /// В каких точках на плоскости пересекается ss с offset из oioa
                            /// </summary>
                            std::unordered_map<HDS_Halfedge_const_handle, Point_2> map__halfedge__offset_points;

                            SS_Offset_HalfEdge_Intersections_Point_2(FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude, /*std::shared_ptr<Ss> _ss,*/ std::unordered_map<HDS_Halfedge_const_handle, Point_2> _map__halfedge__offset_points)
                              : oioa_offset(_oioa_offset), oioa_offset_index(_oioa_offset_index), oioa_altitude(_oioa_altitude), /*ss(_ss),*/ map__halfedge__offset_points(_map__halfedge__offset_points) {
                            }
                          };

                          // Список всех offset у ss_id с информацией какие halfedge от исходного SS пересекаются с этим offset с сортировкой по offset_index (от меньшего к большему)
                          std::map<int /*ss_id*/, std::vector<SS_Offset_HalfEdge_Intersections_Point_2>> map__ss_id__he_intersections_point2;

                          // Обработку контуров надо вести либо по map_adjacent_mesh_faces (если есть отрицательные SS), либо по всем faces от рассчитанных SS (если SS только положительные).
                          {
                            for (auto& [ss_id, oioa_offset_ss_params] : map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_index]) {
                              for (auto& offset_params : oioa_offset_ss_params) {
                                for (auto& he1 : offset_params.vect__SS__HalfEdge__Point_2 /*Какие точки пересечений offset и SS соответствуют halfedge*/) {
                                  if (map__ss_id__he_intersections_point2.find(he1.ss_id) == map__ss_id__he_intersections_point2.end()) {
                                    map__ss_id__he_intersections_point2[he1.ss_id] = std::vector<SS_Offset_HalfEdge_Intersections_Point_2>();
                                  }
                                  map__ss_id__he_intersections_point2[he1.ss_id].push_back(SS_Offset_HalfEdge_Intersections_Point_2(offset_params.oioa1.offset, offset_params.oioa1.offset_index, offset_params.oioa1.altitude, /*he1.ss,*/ he1.map__halfedge__offset_points));
                                }
                              }
                            }
                            // Отсортировать точки пересечения по offset, т.к. только в этом порядке можно однозначно предположить/задать последовательность пересечения списка offset-ов с этим halfedge,
                            // чтобы позже, при обходе контура, выстроить точки пересечения в каждом halfedge по порядку offset_index (с учётом наклона/slope halfedge):
                            for (auto& [ss_id, vect_intersection_points] : map__ss_id__he_intersections_point2) {
                              std::sort(vect_intersection_points.begin(), vect_intersection_points.end(), [](const SS_Offset_HalfEdge_Intersections_Point_2& elem1, const SS_Offset_HalfEdge_Intersections_Point_2& elem2) {
                                bool res = false;
                                if (elem1.oioa_offset != elem2.oioa_offset) {
                                  res = elem1.oioa_offset < elem2.oioa_offset;
                                } else {
                                  // Для равных offset всё равно нужно определить направление обхода. Сначала сравнить altitude:
                                  if (elem1.oioa_altitude != elem2.oioa_altitude) {
                                    res = elem1.oioa_altitude < elem2.oioa_altitude;
                                  } else {
                                    // Индексы никогда не равны друг другу. Будем считать, что при обходе по часовой стрелке контура faces segment больший индекс должен быть первым.
                                    res = elem1.oioa_offset_index > elem2.oioa_offset_index;
                                  }
                                }
                                return res;
                                });
                            }

                            // offset-ы строятся на сегментах контуров SS. Если нет отрицательных SS, то нужно выбрать один face (любой из faces SS находится на контуре), если есть отрицательные SS, то нужно выбрать два смежных face.
#ifdef _DEBUG
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
#endif

                            /// <summary>
                            /// Точка, используемая для точек плана Beveled SS. Будет использоваться для рассчёта Beveled SS. Список этих точек представляет собой "план" по которому уже в 3D будет строится полная версия Beveled SS.
                            /// Beveled SS
                            /// </summary>
                            struct SHAPE_POINT {
                              /// <summary>
                              /// Индекс точки в общем массиве точек. Отрицательные индексы испольуются при рассчёте mesh, чтобы не занимать индексы с неотрицательными значениями,
                              /// которые будут использоваться при рассчёте вершин для точек пересечений he с SS и преобразования проектных точек в реальные точки.
                              /// </summary>
                              int index;

#ifdef _DEBUG
                              /// <summary>
                              /// Если эта точка типа PROJECT_POINT, то при рассчёте BEVEL SS записать в этот вектор индекс точки, в кого эта точка преобразована.
                              /// Примечание: Точка типа PROJECT_POINT не может напрямую использоваться при построении faces, ограниченной двумя offset, т.к. пары offset или цепочки offset могут пересекать эту точку несколько раз.
                              /// При этом сами offset определены однозначно и у них дубликатов быть не может, даже если они совпадают. Также возможны ситуации, что точка никогда не понадобится, если эта точка не будет 
                              /// находится между двумя offsets (хотя в будущем точка может потребоваться для caps-ов). Применяется для справки. Непосредственно для обработки bevel - не нужен.
                              /// <image url="..\code_images\file_0073.png" scale=".2"/>
                              /// <image url="..\code_images\file_0088.png" scale=".5"/>
                              /// </summary>
                              std::vector<std::uint32_t> beveled_indexes; // для тестов. Для рассчётов не требуется
#endif
                              /// <summary>
                              /// Важно для точки типа PROJECT_POINT. Индекс точки из какой проектной точки она получена.
                              /// </summary>
                              int project_point_index;
                              int application_index;

                              /// <summary>
                              /// Сама точка. Только высота по Z у точки рассчитывается не всегда. Высота заранее известна только у точек, лежащих на offset
                              /// (т.к. по исходным данных при задании offset задаётся ещё и altitude).
                              /// Если точка находится на контуре, то её высоту надо будет рассчитать позже. Это будет сделано только после обработки всех
                              /// точек.
                              /// </summary>
                              Point_3 point;

                              /// <summary>
                              /// Рассчитана ли высота у точки? (true - рассчитана, false - не рассчитана). Важно только для точек типа PROJECT_POINT
                              /// </summary>
                              bool is_altitude_calc;

                              /// <summary>
                              /// Является ли точка внешней? Важно только для крайнего контура SS. При построении BEVELED SS важно знать, является ли точка крайней или нет. Если она крайняя,
                              /// то на отрицательных SS она не будет участвовать в построении BEVELED SS или Caps.
                              /// </summary>
                              bool is_border; // Забавно, что в SS используется термин is_boundary, а в mesh - is_border. Т.к. работаю с mesh, то использую is_border.

                              /// <summary>
                              /// Какой знак у SS, на котором находится текущая точка
                              /// </summary>
                              MeshFacePropertyStruct::SS_SIGN ss_sign;
                              int ss_id;
                              int vertex_id;

                              /// <summary>
                              /// Для точек типа PROJECT_POINT задаёт величину time, на значении которой образовалась точка (чтобы потом сравнить offset-ы между которыми находится эта
                              /// точка и определить её высоту). См. функцию construct_lateral_faces.
                              /// </summary>
                              FT event_time;

                              enum POINT_TYPE {
                                
                                /// <summary>
                                /// Применяется для конструктора по умолчанию. В работе не используется.
                                /// </summary>
                                UNDEFINED=-1,
                                
                                /// <summary>
                                /// Применяется к точкам типа PROJECT_POINT после их рассчёта в определённом ss_id.
                                /// </summary>
                                APPLICATED_PROJECT_POINT,

                                /// <summary>
                                /// Точка из контура (не пересекающая edge). По сути - это все точки SS. Одна и таже точка типа PROJECT_POINT может встречаться несколько раз при обработке переходов offset через эту точку при реверсных направлениях пар offset.
                                /// </summary>
                                PROJECT_POINT,

                                /// <summary>
                                /// Точка из offset. Точка пересекает offset вверх (в контуре. В faces это направление должно учитывать относительные положения edges, из которых строится face между 2-мя offset) (slope=+1)
                                /// </summary>
                                UP,

                                /// <summary>
                                /// Точка из offset. Точка пересекает edge вниз (slope=-1)
                                /// </summary>
                                DOWN,
                              };
                              POINT_TYPE point_type;

                              /// Параметры offset по которым строилась эта точка
                              FT  oioa_offset;
                              int oioa_offset_index;
                              FT  oioa_altitude;

                              SHAPE_POINT(
                                int         _index,
                                int         _project_point_index,
                                POINT_TYPE  _point_type,
                                Point_3&    _point,
                                bool        _is_altitude_calc,
                                FT          _oioa_offset,
                                int         _oioa_offset_index,
                                FT          _oioa_altitude,
                                bool        _is_border,
                                int         _ss_id,
                                MeshFacePropertyStruct::SS_SIGN _ss_sign,
                                int         _vertex_id,
                                int         _application_index,
                                FT          _event_time
                              )
                                : 
                                index               (_index),
                                project_point_index (_project_point_index),
                                point_type          (_point_type),
                                point               (_point),
                                is_altitude_calc    (_is_altitude_calc),
                                oioa_offset         (_oioa_offset),
                                oioa_offset_index   (_oioa_offset_index),
                                oioa_altitude       (_oioa_altitude),
                                is_border           (_is_border),
                                ss_id               (_ss_id),
                                ss_sign             (_ss_sign),
                                vertex_id           (_vertex_id),
                                application_index   (_application_index),
                                event_time          (_event_time) {

                              }

                              SHAPE_POINT()
                                :
                                index(-1),
                                project_point_index(-1),
                                point_type(POINT_TYPE::UNDEFINED),
                                point(0.0, 0.0, 0.0),
                                is_altitude_calc(false),
                                oioa_offset(0.0),
                                oioa_offset_index(-1),
                                oioa_altitude(0.0),
                                is_border(false),
                                ss_id(0),
                                ss_sign(MeshFacePropertyStruct::SS_SIGN::UNDEFINED),
                                vertex_id(0),
                                application_index(-1),
                                event_time(0.0) {

                              }
                            };

                            /// <summary>
                            /// Вектор результирующих точек SS-Beveled-проекции
                            /// </summary>
                            struct SHAPE_POINTS {
                              /// <summary>
                              /// Список точек с индексами.
                              /// </summary>
                              //std::unordered_map<int /*point_index*/, SHAPE_POINT> map__point_index__calc_point; // Заменил на std::unordered_map, но это подняло производительность излечения элемента только в 2 раза.
                              std::vector<SHAPE_POINT> vect__point_index__mesh_point;
                              std::vector<SHAPE_POINT> vect__point_index__project_point;

                              /// <summary>
                              /// Внутренний счётчик точек при создании новых точек пересечения, которые уже точно пойдут в финальную фигуру.
                              /// </summary>
                              int ss_intersects_he_points_counter = 0;

                              /// <summary>
                              /// Внутренний счётчик точек при создании точек типа PROJECT_POINT (поэтому Ppoint). И он отрицательный для таких точек. Эти точки
                              /// в будущем не будут использоваться напрямую, а будут участвовать в построении BEVELED SS, потому что могут встречаться в
                              /// конечной фигуре по много раз. (TODO: Сделать скрин лесенки, проходящей нат такой точкой, а потом объёмное изображение
                              /// ситуации). Примечание: сделан отрицательным, чтобы не пересчитывать точки, которые уже точно пойдут в финальную фигуру
                              /// (по счётчику ss_intersects_he_points_counter).
                              /// </summary>
                              int ss_ppoints_counter = -1;

                              /// <summary>
                              /// Сопоставление индексов результирующих точек пересечениям с полурёбрами he во время обхода segment-ов.
                              /// Параметр map__ss_id__he_handle и map__ss_id__vertex_id работают в паре при определении результирующих точек.
                              /// </summary>
                              std::map<std::tuple<int /*ss_id*/, HDS_Halfedge_const_handle, int /*offset_index*/>, int /*уникальный идентификатор точки*/> map__ss_id__he_handle;

                              /// <summary>
                              /// Сопоставление индексов результирующих точек с точками самого SS во время обхода segment-ов.
                              /// Параметр map__ss_id__he_handle и map__ss_id__vertex_id работают в паре при определении результирующих точек.
                              /// </summary>

                              std::map<std::tuple<
                                int /*ss_id*/,
                                int /*vertex_id*/,
                                int, // индекс профиля profile_face_index (раньше mesh был один на все профильные face, а теперь он разделён по отдельности на подобъекты и об этом нужно помнить с помощью этого индекса) <image url="..\code_images\file_0087.png" scale=".1"/>
                                int /*количество использований. -1 - точка типа PROJECT_POINT, которая ещё не использовалась. 0 и больше - уже применялась*/>,
                                int /*уникальный идентификатор точки*/
                              > map__ss_id__vertex_id;

                              // Добавить точку he в финальный mesh
                              int get_index_or_append_he__offset_index(
                                int ss_id,
                                SHAPE_POINT::POINT_TYPE point_type,
                                HDS_Halfedge_const_handle& he_handle,
                                int offset_index,
                                Point_3& _point,
                                bool _is_altitude_calc,
                                FT _oioa_offset,
                                int _oioa_offset_index,
                                FT _oioa_altitude,
                                MeshFacePropertyStruct::SS_SIGN ss_sign
                              ) {
                                const auto& ss_id__he_handle__offset_index = std::tuple(ss_id, he_handle, offset_index);
                                const auto& it = map__ss_id__he_handle.find(ss_id__he_handle__offset_index);
                                int res_counter = -1;
                                if (it == map__ss_id__he_handle.end()) {
                                  SHAPE_POINT mesh_point(
                                    -1/*индекс ещё неизвестен*/,
                                    -1 /* для точки пересечения he не важно. Важно только для точки типа PROJECT_POINT */,
                                    point_type,
                                    _point,
                                    _is_altitude_calc,
                                    _oioa_offset,
                                    _oioa_offset_index,
                                    _oioa_altitude,
                                    false /*точки пересечения точно не могут быть крайними, т.к. они находятся между точками SS*/,
                                    ss_id,
                                    ss_sign,
                                    0 /*это не тип PROJECT_POINT,значение не важно*/,
                                    -1 /*_application_index для точки пересечения he не важен */,
                                    0.0 /*event_time - тут не имеет значения*/
                                  );
                                  vect__point_index__mesh_point.push_back(mesh_point);
                                  res_counter = vect__point_index__mesh_point.size() - 1;
                                  vect__point_index__mesh_point.back().index =
                                    map__ss_id__he_handle[ss_id__he_handle__offset_index] =
                                    res_counter;
                                } else {
                                  auto& [tpl, _res_counter] = (*it);
                                  res_counter = _res_counter;
                                }
                                return res_counter;
                              };

                              /// <summary>
                              /// Добавить точку типа PROJECT_POINT в отрицательную зону проиндексированных точек. Это предварительная точка, которая будет использоваться для рассчёта BEVELED SS
                              /// (это точка не является точкой пересечения he с SS, высота её неизвестна, и будет рассчитана позже по соседним точкам в используемом face. Если точка используется
                              /// в нескольких face, то рассчитывать её нужно только один раз, после этого у неё выставляется признак и если она встретится в другом face, то уже будет известно, 
                              /// что она рассчитана).
                              /// </summary>
                              int get_index_or_append_vertex_id(
                                int   ss_id,
                                MeshFacePropertyStruct::SS_SIGN ss_sign,
                                int   vertex_id,
                                SHAPE_POINT::POINT_TYPE point_type,
                                Point_3& _point,
                                bool _is_altitude_calc,
                                FT   _oioa_offset,
                                int  _oioa_offset_index,
                                FT   _oioa_altitude,
                                bool  is_border,
                                int   application_index,
                                FT   _event_time
                              ) {
                                const auto& ss_id__vertex_id =
                                  std::tuple(
                                    ss_id,
                                    vertex_id,
                                    -1/*при инициализации профильной точки profile_face_index не важен, т.к. эта точка не используется как реальная, а является пока ещё только "проектной" и на её основе позже будет создано столько точек, сколько раз она потребуется при построении.*/,
                                    -1
                                  );
                                const auto& it = map__ss_id__vertex_id.find(ss_id__vertex_id);
                                int res_counter = -1;
                                if (it == map__ss_id__vertex_id.end()) {
                                  // res_counter = ss_ppoints_counter--;
                                  
                                  // Если точка не найдена, то создать её с новым индексом
                                  SHAPE_POINT project_point(
                                      -1/*индекс ещё не известен*/,
                                      -1 /* project_point_index - изначально она ни на кого не ссылается */,
                                      point_type,
                                      _point,
                                      _is_altitude_calc,
                                      _oioa_offset,
                                      _oioa_offset_index,
                                      _oioa_altitude,
                                      is_border,
                                      ss_id,
                                      ss_sign,
                                      vertex_id,
                                      application_index,
                                      _event_time
                                    );
                                  vect__point_index__project_point.push_back(project_point);
                                  res_counter = -vect__point_index__project_point.size() /*1. Индексы проектных точек отрицательные (<0), 2. поэтому самый первый индекс project_point должен быть -1, затем -2, -3 и т.д. Диапазон индексов mesh_point (отрицательные индексы) и project_point (положительные индексы) должен быть неразрывным и непересекающимся*/;
                                  vect__point_index__project_point.back().index =
                                    map__ss_id__vertex_id[ss_id__vertex_id] =
                                    res_counter;
                                } else {
                                  auto& [tpl, _res_counter] = (*it);
                                  res_counter = _res_counter;
                                }
                                return res_counter;
                              };

                              /// <summary>
                              /// Превратить точку типа PROJECT_POINT в реальную точку (эта точка изначально находится в отрицательной зоне [с отрицательным индексом]) - добавить её в массив точек с положительным индексом
                              /// и сразу выполнить рассчёт её высоты, если такая точка не была определена до этого момента.
                              /// </summary>
                              int get_index_or_append__project_point__vertex__application_counter(
                                int point_index,
                                int profile_face_index,
                                int application_index,
                                FT oioa0_altitude,
                                FT oioa1_altitude,
                                FT oioa0_offset,
                                FT oioa1_offset
                              /*, boost::mutex& mtx_*/
                              ) {
                                // Получить проектную точку
                                //auto& ppoint = map__point_index__calc_point[point_index];
                                auto& project_point = vect__point_index__project_point[-(point_index+1 /*индексы для проектных точек отрицательные и всегда начинаеются с -1, поэтому, чтобы попасть в 0, нужно прибавлять 1 */ )];
                                // Проверить, имеется ли у неё указанное применение:
                                const auto& ss_id__vertex_id = std::tuple(project_point.ss_id, project_point.vertex_id, profile_face_index, application_index);
                                int res_counter = -1;
                                //mtx_.lock();
                                const auto& it = map__ss_id__vertex_id.find(ss_id__vertex_id);
                                if (it == map__ss_id__vertex_id.end()) {
                                  FT z_time = (oioa1_altitude - oioa0_altitude) * (project_point.event_time - oioa0_offset) / (oioa1_offset - oioa0_offset) + oioa0_altitude; // Определить высоту точки. Очень простая формула, для которой и выполнялась вся эта программа! )))
                                  Point_3 p3 = Point_3(project_point.point.x(), project_point.point.y(), z_time); // Создать точку для mesh.
                                  // Добавить её в список результирующих точек
                                  SHAPE_POINT mesh_point(
                                    -1 /*индекс ещё не известен*/,
                                    point_index, /*Из какой проектной точки она получена - для информации по время отладки */
                                    SHAPE_POINT::POINT_TYPE::APPLICATED_PROJECT_POINT, /*Признак преобразованной точки*/
                                    p3,
                                    true /*положение по высоте рассчитано*/,
                                    0.0 /*offset - не имеет значения, т.к. это не offset, а PROJECT_POINT*/,
                                    0   /*_oioa_offset_index - не имеет значения, т.к. это не offset, а PROJECT_POINT*/,
                                    z_time   /*_oioa_altitude - высота этой точки рассчитана*/,
                                    project_point.is_border /*не важно, но сохраню*/,
                                    project_point.ss_id,
                                    project_point.ss_sign /*вроде не важно, но сохраню*/,
                                    project_point.vertex_id /*Используется для поиска, что такая точка типа PROJECT_POINT уже применялась/использовалась*/,
                                    application_index, /* Порядковый номер проектной точки в пределах profile_face_index */
                                    project_point.event_time /*Используется для получения высоты точки, если она находится между двумя offset-ами*/
                                  );
#ifdef _DEBUG
                                  mesh_point.project_point_index = point_index;
#endif
                                  // Создать реальную точку из проектной
                                  vect__point_index__mesh_point.push_back(mesh_point);
                                  res_counter = vect__point_index__mesh_point.size()-1;
                                  // Запомнить её индекс
                                  vect__point_index__mesh_point.back().index = 
                                    map__ss_id__vertex_id[ss_id__vertex_id]  =
                                    res_counter;
#ifdef _DEBUG
                                  // Запомнить индексы точек, в которые была преобразована SS точка PROJECT_POINT
                                  project_point.beveled_indexes.push_back(res_counter);

                                  // Запомнить в новой точке индекс, от которого эта точка получена. Для тестов.
                                  vect__point_index__mesh_point.back().beveled_indexes.push_back(project_point.index);
#endif
                                } else {
                                  // Точка найдена, вернуть её индекс
                                  // TODO: добавить картинку с примером как одна такая точка может встретиться несколько раз.
                                  auto& [tpl, _res_counter] = (*it);
                                  res_counter = _res_counter;
                                }
                                //mtx_.unlock();
                                return res_counter;
                              };

                              SHAPE_POINTS(size_t size_calc_points, size_t size_project_points) {
                                vect__point_index__mesh_point   .reserve(size_calc_points);
                                vect__point_index__project_point.reserve(size_project_points);
                              }
                            };

                            // Определить примерное максимальное количество точек в результирующей фигуре.
                            size_t profile_points_size = 0;
                            size_t ss_id__count = 0;
                            for (auto& [ss_id, oioa_offset_ss_params] : map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_index]) {
                              profile_points_size += oioa_offset_ss_params.size();
                            }
                            // Количество уже известных точек для mesh (это точки пересечения he с SS):
                            size_t intersections_points_size = 0;
                            for (auto& [ss_id, vect__he_intersections_point2 /*относится только к одному ss_id*/] : map__ss_id__he_intersections_point2) {
                              // Сколько пересечений имеет этот offset с в текущем SS (по ss_id)
                              for (auto& elem /*один offset в ss_id*/ : vect__he_intersections_point2) {
                                intersections_points_size += elem.map__halfedge__offset_points /*точки пересечения этого offset в ss_id*/.size();
                              }
                              ss_id__count++;
                            }
                            size_t size_project_points = mesh_ss_merged.number_of_vertices(); // Резервирование памяти под проектные точки (они образуются из проектных точек SS). Эмпирическая величина.
                            // Эмпирическая величина. Попытка заранее предсказать количество точек в mesh, чтобы дополнительное выделение память под std::vector для точек mesh по возможности не происходило.
                            // По идее достаточно трудно предсказать количество точек в результирующем mesh и вполне возможна ситуация, когда этого не хватит. Ну, что ж. Значит std::vector
                            // пусть сам разбирается с этой ситуацией (TODO: может удасться в итоге предсказать предположительное количество точек? Надо подумать при случае, чтобы не выделать слишком много памяти.)
                            size_t size_calc_points = (intersections_points_size); // *2; - агрессивный захват памяти. - update: не заметил разницы во времени со значением, когда нет увеличения количества на 2.
                            // Список точек BEVELED SS, образующихся в результате рассчётов сегментов.
                            // Примечание: Beveled SS строится на базе нескольких faces от профиля, поэтому тут находятся точки в сумме от всего рассчёта. Деление на группы по этим faces будет позже.
                            SHAPE_POINTS calc_points(size_calc_points, size_project_points*2);

                            /// <summary>
                            /// Точки, используемые при обработке плана Project Beveled SS. Сами точки SHAPE_POINT не могут использоваться в плане, потому что
                            /// кроме индекса самой точки надо ещё помнить в каком направлении (slope) алгоритм через неё проходит при обходе faces против часовой стрелки.
                            /// Если SHAPE_POINT на месте пересечения находится на границе между двух faces, то в контур она попадает в двух экземплярах
                            /// для обоих соседних faces (соседние - в общем смысле, а не только смежные в сегменте) с учётом slope (типа UP или DOWN).
                            /// Если SHAPE_POINT типа PROJECT_POINT, то она попадает в контур в одном экземпляре и её slope не учитывается.
                            /// Пример 1: <image url="..\code_images\file_0075.png" scale=".4"/>
                            /// Пример 2: <image url="..\code_images\file_0079.png" scale=".2"/>
                            /// </summary>
                            struct CONTOUR_SHAPE_POINT {
                              int point_index;
                              SHAPE_POINT::POINT_TYPE type;

                              CONTOUR_SHAPE_POINT(int _point_index, SHAPE_POINT::POINT_TYPE _type)
                                :point_index(_point_index), type(_type) {

                              }
                            };

                            // Подобрать данные по offset на каждый сегмент всех SS, но учесть, что если нет отрицательных SS, то все сегменты будут принадлежать соответствующим SS,
                            // но при наличии отрицательных SS требуется построить все сегменты только на отрицательных SS и привязать результаты подбора сегментов к ss_face_id с отрицательными SS. 
                            // Примечание: Если все offset сосредоточены только в положительных или только в отрицательных offset, то map__adjacent_mesh_faces будет пустым.
                            // Пример привязки offset к ss_face_id из отрицательного SS: <image url="..\code_images\file_0063.png" scale=".2"/>

                            // Информация о сегментах для offset с привязкой к ss_id:
                            // <image url="..\code_images\file_0065.png" scale=".2"/>
                            std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector<MeshFacePropertyStruct>>>              map__ss_id__mesh_face_id__segment_faces; // Рассчитать segment_faces
                            std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector<CONTOUR_SHAPE_POINT>>>                 map__ss_id__mesh_face_id__segment_contour; // Преобразовать один сегмент, состоящего из одного или двух faces (если есть отрицательные SS, то их точно будет 2 faces,а если есть только положительные SS, то один face)) в набор точек segment_contour
                            std::map<int /*индекс face профиля (порядковый номер)*/, std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>>> map__profile_face_index__ss_id__mesh_face_id__faces_info; // Состав одного сегмента (сразу привязывается к mesh_face_id), распределённый по входным элементам геометрии профиля (в основном по faces профиля. В профиле несколько faces и других элементов).


                            // Отложенные проектные точки, которые надо рассчитать для mesh
                            struct PROJECT_POINT__APPLICATION_COUNTER {
                              int project_point /*Уникальный в пределах всего объекта*/;
                              int application_counter;
                              int profile_face_index; // индекс face_index профиля
                              PROJECT_POINT__APPLICATION_COUNTER(int _project_point, int _profile_face_index, int _application_counter)
                                :project_point(_project_point), profile_face_index(_profile_face_index), application_counter(_application_counter) {

                              }
                            };
                            // Временное хранилище проектных точек во время рассчёта точек сегмента, чтобы не вызывать mutex каждый раз при чтении любой точки сегмента,
                            // что ужасно понижает производительность с учётом того, что проектные точки попадаются реже точек пересечения (есть исключения, но,
                            // скорее это редкость, чем правило, к тому же это исключение только для частоты встречь таких случаев, а не что это приведёт к замедлению
                            // рассчёта - не приведёт)
                            std::vector<PROJECT_POINT__APPLICATION_COUNTER> vect__stack__project_points;

                            {
                              {

                                // Обработать сегменты каждого SS (сегмент - это отрезок у "подножия" SS, граничащий либо с пустотой вокруг фигуры, если нет отрицательных SS либо с отрицательными SS):
                                for (auto& [ss_id, vect__he_intersections_point2 /*относится только к одному ss_id*/] : map__ss_id__he_intersections_point2) {

                                  // Если есть смежные faces, то для расчёта нужно обрабатывать смежные faces по этой таблице.
                                  // Примечание: если в таблице смежных faces есть записи, то все faces во всех SS имеют смежные faces.
                                  // Если таблица смежных faces пустая, то можно обработать текущий ss_id в одиночку.
                                  // Примечание: если в таблице смежных faces нет записей, то вообще не существует смежных faces (это возможно когда все offset сосредоточены в положительных SS
                                  // или когда все offset сосредоточены в отрицательных offset)
                                  if (map__adjacent_mesh_faces.size() > 0) {
                                    // Проверить, есть ли в таблице смежных faces выбрать первый face выбранного ss_id (достаточно проверить только первый face, потому что все mesh_face_id в ключе будут от одного ss_id)
                                    auto& mesh_face_id = map__ss_id__ss_face_id__face_info[ss_id].begin()->second.mesh_face_id; /* выбирать первый элемент через begin() потому что если попадётся внешний отрицательный контур SS, полученный через виртуальный frame, то нумерация начнётся не с 0, а с 4, поэтому выбираем итератором */
                                    if (map__adjacent_mesh_faces.find(mesh_face_id) == map__adjacent_mesh_faces.end()) {
                                      // Не обрабатывать этот SS, т.к. faces этого SS будут обработаны с помощью смежных faces (либо уже были обработаны).
                                      continue;
                                    }
                                  }

                                  // Рассчитать какие faces надо использовать для сегментов:
                                  for (auto& [ss_face_id, mesh_face_id0_info] : map__ss_id__ss_face_id__face_info[ss_id]) {

                                    // Рассчитать, какие faces нужны для построения набора offset на внешнем ребре SS (по одному face или парами смежных face?).
                                    std::vector<MeshFacePropertyStruct> segment_faces;
                                    segment_faces.push_back(mesh_face_id0_info); // Примечание - таблица map__ajacent_faces не гарантирует, что ключём будут только отрицательные SS, в value - +SS. Может быть и наоборот. Поэтому позже будет произведён реверс списка, чтобы на первом месте была ss_face_id из отрицательной SS.
                                    // Запомнить базовый mesh_face к которому будет привязка offset-ов
                                    MeshFacePropertyStruct& target_ss_face_info = mesh_face_id0_info;
                                    // Определить нужно ли обрабатывать только этот ss_face_id как один face или, если есть смежный face, то уже нужно обработать уже двоих вместе со смежным face?
                                    if (map__adjacent_mesh_faces.size() > 0) {
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
                                          HDS_Halfedge_const_handle hds_h = ss_face_handle->halfedge(), hds_h_start = hds_h, hds_h_finish = hds_h;
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
                                              FT  oioa_offset;
                                              int oioa_offset_index;
                                              FT  oioa_altitude;

                                              /// <summary>
                                              /// Запомнить знак SS этой точки
                                              /// </summary>
                                              MeshFacePropertyStruct::SS_SIGN ss_sign;


                                              /// <summary>
                                              /// Точка пересечения с halfedge
                                              /// </summary>
                                              Point_2 point;

                                              POINT_INFO(FT _oioa_offset, int _oioa_offset_index, FT _oioa_altitude, Point_2& _point, MeshFacePropertyStruct::SS_SIGN _ss_sign)
                                                :oioa_offset(_oioa_offset), oioa_offset_index(_oioa_offset_index), oioa_altitude(_oioa_altitude), point(_point), ss_sign(_ss_sign) {
                                              }
                                            };

                                            // Список offset с которыми произошли пересечения
                                            std::vector<POINT_INFO> vect__intersection_point_info;
                                            // В зависимости от уклона ребра SS и знака SS необходимо правильно определить последовательность добавления точек пересечения с ребром.
                                            // При движении "вверх" при положительном +SS индексы добавляются последовательно как они и объявлены клиентом,
                                            // а при движении вниз - в обратном порядке. В отрицательной области -SS надо делать наоборот - добавлять индексы
                                            // в обратной последовательности при движении (-)"наружу"/"вниз" и добавлять в прямой последовательности при
                                            // движении "вверх"/"внутрь" по направлении линии 0-0: <image url="..\code_images\file_0071.png" scale=".1"/>
                                            // Добавлять точки пересечения с he в направлении обхода контура
                                            if (hds_h_slope == CGAL::NEGATIVE && segment_faces_I.ss_sign == MeshFacePropertyStruct::SS_SIGN::POSITIVE /* +SS */
                                              ||
                                              hds_h_slope == CGAL::POSITIVE && segment_faces_I.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE /* -SS */
                                              ) {
                                              // <image url="..\code_images\file_0076.png" scale=".3"/>
                                              for (auto elem_it = map__ss_id__he_intersections_point2[segment_faces_I.ss_id /* ss_id брать не из цикла, а из рассматриваемой segment_faces_I, чтобы искать пересечения offset-ов со своим he */].rbegin(); elem_it != map__ss_id__he_intersections_point2[segment_faces_I.ss_id].rend(); ++elem_it) {
                                                auto& intersection_points = (*elem_it);
                                                auto offset0_he_intersect_iter = intersection_points.map__halfedge__offset_points.find(hds_offset_h);
                                                if (offset0_he_intersect_iter != intersection_points.map__halfedge__offset_points.end()) {
                                                  auto& [HDS_Halfedge_const_handle, point] = *offset0_he_intersect_iter;
                                                  vect__intersection_point_info.push_back(POINT_INFO(intersection_points.oioa_offset, intersection_points.oioa_offset_index, intersection_points.oioa_altitude, point, segment_faces_I.ss_sign));
                                                }
                                              }
                                            } else {
                                              // <image url="..\code_images\file_0077.png" scale=".3"/>
                                              for (auto elem_it = map__ss_id__he_intersections_point2[segment_faces_I.ss_id].begin(); elem_it != map__ss_id__he_intersections_point2[segment_faces_I.ss_id].end(); ++elem_it) {
                                                auto& intersection_points = *elem_it;
                                                auto offset0_he_intersect_iter = intersection_points.map__halfedge__offset_points.find(hds_offset_h);
                                                if (offset0_he_intersect_iter != intersection_points.map__halfedge__offset_points.end()) {
                                                  auto& [HDS_Halfedge_const_handle, point] = *offset0_he_intersect_iter;
                                                  vect__intersection_point_info.push_back(POINT_INFO(intersection_points.oioa_offset, intersection_points.oioa_offset_index, intersection_points.oioa_altitude, point, segment_faces_I.ss_sign));
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
                                                  // <image url="..\code_images\file_0062.png" scale=".5"/> 
                                                  // Нашлось исключение: <image url="..\code_images\file_0072.png" scale=".3"/> (Надо будет проверить может ли оно влиять на рассчёт?)
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
                                                      default:
                                                        break;
                                                    }
                                                  }
                                                  // TODO - так-то их можно сразу обрабатывать на месте обнаружения. Подумать, как сделать это прямо в точке обнаружения выше.
                                                  // Сохранить все точки пересечения в векторе точек контура текущего сегмента.
                                                  for (auto& he_offset : vect__intersection_point_info) {
                                                    Point_3 p3(he_offset.point.x(), he_offset.point.y(), he_offset.oioa_altitude);
                                                    int point_index = calc_points.get_index_or_append_he__offset_index(
                                                      segment_faces_I.ss_id, //ss_id,
                                                      pt,
                                                      hds_offset_h,
                                                      he_offset.oioa_offset_index,
                                                      p3,
                                                      true,
                                                      he_offset.oioa_offset,
                                                      he_offset.oioa_offset_index,
                                                      he_offset.oioa_altitude,
                                                      he_offset.ss_sign
                                                    );
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
                                                // <image url="..\code_images\file_0074.png" scale=".2"/>
                                              } else {
                                                int he_vertex_id = hds_h->vertex()->id();
                                                auto he_vertex_time = hds_h->vertex()->time();
                                                if (segment_faces_I.ss_sign == MeshFacePropertyStruct::SS_SIGN::NEGATIVE) {
                                                  he_vertex_time = -he_vertex_time;
                                                }
                                                bool is_border = hds_h->is_border();
                                                Point_3 p3(hds_h->vertex()->point().x(), hds_h->vertex()->point().y(), 0 /*Z не выставляется (это altitude), будет считаться позже*/);
                                                SHAPE_POINT::POINT_TYPE pt = SHAPE_POINT::POINT_TYPE::PROJECT_POINT;
                                                int point_index = calc_points.get_index_or_append_vertex_id(
                                                  segment_faces_I.ss_id, // Не всегда соответствует ss_id.
                                                  segment_faces_I.ss_sign,
                                                  he_vertex_id,
                                                  pt,
                                                  p3,
                                                  false /*altitude не рассчитана, будет считаться позже*/,
                                                  0.0,
                                                  -1 /*Это точка контура, поэтому offset-а у неё нет*/,
                                                  0.0,
                                                  is_border,
                                                  -1 /*_application_index для точки пересечения he не важен */,
                                                  he_vertex_time
                                                );
                                                vect_contour_shape_point.push_back(CONTOUR_SHAPE_POINT(point_index, pt));
                                              }
                                            }
                                            vect_segment__halfedges_handle.back().push_back(hds_offset_h);
                                            face_slopes.push_back(hds_h->slope()); // Предположительно это указывает направление движения к или от центра SS. update - нет, если представить себе, что SS является объёмной фигурой типа гор с вершинами, то edges являются ребрами на пиках, а slope указывает как смещается это полуребро по вертикали -1 - вниз, 0 - горизонтальное, +1 - вверх. Соответственно противоположное полуребро является инверсией по направлению.
                                            hds_h = hds_h->next();
                                          } while (hds_h != hds_h_finish);

                                          // Проверка, что наклонные slope() на противоположных (opposite()) рёбрах ожидаются противоположных знаков (кроме 0 - инверсия от 0 - всегда 0).
                                          // <image url="..\code_images\file_0061.png" scale=".1"/>
#ifdef _DEBUG
                                          if (verbose) {
                                            printf("\n " VAL2STR(Msg::_0034) " ss_id=% 2u, ss_face_id % 2u, mesh_face_id % 2u, slopes: ", segment_faces_I.ss_id, segment_faces_I.ss_face_id, segment_faces_I.mesh_face_id);
                                            for (auto& slope : face_slopes) {
                                              printf("% 2d ", slope);
                                            }
                                          }
#endif
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

                                // На этот момент получены списки offset-ов на всех сегментах. Все offset распределены по индексам относительно линии 0-0 и знаку SS. Надо создать faces.
                                // faces создаются на основе offset_index

                                // Таймер для определения времени, потраченного на параллельный рассчёт
                                CGAL::Real_timer timer1;
                                timer1.start();
                                if (verbose) {
                                  printf("\n " VAL2STR(Msg::_0076) ". SS 2D Offset. calc faces of segments. ------------------------------------------------");
                                }
                                // Сколько раз проектные точки преобразовывались в точки mesh?
                                int project_point__to__mesh__counter = 0;
                                {

                                  // Параметры курсора для обхода точек
                                  struct COLLECT_CURSOR {
                                    /// <summary>
                                    /// (true) Собирать или (false) пропускать точки тип PROJECT_POINT
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
                                    /// Индекс для второго пересечения, после которого надо добавить новую точку PROJECT_POINT (+ добавить эту точку)
                                    /// </summary>
                                    int offset_index1;

                                    /// <summary>
                                    /// Тип второго пересечения, после которого нужно собирать точки и добавлять их в последний face
                                    /// </summary>
                                    SHAPE_POINT::POINT_TYPE offset_type1;

                                    COLLECT_CURSOR(bool _do_collect_points, bool _do_create_new_faces, bool _do_reverse, int _index_offset0, SHAPE_POINT::POINT_TYPE _offset_type0, int _index_offset1, SHAPE_POINT::POINT_TYPE _offset_type1)
                                      :do_collect_points(_do_collect_points), do_create_new_faces(_do_create_new_faces), do_reverse(_do_reverse), offset_index0(_index_offset0), offset_type0(_offset_type0), offset_index1(_index_offset1), offset_type1(_offset_type1) {

                                    }
                                  };


                                  /// <summary>
                                  /// Параметры offset-ов между которыми производится рассчёт
                                  /// </summary>
                                  struct EDGE_INFO {
                                    /// <summary>
                                    /// Является ли первый offset caps-ом? true - является, false - не является
                                    /// </summary>
                                    bool is_oioa0_offset_caps;

                                    /// <summary>
                                    /// index первого offset
                                    /// </summary>
                                    int oioa0_offset_index;

                                    /// <summary>
                                    /// величина первого offset
                                    /// </summary>
                                    FT  oioa0_offset;

                                    /// <summary>
                                    /// величина первого altitude
                                    /// </summary>
                                    FT  oioa0_altitude;

                                    /// <summary>
                                    /// Является ли второй offset caps-ом? true - является, false - не является
                                    /// </summary>
                                    bool is_oioa1_offset_caps;

                                    /// <summary>
                                    /// index второго offset
                                    /// </summary>
                                    int oioa1_offset_index;

                                    /// <summary>
                                    /// величина второго offset
                                    /// </summary>
                                    FT  oioa1_offset;

                                    /// <summary>
                                    /// величина второго 
                                    /// </summary>
                                    FT  oioa1_altitude;

                                    EDGE_INFO(int _oioa0_offset_index, FT _oioa0_offset, FT _oioa0_altitude, bool _is_oioa0_offset_caps, int _oioa1_offset_index, FT _oioa1_offset, FT _oioa1_altitude, bool _is_oioa1_offset_caps)
                                      : oioa0_offset_index(_oioa0_offset_index), oioa0_offset(_oioa0_offset), oioa0_altitude(_oioa0_altitude), is_oioa0_offset_caps(_is_oioa0_offset_caps), oioa1_offset_index(_oioa1_offset_index), oioa1_offset(_oioa1_offset), oioa1_altitude(_oioa1_altitude), is_oioa1_offset_caps(_is_oioa1_offset_caps) {

                                    }
                                  };

                                  auto& map__ss_id__OIOA_OFFSET_SS_PARAMS = map__object_id__ss_id__OIOA_OFFSET_SS_PARAMS[object_index];
                                  if (map__ss_id__OIOA_OFFSET_SS_PARAMS.size() == 0) {
                                    // Если у объекта не рассчитан ss, то и offsets у него тоже не могут быть посчитаны. mesh у такого объекта будет пустой. Это вполне допустимо,
                                    // по аналогии когда у объекта в результате расчётов не образуется mesh (например, профиль слишком глубоко в shape).
                                  } else {
                                    auto object_original_index = map__object_index__poms[object_index][0 /*Тут может быть один или два ss_id в зависимости от наличия отрицательных offsets, но 0-й будет всегда, поэтому его и берём*/ ].object_original_index;
                                    auto& vect__offset_faces = map__object_index__profile_faces[object_original_index];

                                    // Набор индексов offset-ов;
                                    std::set<int> set__offset_index;
                                    std::vector<OIOA> vect_object_offsets;
                                    std::vector<std::vector<EDGE_INFO>> vect__faces__edges;
                                    {
                                      // Список информации об offsets текущего объекта
                                      for (auto& [ss_id, oioa_offset_ss_params] : map__ss_id__OIOA_OFFSET_SS_PARAMS) {
                                        for (auto& offset_params : oioa_offset_ss_params) {
                                          if (set__offset_index.find(offset_params.oioa1.offset_index) == set__offset_index.end()) {
                                            set__offset_index.insert(offset_params.oioa1.offset_index);
                                            vect_object_offsets.push_back(offset_params.oioa1);
                                          }
                                        }
                                      }

                                      // Отсортировать offset-ы в порядке как они задавались (по индексам):
                                      std::sort(vect_object_offsets.begin(), vect_object_offsets.end(), [](const OIOA& o1, const OIOA& o2) {
                                        bool res = o1.offset_index < o2.offset_index;
                                        return res;
                                        });
                                    }
                                    
                                    std::vector<int> vect_profile_faces__open_modes = map__object_index__profile_faces__open_modes[object_original_index];
                                    {
                                      int vect_object_offsets_size = vect_object_offsets.size();
                                      int face_index = 0;
                                      for (auto& face : vect__offset_faces) {
                                        if (face.size() == 0) {
                                          face_index++;
                                          // Если face не содержит индексов, то пропустить. Возможно, что сокет faces вообще не подключен, тогда
                                          // и не добавлять вектор vect__edges в vect__faces__edges. Если 
                                          continue;
                                        } else {
                                          int face_open_mode = vect_profile_faces__open_modes[face_index++];
                                          std::vector<EDGE_INFO> vect__edges;
                                          face_open_mode = face_open_mode == 1 ? 1 : face_open_mode == 2 ? 2 : 0; // Простенький фильтр на допустимые значения
                                          if (face_open_mode == 2 /*INPAIRS*/) {
                                            // Соеденить парами
                                            for (int I = 1; I <= (int)face.size() - 1; I += 2) {
                                              auto& start_index = face[I - 1];
                                              auto& end_index = face[I];
                                              if (start_index == end_index) {
                                                continue;
                                              }
                                              // Простая проверка на допустимый диапазон:
                                              if (0 <= start_index && 0 <= end_index && start_index <= vect_object_offsets_size - 1 && end_index <= vect_object_offsets_size - 1) {
                                                vect__edges.push_back(
                                                  EDGE_INFO(
                                                    vect_object_offsets[start_index].offset_index,
                                                    vect_object_offsets[start_index].offset,
                                                    vect_object_offsets[start_index].altitude,
                                                    false,
                                                    vect_object_offsets[end_index].offset_index,
                                                    vect_object_offsets[end_index].offset,
                                                    vect_object_offsets[end_index].altitude,
                                                    false));
                                              }
                                            }
                                          } else if (face_open_mode == 1 /*CLOSED*/ || face_open_mode == 0 /*OPENED*/) {
                                            // Если задано замкнуть face, то добавить в начало вектора замыкающий edge:
                                            if (face_open_mode == 1 && face.size() > 2) {
                                              int last_index = face.size() - 1;
                                              int first_index = 0;
                                              if (first_index == last_index) {
                                                // Если Индексы выходят за допустимый диапазон, то не делать соединения краёв в любом случае (TODO: по хорошему надо информировать пользователя ошибкой)
                                              } else {
                                                vect__edges.push_back(EDGE_INFO(
                                                  vect_object_offsets[face[last_index]].offset_index,
                                                  vect_object_offsets[face[last_index]].offset,
                                                  vect_object_offsets[face[last_index]].altitude,
                                                  false,
                                                  vect_object_offsets[face[first_index]].offset_index,
                                                  vect_object_offsets[face[first_index]].offset,
                                                  vect_object_offsets[face[first_index]].altitude,
                                                  false)
                                                );
                                              }
                                            }
                                            // Просканировать остальные индексы вокруг face:
                                            for (int I = 1; I <= (int)face.size() - 1; I++) {
                                              auto& start_index = face[I - 1];
                                              auto& end_index = face[I];
                                              // Проверять, что такие offset-ы есть. Если одного из offset-ов нет, то пропускать такое ребро (edge):
                                              if (0 <= start_index && 0 <= end_index && start_index <= vect_object_offsets_size - 1 && end_index <= vect_object_offsets_size - 1) {
                                                vect__edges.push_back(
                                                  EDGE_INFO(
                                                    vect_object_offsets[start_index].offset_index,
                                                    vect_object_offsets[start_index].offset,
                                                    vect_object_offsets[start_index].altitude,
                                                    false,
                                                    vect_object_offsets[end_index].offset_index,
                                                    vect_object_offsets[end_index].offset,
                                                    vect_object_offsets[end_index].altitude,
                                                    false)
                                                );
                                              }
                                            }
                                          }
                                          vect__faces__edges.push_back(vect__edges);
                                        }
                                      }
                                    }

                                    // Если данные по faces не были переданы (через сокет profile faces <image url="..\code_images\file_0090.png" scale=".3"/>), то соеденить контур по порядку переданных вершин (без замыкания):
                                    if (vect__faces__edges.size() == 0) {
                                      std::vector<EDGE_INFO> vect__edges;
                                      for (int I = (int)vect_object_offsets.size() - 2; I >= 0; I--) {
                                        int firts_index = I;
                                        int last_index = I + 1;
                                        vect__edges.push_back(
                                          EDGE_INFO(
                                            vect_object_offsets[last_index].offset_index,
                                            vect_object_offsets[last_index].offset,
                                            vect_object_offsets[last_index].altitude,
                                            false,
                                            vect_object_offsets[firts_index].offset_index,
                                            vect_object_offsets[firts_index].offset,
                                            vect_object_offsets[firts_index].altitude,
                                            false
                                          )
                                        );
                                      }
                                      vect__faces__edges.push_back(vect__edges);
                                    }


                                    for (int I = 0; I <= (int)vect__faces__edges.size() - 1; I++) {
                                      int profile_face_index = I;
                                      map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index] = std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>>();
                                      auto& vect__edges = vect__faces__edges[I];
                                      if (vect__edges.size() >= 2) {
                                        for (const auto& [ss_id, _] : map__ss_id__mesh_face_id__segment_contour) {
                                          auto& pfi_face_info = map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index];
                                          if (pfi_face_info.find(ss_id) == pfi_face_info.end()) {
                                            pfi_face_info[ss_id] = std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>();
                                          }
                                        }
                                      }
                                    }

                                    // <image url="..\code_images\file_0082.png" scale=".3"/> В некоторых случаях производительность поднимается очень хорошо.
                                    const int threadNumbers = boost::thread::hardware_concurrency();
                                    boost::asio::thread_pool pool1(threadNumbers);
                                    boost::mutex mtx_;


                                    for (const auto& [_ss_id, _map__mesh_face_id__segment_points] : map__ss_id__mesh_face_id__segment_contour) {
                                      // Без этого преобразования возникают проблемы с передачей параметров в многопоточность в clang.
                                      auto& ss_id = _ss_id;;
                                      auto& map__mesh_face_id__segment_points = _map__mesh_face_id__segment_points;
                                      auto& _object_index = object_index;

                                      // Рассчёт одного SS (известно его ss_id). TODO: в случае, когда имеется отрицательный SS, то он очень сильно перетягивает на себя загрузку, потому что
                                      // часто отрицательный SS остаётся единственным, обрабатываемым SS (или их мало в любом случае), т.к. является внешним. Бывает, что к нему подключаются и другие отрицательные SS, 
                                      // Но такое распараллеливание не сильно повышает производительность. // <image url="..\code_images\file_0092.png" scale=".2"/>
                                      boost::asio::post(pool1, [
                                        ss_id,
                                        &map__profile_face_index__ss_id__mesh_face_id__faces_info,
                                        &map__mesh_face_id__segment_points,
                                        &vect_object_offsets,
                                        &vect__faces__edges,
                                        &calc_points,
                                        verbose,
                                        _object_index,
                                        &vect__stack__project_points,
                                        &project_point__to__mesh__counter,
                                        &mtx_
                                      ] {
                                        // Таймер для определения времени рассчёта одного потока.
                                        CGAL::Real_timer timer1;
                                        timer1.start();
                                        //CGAL::Real_timer timer7;
                                        //double get_points7 = 0;
                                        {
                                          {
                                            SHAPE_POINT sp_default;
                                            for (int I = 0; I <= (int)vect__faces__edges.size() - 1; I++) {
                                              int profile_face_index = I;
                                              auto& vect__edges = vect__faces__edges[I];
                                              if (vect__edges.size() >= 1 /*Для расчёта требуется, чтобы хоть одно ребро было*/) {
                                                //mtx_.lock();
                                                //if (map__profile_face_index__ss_id__mesh_face_id__faces_info.find(profile_face_index) == map__profile_face_index__ss_id__mesh_face_id__faces_info.end()) {
                                                //  map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index] = std::map<int /*ss_id*/, std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>>();
                                                //}
                                                //if (map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index].find(ss_id) == map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index].end()) {
                                                //  map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index][ss_id] = std::map<int /*mesh_face_id*/, std::vector/*faces*/<FACE_INFO>>();
                                                //}
                                                //mtx_.unlock();

                                                for (auto& [mesh_face_id, vect_segment_points] : map__mesh_face_id__segment_points) {
                                                  // Считать только если есть минимум пара offset
                                                  // TODO - в будущем учесть, что высоты могут быть одинаковыми и caps-ы могут перекрыться. Пока это не учитывается.

                                                  ////vect_edges.push_back(EDGE_INFO(0, 0.0, true, vect_oioa[0].offset_index, vect_oioa[0].offset_index, true)); // Первый индекс будет с caps-ом - update - пока пропускаю
                                                  //for (int I = 0; I <= (int)vect_object_offsets.size() - 2; I++) {
                                                  //  vect_edges.push_back(EDGE_INFO(vect_object_offsets[I + 0].offset_index, vect_object_offsets[I + 0].offset, vect_object_offsets[I + 0].altitude, false, vect_object_offsets[I + 1].offset_index, vect_object_offsets[I + 1].offset, vect_object_offsets[I + 1].altitude, false));
                                                  //}
                                                  ////vect_edges.push_back(EDGE_INFO(vect_oioa[0].offset_index, vect_oioa[0].offset_index, false, 0, 0.0, true)); // Последний индекс будет с caps-ом - update - пока пропускаю

                                                  std::vector<FACE_INFO> vect_faces /*список faces на текущем сегменте*/; // = std::vector/*faces*/<FACE_INFO>();
                                                  std::vector<PROJECT_POINT__APPLICATION_COUNTER> _vect__stack__project_points; // Тоже самое, только в пределах одного потока. Потом пересчитать индексы faces для глобального параметра vect__stack__project_points 

                                                  // Счётчик сколько раз встретилась проектная точка при рассчёте этого face.
                                                  // Такой подход использует симметричность получения параметра map__point_index__counter[calc_point.index. Этот параметр при
                                                  // проходе контура по offset вычисляется одинаково для одной и той же точки, т.к. эта точка используется одинаковое количество раз
                                                  // во всех смежных faces в этой точке.
                                                  // TODO: Сделать рисунок или анимацию.
                                                  std::map<int, int> map__point_index__counter;

                                                  for (auto& edge_info : vect__edges) {

                                                    auto& _oioa0_offset_index = edge_info.oioa0_offset_index;
                                                    auto& _oioa0_offset = edge_info.oioa0_offset;
                                                    auto& _oioa0_altitude = edge_info.oioa0_altitude;
                                                    auto& _oioa1_offset_index = edge_info.oioa1_offset_index;
                                                    auto& _oioa1_offset = edge_info.oioa1_offset;
                                                    auto& _oioa1_altitude = edge_info.oioa1_altitude;

                                                    // Сначала определить какой сегмент по абсолютному значению ближе к линии 0-0?
                                                    //mtx_.lock();
                                                    auto& elem0 = vect_object_offsets[_oioa0_offset_index];
                                                    auto& elem1 = vect_object_offsets[_oioa1_offset_index];
                                                    auto  elem0_offset = elem0.offset;
                                                    auto  elem0_offset_index = elem0.offset_index;
                                                    auto  elem0_altitude = elem0.altitude;
                                                    auto  elem1_offset = elem1.offset;
                                                    auto  elem1_offset_index = elem1.offset_index;
                                                    auto  elem1_altitude = elem1.altitude;
                                                    //mtx_.unlock();

                                                    // Подсчёт faces всегда выполняется против часовой стрелки от первой точки после старта.
                                                    COLLECT_CURSOR cursor(false, false, false, -1, SHAPE_POINT::POINT_TYPE::PROJECT_POINT, -1, SHAPE_POINT::POINT_TYPE::PROJECT_POINT);
                                                    // Рассчёт зависит от того в каком отношении к 0-0 находятся рассматриваемые offset. Мысленно представляем куда движется стерка обхода и
                                                    // в какой последовательности она должна пересечь offset-ы.
                                                    if (elem0_offset < 0 && elem1_offset < 0) {
                                                      // оба offset ниже 0-0
                                                      if (elem0_offset > elem1_offset) {
                                                        // 0
                                                        // 1
                                                        cursor = COLLECT_CURSOR(false /*collect_points*/, true /*create new face*/, true /*do_reverse*/, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                                      } else if (elem1_offset > elem0_offset) {
                                                        // 1
                                                        // 0
                                                        cursor = COLLECT_CURSOR(false, true, false, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                                      } else {
                                                        /* остальные, когда elem0_offset == elem1_offset */
                                                        bool resc = false;
                                                        // Для равных offset всё равно нужно определить направление обхода. Сначала сравнить altitude:
                                                        // Важно, чтобы эти условия совпадали с условиями сортировки правила пересечения offset-ов с he
                                                        if (elem0_altitude != elem1_altitude) {
                                                          resc = (elem0_altitude < elem1_altitude);
                                                        } else {
                                                          // Индексы никогда не равны друг другу. Будем считать, что при обходе по часовой стрелке контура faces segment больший индекс должен быть первым.
                                                          resc = (elem0_offset_index > elem1_offset_index);
                                                        }

                                                        if (resc == false) {
                                                          // 1
                                                          // 0
                                                          cursor = COLLECT_CURSOR(false, true, true, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                                        } else {
                                                          // 0
                                                          // 1
                                                          cursor = COLLECT_CURSOR(false, true, false, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                                        }
                                                      }
                                                    } else if (elem0_offset >= 0 && elem1_offset >= 0) {
                                                      // оба offset выше 0-0
                                                      if (elem0_offset > elem1_offset) {
                                                        // 0
                                                        // 1
                                                        cursor = COLLECT_CURSOR(false, true, true, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                                      } else if (elem1_offset > elem0_offset) {
                                                        // 1
                                                        // 0
                                                        cursor = COLLECT_CURSOR(false, true, false, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                                      } else {
                                                        /* остальные, когда elem0_offset == elem1_offset */
                                                        // Для равных offset всё равно нужно определить направление обхода. Сначала сравнить altitude:
                                                        // Важно, чтобы эти условия совпадали с условиями сортировки правила пересечения offset-ов с he
                                                        bool res = false;
                                                        if (elem0_altitude != elem1_altitude) {
                                                          res = (elem0_altitude < elem1_altitude);
                                                        } else {
                                                          // Индексы никогда не равны друг другу. Будем считать, что при обходе по часовой стрелке контура faces segment больший индекс должен быть первым.
                                                          res = (elem0_offset_index > elem1_offset_index);
                                                        }
                                                        if (res == false) {
                                                          // 0
                                                          // 1
                                                          cursor = COLLECT_CURSOR(false, true, true, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                                        } else {
                                                          // 1
                                                          // 0
                                                          cursor = COLLECT_CURSOR(false, true, false, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN);
                                                        }
                                                      }
                                                    } else if (elem0_offset < 0 && elem1_offset >= 0) {
                                                      // Добавить новый face и разрешить добавлять в него точки типа PROJECT_POINT сразу
                                                      cursor = COLLECT_CURSOR(true, false, false, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                                      vect_faces.push_back(FACE_INFO(_oioa0_offset_index, _oioa0_offset, _oioa0_altitude, _oioa1_offset_index, _oioa1_offset, _oioa1_altitude, cursor.do_reverse));
                                                    } else if (elem0_offset >= 0 && elem1_offset < 0) {
                                                      // Добавить новый face и разрешить добавлять в него точки типа PROJECT_POINT сразу
                                                      cursor = COLLECT_CURSOR(true, false, true, _oioa0_offset_index, SHAPE_POINT::POINT_TYPE::DOWN, _oioa1_offset_index, SHAPE_POINT::POINT_TYPE::UP);
                                                      vect_faces.push_back(FACE_INFO(_oioa0_offset_index, _oioa0_offset, _oioa0_altitude, _oioa1_offset_index, _oioa1_offset, _oioa1_altitude, cursor.do_reverse));
                                                    } else {
                                                      // TODO: Осталось, когда одна из точек в нуле. Проверить, может быть переделать условие "elem0_offset > 0 && elem1_offset > 0" в "elem0_offset >= 0 && elem1_offset >= 0" ? - update - сработало. Пока можно не делать этот раздел. Оставлю на более подробное тестирование.
                                                    }

                                                    for (auto& contour_point : vect_segment_points) {

                                                      //timer7.start();

                                                      //mtx_.lock();
                                                      // Только что узнал, что std::map не является потокобезопасным не только на запись, но и на чтение.
                                                      // Поэтому все операции чтения/записи при многопоточной работе надо оборачивать mutex-ом.
                                                      // Пример сбоя без mutex: <image url="..\code_images\file_0084.png" scale=".2"/>
                                                      // Этот map используется и при рассчёте/получении нового индекса для проектных точек. Оказалось, что иногда идёт попытка чтения в момент добавления
                                                      // объекта в этот map из другого потока! А один раз даже возникла блокировка и blender повис. Хорошо ещё догадался войти отладчиком, а не сбросить
                                                      // blender. Самое обидное - очень редко вылезает. Много тестов прошло без проблем!
                                                      // Update. Более того, обнаружилось, что операция чтения идёт очень медленно! И именно она занимает больше всего времени - даже больше,
                                                      // чем время рассчёта! Заменил на std::unordered_map, но это подняло производительность излечения элемента только в 2 раза.
                                                      // Update: заменил map на вектор. Отрицательные индексы хранятся в векторе проектных точек

                                                      auto& calc_point = contour_point.point_index >= 0 ? calc_points.vect__point_index__mesh_point[contour_point.point_index] : calc_points.vect__point_index__project_point[-(contour_point.point_index + 1)];

                                                      //mtx_.unlock();
                                                      //timer7.stop();
                                                      //get_points7 += timer7.time();

                                                      if (contour_point.type == SHAPE_POINT::POINT_TYPE::PROJECT_POINT) {
                                                        if (cursor.do_collect_points == true) {
                                                          // TODO: Надо подумать, как избежать такого артифакта при shade_smooth <image url="..\code_images\file_0083.png" scale=".5"/>. Это происходит из-за попадания 0-й точки в отрезок.

                                                          // Определить индекс для точки типа PROJECT_POINT на основе перекрытия. Если такой точки нет, то создать её и вернуть её индекс.
                                                          int application_index = ++map__point_index__counter[calc_point.index]; // Используется свойство map, что если элемента в map ещё не было с этим ключём, то по умолчанию в int создаётся 0.
                                                          _vect__stack__project_points.push_back(PROJECT_POINT__APPLICATION_COUNTER(contour_point.point_index, profile_face_index, application_index));
                                                          //mtx_.lock();
                                                          // Оставлю для истории, какая ошибка появилась в алгоритме с появлением profile_face_index: <image url="..\code_images\file_0085.png" scale=".1"/>
                                                          //int point_index = calc_points.get_index_or_append__project_point__vertex__application_counter(calc_point.index, profile_face_index, application_index, mtx_);
                                                          //mtx_.unlock();
                                                          //vect_faces.back().face_verts_indexes.push_back(point_index);
                                                          vect_faces.back().face_verts_indexes.push_back(-_vect__stack__project_points.size()); // Отрицательные индексы указывают, что в данном месте требуется проектная точка (которую надо рассчитать позже)
                                                        }
                                                        continue;
                                                      }

                                                      if (calc_point.oioa_offset_index == cursor.offset_index0 && contour_point.type == cursor.offset_type0) {
                                                        if (cursor.do_create_new_faces == true) {
                                                          // Создать новый face и запомнить эту точку:
                                                          vect_faces.push_back(FACE_INFO(_oioa0_offset_index, _oioa0_offset, _oioa0_altitude, _oioa1_offset_index, _oioa1_offset, _oioa1_altitude, cursor.do_reverse));
                                                        }
                                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                                        cursor.do_collect_points = true;
                                                      } else if (calc_point.oioa_offset_index == cursor.offset_index0 && contour_point.type != cursor.offset_type0) {
                                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                                        cursor.do_collect_points = false;
                                                      } else if (calc_point.oioa_offset_index == cursor.offset_index1 && contour_point.type == cursor.offset_type1) {
                                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                                        cursor.do_collect_points = true;
                                                      } else if (calc_point.oioa_offset_index == cursor.offset_index1 && contour_point.type != cursor.offset_type1) {
                                                        vect_faces.back().face_verts_indexes.push_back(calc_point.index);
                                                        cursor.do_collect_points = false;
                                                      }
                                                    }
                                                  }
                                                  // Проверить vect_faces. Если там были faces с признаком reverse, то инвертировать порядок индексов:
                                                  for (auto& face_info : vect_faces) {
                                                    if (face_info.do_reverse == false /*Так получилось эмпирически)))*/) {
                                                      std::reverse(face_info.face_verts_indexes.begin(), face_info.face_verts_indexes.end());
                                                    }
                                                  }
                                                  mtx_.lock();
                                                  //// Выполнить рассчёт проектных точек
                                                  //std::vector<int> vect__replace_project_point_index__into__mesh_point_index;
                                                  //for (auto& elem : vect__stack__project_points) {
                                                  //  int point_index = calc_points.get_index_or_append__project_point__vertex__application_counter(elem.project_point, profile_face_index, elem.application_counter, mtx_);
                                                  //  vect__replace_project_point_index__into__mesh_point_index.push_back(point_index);
                                                  //}
                                                  //for (auto& face_info : vect_faces) {
                                                  //  std::vector/*point indexes*/<int> face_verts_indexes;
                                                  //  for (auto& point_index : face_info.face_verts_indexes) {
                                                  //    int new_point_index = point_index;
                                                  //    if (point_index < 0) {
                                                  //      new_point_index = vect__replace_project_point_index__into__mesh_point_index[-(point_index+1/*были только отрицательные индексы, начинались от -1, теперь они стали положительными, начинаясь от 0, поэтому к старому отрицательному индексу надо добавить 1, чтобы индексация вернулась от 0*/)];
                                                  //      project_point__to__mesh__counter++;
                                                  //    }
                                                  //    face_verts_indexes.push_back(new_point_index);
                                                  //  }
                                                  //  face_info.face_verts_indexes = face_verts_indexes;
                                                  //}
                                                  // Пересчитать индексы локальных проектных точек на индексы в глобальном списке:
                                                  for (auto& face_info : vect_faces) {
                                                    std::vector/*point indexes*/<int> face_verts_indexes;
                                                    for (auto& point_index : face_info.face_verts_indexes) {
                                                      int new_point_index = point_index;
                                                      if (point_index < 0) {
                                                        auto& project_point__to__mesh = _vect__stack__project_points[-(point_index + 1)];
                                                        vect__stack__project_points.push_back(project_point__to__mesh);
                                                        new_point_index = -vect__stack__project_points.size(); // Пересчёт индекса локальной project_point в глобальный (позже именно по нему будет пересчёт проектной точки в точку mesh)
                                                        project_point__to__mesh__counter++;
                                                      }
                                                      face_verts_indexes.push_back(new_point_index);
                                                    }
                                                    face_info.face_verts_indexes = face_verts_indexes;
                                                  }

                                                  map__profile_face_index__ss_id__mesh_face_id__faces_info[profile_face_index][ss_id][mesh_face_id] = vect_faces;
                                                  mtx_.unlock();
                                                }
                                              }
                                            }
                                          }
                                        }
                                        timer1.stop();
                                        if (verbose == true) {
                                          printf("\n " VAL2STR(Msg::_0052) ". SS 2D Offset. calc faces of segments. || object_index: % 3d, ss_id=% 3d, calc time: % 10.5f.", _object_index, ss_id, timer1.time());
                                        }
                                      }
                                      );
                                    }
                                    pool1.join();
                                  }
                                }
                                timer1.stop();
                                if (verbose == true) {
                                           printf("\n " VAL2STR(Msg::_0053) ". SS 2D Offset. calc faces of segments.                               general time: % 10.5f, projects to mesh: %4d", timer1.time(), project_point__to__mesh__counter);
                                }
#ifdef _DEBUG
                                //if (verbose) {
                                //	printf("\n ====================== LIST OF CONTOURS of object_index %u: =========================\n", object_index);
                                //	for (auto& [ss_id, map__mesh_face_id__segment_points] : map__ss_id__mesh_face_id__segment_contour) {
                                //		for (auto& [mesh_face_id, vect_segment_points] : map__mesh_face_id__segment_points) {
                                //			printf("\n " VAL2STR(Msg::_0035) " contour: ss_id %u, mesh_face_id %u, length %u:", ss_id, mesh_face_id, (int)vect_segment_points.size());
                                //			int counter = 0;
                                //			for (auto& contour_point : vect_segment_points) {
                                //				auto& calc_point = calc_points.map__point_index__calc_point[contour_point.point_index];

                                //				// вывод параметров точек текущего контура
                                //				printf("\n  % 2.8lf, % 2.8lf, % 2.8lf, type: %s, alt: %s:% 2.8lf, offset_index: % 3d, offset: % 2.8lf; % 3u / % 4d",
                                //					CGAL::to_double(calc_point.point.x()), CGAL::to_double(calc_point.point.y()), CGAL::to_double(calc_point.point.z()),
                                //					contour_point.type == 0 ? "PROJECT_POINT" : contour_point.type == 1 ? "   UP" : " DOWN",
                                //					calc_point.is_altitude_calc == true ? "    set" : "not set",
                                //					CGAL::to_double(calc_point.oioa_altitude), calc_point.oioa_offset_index /*-1 если тип точки PROJECT_POINT, то offset-а у неё нет*/, CGAL::to_double(calc_point.oioa_offset), counter++, contour_point.point_index
                                //				);

                                //				// Вывести список индексов результирующих точек
                                //				if (calc_point.beveled_indexes.size() > 0) {
                                //					// Вывести список в кого преобразована эта точка:
                                //					printf(" => [");
                                //					for (auto& idx : calc_point.beveled_indexes) {
                                //						printf(" % 3u", idx);
                                //					}
                                //					printf("]");
                                //				}
                                //			}
                                //			// Вывести список индексов в которые была преобразована точка, если она типа PROJECT_POINT

                                //			// Вывод только тех точек, которые образуют faces:
                                //			int face_id = 0;
                                //			for (auto& face_info : map__profile_face_index__ss_id__mesh_face_id__faces_info[ss_id][mesh_face_id]) {
                                //				printf("\n      face_id % 3u offsets: (% 2d, % 2d : [% 2.8lf, % 2.8lf; % 2.8lf, % 2.8lf], reversed: %s); len: (% 3d): ", face_id, face_info.oioa0_offset_index, face_info.oioa1_offset_index, CGAL::to_double(face_info.oioa0_offset), CGAL::to_double(face_info.oioa0_altitude), CGAL::to_double(face_info.oioa1_offset), CGAL::to_double(face_info.oioa1_altitude), face_info.do_reverse == true ? "yes" : " no", (int)face_info.face_verts_indexes.size());
                                //				for (auto& vert_index : face_info.face_verts_indexes) {
                                //					printf(" % 3u", vert_index);
                                //				}
                                //				face_id++;
                                //			}
                                //		}
                                //	}
                                //	printf("\n ====================== /LIST OF CONTOURS of object_index %u: =========================\n", object_index);
                                //}
#endif
                              }
                            }

                            {
                              CGAL::Real_timer timer1; // Общее время рассчёта всех SS (вместе с загрузкой данных)
                              timer1.start();
                              // Выполнить рассчёт точек, для которых ещё не определена высота/положение:
                              // update: рассчёт удалось выполнить за один проход благодаря параметру shape_point.event_time.
                              for (auto& [profile_face_index, map__ss_id__mesh_face_id__faces_info] : map__profile_face_index__ss_id__mesh_face_id__faces_info) {
                                for (auto& [ss_id, map__mesh_face_id__faces_info] : map__ss_id__mesh_face_id__faces_info) {
                                  for (auto& [mesh_face_id, vect__faces_info] : map__mesh_face_id__faces_info) {
                                    for (auto& face_info : vect__faces_info) {
                                      // Разбить все точки face на те, которые рассчитаны или заданы и те, которые не заданы
                                      for (auto& point_index : face_info.face_verts_indexes) {
                                        if (point_index < 0) {
                                          // Если индекс отрицательный, то имеем дело с проектной точкой. Нужно преобразовать её в точку mesh.
                                          auto& project_point_info = vect__stack__project_points[-(point_index + 1)];
                                          point_index = calc_points.get_index_or_append__project_point__vertex__application_counter(
                                            project_point_info.project_point,
                                            project_point_info.profile_face_index,
                                            project_point_info.application_counter,
                                            face_info.oioa0_altitude,
                                            face_info.oioa1_altitude,
                                            face_info.oioa0_offset,
                                            face_info.oioa1_offset
                                          );
                                        }
                                        // TODO: Вроде как тут можно отказаться от этого условия и рассчёта, если передавать параметры для
                                        // рассчёта точки сразу в предыдущее условие, чтобы там и считать, если такая точка ещё
                                        // не была создана.
                                        // Update: нет, нельзя. <image url="..\code_images\file_0093.png" scale=".4"/>. TODO: требуется объяснение.
//                                        auto& shape_point = calc_points.vect__point_index__mesh_point[point_index];
//                                        if (shape_point.is_altitude_calc == false) {
//                                          // Формула рассчёта высоты проектных точек, по данным из волнового фронта (shape_point.event_time)
//                                          FT z_time = (face_info.oioa1_altitude - face_info.oioa0_altitude) * (shape_point.event_time - face_info.oioa0_offset) / (face_info.oioa1_offset - face_info.oioa0_offset) + face_info.oioa0_altitude;
//#ifdef _DEBUG
//                                          double oioa0_altitude = CGAL::to_double(face_info.oioa0_altitude);
//                                          double oioa1_altitude = CGAL::to_double(face_info.oioa1_altitude);
//                                          double z_time_d = CGAL::to_double(z_time);
//#endif
//                                          //calc_points.vect__point_index__mesh_point[shape_point.index].is_altitude_calc = true;
//                                          //calc_points.vect__point_index__mesh_point[shape_point.index].point = Point_3(shape_point.point.x(), shape_point.point.y(), z_time);
//                                          //calc_points.vect__point_index__mesh_point[shape_point.index].oioa_altitude = z_time;
//                                          shape_point.is_altitude_calc = true;
//                                          shape_point.point = Point_3(shape_point.point.x(), shape_point.point.y(), z_time);
//                                          shape_point.oioa_altitude = z_time;
//                                        }
//#ifdef _DEBUG
//                                        double oioa0_altitude = CGAL::to_double(face_info.oioa0_altitude);
//                                        double oioa1_altitude = CGAL::to_double(face_info.oioa1_altitude);
//                                        double oioa_altitude = CGAL::to_double(calc_points.vect__point_index__mesh_point[shape_point.index].oioa_altitude);
//                                        int x = 0;
//#endif
                                      }
                                    }
                                  }
                                }
                              }
                              timer1.stop();
                              if (verbose == true) {
                                printf("\n " VAL2STR(Msg::_0051) ". SS 2D Offset. Calculate altitudes of SS points of object_index: % 2d, build time: % 10.5f", object_index, timer1.time());
                              }
                            }

                            {
                              CGAL::Real_timer timer1; // Общее время рассчёта всех SS (вместе с загрузкой данных)
                              timer1.start();

                              std::vector<SM_Face_index> vect_faces_indexes;
                              map__object_index__mesh_ss_merged[object_index] = std::vector<Mesh>();

                              //if (results_join_mode != 0) {
                              //  // Не применять этот метод при других join_mode.
                              //  bevel_more_split = false;
                              //}

                              std::vector<Mesh> vect_mesh;
                              if (results_join_mode != 0) {
                                // Для results_join_mode!='SPLIT' всё собирать только в один Mesh
                                vect_mesh.push_back(Mesh());
                              }
                              std::map<std::tuple<int, int, int> /*point_index, ss_id, _application_index*/, Mesh::vertex_index> map__beveled_points_vertex_indexes;
                              for (auto& [profile_face_index, map__ss_id__mesh_face_id__faces_info] : map__profile_face_index__ss_id__mesh_face_id__faces_info) {
                                if(results_join_mode==0 && bevel_more_split == false /*Применять bevel_more_split только при results_join_mode==0*/) {
                                  // Для results_join_mode=='SPLIT', но дробление только по profile_face_index, то на каждый profile_face_index создавать по Mesh.
                                  vect_mesh.push_back(Mesh());
                                  map__beveled_points_vertex_indexes.clear();
                                }
                                for (auto& [ss_id, map__mesh_face_id__faces_info] : map__ss_id__mesh_face_id__faces_info) {
                                  if(results_join_mode==0 && bevel_more_split == true /*Применять bevel_more_split только при results_join_mode==0*/) {
                                    // Для results_join_mode=='SPLIT', но дробление кроме как по profile_face_index, ещё и по ss_id то на каждый profile_face_index/ss_id создавать по Mesh.
                                    vect_mesh.push_back(Mesh());
                                    map__beveled_points_vertex_indexes.clear();
                                  }

                                  // Собрать mesh для Bevel SS:
                                  Mesh& beveled_mesh = vect_mesh.back();
                                  for (auto& [mesh_face_id, vect__faces_info] : map__mesh_face_id__faces_info) {
                                    // Работа с данными одной из face в профиле

                                    for (auto& face_info : vect__faces_info) {
                                      std::vector<Mesh::vertex_index> vect_face_indexes;
                                      for (auto& point_index : face_info.face_verts_indexes) {
                                        SM_Vertex_index point_id;
                                        int _point_index = point_index;
                                        int _application_index = 0;
                                        int _profile_face_index = profile_face_index;
                                        auto& shape_point = calc_points.vect__point_index__mesh_point[point_index];
                                        if (shape_point.point_type == SHAPE_POINT::POINT_TYPE::APPLICATED_PROJECT_POINT) {
                                          // Такая точка может попасться несколько раз в зависимости от _application_index
                                          if (results_join_mode==1 || results_join_mode==2) {
                                            _point_index = shape_point.project_point_index;
                                            _application_index = shape_point.application_index;
                                          } else {
                                            // _application_index остаётся как был, т.к. не меняется point_index, который в пределах своего ss_id однозначно указывает на любую точку (не только в точке UP/DOWN).
                                            _point_index = point_index;
                                            _application_index = shape_point.application_index;
                                          }
                                        } else {
                                          // Для остальных типов точек (UP/DOWN), которые являются точками пересечения he и SS _application_index всегда=1, т.к. такая точка никогда не граничит с другим ss_id,
                                          // т.к. "двигается" вдоль контура своего ss_id и поэтому встречается только один раз внутри ss_id и не на границе ss_id.
                                          _point_index = point_index;
                                          _application_index = 1;
                                        }
                                        const auto& beveled_point_index = std::tuple(_point_index, _profile_face_index, _application_index);
                                        const auto& it_beveled_points = map__beveled_points_vertex_indexes.find(beveled_point_index);
                                        if (it_beveled_points == map__beveled_points_vertex_indexes.end()) {
#ifdef _DEBUG
                                          Point_3 p_test = Point_3(shape_point.point);
                                          double x = CGAL::to_double(p_test.x());
                                          double y = CGAL::to_double(p_test.y());
                                          double z = CGAL::to_double(p_test.z());
#endif
                                          point_id = beveled_mesh.add_vertex(shape_point.point);
                                          map__beveled_points_vertex_indexes[beveled_point_index] = point_id;
#ifdef _DEBUG
                                          auto point_id_idx = point_id.idx();
                                          shape_point.beveled_indexes.push_back(point_id_idx);
#endif
                                        } else {
                                          auto& [_point_index, _point_id] = *it_beveled_points;
                                          point_id = _point_id;
                                        }
                                        vect_face_indexes.push_back(point_id);
                                      }
                                      const auto& face_id = beveled_mesh.add_face(vect_face_indexes);

                                      if (face_id == Mesh::null_face()) {
                                        if (verbose == true) {
                                          printf("\nERROR " VAL2STR(Msg::_0041) " Failed to add face with verts (ss_id: %3u, mesh_face_id: %3u): [", ss_id, mesh_face_id);
                                          for (auto& vert_index : vect_face_indexes) {
                                            printf(" %3u", vert_index.idx());
                                          }
                                          printf("]");
                                        }
                                      } else {
                                        vect_faces_indexes.push_back(face_id);
                                      }
                                    }
                                  }
                                }
                              }
                              for (auto& mesh : vect_mesh) {
                                map__object_index__mesh_ss_merged[object_index].push_back(mesh);
                              }
                              timer1.stop();
                              if (verbose == true) {
                                printf("\n " VAL2STR(Msg::_0050) ". SS 2D Offset. Build result mesh of object_index: % 2d, build time: % 10.5f", object_index, timer1.time());
                              }
                            }
                          }
                          if (verbose) {
                            printf("\n " VAL2STR(Msg::_0033) ". SS 2D Offset. object_index %u, ss timer_bevel %.5g, %.5g", object_index, timer_bevel.time(), t1);
                          }
                        }
                      }
                    } // for (auto& [object_index, mesh_ss_joined] : map__object_index__mesh_ss_joined)
                  }

                  timer1.stop();
                  if (verbose) {
                    printf("\n " VAL2STR(Msg::_0049) ". SS 2D Offset. ss faces ajacent timer %.5g", timer1.time());
                  }

                  int split_object_index = 0; // для join_mode==SPLIT индекс нового объекта назначается группе контуров с одним index_offset (а не каждому контуру, чтобы количество объектов в режиме Edges и Faces совпадало)
                  int  keep_object_index = 0; // Для join_mode==KEEP  индексы результирующих объектов не меняются
                  int merge_object_index = 0; // Для join_mode==MERGE индекс результирующего объекта один - 0, т.к. объект будет единственным.
                  Mesh mesh_join_all_objects; // При join_mode==MERGE выполнить join всех mesh

                  // Сгруппировать или разгруппировать объекты по параметру results_join_mode
                  for (auto& [object_index, vect__meshes] : map__object_index__mesh_ss_merged) {
                    if (results_join_mode == 0) {
                      for (auto& mesh : vect__meshes) {
                        if (map_join_mesh.find(split_object_index) == map_join_mesh.end()) {
                          map_join_mesh[split_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                        }
                        OIOA_OFFSET_SS_PARAMS offset_info(OIOA(split_object_index, map__object_index__poms[object_index][0].object_original_index, 0, 0.0, 0.0), std::vector<Polygon_with_holes_2>(), std::vector<SS__HalfEdge__Point_2>());
                        offset_info.mesh_ss_02_merged = mesh;
                        map_join_mesh[split_object_index].push_back(offset_info);
                        // Запомнить привязку индекса нового объекта к индексу оригинального объекта
                        map__new_object_index__original_object_index[split_object_index] = map__object_index__poms[object_index][0].object_original_index;
                        split_object_index++;
                      }
                    } else if (results_join_mode == 1) {
                      // results_join_mode уже учтён, поэтому в векторе должен находится один объект, поэтому merge делать не нужно
                      keep_object_index = object_index;
                      map__new_object_index__original_object_index[object_index] = object_index;
                      Mesh mesh_joined;
                      for (auto& mesh : vect__meshes) {
                        mesh_joined.join(mesh);
                      }
                      map_join_mesh[keep_object_index] = std::vector<OIOA_OFFSET_SS_PARAMS>();
                      OIOA_OFFSET_SS_PARAMS offset_info(OIOA(keep_object_index, map__object_index__poms[keep_object_index][0].object_original_index, 0, 0.0, 0.0), std::vector<Polygon_with_holes_2>(), std::vector<SS__HalfEdge__Point_2>());
                      // Запомнить привязку индекса нового объекта к индексу оригинального объекта
                      map__new_object_index__original_object_index[keep_object_index] = map__object_index__poms[keep_object_index][0].object_original_index;
                      offset_info.mesh_ss_02_merged = mesh_joined;
                      map_join_mesh[keep_object_index].push_back(offset_info);
                    } else if (results_join_mode == 2) {
                      // results_join_mode уже учтён, поэтому в векторе должен находится один объект, поэтому merge делать не нужно
                      for (auto& mesh : vect__meshes) {
                        mesh_join_all_objects.join(mesh);
                      }
                    }
                  }

                  if (results_join_mode == 2) {

                    Mesh mesh_merged_all_objects; // При join_mode==MERGE выполнить join всех mesh
                    OIOA oioa1(merge_object_index, map__object_index__poms[merge_object_index][0].object_original_index, 0, 0.0, 0.0);
                    OIOA_OFFSET_SS_PARAMS offset_info(oioa1, std::vector<Polygon_with_holes_2>(), std::vector<SS__HalfEdge__Point_2>());
                    // Параметры для получения результатов merge:
                    std::vector<Point_3>                  result_ss_mesh_merged_points;
                    std::vector<std::vector<std::size_t>> result_ss_mesh_merged_polygons;
                    {
                      CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(mesh_join_all_objects, result_ss_mesh_merged_points, result_ss_mesh_merged_polygons);
                      CGAL::Polygon_mesh_processing::merge_duplicate_points_in_polygon_soup(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons);
                      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(result_ss_mesh_merged_points, result_ss_mesh_merged_polygons, mesh_merged_all_objects);
                    }
                    offset_info.mesh_ss_02_merged = mesh_merged_all_objects;
                    map_join_mesh[merge_object_index].push_back(offset_info);
                    // Запомнить привязку индекса нового объекта к индексу оригинального объекта
                    map__new_object_index__original_object_index[merge_object_index] = map__object_index__poms[merge_object_index][0].object_original_index;
                  }
                } // else if (result_type == ResType::STRAIGHT_SKELETON)
              } // Конец рассчёта map_join_mesh. В настоящий момент известен количественный состав результата.

              mesh_data->nn_objects = map_join_mesh.size();
              if (map_join_mesh.size() == 0) {
              }else{
                CGAL::Real_timer timer1;
                timer1.start();
                if (verbose) {
                  printf("\n " VAL2STR(Msg::_0044) ". SS 2D Offset. Building result data. Count of objects: % 3zd", map_join_mesh.size());
                }
                mesh_data->nn_objects_indexes           = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
                mesh_data->nn_objects_original_indexes  = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
                mesh_data->nn_offsets_counts            = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
                mesh_data->nn_verts                     = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
                mesh_data->nn_edges                     = (int*)malloc(sizeof(int) * mesh_data->nn_objects);
                mesh_data->nn_faces                     = (int*)malloc(sizeof(int) * mesh_data->nn_objects);

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

                for (auto& [object_index, vect_pwhs_with_offset_index] : map_join_mesh) {
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
                  //} else if (result_type == ResType::BEVEL) {
                  } else if (result_type == ResType::BEVEL || result_type == ResType::STRAIGHT_SKELETON) {
                    // TODO: На этом уровне данные от Extrude должны быть подготовлены. Тут их нужно получить, только "извлечь" для отправки.
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
                      for (const auto& ei : pwh_altitude.mesh_ss_02_merged.edges()) {
                        const auto& vi0 = pwh_altitude.mesh_ss_02_merged.vertex(ei, 0);
                        const auto& vi1 = pwh_altitude.mesh_ss_02_merged.vertex(ei, 1);
                        vect_res_object1_edges.push_back({ vi0, vi1 });
                      }

                      for (const auto& fi : pwh_altitude.mesh_ss_02_merged.faces() ) {
                        std::vector<unsigned int> vect_face_vertex_id;
                        const auto& hf = pwh_altitude.mesh_ss_02_merged.halfedge(fi);
                        for (auto& hi : halfedges_around_face(hf, pwh_altitude.mesh_ss_02_merged)) {
                          const auto& vi = target(hi, pwh_altitude.mesh_ss_02_merged);
                          vect_face_vertex_id.push_back( vi.idx() );
                        }
                        vect_res_object1_faces.push_back(vect_face_vertex_id);
                      }
                    }
                  }

                  if (verbose) {
                    printf("\n " VAL2STR(Msg::_0064) ". SS 2D Offset. object_index: % 3d, vertices: % 9zd, edges: % 9zd, faces: % 9zd", object_index, vect_res_object1_verts.size(), vect_res_object1_edges.size(), vect_res_object1_faces.size());
                  }
                  mesh_data->nn_objects_indexes         [object_index] = object_index; // Индекс объекта можно взять из первого результата
                  mesh_data->nn_objects_original_indexes[object_index] = map__new_object_index__original_object_index[object_index];
                  mesh_data->nn_offsets_counts          [object_index] = set_object1_offsets_indexes.size();
                  mesh_data->nn_verts                   [object_index] = vect_res_object1_verts.size();
                  mesh_data->nn_edges                   [object_index] = vect_res_object1_edges.size();
                  mesh_data->nn_faces                   [object_index] = vect_res_object1_faces.size();

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

                if (verbose) {
                  printf("\n " VAL2STR(Msg::_0065) ". SS 2D Offset. --------------------------------------------------------------------------");
                  printf("\n " VAL2STR(Msg::_0065) ". SS 2D Offset. Finally mesh:      vertices: % 9zd, edges: % 9zd, faces: % 9zd", vect_objects_verts.size(), vect_objects_edges.size(), vect_objects_faces.size());
                }

                // записать результаты vertices, edges и faces
                // Временно отключу, т.к. вернуть надо столько объектов, сколько и поступило.
                //if (vect_objects_verts.size() == 0) {
                //  mesh_data->nn_offsets_indexes = NULL;
                //  mesh_data->vertices = NULL;
                //  mesh_data->edges = NULL;
                //  mesh_data->faces = NULL;
                //  mesh_data->nn_faces_indexes_counters = NULL;
                //} else 
                {
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
                    }
                  }
                }

                timer1.stop();
                if (verbose == true) {
                  printf("\n " VAL2STR(Msg::_0066) ". SS 2D Offset. Output data time: %.5g", timer1.time());
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
                      res_errors.push_back(ObjectError(I, contour1_points, VAL2STR(Msg::_0009)". Polygon contains less 3 points."));
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
                      res_errors.push_back(ObjectError(I, plane1[0].c1, VAL2STR(Msg::_0010)". Boundary polygon is not simple."));
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
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Msg::_0010)". Hole is not simple."));
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
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Msg::_0023)". Hole is not inside boundary."));
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
                          res_errors.push_back(ObjectError(I, c1, VAL2STR(Msg::_0024)". Hole intersect with hole."));
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
                  //  res_errors.push_back(ObjectError(object1.object_index, *ppwh.get(), VAL2STR(Msg::_0022) ". Holes of the Polygon intersect amongst themselves or with outer boundary ()"));
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
                        printf("\n" VAL2STR(Msg::_0015) ". Polygons (%u, %u) INTERSECTS. Both marked as error", I, IJ);
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
                  res_errors.push_back(ObjectError(object1.object_index, *object1.vect_polygons_with_holes[I].get(), VAL2STR(Msg::_0011)". Polygons intersects"));
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
                          printf("\n" VAL2STR(Msg::_0027) ". SS Extrude: Thread: %i, Object: %i, ss_valid - %d, contours vertices: %u, extrude time - %.5g sec.", threads_counts, object1.object_index, object1.ss_extrude_is_valid, vertices_count, timer1.time());
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
          if (md->nn_objects_original_indexes != NULL) {
            free(md->nn_objects_original_indexes);
            md->nn_objects_original_indexes = NULL;
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