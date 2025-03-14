import ctypes
import numpy as np
import sys
import traceback
import itertools
from time import time

#print("== pySVCGAL.py start =================================================")

from pySVCGAL.clib import load_library
SVCGAL_clib = load_library.load_library()

# Результаты рассчёта в виде Mesh
class MESH_DATA2(ctypes.Structure):
    _fields_ = [
        ("has_error", ctypes.c_bool, ),
        ("str_error", ctypes.c_char_p, ),

        ("polygon_id", ctypes.c_int),
        
        ("nn_objects", ctypes.c_int),
        ("nn_objects_indexes", ctypes.POINTER( ctypes.c_int), ),
        ("nn_objects_original_indexes", ctypes.POINTER( ctypes.c_int), ),
        ("nn_offsets_counts" , ctypes.POINTER( ctypes.c_int), ),
        ("nn_offsets_indexes", ctypes.POINTER( ctypes.c_int), ),

        ("nn_verts", ctypes.POINTER( (ctypes.c_int)), ),
        ("nn_edges", ctypes.POINTER( (ctypes.c_int)), ),
        ("nn_faces", ctypes.POINTER( (ctypes.c_int)), ),
        ("nn_faces_indexes_counters", ctypes.POINTER( (ctypes.c_int)), ), # Счётчик индексов вершин для faces (количество счётчиков равно количеству nn_faces). Этот параметр появился, потому что количество вершин в faces стало переменным.
        ("vertices", ctypes.POINTER( (ctypes.c_float)*3), ), # Координаты вершин - тройной float
        ("edges"   , ctypes.POINTER( (ctypes.c_int)*2), ), # Индексы рёбер (всегда по 2)
        ("faces"   , ctypes.POINTER( (ctypes.c_int)  ), ), # Индексы вершин в гранях
        #("faces"   , ctypes.POINTER( (ctypes.c_int)*3), ),

        ("nn_source_objects_count"   , ctypes.c_int                          ), # Количество исходных объектов для чтения и выдачи информации по ошибкам по исходных объектам. Количество результирующих mesh может отличаться от исходных объектов по параметру results_join_mode
        ("nn_source_objects_indexes" , ctypes.POINTER( (ctypes.c_int   )  ), ), # индексы исходных объектов
        ("nn_errors_count"           , ctypes.POINTER( (ctypes.c_int   )  ), ),
        ("nn_description1_per_error1", ctypes.POINTER( (ctypes.c_char_p)  ), ),
        ("nn_contours_per_error1"    , ctypes.POINTER( (ctypes.c_int   )  ), ),
        ("nn_vertices_per_contour1"  , ctypes.POINTER( (ctypes.c_int   )  ), ),
        ("vertices_of_errors"        , ctypes.POINTER( (ctypes.c_float )*3), ),

    ]

def pySVCGAL_straight_skeleton_2d_offset(data):
    """
    Documentation
    """   
    time01 = time()

    # Оставить на будущее для дальнейших работ по увеличению производительности:
    # fill_vertices = SVCGAL_clib.fill_vertices
    # fill_vertices.argtypes = [
    #     ctypes.POINTER((ctypes.c_float)*3 ),  # src
    #     np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),  # target
    #     ctypes.c_int,                         # count
    # ]
    # fill_vertices.restype=None

    # https://docs.python.org/3/library/ctypes.html
    straight_skeleton_2d_offset = SVCGAL_clib.straight_skeleton_2d_offset
    straight_skeleton_2d_offset.argtypes = [
        ctypes.c_int,                       # in_count_of_objects
        ctypes.POINTER((ctypes.c_int  ) ),  # in_shapes_mode
        ctypes.POINTER((ctypes.c_int  ) ),  # in_count_of_offsets (chained)
        ctypes.POINTER((ctypes.c_float) ),  # in_offsets 
        ctypes.POINTER((ctypes.c_int  ) ),  # in_count_of_altitudes (chained)
        ctypes.POINTER((ctypes.c_float) ),  # in_altitudes (chained)
        ctypes.POINTER((ctypes.c_int  ) ),  # in__profile_faces__counters
        ctypes.POINTER((ctypes.c_int  ) ),  # in__profile_faces__indexes_counters
        ctypes.POINTER((ctypes.c_int  ) ),  # in__profile_faces__close_mode
        ctypes.POINTER((ctypes.c_int  ) ),  # in__profile_faces__indexes_array - (array)
        ctypes.POINTER((ctypes.c_int  ) ),  # in_count_of_planes (chained)
        ctypes.POINTER((ctypes.c_int  ) ),  # in_contours_in_planes (chained)

        ctypes.POINTER((ctypes.c_int  ) ),  # in_vertices_in_contours - count of vertices in contours (chained)
        ctypes.POINTER((ctypes.c_float)*3), # in_vertices - (array)
        ctypes.c_bool,                      # only_tests_for_valid
        ctypes.c_int,                       # result type: 0 - contours, 1 - Mesh
        ctypes.c_bool,                      # force_z_zero
        ctypes.c_int,                       # source_objects_join_mode 0 - split, 1 - keep, 2 - merge all meshes
        ctypes.c_int,                       # results_join_mode 0 - split, 1 - keep, 2 - merge all meshes
        ctypes.c_bool,                      # verbose (additional info output to a console)
        ctypes.c_bool,                      # use_cache_of_straight_skeleton
        ctypes.c_bool,                      # bevel_more_split - if negative offset or several contours in source shape
    ]
    straight_skeleton_2d_offset.restype = ctypes.POINTER(MESH_DATA2)

    free_mem2 = SVCGAL_clib.free_mem2
    free_mem2.argtypes = [
        ctypes.POINTER(MESH_DATA2)
    ]

    force_z_zero                   = data['force_z_zero']
    res_type                       = data['res_type']
    only_tests_for_valid           = data['only_tests_for_valid']
    source_objects_join_mode       = data['source_objects_join_mode']
    results_join_mode              = data['results_join_mode']
    verbose                        = data['verbose']
    use_cache_of_straight_skeleton = data['use_cache_of_straight_skeleton']
    bevel_more_split               = data['bevel_more_split']

    offsets_counters                = []
    offsets_array                   = []
    altitudes_array                 = []
    matrix_array                    = []
    profile_faces__counters         = [] # Счётчик количества faces на один объект (количество этого параметра равно количеству объектов)
    profile_faces__indexes_counters = [] # Счётчик количества индексов у каждого face
    profile_faces__indexes_array    = [] # Индексы каждого face
    profile_faces__close_modes_array = [] # режим открытия/закрытия каждого профильного face (0-закрыто, 1-открыто)
    planes_in_object                = [] # Сколько в Объекте planes [1-элемент на объект]
    contours_in_planes_counter      = [] # Сколько в каждом плане контуров
    vertices_in_contour_counters    = []
    vertices_list                   = []
    shapes_mode                     = [] # Количество соответствует количеству объектов

    for I, object1 in enumerate(data['objects']):
        offsets_counters.append(len(object1['offsets']))
        offsets_array   .extend(object1['offsets'])
        altitudes_array .extend(object1['altitudes'])
        matrix_array    .append(object1['matrix'])

        # Собрать данные по faces, из которых состоит профиль для SS
        profile_faces__counters.append(len(object1['profile_faces__indexes']))
        for face1 in object1['profile_faces__indexes']:
            profile_faces__indexes_counters.append(len(face1))
            profile_faces__indexes_array.extend(face1)
            pass
        profile_faces__close_modes = object1['profile_faces__close_modes']
        profile_faces__close_modes_array.extend(profile_faces__close_modes)
        planes_in_object.append(len(object1['planes'])) # количество планов текущего объекта
        for plane1 in object1['planes']:
            contours_in_planes_counter.append(len(plane1)) # Количество контуров по каждому плану
            for contour1 in plane1:
                vertices_in_contour_counters.append(len(contour1))
                vertices_list.extend(contour1)
            pass
        shapes_mode.append(object1['shape_mode'])
        pass
    pass

    params = [
        len(data['objects']),            #in_count_of_objects
        shapes_mode,                     # 0-ORIGINAL_BOUNDARIES, 1-EXCLUDE_HOLES, 2-INVERT_HOLES
        offsets_counters,                # in_count_of_offsets - array (len is count of objects)
        offsets_array,
        offsets_counters,                # in_count_of_altitudes - array (len is count of objects)
        altitudes_array,                 # in_altitudes
        profile_faces__counters,         # 
        profile_faces__indexes_counters, # in_offset_edges
        profile_faces__indexes_array,    #
        profile_faces__close_modes_array, #
        planes_in_object,                # in_count_of_planes
        contours_in_planes_counter,      # in_contours_in_planes
        vertices_in_contour_counters,    # in_vertices_in_contours
        vertices_list,                   # in_vertices
        only_tests_for_valid,
        res_type,
        force_z_zero,
        source_objects_join_mode,
        results_join_mode,
        verbose
    ]

    # https://docs.python.org/3/library/ctypes.html
    ctypes_in_count_of_objects      = len(data['objects'])
    ctypes_in_shapes_mode           = ctypes.ARRAY( len(shapes_mode)                 ,  ctypes.c_int     )(*[(ctypes.c_int    )( s1) for s1 in shapes_mode ])
    ctypes_in_count_of_offsets      = ctypes.ARRAY( len(offsets_counters)            ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in offsets_counters ])
    ctypes_in_offsets               = ctypes.ARRAY( len(offsets_array)               ,  ctypes.c_float   )(*[(ctypes.c_float  )( o1) for o1 in offsets_array    ])
    ctypes_in_count_of_altitudes    = ctypes.ARRAY( len(offsets_counters)            ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in offsets_counters ])
    ctypes_in_altitudes             = ctypes.ARRAY( len(altitudes_array)             ,  ctypes.c_float   )(*[(ctypes.c_float  )( a1) for a1 in altitudes_array  ])
    
    ctypes_in__profile_faces__counters              = ctypes.ARRAY( len(profile_faces__counters)         ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in profile_faces__counters ])
    ctypes_in__profile_faces__indexes_counters      = ctypes.ARRAY( len(profile_faces__indexes_counters) ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in profile_faces__indexes_counters ])
    ctypes_in__profile_faces__indexes_array         = ctypes.ARRAY( len(profile_faces__indexes_array)    ,  ctypes.c_int     )(*[(ctypes.c_int    )( f1) for f1 in profile_faces__indexes_array ])
    ctypes_in__profile_faces__close_modes_array     = ctypes.ARRAY( len(profile_faces__close_modes_array),  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in profile_faces__close_modes_array ])
    
    #ctypes_in_count_of_offset_edges = ctypes.ARRAY( len(offset_edges_counters)       ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in offset_edges_counters ])
    #ctypes_in_offset_edges          = ctypes.ARRAY( len(offset_edges_array)          ,  ctypes.c_int*2   )(*[(ctypes.c_int*2  )(*e1) for e1 in offset_edges_array  ])
    
    ctypes_in_count_of_planes       = ctypes.ARRAY( len(planes_in_object)            ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in planes_in_object ])
    ctypes_in_contours_in_planes    = ctypes.ARRAY( len(contours_in_planes_counter)  ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in contours_in_planes_counter ])
    ctypes_in_vertices_in_contours  = ctypes.ARRAY( len(vertices_in_contour_counters),  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in vertices_in_contour_counters ])
    ctypes_in_vertices              = ctypes.ARRAY( len(vertices_list)               ,  ctypes.c_float*3 )(*[(ctypes.c_float*3)(*v1) for v1 in vertices_list ])
    
    md = None
    new_mesh = {
        'objects': [], # Массив объектов
        'source_objects_errors':[], # массив ошибок по исходным объектам
        'vertices'  : [],
        'edges'     : [],
        'faces'     : [],
    }

    time01 = time()-time01
    if verbose==True:
        print(f'\nPrepare data to ctypes before call library: {time01} ms')

    try:
        timer_process_external_call = time()
        md = straight_skeleton_2d_offset(
            ctypes_in_count_of_objects,
            ctypes_in_shapes_mode,
            ctypes_in_count_of_offsets,
            ctypes_in_offsets,
            ctypes_in_count_of_altitudes,
            ctypes_in_altitudes,
            ctypes_in__profile_faces__counters,
            ctypes_in__profile_faces__indexes_counters,
            ctypes_in__profile_faces__indexes_array,
            ctypes_in__profile_faces__close_modes_array,
            #ctypes_in_count_of_offset_edges,
            #ctypes_in_offset_edges,
            ctypes_in_count_of_planes,
            ctypes_in_contours_in_planes,
            ctypes_in_vertices_in_contours,
            ctypes_in_vertices,
            only_tests_for_valid,
            res_type,
            force_z_zero,
            source_objects_join_mode,
            results_join_mode,
            verbose,
            use_cache_of_straight_skeleton,
            bevel_more_split,
        )
        if verbose==True:
            timer_process_external_call = time()-timer_process_external_call
            print(f'\nexternal process finished in {timer_process_external_call} ms')
        ################ Get Results ################
        mdc = md.contents

        str_error = None
        if(mdc.str_error is not None):
            str_error = mdc.str_error.decode("ascii")

        # Собрать результаты, чтобы хоть что-то показать пользователю.
        new_vertices = []
        new_edges    = []
        new_faces    = []
        
        time01 = time()
        if mdc.nn_objects>0:
            verts_I0   = 0
            edges_I0   = 0
            faces_I0   = 0
            faces_indexes_I0 = 0
            offsets_I0 = 0

            error_I0                        = 0
            contours_per_error_cursor       = 0
            nn_vertices_per_contour1_cursor = 0
            vertices_of_errors_cursor       = 0
            offsets_indexes_cursor          = 0

            mdc_vertices = mdc.vertices
            mdc_edges = mdc.edges
            mdc_faces = mdc.faces
            mdc_nn_faces = mdc.nn_faces
            mdc_nn_faces_indexes_counters = mdc.nn_faces_indexes_counters

            #summ_vertices_count = 0 # Для тестирования загрузки данных в numpy array

            for I in range(mdc.nn_objects):
                object_index            = mdc.nn_objects_indexes[I]
                object_original_index   = mdc.nn_objects_original_indexes[I]
                if verbose==True:
                    ## for developer:
                    # print("")
                    # print(f'object index: {object_index}, ', end="")
                    pass

                #print(f"Loading data for object {I}")
                #### Extract New Vertices #### 
                vertices_count = mdc.nn_verts[I]
                if verbose==True:
                    ## for developer:
                    #print(f'verts: [pos: {verts_I0}, count: {vertices_count}], ', end="")
                    pass
                #summ_vertices_count += vertices_count # Для тестирования загрузки данных через numpy array
                new_vertices1 = [ tuple(mdc_vertices[verts_I0+i]) for i in range(vertices_count)]
                verts_I0 += vertices_count

                #### Extract New Edges ####
                edges_count = mdc.nn_edges[I]
                if verbose==True:
                    ## for developer:
                    #print(f'edges: [pos: {edges_I0}, count: {edges_count}], ', end="")
                    pass
                new_edges1 = [ tuple(mdc_edges[edges_I0+i]) for i in range(edges_count)]
                edges_I0 += edges_count

                #### Extract New Faces #### 
                # extract face indexes counters:
                faces_count = mdc_nn_faces[I]
                faces_indexes_counters = [ mdc_nn_faces_indexes_counters[faces_I0+i] for i in range(faces_count)]
                if verbose==True:
                    ## for developer:
                    #print(f'faces: [pos: {faces_indexes_I0}, count: {faces_indexes_counters}], ', end="")
                    pass
                faces_indexes_pos = 0
                new_faces1 = []
                new_faces1_append = new_faces1.append
                for counter in faces_indexes_counters:
                    face1 = []
                    face1_append = face1.append
                    for i in range(counter):
                        face1_append(mdc_faces[faces_indexes_I0+faces_indexes_pos])
                        faces_indexes_pos+=1
                        pass
                    new_faces1_append(tuple(face1))
                faces_I0 += faces_count
                faces_indexes_I0 += faces_indexes_pos

                #### Extract indexes of object
                offsets_counts = mdc.nn_offsets_counts[I]
                offsets_indexes = [ mdc.nn_offsets_indexes[offsets_I0+i] for i in range(offsets_counts)]
                offsets_I0 += offsets_counts
                
                new_mesh['objects'].append({
                    'object_index'              : object_index,
                    'object_original_index'     : object_original_index,
                    'offsets_indexes'           : offsets_indexes,
                    'vertices'                  : new_vertices1,
                    'edges'                     : new_edges1,
                    'faces'                     : new_faces1,
                })

                pass
            pass

        
            # if summ_vertices_count>0:
            #     # Тестирование загрузки vertices в numpy массив. Оставить на будущее для дальнейших работ по увеличению производительности:
            #     np_vertices = np.zeros((summ_vertices_count, 3), dtype=np.float32)
            #     time_10 = time()
            #     fill_vertices(mdc.vertices, np_vertices, summ_vertices_count)
            #     time_10 = time()-time_10
            # pass

        time01 = time()-time01
        if verbose==True:
            print(f'\nOffset Straight Skeleton. Load data from external call: {time01} ms')

        # Извлечь ошибки по исходным объектам
        if mdc.nn_source_objects_count>0:
            error_I0 = 0
            contours_per_error_cursor = 0
            nn_vertices_per_contour1_cursor = 0
            vertices_of_errors_cursor = 0
            for I in range(mdc.nn_source_objects_count):
                object_index = mdc.nn_source_objects_indexes[I]
                ### Извлечь количество ошибок в объекте
                nn_errors_count = mdc.nn_errors_count[I]
                description1_per_error1 = []
                nn_contours_per_error1 = []
                nn_vertices_per_contour1 = []
                vertices_of_errors = []
                for I in range(nn_errors_count):
                    error_description = mdc.nn_description1_per_error1[error_I0+I]
                    nn_contours_per_error1_I = mdc.nn_contours_per_error1[error_I0+I]

                    description1_per_error1.append(error_description)
                    nn_contours_per_error1.append(nn_contours_per_error1_I)
                    vertices_of_errors1 = []
                    for IJ in range(nn_contours_per_error1_I):
                        nn_vertices_per_contour1_IJ =  mdc.nn_vertices_per_contour1[nn_vertices_per_contour1_cursor+IJ]
                        nn_vertices_per_contour1.append(nn_vertices_per_contour1_IJ)
                        for IJK in range(nn_vertices_per_contour1_IJ):
                            vertices_of_errors1.append(mdc.vertices_of_errors[vertices_of_errors_cursor+IJK][:])
                            pass
                        vertices_of_errors_cursor+=nn_vertices_per_contour1_IJ
                        pass
                    vertices_of_errors.extend(vertices_of_errors1)
                    nn_vertices_per_contour1_cursor+=nn_contours_per_error1_I
                        
                    pass
                error_I0+=nn_errors_count
                
                new_mesh['source_objects_errors'].append({
                    'object_index'              : object_index,
                    'vertices_of_errors'        : vertices_of_errors,
                    'descriptions_per_errors'   : description1_per_error1,
                })

                #print(f"Finish loading data for object {I}")

        #print("free_mem2 start")
        free_mem2(mdc)
        #print("free_mem2 finish")

        return new_mesh

    except Exception as _ex:
        str_error = ""
        if md is not None:
            mdc = md.contents
            mdc_str_error = ""
            if mdc.str_error is not None:
                mdc_str_error = mdc.str_error.decode("ascii")
                str_error = f"Unexpected exception while calculating data: {mdc_str_error}."
            else:
                str_error = f"Unexpected exception while calculating data. No internal addition info."
            free_mem2(mdc)
        else:
            str_error = "General unexpected exception."
        str_error_summ = f"{str_error} {_ex}."
        print(str_error_summ)
        raise _ex
        #return new_mesh


def pySVCGAL_straight_skeleton_2d_extrude(data):
    """
    Documentation
    """   

    # https://docs.python.org/3/library/ctypes.html
    straight_skeleton_2d_extrude = SVCGAL_clib.straight_skeleton_2d_extrude
    straight_skeleton_2d_extrude.argtypes = [
        ctypes.c_int,                       # in_count_of_objects
        ctypes.POINTER((ctypes.c_int  ) ),  # in_shapes_mode
        ctypes.POINTER((ctypes.c_bool ) ),  # in_restrict_height
        ctypes.POINTER((ctypes.c_double) ), # in_height
        ctypes.POINTER((ctypes.c_double) ), # in_angles - Должны точно соответствовать параметрам контуров in_count_of_planes и in_contours_in_planes
        ctypes.POINTER((ctypes.c_int  ) ),  # in_count_of_planes (chained)
        ctypes.POINTER((ctypes.c_int  ) ),  # in_contours_in_planes (chained)

        ctypes.POINTER((ctypes.c_int  ) ),  # in_vertices_in_contours - count of vertices in contours (chained)
        ctypes.POINTER((ctypes.c_float)*3), # in_vertices - (array)
        ctypes.c_bool,                      # only_tests_for_valid
        ctypes.c_bool,                      # force_z_zero
        ctypes.c_int,                       # join_mode 0 - split, 1 - keep, 2 - merge all meshes
        ctypes.c_bool,                      # verbose (additional info output to a console)
    ]
    straight_skeleton_2d_extrude.restype = ctypes.POINTER(MESH_DATA2)

    free_mem2 = SVCGAL_clib.free_mem2
    free_mem2.argtypes = [
        ctypes.POINTER(MESH_DATA2)
    ]

    force_z_zero         = data['force_z_zero']
    only_tests_for_valid = data['only_tests_for_valid']
    join_mode            = data['join_mode']
    verbose              = data['verbose']

    planes_in_object = []   # Сколько в Объекте planes [1-элемент на объект]
    angles_in_object = []   # Сколько в Объекте planes [1-элемент на объект]
    contours_in_planes_counter = [] # Сколько в каждом плане контуров
    vertices_in_contour_counters = []
    vertices_list = []
    angles_list = []
    shapes_mode_list = [] # Количество соответствует количеству объектов
    restrict_height_list = []
    height_list = []

    for I, object1 in enumerate(data['objects']):
        shapes_mode_list    .append(object1['shape_mode'])
        restrict_height_list.append(object1['restrict_height'])
        height_list         .append(object1['height'])

        planes_in_object.append(len(object1['planes'])) # количество планов текущего объекта
        for plane1 in object1['planes']:
            contours_in_planes_counter.append(len(plane1)) # Количество контуров по каждому плану
            for contour1 in plane1:
                vertices_in_contour_counters.append(len(contour1))
                vertices_list.extend(contour1)
            pass
        pass
        for angles1 in object1['angles']:
            for angles_in_contour1 in angles1:
                angles_list.extend(angles_in_contour1)
            pass
        pass
    pass

    params = [
        len(data['objects']),           #in_count_of_objects
        shapes_mode_list,               # 0-ORIGINAL_BOUNDARIES, 1-EXCLUDE_HOLES, 2-INVERT_HOLES
        restrict_height_list,
        height_list,
        angles_list,
        planes_in_object,               # in_count_of_planes
        contours_in_planes_counter,     # in_contours_in_planes
        vertices_in_contour_counters,   # in_vertices_in_contours
        vertices_list,                  # in_vertices
        only_tests_for_valid,
        force_z_zero,
        join_mode,
        verbose
    ]

    # https://docs.python.org/3/library/ctypes.html
    ctypes_in_count_of_objects     = len(data['objects'])
    ctypes_in_shapes_mode_list     = ctypes.ARRAY( len(shapes_mode_list)            ,  ctypes.c_int     )(*[(ctypes.c_int    )( s1) for s1 in shapes_mode_list ])
    ctypes_in_restrict_height_list = ctypes.ARRAY( len(restrict_height_list)        ,  ctypes.c_bool    )(*[(ctypes.c_bool   )( s1) for s1 in restrict_height_list ])
    ctypes_in_height_list          = ctypes.ARRAY( len(height_list)                 ,  ctypes.c_double    )(*[(ctypes.c_double )( s1) for s1 in height_list ])
    ctypes_in_angles_list          = ctypes.ARRAY( len(angles_list)                 ,  ctypes.c_double    )(*[(ctypes.c_double )( s1) for s1 in angles_list ])
    ctypes_in_count_of_planes      = ctypes.ARRAY( len(planes_in_object)            ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in planes_in_object ])
    ctypes_in_contours_in_planes   = ctypes.ARRAY( len(contours_in_planes_counter)  ,  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in contours_in_planes_counter ])
    ctypes_in_vertices_in_contours = ctypes.ARRAY( len(vertices_in_contour_counters),  ctypes.c_int     )(*[(ctypes.c_int    )( c1) for c1 in vertices_in_contour_counters ])
    ctypes_in_vertices             = ctypes.ARRAY( len(vertices_list)               ,  ctypes.c_float*3 )(*[(ctypes.c_float*3)(*v1) for v1 in vertices_list ])

    md = None
    new_mesh = {
        'objects': [], # Массив объектов
        'source_objects_errors':[], # массив ошибок по исходным объектам
        'vertices'  : [],
        'edges'     : [],
        'faces'     : [],
    }

    try:
        md = straight_skeleton_2d_extrude(
            ctypes_in_count_of_objects,
            ctypes_in_shapes_mode_list,
            ctypes_in_restrict_height_list,
            ctypes_in_height_list,
            ctypes_in_angles_list,
            ctypes_in_count_of_planes,
            ctypes_in_contours_in_planes,
            ctypes_in_vertices_in_contours,
            ctypes_in_vertices,
            only_tests_for_valid,
            force_z_zero,
            join_mode,
            verbose
        )
        ################ Get Results ################
        mdc = md.contents

        str_error = None
        if(mdc.str_error is not None):
            str_error = mdc.str_error.decode("ascii")

        # Собрать результаты, чтобы хоть что-то показать пользователю.
        new_vertices = []
        new_edges    = []
        new_faces    = []
        
        if mdc.nn_objects>0:
            verts_I0   = 0
            edges_I0   = 0
            faces_I0    = 0
            faces_indexes_I0 = 0
            offsets_I0 = 0

            for I in range(mdc.nn_objects):
                object_index = mdc.nn_objects_indexes[I]

                #### Extract New Vertices #### 
                vertices_count = mdc.nn_verts[I]
                new_vertices1 = [ tuple(mdc.vertices[verts_I0+i]) for i in range(vertices_count)]
                verts_I0 += vertices_count

                #### Extract New Edges ####
                edges_count = mdc.nn_edges[I]
                new_edges1 = [ tuple(mdc.edges[edges_I0+i]) for i in range(edges_count)]
                edges_I0 += edges_count

                #### Extract New Faces #### 
                # faces_count = mdc.nn_faces[I]
                # new_faces1 = [ tuple(mdc.faces[faces_I0+i]) for i in range(faces_count)]
                # faces_I0 += faces_count

                #### Extract New Faces #### 
                # extract face indexes counters:
                faces_count = mdc.nn_faces[I]
                faces_indexes_counters = [ mdc.nn_faces_indexes_counters[faces_I0+i] for i in range(faces_count)]
                faces_indexes_pos = 0
                new_faces1 = []
                for counter in faces_indexes_counters:
                    face1 = []
                    face1_append = face1.append
                    for i in range(counter):
                        face1_append(mdc.faces[faces_indexes_I0+faces_indexes_pos])
                        faces_indexes_pos+=1
                        pass
                    new_faces1.append(tuple(face1))
                faces_I0 += faces_count
                faces_indexes_I0 += faces_indexes_pos


                new_mesh['objects'].append({
                    'object_index'              : object_index,
                    'vertices'                  : new_vertices1,
                    'edges'                     : new_edges1,
                    'faces'                     : new_faces1,
                })

        # Извлечь ошибки по исходным объектам
        if mdc.nn_source_objects_count>0:
            error_I0 = 0
            contours_per_error_cursor = 0
            nn_vertices_per_contour1_cursor = 0
            vertices_of_errors_cursor = 0
            for I in range(mdc.nn_source_objects_count):
                object_index = mdc.nn_source_objects_indexes[I]
                ### Извлечь количество ошибок в объекте
                nn_errors_count = mdc.nn_errors_count[I]
                description1_per_error1 = []
                nn_contours_per_error1 = []
                nn_vertices_per_contour1 = []
                vertices_of_errors = []
                for I in range(nn_errors_count):
                    error_description = mdc.nn_description1_per_error1[error_I0+I]
                    nn_contours_per_error1_I = mdc.nn_contours_per_error1[error_I0+I]

                    description1_per_error1.append(error_description)
                    nn_contours_per_error1.append(nn_contours_per_error1_I)
                    vertices_of_errors1 = []
                    for IJ in range(nn_contours_per_error1_I):
                        nn_vertices_per_contour1_IJ =  mdc.nn_vertices_per_contour1[nn_vertices_per_contour1_cursor+IJ]
                        nn_vertices_per_contour1.append(nn_vertices_per_contour1_IJ)
                        for IJK in range(nn_vertices_per_contour1_IJ):
                            vertices_of_errors1.append(mdc.vertices_of_errors[vertices_of_errors_cursor+IJK][:])
                            pass
                        vertices_of_errors_cursor+=nn_vertices_per_contour1_IJ
                        pass
                    vertices_of_errors.extend(vertices_of_errors1)
                    nn_vertices_per_contour1_cursor+=nn_contours_per_error1_I
                        
                    pass
                error_I0+=nn_errors_count
                
                new_mesh['source_objects_errors'].append({
                    'object_index'              : object_index,
                    'vertices_of_errors'        : vertices_of_errors,
                    'descriptions_per_errors'   : description1_per_error1,
                })

                #print(f"Finish loading data for object {I}")

        #print("free_mem2 start")
        free_mem2(mdc)
        #print("free_mem2 finish")

        return new_mesh

    except Exception as _ex:
        #print(_ex)
        str_error = ""
        if md is not None:
            mdc = md.contents
            mdc_str_error = ""
            if mdc.str_error is not None:
                mdc_str_error = mdc.str_error.decode("ascii")
                str_error = f"Unexpected exception while calculating data: {mdc_str_error}."
            else:
                str_error = f"Unexpected exception while calculating data. No internal addition info."
            free_mem2(mdc)
        else:
            str_error = "General unexpected exception."
        str_error_summ = f"{str_error} {_ex}."
        print(str_error_summ)
        raise _ex
        #return new_mesh

if __name__ == "__main__":
    pass