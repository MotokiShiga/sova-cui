#ifdef _WIN32

#define HAVE_BOOLEAN

#include <windows.h> /* required for all Windows applications */
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct subgrid_cell {
    int num_atoms;
    int *atom_positions;
    int num_domains;
    int *domain_points;
    int *domain_indices;
} subgrid_cell_t;

typedef struct subgrid {
    subgrid_cell_t *a;
    int cubesize;
    int ncells;
    int dimensions[3];
    int strides[3];
} subgrid_t;

#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
#define INDEXDISCGRID(i,j,k) ((int64_t)(i)*discgrid_strides[0]+(j)*discgrid_strides[1]+(k)*discgrid_strides[2])
DLLEXPORT void atomstogrid(
        int64_t *grid, int dimensions[3], int strides[3],
        int natoms, int *atom_positions, int *radii_indices,
        int nradii, int *radii,
        int ntranslations, int *translations,
        char *discretization_grid, int discgrid_strides[3]);
#undef INDEXGRID
#undef INDEXDISCGRID

DLLEXPORT subgrid_t *subgrid_create(int cubesize, int grid_dimensions[3]);
DLLEXPORT void subgrid_destroy(subgrid_t *sg);
DLLEXPORT int subgrid_index(subgrid_t *sg, int *pos);
DLLEXPORT void subgrid_add_atoms(subgrid_t *sg,
        int natoms, int *atom_positions,
        int ntranslations, int *translations);
DLLEXPORT void subgrid_add_domains(subgrid_t *sg,
        int npoints, int *domain_indices, int *domain_points,
        int ntranslations, int *translations);

#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
#define INDEXDISCGRID(i,j,k) ((int64_t)(i)*discgrid_strides[0]+(j)*discgrid_strides[1]+(k)*discgrid_strides[2])
DLLEXPORT void mark_cavities(int64_t *grid, int64_t *domain_grid, int dimensions[3], int strides[3],
        char *discretization_grid, int discgrid_strides[3],
        subgrid_t *sg, int use_surface_points);
#undef INDEXGRID
#undef INDEXDISCGRID

#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
#define INDEXDISCGRID(i,j,k) ((int64_t)(i)*discgrid_strides[0]+(j)*discgrid_strides[1]+(k)*discgrid_strides[2])
DLLEXPORT int cavity_triangles(
        int64_t *cavity_grid,
        int dimensions[3],
        int strides[3],
        int ncavity_indices,
        int *cavity_indices,
        int isolevel,
        float step[3],
        float offset[3],
        int8_t *discretization_grid,
        int discgrid_strides[3],
        float **vertices,
        float **normals,
        float *surface_area);
#undef INDEXGRID
#undef INDEXDISCGRID

DLLEXPORT void free_float_p(float *p);

#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
DLLEXPORT void cavity_intersections(
        int64_t *grid,
        int dimensions[3],
        int strides[3],
        int num_domains,
        int8_t *intersection_table);
#undef INDEXGRID

#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
DLLEXPORT void mark_translation_vectors(int8_t *grid,
                              int dimensions[3],
                              int strides[3],
                              int ntranslations,
                              int *translations);
#undef INDEXGRID

#ifdef __cplusplus
}
#endif