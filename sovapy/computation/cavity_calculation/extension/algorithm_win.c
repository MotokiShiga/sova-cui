#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "algorithm_win.h"

#define ABS(x) ((x) < 0 ? -(x) : (x))
#define INDEX(x, y, z) ((x)*mcdata.stride[0] + (y)*mcdata.stride[1] + (z)*mcdata.stride[2])
#define IDX2D(y, z) ((y)*mcdata.dim[2] + (z))
#define SQUARE(x) ((x)*(x))
#define CLIP(x,a,b) ((x)<(a)?(a):((x)>(b)?(b):(x)))
#define GR3_MC_DTYPE unsigned short

/* speedup does not grow much with a high number of threads */
#define THREADLIMIT 16


#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
#define INDEXDISCGRID(i,j,k) ((int64_t)(i)*discgrid_strides[0]+(j)*discgrid_strides[1]+(k)*discgrid_strides[2])

/**
 * Mark spheres around atoms on the grid.
 * For each discretized atom and its equivalents in adjacent cells:
 * for each cell of the discretized sphere around them:
 * find the grid cells which are inside the cutoff radius.
 * For each of this grid cells:
 * check if the cell is inside the volume,
 * check if this atom is the closest to this cell
 * and write the atom index (+1) into it.
 */
void atomstogrid(
        int64_t *grid, int dimensions[3], int strides[3],
        int natoms, int *atom_positions, int *radii_indices,
        int nradii, int *radii,
        int ntranslations, int *translations,
        char *discretization_grid, int discgrid_strides[3])
{
    int i, j, k;
    int radius;
    int cubesize;
    int atompos[3];
    int transpos[3];
    int sphereindex[3];
    int gridpos[3];
    int grid_index;
    int grid_value;
    int this_squared_distance;
    int other_squared_distance;
    int other_atompos[3];
    int other_transpos[3];

    (void) nradii;

    for (i = 0; i < natoms; i++) {
        radius = radii[radii_indices[i]];
        cubesize = 2 * radius + 1;
        atompos[0] = atom_positions[i * 3 + 0];
        atompos[1] = atom_positions[i * 3 + 1];
        atompos[2] = atom_positions[i * 3 + 2];
        for (j = 0; j < ntranslations; j++) {
            transpos[0] = atompos[0] + translations[j * 3 + 0];
            transpos[1] = atompos[1] + translations[j * 3 + 1];
            transpos[2] = atompos[2] + translations[j * 3 + 2];
            if (transpos[0] + radius < 0 || transpos[0] - radius >= dimensions[0]
                    || transpos[1] + radius < 0 || transpos[1] - radius >= dimensions[1]
                    || transpos[2] + radius < 0 || transpos[2] - radius >= dimensions[2]) {
                /* entire cube is outside */
                continue;
            }
            for (sphereindex[0] = 0; sphereindex[0] < cubesize; sphereindex[0]++) {
                gridpos[0] = transpos[0] + sphereindex[0] - radius;
                if (gridpos[0] < 0 || gridpos[0] >= dimensions[0]) {
                    continue;
                }
                for (sphereindex[1] = 0; sphereindex[1] < cubesize; sphereindex[1]++) {
                    gridpos[1] = transpos[1] + sphereindex[1] - radius;
                    if (gridpos[1] < 0 || gridpos[1] >= dimensions[1]) {
                        continue;
                    }
                    for (sphereindex[2] = 0; sphereindex[2] < cubesize; sphereindex[2]++) {
                        gridpos[2] = transpos[2] + sphereindex[2] - radius;
                        if (gridpos[2] < 0 || gridpos[2] >= dimensions[2]) {
                            continue;
                        }
                        if (SQUARE(sphereindex[0] - radius)
                                + SQUARE(sphereindex[1] - radius)
                                + SQUARE(sphereindex[2] - radius)
                                <= SQUARE(radius) 
                                && discretization_grid[INDEXDISCGRID(gridpos[0], gridpos[1], gridpos[2])] == 0) {
                            grid_index = INDEXGRID(gridpos[0], gridpos[1], gridpos[2]);
                            grid_value = grid[grid_index];
                            /* check if it is the closest atom */
                            if (grid_value == 0) {
                                grid[grid_index] = i + 1;
                            } else {
                                this_squared_distance = SQUARE(transpos[0] - gridpos[0])
                                        + SQUARE(transpos[1] - gridpos[1])
                                        + SQUARE(transpos[2] - gridpos[2]);
                                other_atompos[0] = atom_positions[3 * (grid_value - 1) + 0];
                                other_atompos[1] = atom_positions[3 * (grid_value - 1) + 1];
                                other_atompos[2] = atom_positions[3 * (grid_value - 1) + 2];
                                for (k = 0; k < ntranslations; k++) {
                                    other_transpos[0] = other_atompos[0] + translations[k * 3 + 0];
                                    other_transpos[1] = other_atompos[1] + translations[k * 3 + 1];
                                    other_transpos[2] = other_atompos[2] + translations[k * 3 + 2];
                                    other_squared_distance = SQUARE(other_transpos[0] - gridpos[0])
                                            + SQUARE(other_transpos[1] - gridpos[1])
                                            + SQUARE(other_transpos[2] - gridpos[2]);
                                    if (other_squared_distance <= this_squared_distance) {
                                        break;
                                    }
                                }
                                if (this_squared_distance < other_squared_distance) {
                                    grid[grid_index] = i + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#undef INDEXGRID
#undef INDEXDISCGRID


/**
 * Routines to work with subgrids
 */

subgrid_t *subgrid_create(int cubesize, int grid_dimensions[3])
{
    subgrid_t *sg;
    int k;

    sg = malloc(sizeof(subgrid_t));
    sg->cubesize = cubesize;
    for (k = 0; k < 3; k++) {
        sg->dimensions[k] = (int) ceil((double) grid_dimensions[k] / cubesize) + 4;
    }
    sg->ncells = sg->dimensions[0] * sg->dimensions[1] * sg->dimensions[2];
    sg->strides[0] = sg->dimensions[1] * sg->dimensions[2];
    sg->strides[1] = sg->dimensions[2];
    sg->strides[2] = 1;

    sg->a = calloc(sg->ncells, sizeof(subgrid_cell_t));

    return sg;
}

void subgrid_destroy(subgrid_t *sg)
{
    int i;

    for (i = 0; i < sg->ncells; i++) {
        free(sg->a[i].atom_positions);
        free(sg->a[i].domain_points);
        free(sg->a[i].domain_indices);
    }
    free(sg->a);
    free(sg);
}

static
int subgrid_index(subgrid_t *sg, int *pos)
{
    int k;
    int z;
    int index;

    index = 0;
    for (k = 0; k < 3; k++) {
        /* python's floor division */
        z = (int) (floor((double) pos[k] / sg->cubesize)) + 2;
        index += CLIP(z, 0, sg->dimensions[k] - 1) * sg->strides[k];
    }
    return index;
}

void subgrid_add_atoms(subgrid_t *sg,
        int natoms, int *atom_positions,
        int ntranslations, int *translations)
{
    int i, j;
    int real_pos[3];
    subgrid_cell_t *cell;

    for (i = 0; i < natoms; i++) {
        for (j = 0; j < ntranslations; j++) {
            real_pos[0] = atom_positions[i * 3 + 0] + translations[j * 3 + 0];
            real_pos[1] = atom_positions[i * 3 + 1] + translations[j * 3 + 1];
            real_pos[2] = atom_positions[i * 3 + 2] + translations[j * 3 + 2];
            cell = sg->a + subgrid_index(sg, real_pos);
            cell->atom_positions = realloc(cell->atom_positions,
                    (cell->num_atoms + 1) * 3 * sizeof(int));
            cell->atom_positions[3 * cell->num_atoms + 0] = real_pos[0];
            cell->atom_positions[3 * cell->num_atoms + 1] = real_pos[1];
            cell->atom_positions[3 * cell->num_atoms + 2] = real_pos[2];
            cell->num_atoms++;
        }
    }
}

void subgrid_add_domains(subgrid_t *sg,
        int npoints, int *domain_indices, int *domain_points,
        int ntranslations, int *translations)
{
    int i, j;
    int real_pos[3];
    subgrid_cell_t *cell;

    for (i = 0; i < npoints; i++) {
        for (j = 0; j < ntranslations; j++) {
            real_pos[0] = domain_points[i * 3 + 0] + translations[j * 3 + 0];
            real_pos[1] = domain_points[i * 3 + 1] + translations[j * 3 + 1];
            real_pos[2] = domain_points[i * 3 + 2] + translations[j * 3 + 2];
            cell = sg->a + subgrid_index(sg, real_pos);
            cell->domain_indices = realloc(cell->domain_indices,
                    (cell->num_domains + 1) * sizeof(int));
            cell->domain_indices[cell->num_domains] = domain_indices[i];
            cell->domain_points = realloc(cell->domain_points,
                    (cell->num_domains + 1) * 3 * sizeof(int));
            cell->domain_points[3 * cell->num_domains + 0] = real_pos[0];
            cell->domain_points[3 * cell->num_domains + 1] = real_pos[1];
            cell->domain_points[3 * cell->num_domains + 2] = real_pos[2];
            cell->num_domains++;
        }
    }
}


#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
#define INDEXDISCGRID(i,j,k) ((int64_t)(i)*discgrid_strides[0]+(j)*discgrid_strides[1]+(k)*discgrid_strides[2])

/**
 * For each cell, determine if it is closer to a cavity domain than
 * to an atom center. If so, mark the cell in the grid.
 */
void mark_cavities(int64_t *grid, int64_t *domain_grid, int dimensions[3], int strides[3],
        char *discretization_grid, int discgrid_strides[3],
        subgrid_t *sg, int use_surface_points)
{
    int pos[3];
    int grid_index;
    int grid_value;
    int sg_index;
    int min_squared_atom_distance;
    int squared_atom_distance;
    int neigh[3];
    int neigh_index;
    subgrid_cell_t *cell;
    int i;
    int breaknext;
    int squared_domain_distance;

    for (pos[0] = 0; pos[0] < dimensions[0]; pos[0]++) {
        for (pos[1] = 0; pos[1] < dimensions[1]; pos[1]++) {
            for (pos[2] = 0; pos[2] < dimensions[2]; pos[2]++) {
                grid_index = INDEXGRID(pos[0], pos[1], pos[2]);
                if (use_surface_points) {
                    grid_value = domain_grid[grid_index];
                    if (grid_value == 0) {
                        /* outside the volume */
                        grid[grid_index] = 0;
                        continue;
                    } else if (grid_value < 0) {
                        /* cavity domain (stored as: -index-1), therefore guaranteed to be in a cavity */
                        grid[grid_index] = grid_value;
                        continue;
                    } else {
                        grid[grid_index] = 0;
                    }
                } else {
                    if (discretization_grid[INDEXDISCGRID(pos[0], pos[1], pos[2])] != 0) {
                        continue;
                    }
                }
                /* step 5 */
                min_squared_atom_distance = INT_MAX;
                sg_index = subgrid_index(sg, pos);
                for (neigh[0] = -1; neigh[0] <= 1; neigh[0]++) {
                    for (neigh[1] = -1; neigh[1] <= 1; neigh[1]++) {
                        for (neigh[2] = -1; neigh[2] <= 1; neigh[2]++) {
                            neigh_index = sg_index + neigh[0] * sg->strides[0]
                                    + neigh[1] * sg->strides[1]
                                    + neigh[2] * sg->strides[2];
                            cell = sg->a + neigh_index;
                            for (i = 0; i < cell->num_atoms; i++) {
                                squared_atom_distance = 
                                        SQUARE(cell->atom_positions[i * 3 + 0] - pos[0])
                                        + SQUARE(cell->atom_positions[i * 3 + 1] - pos[1])
                                        + SQUARE(cell->atom_positions[i * 3 + 2] - pos[2]);
                                if (squared_atom_distance < min_squared_atom_distance) {
                                    min_squared_atom_distance = squared_atom_distance;
                                }
                            }
                        }
                    }
                }
                breaknext = 0;
                for (neigh[0] = -1; neigh[0] <= 1; neigh[0]++) {
                    for (neigh[1] = -1; neigh[1] <= 1; neigh[1]++) {
                        for (neigh[2] = -1; neigh[2] <= 1; neigh[2]++) {
                            neigh_index = sg_index + neigh[0] * sg->strides[0]
                                    + neigh[1] * sg->strides[1]
                                    + neigh[2] * sg->strides[2];
                            cell = sg->a + neigh_index;
                            for (i = 0; i < cell->num_domains; i++) {
                                squared_domain_distance = 
                                        SQUARE(cell->domain_points[i * 3 + 0] - pos[0])
                                        + SQUARE(cell->domain_points[i * 3 + 1] - pos[1])
                                        + SQUARE(cell->domain_points[i * 3 + 2] - pos[2]);
                                if (squared_domain_distance < min_squared_atom_distance) {
                                    grid[grid_index] = -cell->domain_indices[i] - 1;
                                    breaknext = 1;
                                    break; /* i */
                                }
                            }
                            if (breaknext) {
                                break; /* neigh[2] */
                            }
                        }
                        if (breaknext) {
                            break; /* neigh[1] */
                        }
                    }
                    if (breaknext) {
                        break; /* neigh[0] */
                    }
                }
            }
        }
    }
}
#undef INDEXGRID
#undef INDEXDISCGRID


/* calculate the gradient via difference qoutient */
static gr3_coord_t getgrad(mcdata_t mcdata, int x, int y, int z)
{
  int v[3];
  int neigh[3][2];
  int i;
  gr3_coord_t n;

  v[0] = x;
  v[1] = y;
  v[2] = z;

  for (i = 0; i < 3; i++)
    {
      if (v[i] > 0)
        neigh[i][0] = v[i] - 1;
      else
        neigh[i][0] = v[i];
      if (v[i] < mcdata.dim[i] - 1)
        neigh[i][1] = v[i] + 1;
      else
        neigh[i][1] = v[i];
    }
  n.x = (float)(mcdata.data[INDEX(neigh[0][1], y, z)] - mcdata.data[INDEX(neigh[0][0], y, z)]) /
        (neigh[0][1] - neigh[0][0]) / mcdata.step[0];
  n.y = (float)(mcdata.data[INDEX(x, neigh[1][1], z)] - mcdata.data[INDEX(x, neigh[1][0], z)]) /
        (neigh[1][1] - neigh[1][0]) / mcdata.step[1];
  n.z = (float)(mcdata.data[INDEX(x, y, neigh[2][1])] - mcdata.data[INDEX(x, y, neigh[2][0])]) /
        (neigh[2][1] - neigh[2][0]) / mcdata.step[2];

  return n;
}

/* interpolate points and calulate normals */
static void interpolate(mcdata_t mcdata, int px, int py, int pz, GR3_MC_DTYPE v1, int qx, int qy, int qz,
                        GR3_MC_DTYPE v2, gr3_coord_t *p, gr3_coord_t *n)
{
  double mu;
  gr3_coord_t n1, n2;
  double norm;

  if (ABS(mcdata.isolevel - v1) < 0.00001)
    mu = 0.0;
  else if (ABS(mcdata.isolevel - v2) < 0.00001)
    mu = 1.0;
  else if (ABS(v1 - v2) < 0.00001)
    mu = 0.5;
  else
    mu = 1.0 * (mcdata.isolevel - v1) / (v2 - v1);

  p->x = (px + mu * (qx - px)) * mcdata.step[0] + mcdata.offset[0];
  p->y = (py + mu * (qy - py)) * mcdata.step[1] + mcdata.offset[1];
  p->z = (pz + mu * (qz - pz)) * mcdata.step[2] + mcdata.offset[2];

  n1 = getgrad(mcdata, px, py, pz);
  n2 = getgrad(mcdata, qx, qy, qz);
  n->x = -(n1.x + mu * (n2.x - n1.x));
  n->y = -(n1.y + mu * (n2.y - n1.y));
  n->z = -(n1.z + mu * (n2.z - n1.z));

  norm = sqrt(n->x * n->x + n->y * n->y + n->z * n->z);
  if (norm > 0.0)
    {
      n->x /= norm;
      n->y /= norm;
      n->z /= norm;
    }
}

/*!
 * marching cubes algorithm for one x-layer.
 * created vertices are cached between calls using vindex.
 * vindex associates the intersected edge with the vertex index.
 * the edge is identified by its location (low, high), direction (x, y, z)
 * and coordinates (py, pz) of its starting point.
 * direction and location are the first index:
 * (x, y_low, z_low, y_high, z_high) (see mc_edgeprop)
 * second index is py * mcdata.dim[1] + pz.
 * py and pz are the coordinates of the lower one of both edge vertices
 */
static void layer(mcdata_t mcdata, int x, int **vindex, unsigned int *num_vertices, gr3_coord_t **vertices,
                  gr3_coord_t **normals, unsigned int *vertcapacity, unsigned int *num_faces, unsigned int **indices,
                  unsigned int *facecapacity)
{
  int i, j;
  int y, z;
  int cubeindex;
  GR3_MC_DTYPE cubeval[8]; /* also cache between adjacent cubes */

  for (y = 0; y < mcdata.dim[1] - 1; y++)
    {
      /* init z-cache */
      for (i = 0; i < 4; i++)
        {
          int zi = mc_zvertices[0][i];

          cubeval[mc_zvertices[1][i]] =
              mcdata.data[INDEX(x + mc_cubeverts[zi][0], y + mc_cubeverts[zi][1], 0 + mc_cubeverts[zi][2])];
        }
      for (z = 0; z < mcdata.dim[2] - 1; z++)
        {
          cubeindex = 0;
          /* shift old values (z-cache) */
          for (i = 0; i < 4; i++)
            {
              int zi = mc_zvertices[0][i];

              cubeval[zi] = cubeval[mc_zvertices[1][i]];
              if (cubeval[zi] < mcdata.isolevel)
                {
                  cubeindex |= 1 << zi;
                }
            }
          /* read new cube values */
          for (i = 0; i < 4; i++)
            {
              int zi = mc_zvertices[1][i];

              cubeval[zi] =
                  mcdata.data[INDEX(x + mc_cubeverts[zi][0], y + mc_cubeverts[zi][1], z + mc_cubeverts[zi][2])];
              if (cubeval[zi] < mcdata.isolevel)
                {
                  cubeindex |= 1 << zi;
                }
            }
          if (cubeindex != 0 && cubeindex != 255)
            {
              /* create triangles */
              for (i = 0; i < mc_tricount[cubeindex]; i++)
                {
                  if (*facecapacity <= *num_faces)
                    {
                      (*facecapacity) = (unsigned int)(*num_faces * 1.5) + 50;
                      *indices = realloc(*indices, (*facecapacity) * 3 * sizeof(int));
                    }
                  /* create triangle vertices */
                  for (j = 0; j < 3; j++)
                    {
                      int trival = mc_tritable[cubeindex][i * 3 + j];
                      const int *edge = mc_cubeedges[trival];
                      int dir = mc_edgeprop[trival];
                      int px = x + mc_cubeverts[edge[0]][0];
                      int py = y + mc_cubeverts[edge[0]][1];
                      int pz = z + mc_cubeverts[edge[0]][2];
                      /* lookup if vertex already exists */
                      int node = vindex[dir][IDX2D(py, pz)];
                      if (node < 0)
                        {
                          /* it does not, create it */
                          GR3_MC_DTYPE v1 = cubeval[edge[0]];
                          GR3_MC_DTYPE v2 = cubeval[edge[1]];
                          if (*vertcapacity <= *num_vertices)
                            {
                              (*vertcapacity) = (unsigned int)(*num_vertices * 1.5) + 50;
                              *vertices = realloc(*vertices, (*vertcapacity) * sizeof(gr3_coord_t));
                              *normals = realloc(*normals, (*vertcapacity) * sizeof(gr3_coord_t));
                            }
                          node = *num_vertices;
                          interpolate(mcdata, px, py, pz, v1, x + mc_cubeverts[edge[1]][0],
                                      y + mc_cubeverts[edge[1]][1], z + mc_cubeverts[edge[1]][2], v2, *vertices + node,
                                      *normals + node);
                          vindex[dir][IDX2D(py, pz)] = node;
                          (*num_vertices)++;
                        }
                      /* add vertex index to the element array */
                      (*indices)[*num_faces * 3 + j] = node;
                    }
                  (*num_faces)++;
                }
            }
        }
    }
}

/*!
 * handle consecutive calls to layer
 */
static void layerblock(mcdata_t mcdata, int from, int to, unsigned int *num_vertices, gr3_coord_t **vertices,
                       gr3_coord_t **normals, unsigned int *num_faces, unsigned int **faces)
{
  int x;
  int y;
  int z;
  unsigned int vertcapacity;
  unsigned int facecapacity;
  /* cache for the vertex indices of the x-layer
   * [x, y_bot, z_bot, y_top, z_top] */
  int *vindex[5], *ntmp;

  *num_vertices = 0;
  vertcapacity = 0;
  *vertices = NULL;
  *normals = NULL;
  *num_faces = 0;
  facecapacity = 0;
  *faces = NULL;

  vindex[0] = malloc(5 * mcdata.dim[1] * mcdata.dim[2] * sizeof(int));
  vindex[1] = vindex[0] + mcdata.dim[1] * mcdata.dim[2];
  vindex[2] = vindex[0] + 2 * mcdata.dim[1] * mcdata.dim[2];
  vindex[3] = vindex[0] + 3 * mcdata.dim[1] * mcdata.dim[2];
  vindex[4] = vindex[0] + 4 * mcdata.dim[1] * mcdata.dim[2];

  for (y = 0; y < mcdata.dim[1]; y++)
    {
      for (z = 0; z < mcdata.dim[2]; z++)
        {
          vindex[0][IDX2D(y, z)] = -1;
          vindex[3][IDX2D(y, z)] = -1;
          vindex[4][IDX2D(y, z)] = -1;
        }
    }
  /*
   * iterate layer-by-layer through the data
   * create an indexed mesh
   * indices are cached in vindex[direction of the edge][y, z]
   * the top cache becomes the bottom in the next iterarion
   */
  for (x = from; x < to; x++)
    {
      ntmp = vindex[1];
      vindex[1] = vindex[3];
      vindex[3] = ntmp;
      ntmp = vindex[2];
      vindex[2] = vindex[4];
      vindex[4] = ntmp;
      for (y = 0; y < mcdata.dim[1]; y++)
        {
          for (z = 0; z < mcdata.dim[2]; z++)
            {
              vindex[0][IDX2D(y, z)] = -1;
              vindex[3][IDX2D(y, z)] = -1;
              vindex[4][IDX2D(y, z)] = -1;
            }
        }
      layer(mcdata, x, vindex, num_vertices, vertices, normals, &vertcapacity, num_faces, faces, &facecapacity);
    }
  free(vindex[0]);
}

/*!
 * Create an isosurface (as indexed mesh) from voxel data
 * with the marching cubes algorithm.
 * This function manages the parallelization:
 * Divide the data into blocks along the x-axis. Allocate memory,
 * call layerblock and merge the individual meshes into a single one.
 *
 * \param [in]  data          the volume (voxel) data
 * \param [in]  isolevel      value where the isosurface will be extracted
 * \param [in]  dim_x         number of elements in x-direction
 * \param [in]  dim_y         number of elements in y-direction
 * \param [in]  dim_z         number of elements in z-direction
 * \param [in]  stride_x      number of elements to step when traversing
 *                            the data in x-direction
 * \param [in]  stride_y      number of elements to step when traversing
 *                            the data in y-direction
 * \param [in]  stride_z      number of elements to step when traversing
 *                            the data in z-direction
 * \param [in]  step_x        distance between the voxels in x-direction
 * \param [in]  step_y        distance between the voxels in y-direction
 * \param [in]  step_z        distance between the voxels in z-direction
 * \param [in]  offset_x      coordinate origin
 * \param [in]  offset_y      coordinate origin
 * \param [in]  offset_z      coordinate origin
 * \param [out] num_vertices  number of vertices created
 * \param [out] vertices      array of vertex coordinates
 * \param [out] normals       array of vertex normal vectors
 * \param [out] num_indices   number of indices created
 *                            (3 times the number of triangles)
 * \param [out] indices       array of vertex indices that make the triangles
 */
void gr3_triangulateindexed(const GR3_MC_DTYPE *data, GR3_MC_DTYPE isolevel, unsigned int dim_x,
                                   unsigned int dim_y, unsigned int dim_z, unsigned int stride_x, unsigned int stride_y,
                                   unsigned int stride_z, double step_x, double step_y, double step_z, double offset_x,
                                   double offset_y, double offset_z, unsigned int *num_vertices, gr3_coord_t **vertices,
                                   gr3_coord_t **normals, unsigned int *num_indices, unsigned int **indices)
{
  int num_threads;
  unsigned int num_faces;
  unsigned int *num_t_vertices, *num_t_faces, **t_faces;
  gr3_coord_t **t_vertices, **t_normals;
  unsigned int *vertblock, *faceblock;
  mcdata_t mcdata;
#if defined(_OPENMP) && defined(THREADLIMIT)
  int max_threads;

  max_threads = omp_get_max_threads();
  if (max_threads > THREADLIMIT) omp_set_num_threads(THREADLIMIT);
#endif

  if (stride_x == 0) stride_x = dim_z * dim_y;
  if (stride_y == 0) stride_y = dim_z;
  if (stride_z == 0) stride_z = 1;

  mcdata.data = data;
  mcdata.isolevel = isolevel;
  mcdata.dim[0] = dim_x;
  mcdata.dim[1] = dim_y;
  mcdata.dim[2] = dim_z;
  mcdata.stride[0] = stride_x;
  mcdata.stride[1] = stride_y;
  mcdata.stride[2] = stride_z;
  mcdata.step[0] = step_x;
  mcdata.step[1] = step_y;
  mcdata.step[2] = step_z;
  mcdata.offset[0] = offset_x;
  mcdata.offset[1] = offset_y;
  mcdata.offset[2] = offset_z;

  *num_vertices = 0;
  *vertices = NULL;
  *normals = NULL;
  *num_indices = 0;
  *indices = NULL;

#ifdef _OPENMP
#pragma omp parallel default(none)                                                                                 \
    shared(num_threads, num_t_vertices, t_vertices, t_normals, num_t_faces, t_faces, mcdata, vertblock, faceblock, \
           num_vertices, num_faces, vertices, normals, indices)
#endif
  {
    int thread_id;
    unsigned int from, to;
    unsigned int i;
#ifdef _OPENMP
#pragma omp single
#endif
    {
      /* allocate temporary memory for each thread */
#ifdef _OPENMP
      num_threads = omp_get_num_threads();
#else
      num_threads = 1;
#endif
      num_t_vertices = malloc(num_threads * sizeof(unsigned int));
      t_vertices = malloc(num_threads * sizeof(gr3_coord_t *));
      t_normals = malloc(num_threads * sizeof(gr3_coord_t *));
      num_t_faces = malloc(num_threads * sizeof(unsigned int));
      t_faces = malloc(num_threads * sizeof(unsigned int *));
    }
    /* create a mesh per thread */
#ifdef _OPENMP
    thread_id = omp_get_thread_num();
#else
    thread_id = 0;
#endif
    from = thread_id * (mcdata.dim[0] - 1) / num_threads;
    to = (thread_id + 1) * (mcdata.dim[0] - 1) / num_threads;
    num_t_vertices[thread_id] = 0;
    t_vertices[thread_id] = NULL;
    t_normals[thread_id] = NULL;
    num_t_faces[thread_id] = 0;
    t_faces[thread_id] = NULL;
    layerblock(mcdata, from, to, num_t_vertices + thread_id, t_vertices + thread_id, t_normals + thread_id,
               num_t_faces + thread_id, t_faces + thread_id);
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      /* calculate beginning indices of thread blocks */
      vertblock = malloc((num_threads + 1) * sizeof(unsigned int));
      vertblock[0] = 0;
      faceblock = malloc((num_threads + 1) * sizeof(unsigned int));
      faceblock[0] = 0;
      for (i = 0; i < (unsigned int)num_threads; i++)
        {
          vertblock[i + 1] = vertblock[i] + num_t_vertices[i];
          faceblock[i + 1] = faceblock[i] + num_t_faces[i];
        }
      *num_vertices = vertblock[num_threads];
      num_faces = faceblock[num_threads];
      *vertices = realloc(*vertices, *num_vertices * sizeof(gr3_coord_t));
      *normals = realloc(*normals, *num_vertices * sizeof(gr3_coord_t));
      *indices = realloc(*indices, num_faces * 3 * sizeof(unsigned int));
    }
    /* copy thread meshes into the arrays */
    memmove(*vertices + vertblock[thread_id], t_vertices[thread_id], num_t_vertices[thread_id] * sizeof(gr3_coord_t));
    memmove(*normals + vertblock[thread_id], t_normals[thread_id], num_t_vertices[thread_id] * sizeof(gr3_coord_t));
    /* translate thread indices to global indices */
    for (i = 0; i < num_t_faces[thread_id]; i++)
      {
        (*indices)[(faceblock[thread_id] + i) * 3 + 0] = t_faces[thread_id][i * 3 + 0] + vertblock[thread_id];
        (*indices)[(faceblock[thread_id] + i) * 3 + 1] = t_faces[thread_id][i * 3 + 1] + vertblock[thread_id];
        (*indices)[(faceblock[thread_id] + i) * 3 + 2] = t_faces[thread_id][i * 3 + 2] + vertblock[thread_id];
      }
    free(t_vertices[thread_id]);
    free(t_normals[thread_id]);
    free(t_faces[thread_id]);
  }
  free(faceblock);
  free(vertblock);
  free(t_faces);
  free(num_t_faces);
  free(t_normals);
  free(t_vertices);
  free(num_t_vertices);
  *num_indices = num_faces * 3;
#if defined(_OPENMP) && defined(THREADLIMIT)
  omp_set_num_threads(max_threads);
#endif
}

/*!
 * Create an isosurface (as mesh) from voxel data
 * with the marching cubes algorithm.
 * This function calls gr3_triangulateindexed and copies the values.
 *
 * \param [in]  data          the volume (voxel) data
 * \param [in]  isolevel      value where the isosurface will be extracted
 * \param [in]  dim_x         number of elements in x-direction
 * \param [in]  dim_y         number of elements in y-direction
 * \param [in]  dim_z         number of elements in z-direction
 * \param [in]  stride_x      number of elements to step when traversing
 *                           the data in x-direction
 * \param [in]  stride_y      number of elements to step when traversing
 *                           the data in y-direction
 * \param [in]  stride_z      number of elements to step when traversing
 *                           the data in z-direction
 * \param [in]  step_x        distance between the voxels in x-direction
 * \param [in]  step_y        distance between the voxels in y-direction
 * \param [in]  step_z        distance between the voxels in z-direction
 * \param [in]  offset_x      coordinate origin
 * \param [in]  offset_y      coordinate origin
 * \param [in]  offset_z      coordinate origin
 * \param [out] triangles_p   array of triangle data
 *
 * \returns the number of triangles created
 */
unsigned int gr3_triangulate(const GR3_MC_DTYPE *data, GR3_MC_DTYPE isolevel, unsigned int dim_x,
                                    unsigned int dim_y, unsigned int dim_z, unsigned int stride_x,
                                    unsigned int stride_y, unsigned int stride_z, double step_x, double step_y,
                                    double step_z, double offset_x, double offset_y, double offset_z,
                                    gr3_triangle_t **triangles_p)
{
  unsigned int num_vertices;
  gr3_coord_t *vertices, *normals;
  unsigned int num_indices;
  unsigned int *indices;
  unsigned int i, j;
#if defined(_OPENMP) && defined(THREADLIMIT)
  int max_threads;

  max_threads = omp_get_max_threads();
  if (max_threads > THREADLIMIT) omp_set_num_threads(THREADLIMIT);
#endif

  gr3_triangulateindexed(data, isolevel, dim_x, dim_y, dim_z, stride_x, stride_y, stride_z, step_x, step_y, step_z,
                         offset_x, offset_y, offset_z, &num_vertices, &vertices, &normals, &num_indices, &indices);

  *triangles_p = malloc(num_indices / 3 * sizeof(gr3_triangle_t));
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j) shared(num_indices, triangles_p, indices, vertices, normals)
#endif
  for (i = 0; i < num_indices / 3; i++)
    {
      for (j = 0; j < 3; j++)
        {
          (*triangles_p)[i].vertex[j] = vertices[indices[i * 3 + j]];
          (*triangles_p)[i].normal[j] = normals[indices[i * 3 + j]];
        }
    }
  free(vertices);
  free(normals);
  free(indices);

#if defined(_OPENMP) && defined(THREADLIMIT)
  omp_set_num_threads(max_threads);
#endif
  return num_indices / 3;
}




#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
#define INDEXDISCGRID(i,j,k) ((int64_t)(i)*discgrid_strides[0]+(j)*discgrid_strides[1]+(k)*discgrid_strides[2])
int cavity_triangles(
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
        float *surface_area)
{
    uint16_t *counts;
    int pos[3];
    int gridindex;
    int gridval;
    int i, j, k;
    int is_cavity;
    int neigh[3];
    int neighindex;
    int bbox[2][3] = {{-1, -1, -1}, {-1, -1, -1}};
    int ntriangles;
    float *triangles_p;
    float *continuous_vertices;
    float *continuous_normals;
    double area;
    int any_outside;
    float *vertex_p;
    float *normal_p;
    int disc_pos[3];
    double a[3], b[3];
    double cross[3];

    counts = calloc(dimensions[0] * dimensions[1] * dimensions[2],
            sizeof(uint16_t));
    for (pos[0] = 1; pos[0] < dimensions[0] - 1; pos[0]++) {
        for (pos[1] = 1; pos[1] < dimensions[1] - 1; pos[1]++) {
            for (pos[2] = 1; pos[2] < dimensions[2] - 1; pos[2]++) {
                gridindex = INDEXGRID(pos[0], pos[1], pos[2]);
                counts[gridindex] += 100;
                gridval = cavity_grid[gridindex];
                is_cavity = 0;
                for (i = 0; i < ncavity_indices; i++) {
                    if (gridval == -cavity_indices[i] - 1) {
                        is_cavity = 1;
                        break;
                    }
                }
                if (!is_cavity) {
                    continue;
                }
                for (neigh[0] = -1; neigh[0] <= 1; neigh[0]++) {
                    for (neigh[1] = -1; neigh[1] <= 1; neigh[1]++) {
                        for (neigh[2] = -1; neigh[2] <= 1; neigh[2]++) {
                            neighindex = gridindex + INDEXGRID(
                                    neigh[0], neigh[1], neigh[2]);
                            counts[neighindex]++;
                        }
                    }
                }
                for (i = 0; i < 3; i++) {
                    if (bbox[0][i] == -1 ||
                            bbox[0][i] > pos[i] - 1) {
                        bbox[0][i] = pos[i] - 1;
                    }
                    if (bbox[1][i] == -1 ||
                            bbox[1][i] < pos[i] + 1) {
                        bbox[1][i] = pos[i] + 1;
                    }
                }
            }
        }
    }
    for (i = 0; i < 3; i++) {
        if (bbox[0][i] >= 1) {
            bbox[0][i]--;
        }
        if (bbox[1][i] < dimensions[i] - 1) {
            bbox[1][i]++;
        }
    }
    
    ntriangles = gr3_triangulate(
            counts + INDEXGRID(bbox[0][0], bbox[0][1], bbox[0][2]),
            100 + isolevel,
            bbox[1][0] - bbox[0][0] + 1,
            bbox[1][1] - bbox[0][1] + 1,
            bbox[1][2] - bbox[0][2] + 1,
            strides[0], strides[1], strides[2],
            1.0, 1.0, 1.0,
            bbox[0][0], bbox[0][1], bbox[0][2],
            (gr3_triangle_t **) &triangles_p);
    free(counts);

    continuous_vertices = malloc(ntriangles * 3 * 3 * sizeof(float));
    continuous_normals = malloc(ntriangles * 3 * 3 * sizeof(float));
    area = 0.0;
    for (i = 0; i < ntriangles; i++) {
        any_outside = 0;
        for (j = 0; j < 3; j++) {
            vertex_p = triangles_p + (i * 3 * 2 + j) * 3;
            normal_p = vertex_p + 3 * 3;
            for (k = 0; k < 3; k++) {
                disc_pos[k] = floor(vertex_p[k] + 0.5);
                continuous_vertices[(i * 3 + j) * 3 + k] = 
                        vertex_p[k] * step[k] + offset[k];
                continuous_normals[(i * 3 + j) * 3 + k] = 
                        normal_p[k] / step[k];
            }
            if (discretization_grid[INDEXDISCGRID(
                    disc_pos[0], disc_pos[1], disc_pos[2])] != 0) {
                any_outside = 1;
            }
        }
        if (!any_outside) {
            for (k = 0; k < 3; k++) {
                a[k] = continuous_vertices[(i * 3 + 1) * 3 + k]
                        - continuous_vertices[(i * 3 + 0) * 3 + k];
                b[k] = continuous_vertices[(i * 3 + 2) * 3 + k]
                        - continuous_vertices[(i * 3 + 0) * 3 + k];
            }
            cross[0] = a[1] * b[2] - a[2] * b[1];
            cross[1] = a[2] * b[0] - a[0] * b[2];
            cross[2] = a[0] * b[1] - a[1] * b[0];
            area += 0.5 * sqrt(cross[0] * cross[0]
                    + cross[1] * cross[1] + cross[2] * cross[2]);
        }
    }
    free(triangles_p);

    *vertices = continuous_vertices;
    *normals = continuous_normals;
    *surface_area = area;
    return ntriangles;
}
#undef INDEXGRID
#undef INDEXDISCGRID


void free_float_p(float *p)
{
    free(p);
}


#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
void cavity_intersections(
        int64_t *grid,
        int dimensions[3],
        int strides[3],
        int num_domains,
        int8_t *intersection_table)
{
    int pos[3];
    int i;
    int neigh[3];
    int gridindex, neighindex;
    int64_t domain1, domain2;
    int offsets[13][3] = {{-1, -1, -1}, {-1, -1, 0}, {-1, -1, 1}, {-1, 0, -1}, {-1, 0, 0}, {-1, 0, 1}, {-1, 1, -1}, {-1, 1, 0}, {-1, 1, 1}, {0, -1, -1}, {0, -1, 0}, {0, -1, 1}, {0, 0, -1}};

    for (pos[0] = 1; pos[0] < dimensions[0] - 1; pos[0]++) {
        for (pos[1] = 1; pos[1] < dimensions[1] - 1; pos[1]++) {
            for (pos[2] = 1; pos[2] < dimensions[2] - 1; pos[2]++) {
                gridindex = INDEXGRID(pos[0], pos[1], pos[2]);
                domain1 = -grid[gridindex] - 1;
                if (domain1 != -1) {
                    for (i = 0; i < 13; i++) {
                        neigh[0] = offsets[i][0];
                        neigh[1] = offsets[i][1];
                        neigh[2] = offsets[i][2];
                        neighindex = gridindex + INDEXGRID(
                                neigh[0], neigh[1], neigh[2]);
                        domain2 = -grid[neighindex] - 1;
                        if (domain2 != -1) {
                            intersection_table[domain1 * num_domains + domain2] = 1;
                            intersection_table[domain2 * num_domains + domain1] = 1;
                        } /* if domain2 */
                    } /* for i */
                } /* if domain1 */
            } /* for pos[2] */
        } /* for pos[1] */
    } /* for pos[0] */
}
#undef INDEXGRID


#define INDEXGRID(i,j,k) ((int64_t)(i)*strides[0]+(j)*strides[1]+(k)*strides[2])
/**
 * Take a discretization grid, where cells inside the volume are 0
 * and cells outside are 1. For each outside cell, find the translation
 * vector that leads back inside the volume and set the cells value to
 * -(index + 1). Also make sure that no cell
 * inside has an equivalent (i.e. reachable through a translation vector)
 * cell inside; and that every cell outside has an equivalent cell inside.
 */
void mark_translation_vectors(int8_t *grid,
                              int dimensions[3],
                              int strides[3],
                              int ntranslations,
                              int *translations)
{
    int pos[3];
    int grid_index;
    int grid_value;
    int i, j;
    int *trans_pos;
    int *trans_valid;
    int trans_index;
    int center_dist, min_center_dist;

    trans_pos = malloc(ntranslations * 3 * sizeof(int));
    trans_valid = malloc(ntranslations * sizeof(int));

    for (pos[0] = 0; pos[0] < dimensions[0]; pos[0]++) {
        for (pos[1] = 0; pos[1] < dimensions[1]; pos[1]++) {
            for (pos[2] = 0; pos[2] < dimensions[2]; pos[2]++) {
                grid_index = INDEXGRID(pos[0], pos[1], pos[2]);
                grid_value = grid[grid_index];
                if (grid_value != 0) {
                    continue;
                }
                for (i = 0; i < ntranslations; i++) {
                    trans_valid[i] = 1;
                    for (j = 0; j < 3; j++) {
                        int tp = pos[j] + translations[i * 3 + j];
                        trans_pos[i * 3 + j] = tp;
                        trans_valid[i] &= tp >= 0 && tp < dimensions[j];
                    }
                    if (trans_valid[i]) {
                        grid[INDEXGRID(trans_pos[i * 3 + 0],
                                       trans_pos[i * 3 + 1],
                                       trans_pos[i * 3 + 2])] = 1;
                    }
                }
            }
        }
    }
    for (pos[0] = 0; pos[0] < dimensions[0]; pos[0]++) {
        for (pos[1] = 0; pos[1] < dimensions[1]; pos[1]++) {
            for (pos[2] = 0; pos[2] < dimensions[2]; pos[2]++) {
                grid_index = INDEXGRID(pos[0], pos[1], pos[2]);
                grid_value = grid[grid_index];
                if (grid_value != 1) {
                    continue;
                }
                for (i = 0; i < ntranslations; i++) {
                    trans_valid[i] = 1;
                    for (j = 0; j < 3; j++) {
                        int tp = pos[j] + translations[i * 3 + j];
                        trans_pos[i * 3 + j] = tp;
                        trans_valid[i] &= tp >= 0 && tp < dimensions[j];
                    }
                }
                trans_index = -1;
                for (i = 0; i < ntranslations; i++) {
                    if (trans_valid[i]) {
                        if (grid[INDEXGRID(trans_pos[i * 3 + 0],
                                           trans_pos[i * 3 + 1],
                                           trans_pos[i * 3 + 2])] == 0) {
                            trans_index = i;
                            break;
                        }
                    }
                }
                if (trans_index != -1) {
                    grid[grid_index] = -trans_index - 1;
                } else {
                    min_center_dist = 0;
                    for (j = 0; j < 3; j++) {
                        min_center_dist += SQUARE(pos[j] - dimensions[j] / 2);
                    }
                    for (i = 0; i < ntranslations; i++) {
                        center_dist = 0;
                        for (j = 0; j < 3; j++) {
                            center_dist += SQUARE(trans_pos[i * 3 + j] - dimensions[j] / 2);
                        }
                        if (center_dist < min_center_dist) {
                            trans_index = i;
                            min_center_dist = center_dist;
                        }
                    }
                    if (trans_index != -1) {
                        grid[INDEXGRID(trans_pos[trans_index * 3 + 0],
                                       trans_pos[trans_index * 3 + 1],
                                       trans_pos[trans_index * 3 + 2])] = 0;
                    }
                    /* trans_index == -1: grid[grid_index] = 0 */
                    grid[grid_index] = -trans_index - 1;
                }
            }
        }
    }
    free(trans_pos);
    free(trans_valid);
}
#undef INDEXGRID

