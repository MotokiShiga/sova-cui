#ifdef __APPLE__
#include <Python.h>
#else
#include <Python.h>
#endif
#include <stdio.h>
#include <limits.h>

typedef struct {
  double x, y, z;
} point_t;

static PyObject* calculate_atomsdata_list(PyObject* self, PyObject* args) 
{
    PyObject* atoms_data, *atom, *positions, *ipos, * jpos;
    int atom_count;
    int i, j;
    point_t* atom_pos;
    double xi, yi, zi, xj, yj, zj;
    PyObject* atoms_it;
    PyObject* atoms_list;

    if (!PyArg_ParseTuple(args, "O", &positions))
        return NULL;
    
    atom_count = PyList_Size(positions);
    atoms_it = PyObject_GetIter(positions);
    atom_pos = malloc(atom_count * sizeof(point_t));

    i = 0;
    while ((atom = PyIter_Next(atoms_it))) {
        PyArg_Parse(atom, "(ddd)", &atom_pos[i].x, &atom_pos[i].y, &atom_pos[i].z);
        ++i;
        Py_DECREF(atom);
    }
    Py_DECREF(atoms_it);

    atoms_list = PyList_New(0);

    for (i = 0; i < atom_count; ++i)
    {
        for (j = atom_count - 1; j > i; --j)
        {
            if (atom_pos[i].x != atom_pos[j].x) continue;
            if (atom_pos[i].y != atom_pos[j].y) continue;
            if (atom_pos[i].z != atom_pos[j].z) continue;

            PyList_Append(atoms_list, Py_BuildValue("(ii)", i, j));

            //printf("%d %d %f %f %f\n", i, j, atom_pos[j].x, atom_pos[j].y, atom_pos[j].z);
        }
    }

    return atoms_list;

    /*
    for (i = 0; i < atoms_count; ++i)
    {
        ipos = PySequence_GetItem(positions, i);
        PyArg_Parse(ipos, "(ddd)", &xi, &yi, &zi);
        printf("%d\n", i);

        for (j = atoms_count-1; j > i; --j)
        {
            jpos = PySequence_GetItem(positions, j);
            PyArg_Parse(jpos, "(ddd)", &xj, &yj, &zj);
            
            if (xi == xj && yi == yj && zi == zj)
            {
                printf("%d %d\n", i, j);
            }
            //atom_count = PySequence_Length(atom);
            Py_DECREF(jpos);
        }
        Py_DECREF(ipos);
    } 
    Py_DECREF(positions);
    */
    //Py_RETURN_NONE;
}

//removeindices

static PyObject* search_removeindices(PyObject* self, PyObject* args)
{
    PyObject* atoms_data, * atom, * positions, * id, * lp;
    PyObject* start_idx, * length;
    int atoms_count, positions_count;
    int i, j, k, l, si, sk;
    point_t* atom_pos;
    double xi, yi, zi, xj, yj, zj;
    int* idx, * len;
    PyObject* atoms_it, * idx_it, * len_it;
    PyObject* atoms_list;

    if (!PyArg_ParseTuple(args, "OOO", &positions, &start_idx, &length))
        return NULL;
    
    atoms_count = PyList_Size(start_idx); 

    positions_count = PyList_Size(positions);
    atoms_it = PyObject_GetIter(positions);
    atom_pos = malloc(positions_count * sizeof(point_t));
    idx = malloc(atoms_count * sizeof(int));
    idx_it = PyObject_GetIter(start_idx);
    len = malloc(atoms_count * sizeof(int));
    len_it = PyObject_GetIter(length);
    
    i = 0;
    while ((atom = PyIter_Next(atoms_it))) {
        PyArg_Parse(atom, "(ddd)", &atom_pos[i].x, &atom_pos[i].y, &atom_pos[i].z);        
        ++i;
        Py_DECREF(atom);
    }
    Py_DECREF(atoms_it);
    
    i = 0;
    while ((id = PyIter_Next(idx_it))) {
        PyArg_Parse(id, "i", &idx[i]);
        ++i;
        Py_DECREF(id);
    }
    Py_DECREF(idx_it);
    
    i = 0;
    while ((lp = PyIter_Next(len_it))) {
        PyArg_Parse(lp, "i", &len[i]);
        ++i;
        Py_DECREF(lp);
    }
    Py_DECREF(len_it);
    
    atoms_list = PyList_New(0);

    for (i = 0; i < atoms_count; ++i)
    {
        si = idx[i];
        for (j = 0; j < len[i]; ++j)
        {
            for (k = i + 1; k < atoms_count; ++k)
            {
                sk = idx[k];
                for (l = 0; l < len[k]; ++l)
                {
                    if (atom_pos[si + j].x != atom_pos[sk + l].x) continue;
                    if (atom_pos[si + j].y != atom_pos[sk + l].y) continue;
                    if (atom_pos[si + j].z != atom_pos[sk + l].z) continue;

                    PyList_Append(atoms_list, Py_BuildValue("(ii)", k, l));
                }
            }
        }
    }
    return atoms_list;    

    //Py_RETURN_NONE;
}

static PyMethodDef calculate_atomsdataMethods[] = {
    {"calculate_atomsdata_list", calculate_atomsdata_list, METH_VARARGS, "Calculates atoms data"},
    {"search_removeindices", search_removeindices, METH_VARARGS, "Search remove indices"},
    {NULL, NULL, 0, NULL}        /* end marker */
};

static struct PyModuleDef calculate_atomsdataModule =
{
    PyModuleDef_HEAD_INIT,
    "calculate_atomsdata",          /* name of module */
    "",                         /* module documentation, may be NULL */
    -1,                         /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    calculate_atomsdataMethods
};

PyMODINIT_FUNC
PyInit_calculate_atomsdata(void)
{
    (void)PyModule_Create(&calculate_atomsdataModule);
}

