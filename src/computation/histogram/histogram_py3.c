#include <Python.h>
#include <stdio.h>
#include <math.h>

typedef struct {
  double x, y, z;
} point_t;

double sign(double x) 
{
  return (x > 0) - (x < 0);
}

static PyObject *calc_histogram(PyObject *self, PyObject *args)
{
  PyObject *atoms, *atom, *m, *mat;
  point_t *atom_pos, *metric;
  
  int i, j, ntypes, npar, nr, ic, ig;	
  int *ni, *ncum;
  int type1, type2, itype;
  int atom_count, mat_count;
  double x, y, z;
  
  PyObject *histogram;
  PyObject *atoms_it, *mat_it, *ni_it;
  PyObject *item, *_ni;
  double d, dr;
  int truncated;
  
  int *hist;
  
  if (!PyArg_ParseTuple(args, "OOOddi", &atoms, &mat, &_ni, &d, &dr, &truncated))
    return NULL;
  
  atom_count = PyList_Size(atoms);
  atom_pos = malloc(atom_count * sizeof(point_t));
  atoms_it = PyObject_GetIter(atoms);
  
  mat_count = PyList_Size(mat);
  metric = malloc(mat_count * sizeof(point_t));
  mat_it = PyObject_GetIter(mat);
  
  ntypes = PyList_Size(_ni);
  npar = (int)(ntypes*(ntypes + 1) / 2);
  nr = (int)(d / dr) + 1;
  
  ni = malloc(ntypes * sizeof(int));
  ncum = malloc(ntypes * sizeof(int));

  hist = malloc(nr*npar * sizeof(int));
  for (int i = 0; i<nr*npar; i++) 
  {
    hist[i] = 0;
  }
  
  for (i = 0; i < ntypes; i++) 
  {
    ni[i] = PyLong_AsSsize_t(PyList_GetItem(_ni, (Py_ssize_t)i));
  }

  ncum[0] = ni[0];
  for (itype = 1; itype < ntypes; itype++)
  {
    ncum[itype] = ncum[itype - 1] + ni[itype];
  }
  
  i = 0;
  while ((m = PyIter_Next(mat_it)))
  {
    PyArg_Parse(m, "(ddd)", &metric[i].x, &metric[i].y, &metric[i].z);
    //printf("%d %f %f %f\n", i, metric[i].x, metric[i].y, metric[i].z);
    Py_DECREF(m);
    ++i;
  }
  Py_DECREF(mat_it);

  i = 0;
  while ((atom = PyIter_Next(atoms_it)))
  {
    PyArg_Parse(atom, "(ddd)", &atom_pos[i].x, &atom_pos[i].y, &atom_pos[i].z);
    //printf("%d %f %f %f\n", i, atom_pos[i].x, atom_pos[i].y, atom_pos[i].z);
    Py_DECREF(atom);
    ++i;
  }
  Py_DECREF(atoms_it);

  type2 = 0;	
  for (i = 0; i < atom_count; i++)
  {
    if (i > ncum[type2] - 1)
    {
      type2++;
    }
    type1 = 0;
    
    for (j = 0; j < i; j++)
    {
      if (j > ncum[type1] - 1)
      {
	type1++;
      }
      ic = (int)(type1*(2 * ntypes - type1 - 1) / 2 + type2);

      x = atom_pos[i].x - atom_pos[j].x + 3.0;
      y = atom_pos[i].y - atom_pos[j].y + 3.0;
      z = atom_pos[i].z - atom_pos[j].z + 3.0;

      //printf("%f, %f, %f\n", atom_pos[i].x, atom_pos[i].y, atom_pos[i].z);
      //printf("%f, %f, %f\n", atom_pos[j].x, atom_pos[j].y, atom_pos[j].z);      
      //exit(1);
      
      x = 2.0*(x / 2.0 - (int)(x / 2.0)) - 1.0;
      y = 2.0*(y / 2.0 - (int)(y / 2.0)) - 1.0;
      z = 2.0*(z / 2.0 - (int)(z / 2.0)) - 1.0;
      //printf("%f, %f, %f\n", x, y, z);      
      //exit(1);

      if (truncated == 1 && fabs(x) + fabs(y) + fabs(z) > 1.5)
      {
	x = x - sign(x);
	y = y - sign(y);
	z = z - sign(z);
	//printf("ggg\n");
      }
      //printf("%f, %f, %f\n", x, y, z);
      
      d = metric[0].x*x*x + metric[1].y * y*y + metric[2].z * z*z
      	+ 2.0*(metric[0].y * x*y + metric[0].z * x*z + metric[1].z * y*z);

      //printf("%f\n", d);
      //printf("%f %f %f\n", metric[0].x, metric[1].y, metric[2].z);      
      //printf("%f %f %f\n", metric[0].y, metric[0].z, metric[1].z);
      //exit(1);

      d = sqrt(d);
      
      ig = (int)(round(d / dr));
      
      if (ig < nr) 
      {
	//if(i==1) 
	  {
	    //printf("%d, %d, %d\n", i, ig, ic);
	    //exit(1);
	  }
	hist[ig*npar+ic]++;
      }
    }
  }
  
  /*
  for(i=0; i<nr; i++)
  {
    printf("%f", i*dr); 
    for(j=0; j<npar; j++)
    {
      printf(" %d", hist[i*npar+j]);
    }
	printf("\n");
  }
  */
  histogram = PyList_New(nr);
  for (i = 0; i < nr; i++)
  {
    item = PyList_New(npar);
    for (int j = 0; j < npar; j++) 
    {
      PyList_SetItem(item, j, Py_BuildValue("i", hist[i*npar+j]));
      //PyList_Append(item, Py_BuildValue("(ddd)", i, i, i));
      //PyList_Append(centers_list, item);
    }
    PyList_SetItem(histogram, i, item);
  }

  return histogram;
}

static PyMethodDef histgramMethods[] = {
    {"calc_histogram", calc_histogram, METH_VARARGS, "Calculates histogram"},
    {NULL, NULL, 0, NULL}        /* end marker */
};

// myModule definition struct
static struct PyModuleDef calculate_histogram_def = {
    PyModuleDef_HEAD_INIT,
    "calc_histogram",
    "Python3 C API Module",
    -1,
    histgramMethods
};

PyMODINIT_FUNC PyInit_histogram(void)
{
    return PyModule_Create(&calculate_histogram_def);
}
