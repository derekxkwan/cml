#include "m_pd.h"
#include <math.h>
#include <string.h>

#ifndef CMLMAX
#define CMLMAX 240000
#endif

#ifndef T_PI
#define T_PI 6.283185307179586
#endif


#define CML_K 0
#define CML_W 0
#define CML_EPS 0.5
#define CML_N 100

static t_class *cmlsynth_tilde_class;

typedef struct _cmlsynth_tilde {
	t_object x_obj;
	int x_sr; //samplerate
	double x_lat[CMLMAX];
  double x_prev;
  int x_latn; //lattice size
  unsigned int x_max; // max lattice size

  double x_eps;
  double x_omg;
  double x_k;
  int x_phase; //pos in buffer

  t_symbol * x_tabname;
  t_word * x_tabvec; //pointer to table contentes
  int x_usetable; //if table loaded
  int x_tabn; //size of loaded table

  int x_npts; //actual read size
  
  t_inlet *x_latnlet;
  t_inlet *x_epslet;
  t_inlet *x_wlet;
  t_inlet *x_klet;
  t_outlet *x_outlet;
} t_cmlsynth_tilde;


//array reading stuff


t_word *cmlsynth_tilde_findarray(t_cmlsynth_tilde *x, t_symbol * name, int *bufsize, int indsp, int *success, int complain){
//in dsp = used in dsp, 
  
    if (name && name != &s_){
        //basically trying to find an array by looking up its symbol
	t_garray *arraypointer = (t_garray *)pd_findbyclass(name, garray_class);
	if (arraypointer){
	    int bufsz;
	    t_word *vec;

            //checks if the array is floats, calls garray_npoints to set bufsz, garray_vec to get contents
	    if (garray_getfloatwords(arraypointer, &bufsz, &vec)){
   	        //b->b_len = garray_npoints(ap);

                //used to set flag of garray struct saying we're gonna use this in dsp
		if (indsp){
                    garray_usedindsp(arraypointer);
                };
		if (bufsize){
                    *bufsize = bufsz;
                };
		*success = 1;
		return (vec);
	    }
	    else{
	      pd_error(x,  /* always complain */
			"bad template of array '%s'", name->s_name);
	      *success = 0;
	    };
        }
	else{
            if(complain){
	        pd_error(x, "no such array '%s'", name->s_name);
		*success = 0;
            };
	};
    }
    return (0);
}

static int cmlsynth_tilde_loadarray(t_cmlsynth_tilde *x, t_symbol * name, int complain)
{
  int success = 0;
  x->x_tabvec = cmlsynth_tilde_findarray(x, name, &x->x_tabn, 1, &success, complain);
  return success;
}
static double cmlsynth_map_ipt(double ipt)
{ // [-1, 1] -> [0,1];

  return (ipt + 1) * 0.5;

}

static void cmlsynth_tilde_copyarray(t_cmlsynth_tilde *x)
{

  int i, readmax = x->x_tabn;

  if(readmax > CMLMAX) readmax = CMLMAX;
  if(x->x_usetable == 1)
    {
      for(i=0; i< readmax; i++)
	{
	  x->x_lat[i] = cmlsynth_map_ipt(x->x_tabvec[i].w_float);
	};
    };
}

static double cmlsynth_tilde_noise(void){//range 0 to 1
	double retval;
	static unsigned long int exprandseed = 1997333137;
	exprandseed = exprandseed * 2891336453 + 1500450271;
	exprandseed = exprandseed % 4294967296;
	exprandseed &= 0x7fffffff;
	retval = ((double)exprandseed)/(4294967295.f*0.5f); //biggest num in 32-bit div 2
	return retval;
}

static void cmlsynth_tilde_setbounds(t_cmlsynth_tilde *x)
{
  int tabsize = 0, maxsize = x->x_max;
  if(x->x_usetable == 0)
    {
      tabsize = x->x_latn;
    }
  else{
    //use smaller of specified lattice size and table size if using table (bounded by x->x_max)
    if(x->x_tabn > x->x_latn) tabsize = x->x_latn;
    else tabsize = x->x_tabn;

  };

  if(tabsize > maxsize) tabsize = maxsize;
  x->x_npts = tabsize;
}

static void cmlsynth_tilde_fillnoise(t_cmlsynth_tilde *x)
{
  int i;
  double randval;
  for(i=0; i < x->x_max; i++)
    {
      randval = cmlsynth_tilde_noise();
      x->x_lat[i] = randval;
      //post("%f", randval);
    };
}

static void cmlsynth_tilde_initbuf(t_cmlsynth_tilde *x)
{
  if(x->x_usetable == 0)
      cmlsynth_tilde_fillnoise(x);
  else cmlsynth_tilde_copyarray(x);

  cmlsynth_tilde_setbounds(x);
  
}

static double cmlsynth_tilde_map(t_cmlsynth_tilde *x, double ipt)
{
  double res;
  res = ipt + x->x_omg - (x->x_k * sin(T_PI * ipt));
  return res - floor(res);
}

static double cmlsynth_tilde_couple(t_cmlsynth_tilde * x, double prev, double cur, double fut)
{
  double eps = x->x_eps;

  return ((1.f - eps)*cmlsynth_tilde_map(x, cur)) + (0.5f * eps *( cmlsynth_tilde_map(x, prev)+ cmlsynth_tilde_map(x, fut)));
}


static double cmlsynth_tilde_calcsample(t_cmlsynth_tilde *x, int phase)
{
  double prev, cur, fut, res;
  

  prev = x->x_prev;
  cur = x->x_lat[phase];
  fut = x->x_lat[(phase + 1) % x->x_npts];
  res = cmlsynth_tilde_couple(x, prev, cur, fut);

  return res;
  
}


static void cmlsynth_tilde_updatebuf(t_cmlsynth_tilde *x, double cur, int phase)
{
  x->x_prev = x->x_lat[phase];
  x->x_lat[phase] = cur;

}

static void cmlsynth_tilde_phasereset(t_cmlsynth_tilde *x)
{

  x->x_phase = 0;
  x->x_prev = x->x_lat[x->x_npts - 1];
}

//calcsample[i] -> out[i] rptrptrpt
static t_int *cmlsynth_tilde_perform(t_int *w){
		 t_cmlsynth_tilde *x = (t_cmlsynth_tilde *)(w[1]);
	t_float *out = (t_float *)(w[2]);
	int block = (int)(w[3]);
	int buflen = x->x_npts;
	int i, phase = x->x_phase;
	double curout;
       
	for(i=0; i< block; i++){
	  //printf("phase: %d\n", phase);
	  curout = cmlsynth_tilde_calcsample(x, phase); // [0,1]
	  cmlsynth_tilde_updatebuf(x,curout, phase);
	  out[i] = ((2.f*curout)-1.f); //[0,1] -> [0,2] -> [-1, 1]
	  //out[i] = curout;
	  phase++;

	  if(phase >= buflen)
	    {
	      cmlsynth_tilde_phasereset(x);
	      phase = x->x_phase;
	    };
	    
	  
	};

	x->x_phase = phase;
	
	return(w+4);
}

static void cmlsynth_tilde_setn(t_cmlsynth_tilde * x, t_float n)
{
  unsigned int latmax = x->x_max;
  int arg = (int)n;

  if(arg < 1) arg = 1;
  if(arg > latmax) arg = latmax;

  x->x_latn = arg;

  cmlsynth_tilde_setbounds(x);
}

static void cmlsynth_tilde_seteps(t_cmlsynth_tilde *x, t_float eps)
{
  if(eps < 0) eps = 0.f;
  if(eps > 1) eps = 1.f;

  x->x_eps = eps;
}

static void cmlsynth_tilde_setw(t_cmlsynth_tilde *x, t_float omg)
{
  if(omg < 0) omg = 0.f;
  if(omg > 1) omg = 1.f;

  x->x_omg = omg;
}

static void cmlsynth_tilde_setk(t_cmlsynth_tilde *x, t_float k)
{
  if(k < 0) k = 0;
  if(k > 1) k = 1;

  x->x_k = k;
}


static void cmlsynth_tilde_setmax(t_cmlsynth_tilde *x, unsigned int curmax)
{
  if(curmax < 1) curmax = 1;
  if(curmax > CMLMAX) curmax = CMLMAX;
  x->x_max = curmax;
}

static void cmlsynth_tilde_doreset(t_cmlsynth_tilde *x, int complain)
{

  int success = 0;
  if(x->x_usetable > 0) success = cmlsynth_tilde_loadarray(x, x->x_tabname, complain);
  x->x_usetable = success;
  cmlsynth_tilde_initbuf(x);
  cmlsynth_tilde_phasereset(x);
}

static void cmlsynth_tilde_reset(t_cmlsynth_tilde *x)
{
  cmlsynth_tilde_doreset(x, 1);
}


static void cmlsynth_tilde_read(t_cmlsynth_tilde *x, t_symbol * s, int argc, t_atom * argv)
{
  int curusetable = x->x_usetable, samename = 0;
  if(argc >= 1)
    {
      if(argv[0].a_type == A_SYMBOL)
	{
	  int success = 0;
	  t_symbol * name = atom_getsymbolarg(0, argc, argv);
	  success = cmlsynth_tilde_loadarray(x, name, 1);
	  if(success == 1)
	    {
	      samename = strcmp(name->s_name, x->x_tabname->s_name) == 0;
	      x->x_usetable = 1;
	      post("cmlsynth~: using %s", name->s_name);
	      if(samename == 0) x->x_tabname = name;
	    };
	}
      else{
	post("cmlsynth~: using noise buffer");
	x->x_usetable = 0;
      };
    }
  else{
    post("cmlsynth~: using noise buffer");
    x->x_usetable = 0;
  };


  if(curusetable != x->x_usetable || !samename)
    {
      cmlsynth_tilde_initbuf(x);
      cmlsynth_tilde_phasereset(x);
    };
}

static void *cmlsynth_tilde_new(t_symbol *s, int argc, t_atom *argv){
  //max lattice, lattice size, eps, omega, k
  t_symbol * name = &s_;
  t_cmlsynth_tilde *x = (t_cmlsynth_tilde *)pd_new(cmlsynth_tilde_class);
	int i = 0;

	int max_latn = CMLMAX;
	int n = CML_N;
	double eps = CML_EPS; //[0,1]
	double omg = CML_W; //[0,1]
	double k = CML_K; //[0,1]
	

	x->x_sr = sys_getsr();
	x->x_usetable = 0;

	while(i < argc)
	  {
	    if(argv[i].a_type == A_FLOAT)
	      {
		double cur_arg = atom_getfloatarg(i, argc, argv);
		//printf("%f\n", cur_arg);
		switch(i)
		  {
		  case 0:
		    max_latn = (int)cur_arg;
		    break;
		  case 1:
		    n = (int)cur_arg;
		    break;
		  case 2:
		    eps = cur_arg;
		    break;
		  case 3:
		    omg = cur_arg;
		    break;
		  case 4:
		    k = cur_arg;
		    break;
		  default:
		    break;
		    
		    
		  };
	      }
	    else if(argv[i].a_type == A_SYMBOL){
	      t_symbol * cur_sym = atom_getsymbolarg(i, argc, argv);
	      if((i+1) < argc)
		{ //if there are at least two args left
		  if(argv[i+1].a_type == A_SYMBOL && strcmp(cur_sym->s_name, "-array") == 0)
		  { //second arg must also be a symbol
		    t_symbol * arr_name = atom_getsymbolarg(i+1, argc, argv);
		    name = arr_name;
		    x->x_usetable = 1;
		    i++;
		  }
		else goto errstate;
	      }
		else goto errstate;
		

	    }

	    else goto errstate;

	    i++;
	  };

	if(name != &s_) x->x_tabname = name;
	cmlsynth_tilde_setmax(x, max_latn);
        cmlsynth_tilde_setn(x, n);
	cmlsynth_tilde_seteps(x, eps);
	cmlsynth_tilde_setw(x, omg);
	cmlsynth_tilde_setk(x, k);
	cmlsynth_tilde_doreset(x, 0);
	x->x_latnlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("n"));
	x->x_epslet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("eps"));
	x->x_wlet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("omega"));
	x->x_klet = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("k"));



	x->x_outlet = outlet_new(&x->x_obj, gensym("signal"));
	return (x);
	errstate:
		pd_error(x, "cmlsynth~: improper args");
		return NULL;
}


static void cmlsynth_tilde_dsp(t_cmlsynth_tilde *x, t_signal **sp){

	//(freq*tablelen)/(samplerate) = array values per sample to advance
	// divide by tablelen to map to 0 to 1 range,..freq/samplerate
  // printf("dsp\n");
    int sr = sp[0]->s_sr; //amount to change phase for freq 1
    cmlsynth_tilde_doreset(x, 1);

	if(sr != x->x_sr){
	  printf("n: %d\n", x->x_latn);
		x->x_sr = sr;

	};
	    dsp_add(cmlsynth_tilde_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);

}

static void *cmlsynth_tilde_free(t_cmlsynth_tilde *x){
	inlet_free(x->x_latnlet);
	inlet_free(x->x_epslet);
	inlet_free(x->x_wlet);
	inlet_free(x->x_klet);

	outlet_free(x->x_outlet);
	
	return (void *)x;
}

void cmlsynth_tilde_setup(void){
	cmlsynth_tilde_class = class_new(gensym("cmlsynth~"), (t_newmethod)cmlsynth_tilde_new, 0,
			sizeof(t_cmlsynth_tilde), CLASS_DEFAULT, A_GIMME, 0);

	   class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_reset, gensym("reset"), 0);

	   class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_read, gensym("read"), A_GIMME, 0);

	class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_dsp, gensym("dsp"), A_CANT, 0);
   class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_setn, gensym("n"), A_FLOAT, 0);
   class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_seteps, gensym("eps"), A_FLOAT, 0);
   class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_setw, gensym("omega"), A_FLOAT, 0);
   class_addmethod(cmlsynth_tilde_class, (t_method)cmlsynth_tilde_setk, gensym("k"), A_FLOAT, 0);

}
