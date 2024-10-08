/**
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
#include "df1b2fun.h"
#include "adrndeff.h"
#include "admodel.h"

double better_rand(long int&);

double function_minimizer::get_monte_carlo_value([[maybe_unused]] int nvar,
  const independent_variables& x)
{
  //initial_params::xinit(x);
  double f=0.0;
  if (mcmc2_flag==0 && lapprox)
  {
    dvector g=(*lapprox)(x,f,this);
  }
  else
  {
    dvariable vf=0.0;
    vf=initial_params::reset(dvar_vector(x));
    *objective_function_value::pobjfun=0.0;
    userfunction();
    vf+=*objective_function_value::pobjfun;
    f=value(vf);
  }
  return f;
}

#if !defined(USE_ADPVM)
double function_minimizer::pvm_master_get_monte_carlo_value([[maybe_unused]] int nvar,
  [[maybe_unused]] const dvector& x){return 0.0;}
#else
double function_minimizer::pvm_master_get_monte_carlo_value([[maybe_unused]] int nvar,
  [[maybe_unused]] const dvector& x)
{
  double f=0.0;
  {
    dvariable vf=0.0;
    dvar_vector vx=dvar_vector(x);
    vf=initial_params::reset(vx);
    send_x_to_slaves(vx);
    *objective_function_value::pobjfun=0.0;
    userfunction();
    vf+=*objective_function_value::pobjfun;
    dvar_matrix fslave=get_f_from_slaves();
    vf+=sum(fslave);
    f=value(vf);
  }
  return f;
}
#endif

#if !defined(USE_ADPVM)
void function_minimizer::pvm_slave_get_monte_carlo_value([[maybe_unused]] int nvar){}
#else
void function_minimizer::pvm_slave_get_monte_carlo_value([[maybe_unused]] int nvar)
{
  double f=0.0;
  {
    dvariable vf=0.0;
    dvar_vector vx=get_x_from_master();
    vf=initial_params::reset(vx);
    *objective_function_value::pobjfun=0.0;
    userfunction();
    vf+=*objective_function_value::pobjfun;
    send_int_to_master(67);
    send_f_to_master(vf);
  }
}
#endif

double function_minimizer::get_monte_carlo_value(int nvar,
  const independent_variables& x,dvector& g)
{
  //initial_params::xinit(x);
  double f=0.0;
  if (mcmc2_flag==0 && lapprox)
  {
    g=(*lapprox)(x,f,this);
  }
  else
  {
    dvariable vf=0.0;
    vf=initial_params::reset(dvar_vector(x));
    *objective_function_value::pobjfun=0.0;
    userfunction();
    vf+=*objective_function_value::pobjfun;
    f=value(vf);
    gradcalc(nvar,g);
  }
  return f;
}
