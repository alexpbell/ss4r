/*
 * $Id$
 *
 * Author: David Fournier
 * Copyright (c) 2008-2012 Regents of the University of California
 */
#include "admodel.h"

void ad_make_code_reentrant(void)
{
  initial_params::num_initial_params=0;
  initial_params::mc_phase=0;
  initial_params::mceval_phase=0;
  initial_params::sd_phase=0;
  likeprof_params::num_likeprof_params=0;
  stddev_params::num_stddev_params=0;
  if (ad_comm::global_datafile)
  {
    delete ad_comm::global_datafile;
    ad_comm::global_datafile=NULL;
  }

  if (ad_comm::global_parfile)
  {
    delete ad_comm::global_parfile;
    ad_comm::global_parfile=NULL;
  }

  if (ad_comm::global_bparfile)
  {
    delete ad_comm::global_bparfile;
    ad_comm::global_bparfile=NULL;
  }
}

void parse_dll_options(char *pname, const int& _nopt, char *sp_options,
  char *** pargv)
{
  int& nopt=(int&) _nopt;
  nopt=1;
  char ** a=*pargv;
  char * p=NULL;
  if (strlen(sp_options))
  {
    size_t string_size = (strlen(sp_options)+5)*sizeof(char);
    //int string_size= (strlen(sp_options)+1)*sizeof(char);
    p = (char *) malloc(string_size);
    char * ptmp=NULL;
    strcpy(p,sp_options);
    ptmp = strtok(p," ");
    if (ptmp) nopt++;
    //ssbul_l(ptmp);
    do
    {
      ptmp = strtok(NULL, " ");
      if (ptmp) nopt++;
    }
    while(ptmp);

    //a=(char **)malloc((nopt+5)*sizeof(char *));
    //for (int k=0;k<=nopt;k++) a[k]=NULL;
    nopt=1;
    strcpy(p,sp_options);
    ptmp = strtok(p," ");
    if (ptmp)
    {
      //a[nopt]=(char *)malloc((strlen(ptmp)+5)*sizeof(char));
      strcpy(a[nopt++],ptmp);
    }

    do
    {
      ptmp = strtok(NULL, " ");
      if (ptmp)
      {
        //a[nopt]=(char *)malloc((strlen(ptmp)+5)*sizeof(char));
        strcpy(a[nopt++],ptmp);
      }
    }
    while(ptmp);

    if (p)
    {
      free(p);
      p=NULL;
    }
  }
/*
  else
  {
    //a=(char **)malloc((nopt+5)*sizeof(char *));
    //a[nopt]=NULL;
  }
*/
  //a[0]=(char *)malloc((strlen(pname)+10)*sizeof(char));
  strcpy(a[0],pname);
#ifdef _WIN32
  strcat(a[0],".exe");
#endif
  //return a;
}

char** no_dll_options(char *pname, const int& _nopt)
{
  int& nopt=(int&) _nopt;
  nopt=1;
  unsigned int nopt2 = 2;
  char** a = (char**)malloc(nopt2 * sizeof(char*));
  if (a)
  {
    a[nopt] = NULL;
    a[0] = (char*)malloc((strlen(pname)+5)*sizeof(char));
    if (a[0])
    {
      strcpy(a[0],pname);
#ifdef _WIN32
      strcat(a[0],".exe");
#endif
    }
  }
  return a;
}

void cleanup_argv(int nopt,char *** pa)
{
  if (pa)
  {
    char **a=*pa;
    if (a)
    {
      for (int i=0;i<=nopt;i++)
      {
        if (a[i])
        {
          free(a[i]);
        }
      }
      free(a);
      a=NULL;
    }
  }
}

#if defined(_WIN32)
  #include <windows.h>
#endif

  #if !defined(_MSC_VER)
void get_sp_printf(void)
{
/*
  HINSTANCE h=LoadLibrary("sqpe.dll");
  if(h)
    ad_printf= (fptr) GetProcAddress(h,"S_newio_printf");
*/
}
  #endif

void do_dll_housekeeping(int argc,char ** argv)
{
#if !defined(_MSC_VER)
  int on = option_match(argc,argv,"-sp");
  if (on > -1)
  {
  #if defined(_WIN32)
    get_sp_printf();
  #endif
    ad_exit=spdll_exit;
  }
  else if ( (on=option_match(argc,argv,"-spexit"))>-1)
  {
    ad_exit=spdll_exit;
  }
  else
  {
    ad_exit=exit;
  }
#else
  ad_exit=exit;

#endif
}

/*
void ssbul_l(char * ptmp){;}
char ** parse_dll_options(char * pname, const int& _nopt,char * sp_options)
{
  int& nopt=(int&) _nopt;
  nopt=1;
  char ** a;

  a=(char**)malloc(61);
  for (int i=0;i<60;i++)
  {
    a[i]=(char*)malloc(101*sizeof(char*));
  }

  char * p=NULL;
  if (strlen(sp_options))
  {
    int string_size= (strlen(sp_options)+5)*sizeof(char);
    //int string_size= (strlen(sp_options)+1)*sizeof(char);
    p = (char *) malloc(string_size);
    char * ptmp=NULL;
    strcpy(p,sp_options);
    ptmp = strtok(p," ");
    if (ptmp) nopt++;
    ssbul_l(ptmp);
    do
    {
      ptmp=strtok('\0'," ");
      if (ptmp) nopt++;
    }
    while(ptmp);

    //a=(char **)malloc((nopt+5)*sizeof(char *));
    //for (int k=0;k<=nopt;k++) a[k]=NULL;
    nopt=1;
    strcpy(p,sp_options);
    ptmp = strtok(p," ");
    if (ptmp)
    {
      //a[nopt]=(char *)malloc((strlen(ptmp)+5)*sizeof(char));
      strncpy(a[nopt++],ptmp,99);
    }
    do
    {
      ptmp=strtok('\0'," ");
      if (ptmp)
      {
        //a[nopt]=(char *)malloc((strlen(ptmp)+5)*sizeof(char));
        strncpy(a[nopt++],ptmp,99);
      }
      if (nopt>59)
      {
        cerr << "too many command line options in parse command line options" << endl;
        ad_exit(1);
      }
    }
    while(ptmp);
    delete p;
    p=NULL;
  }
  else
  {
    //a=(char **)malloc((nopt+5)*sizeof(char *));
    //a[nopt]=NULL;
  }
  //a[0]=(char *)malloc((strlen(pname)+10)*sizeof(char));
  strncpy(a[0],pname,99);
  strncat(a[0],".exe",99-strlen(a[0]));
  return a;
}
*/
/*
class string_parser
{
  char * str;
  char * nt;
  char * ct;
  string_parser(char * s,int n);
};

string_parser::string_parser(char * s,int n)
{
  if (strlen(s)>n)
  {
    cerr "string too long in string_parser" << endl;
    ad_exit(1);
  }

  str=new char[strlen(s+20];
  strcpy(str,s);

  ct=str;

  while (ct==' ')
  {
    ct++;
  }

  nt=ct;

  int ii=0;
  while (nt != ' ' && nt != '\0')
  {
    nt++;
    ii++;
  }
  nlen=ii;

  while (nt==' ')
  {
    nt++;
  }
}

char * string_parser::get_next_option(int & n)
{
  n=nlen;
  char * rs=ct;

  if (nlen!=0)
  {
    ct=nt;
    int ii=0;
    while (nt != ' ' && nt != '\0')
    {
      nt++;
      ii++;
    }
    nlen=ii;

    while (nt==' ')
    {
      nt++;
    }
  }
  return rs;
}

char** parse_dll_options(char *pname, const int& _nopt, char *sp_options)
{
  int& nopt=(int&) _nopt;
  nopt=1;
  char ** a;

  a=(char**)malloc(61);
  for (int i=0;i<60;i++)
  {
    a[i]=(char*)malloc(101*sizeof(char*));
  }
  do
  {
    s=get_next_option(sp,n);
    if(n==0) break;
    strncpy(a[ii++],sp,100);
  }
  while(1);



  strncpy(a[0],pname,99);
  strncat(a[0],".exe",99-strlen(a[0]));
  strncpy(a[1],"-maxfn",99);
  nopt=2;
  return a;
}
*/
class string_parser
{
  char * str;
  char * nt;
  char * ct;
  size_t nlen;
public:
  string_parser(char* s, const size_t n);
  ~string_parser();

  char* get_next_option(size_t& n);
};

string_parser::string_parser(char * s, const size_t n)
{
  if (strlen(s) > n)
  {
    cerr << "string too long in string_parser" << endl;
    ad_exit(1);
  }

  str=new char[strlen(s)+20];
  strcpy(str,s);

  ct=str;

  while (*ct==' ')
  {
    ct++;
  }

  nt=ct;

  size_t ii=0;
  while (*nt != ' ' && *nt != '\0')
  {
    nt++;
    ii++;
  }
  nlen=ii;

  while (*nt==' ')
  {
    nt++;
  }
}
/**
Destructor
*/
string_parser::~string_parser()
{
  if (str != NULL)
  {
    delete [] str;
    str = NULL;
  }
  nt = NULL;
  ct = NULL;
  nlen = 0;
}

char* string_parser::get_next_option(size_t& n)
{
  n=nlen;

  char * rt=ct;
  if (nlen!=0)
  {
    ct=nt;
    size_t ii=0;
    while (*nt != ' ' && *nt != '\0')
    {
      nt++;
      ii++;
    }
    nlen=ii;

    while (*nt==' ')
    {
      nt++;
    }
  }
  return rt;
}

typedef char * chararray;

void davesnobullshitstrncpy(char * a,const char * b, const size_t n)
{
  strncpy(a,b,n);
  a[n]=0;
}

char** parse_dll_options(char *pname, const int& _argc, char *s)
{
  const int maxargs=200;
  char** a = (char**)malloc(sizeof(char*)*(maxargs+2));
  if (a)
  {
    int& argc=(int&)(_argc);
    for (int i=0;i<=maxargs;i++)
    {
      a[i]=0;
    }
    string_parser sp(s,1000);

    size_t n=0;
    int ii=1;
    do
    {
      char * ss=sp.get_next_option(n);
      if (n>0)
      {
        if (ii>maxargs)
        {
          cerr << "maximum number of command line arguments exceeded"
               << endl;
        }
        a[ii]=(char*)malloc(sizeof(char)*(n+2));
        davesnobullshitstrncpy(a[ii],ss,n);
        cout << a[ii] << endl;
        ii++;
      }
    } while (n > 0);
    a[0]=(char*)malloc(sizeof(char)*(strlen(pname)+6));
    if (a[0])
    {
      strcpy(a[0],pname);
#ifdef _WIN32
      strcat(a[0],".exe");
#endif
      argc=ii;
    }
  }
  return a;
}
