//Filename: findFactor.cpp
//This program finds a non-trivial factor of n or proves that n is prime

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;

void findFactorLoop (mpz_t x, mpz_t n,mpz_t sqrtN,mpz_t H, mpz_t factor);
void rationalApprox(mpz_t m,mpz_t n,mpz_t H,mpz_t & A, mpz_t & B);
void trialDivision (mpz_t n, mpz_t v,mpz_t & factor);

int main()
{
  mpz_t n, sqrtN,factor;
  mpz_init(factor);
  mpz_init_set_str (n,"32193549881843477",10);
  //179424673*179424691
  mpz_init (sqrtN);
  mpz_sqrt(sqrtN, n);

  //trialDivision
  trialDivision (n, sqrtN,factor);
  if(mpz_cmp_ui(factor,0) != 0){cout<< "The non-trivial factor of "<< n <<" is "<< factor<<endl;}
  else{cout<< n << " is prime"<<endl;}
  mpz_clear(n);
  mpz_clear(sqrtN);
  mpz_clear(factor);

  return 0;
}


void trialDivision (mpz_t n, mpz_t v,mpz_t &factor)
{
  mpz_t i,a,b,c,d,e,f,g,h;
  mpz_init (i);
  mpz_init (a);
  mpz_init (b);
  mpz_init (c); 
  mpz_init (d);
  mpz_init (e);
  mpz_init (f);
  mpz_init (g);
  mpz_init (h);
  //when n < 30
  if (mpz_cmp_ui(n,30)<0)
    {
      if ( mpz_divisible_ui_p(n,2) != 0 && mpz_cmp_ui (n,2)!=0 && mpz_cmp_ui (v,2) >= 0)
	    {mpz_set_ui (factor,2);return;}
      else if (mpz_divisible_ui_p(n,3) != 0 && mpz_cmp_ui (n,3)!=0 && mpz_cmp_ui (v,3) >= 0)
	    {mpz_set_ui (factor,2);return;}
      else if(mpz_divisible_ui_p(n,5) != 0 && mpz_cmp_ui (n,5)!=0 && mpz_cmp_ui (v,5) >= 0)
	    {mpz_set_ui (factor,2);return;}
      else{return;}
    }
  //when n > 30
  else
    {//check whether n can be divided by 2,3,5 
      if (mpz_divisible_ui_p(n,2) != 0)
	    {if (mpz_cmp_ui (v,2) >= 0){mpz_set_ui (factor,2);}return;}
      else if (mpz_divisible_ui_p(n,3) != 0)
	    {if (mpz_cmp_ui (v,3) >= 0){mpz_set_ui (factor,3);}return;}
      else if(mpz_divisible_ui_p(n,5) != 0)
	    {if (mpz_cmp_ui (v,5) >= 0){mpz_set_ui (factor,5);}return;}
      //if n cannot be divided by 2,3,5, then n cannot be divided by all the multiples of 2,3,5
      else
    	{
    	  for(mpz_set_ui(i,0);mpz_cmp(i,v)<=0;mpz_add_ui(i,i,30))
    	    {
    	      if (mpz_cmp_ui(i,30)>=0){mpz_add_ui(a,i,1);}
    	      else{mpz_add_ui(a,i,2);}
    	      mpz_add_ui(b,i,7);
    	      mpz_add_ui(c,i,11);
    	      mpz_add_ui(d,i,13);
    	      mpz_add_ui(e,i,17);
    	      mpz_add_ui(f,i,19);
    	      mpz_add_ui(g,i,23);
    	      mpz_add_ui(h,i,29);
    	      if (mpz_divisible_p(n,a)!= 0)
        		{   
        		  if(mpz_cmp(a,v)<=0){mpz_set(factor,a);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,b)!= 0)
        		{
        		  if(mpz_cmp(b,v)<=0){mpz_set(factor,b);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,c)!= 0)
        		{
        		  if(mpz_cmp(c,v)<=0){mpz_set(factor,c);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,d)!= 0)
        		{
        		  if(mpz_cmp(d,v)<=0){mpz_set(factor,d);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,e)!= 0)
        		{
        		  if(mpz_cmp(e,v)<=0){mpz_set(factor,e);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,f)!= 0)
        		{
        		  if(mpz_cmp(f,v)<=0){mpz_set(factor,f);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,g)!= 0)
        		{
        		  if(mpz_cmp(g,v)<=0){mpz_set(factor,g);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	      else if (mpz_divisible_p(n,h)!= 0)
        		{
        		  if(mpz_cmp(h,v)<=0){mpz_set(factor,h);}
        		  else{mpz_set_ui(factor,0);}
        		  return;
        		}
    	    }
    	  mpz_set_ui(factor,0);
    	  return;
    	}  
    }
}


