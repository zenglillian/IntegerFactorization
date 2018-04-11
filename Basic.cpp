//Filename: Basic.cpp
//Created by Lillian Zeng
//a basic algorithm for integer factorization derived by Dr. Hiary.
/*This program should be link against gmp library, for example
 g++ Basic.cpp -lgmp -lgmpxx */

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
  mpz_t n, sqrtN,x0, x, H, factor;
  mpz_init(factor);
  mpz_init_set_str (n,"9797",10);
  //x0 = min(ceil(pow(17*n),1/3.0)),sqrt(n))
  mpz_init (x0);
  mpz_mul_ui (x0,n,17);
  mpz_root (x0, x0,3);
  mpz_init (sqrtN);
  mpz_sqrt(sqrtN, n);
  if(mpz_cmp(x0, sqrtN)>0){mpz_set(x0,sqrtN);}
  //x = x0 + 2
  mpz_init (x);
  mpz_add_ui (x,x0,2);
  //H = 1
  mpz_init_set_str (H,"1",10);
  //trialDivision
  trialDivision (n,x0,factor);
  if(mpz_cmp_ui(factor,0) != 0){gmp_printf("The non-trivial factor of %Zd is %Zd \n",n,factor);}
  else
    {
      findFactorLoop (x,n,sqrtN,H,factor);
      if(mpz_cmp_ui(factor,0) == 0){gmp_printf("%Zd is prime \n",n);}
      else{gmp_printf("The non-trivial factor of %Zd is %Zd \n",n,factor);}
    }
  mpz_clear(n);
  mpz_clear(sqrtN);
  mpz_clear(x0);
  mpz_clear(x);
  mpz_clear(H);
  mpz_clear(factor);
  
  return 0;
}


void findFactorLoop (mpz_t x, mpz_t n,mpz_t sqrtN, mpz_t H, mpz_t factor)
{
  mpz_t t;
  mpz_t b,q,a,s,h,n17,XsubH,delta,d,Hm4,Xsq,quad_a,quad_b,quad_c,quad_b2,quad_2a,quad_4ac,test,Xq,Xcubed;
  mpz_init(t); 
  mpz_init(b);
  mpz_init(q);
  mpz_init(a);
  mpz_init(n17); 
  mpz_init (s); 
  mpz_init (delta);
  mpz_init (Hm4);
  mpz_init (Xq);
  mpz_init(XsubH); 
  mpz_init (h);
  mpz_init (d);
  mpz_init(Xsq);
  mpz_init(quad_a);
  mpz_init(quad_b);
  mpz_init(quad_c);
  mpz_init(quad_b2);
  mpz_init(quad_2a);
  mpz_init(quad_4ac);
  mpz_init(test);
  mpz_init(Xcubed);
  mpz_mul_ui(n17,n,17);
  mpz_sub(XsubH,x,H);
  while (mpz_cmp(XsubH,sqrtN)<=0)
    {
      //find the convergent with the largest qj<=4H
      mpz_mul_ui(Hm4,H,4);  //Hm4
      mpz_pow_ui(Xq,x,2);  //Xsq
      rationalApprox(n,Xq,Hm4, b,q);
      //a = q*n/x;
      mpz_mul(a,q,n);
      mpz_fdiv_q (a,a,x);
      //solve
      //delta = (b*x-a)*(b*x-a)-4*b*(q*n-a*x);
      mpz_set(quad_a,b);//quad_a = b
      mpz_mul(quad_b,b,x); 
      mpz_sub (quad_b,quad_b,a);//quad_b=bx-a
      mpz_pow_ui(quad_b2,quad_b,2);//quad_b2
      mpz_mul(quad_c,q,n); 
      mpz_mul(t,a,x); //t
      mpz_sub(quad_c,quad_c,t); //quad_c = qn-ax
      mpz_mul(quad_4ac,quad_c,quad_a); 
      mpz_mul_ui(quad_4ac,quad_4ac,4); //quad_4ac
      mpz_sub(delta,quad_b2,quad_4ac);//delta = quad_b2-quad_4ac
      //check int of sqrt(delta),delta>=0 and quad_a!=0
      if (mpz_cmp_ui(delta,0)>=0 && mpz_cmp_ui(quad_a,0)!=0)
	{ 
	  //s = sqrt(delta);
	  mpz_sqrt(s,delta); 
	  mpz_pow_ui(test,s,2); //test whether sqrt(delta) is integer
	  if(mpz_cmp(test,delta)==0)
	    {
	      //t1 = 2*b;
	      mpz_mul_ui(quad_2a,quad_a,2); //quad_2a
	      //h
	      mpz_neg(h,quad_b);
	      mpz_sub(h,h,s);
	      if(mpz_divisible_p(h,quad_2a)!=0)
		{
		  mpz_fdiv_q(h,h,quad_2a);
		  mpz_add(h,h,x);
		  //check if x+h divides n
		  if (mpz_divisible_p(n,h) != 0)
		    {
		      mpz_set (factor,h);
		      mpz_clear(b);
		      mpz_clear(q);
		      mpz_clear(a);
		      mpz_clear(n17);
		      mpz_clear (s);
		      mpz_clear (delta);
		      mpz_clear (XsubH);
		      mpz_clear (h);
		      mpz_clear(Xsq);
		      mpz_clear(quad_a);
		      mpz_clear(quad_b);
		      mpz_clear(quad_c);
		      mpz_clear(quad_b2);
		      mpz_clear(quad_2a);
		      mpz_clear(quad_4ac);
		      mpz_clear(t);
		      mpz_clear(test);
		      return;
		    }
		}
	      
	      //h ;
	      mpz_add(h,h,s);
	      mpz_add(h,h,s);
	      if(mpz_divisible_p(h,quad_2a)!=0)
		{
		  mpz_fdiv_q(h,h,quad_2a);
		  mpz_add(h,h,x);
		  //check if x+h divides n
		  if (mpz_divisible_p(n,h) != 0)
		    {
		      mpz_set (factor,h);
		      mpz_clear(b);
		      mpz_clear(q);
		      mpz_clear(a);
		      mpz_clear(h);
		      mpz_clear(n17);
		      mpz_clear (s);
		      mpz_clear (delta);
		      mpz_clear (XsubH);
		      mpz_clear (h);
		      mpz_clear(Xsq);
		      mpz_clear(quad_a);
		      mpz_clear(quad_b);
		      mpz_clear(quad_c);
		      mpz_clear(quad_b2);
		      mpz_clear(quad_2a);
		      mpz_clear(quad_4ac);
		      mpz_clear(t);
		      mpz_clear(test);
		      return;
		    }
		}
	    }
	}
      //increment
      //H = floor (pow(n17, (-1/3.0))*x); 
      mpz_pow_ui (Xcubed,x,3); //Xcubed
      mpz_fdiv_q(H, Xcubed,n17);
      mpz_root (H, H,3);
      // x = x + 2 * H + 1;
      mpz_add(x,x,H);
      mpz_add(x,x,H);
      mpz_add_ui(x,x,1);
      mpz_sub(XsubH,x,H);
    }
  mpz_set_ui (factor,0);
  mpz_clear(b);
  mpz_clear(q);
  mpz_clear(a);
  mpz_clear(n17);
  mpz_clear (s);
  mpz_clear (delta);
  mpz_clear (XsubH);
  mpz_clear (h);
  mpz_clear (d);
  mpz_clear(Xsq);
  mpz_clear(quad_a);
  mpz_clear(quad_b);
  mpz_clear(quad_c);
  mpz_clear(quad_b2);
  mpz_clear(quad_2a);
  mpz_clear(quad_4ac);
  mpz_clear(t);
  mpz_clear(test);
  mpz_clear(Xq);
  mpz_clear(Xcubed);
  
  return;
}

void rationalApprox( mpz_t m,mpz_t n,mpz_t H, mpz_t & A, mpz_t & B)
{
  
  mpz_t A_old, A_new, B_old, B_new, a,temp, tempA, tempB,mm,nn;
  mpz_init_set (mm,m);
  mpz_init_set (nn,n);
  mpz_init (A_old);
  mpz_init_set_ui (A_new,1);  
  mpz_init_set_ui (B_old,1);
  mpz_init (B_new); 
  mpz_init (a);
  mpz_init (temp);
  mpz_init (tempA);
  mpz_init (tempB);
  while(mpz_cmp (B_new,H)<=0 && mpz_cmp(B_new,n)<0)
    {
      //a = m/n;
      mpz_fdiv_q (a,mm,nn);
      //temp = m;
      mpz_set (temp, mm);
      //m = n;
      mpz_set (mm, nn);
      //n = temp % n;
      mpz_mod (nn,temp,nn);
      //tempA = A_new;
      mpz_set (tempA, A_new);
      //tempB = B_new;
      mpz_set (tempB, B_new);
      //A_new = a*A_new+A_old;
      mpz_mul(A_new, a, A_new);
      mpz_add(A_new, A_old, A_new);
      //B_new = a*B_new+B_old;
      mpz_mul(B_new, a, B_new);
      mpz_add(B_new, B_old, B_new);
      //A_old=tempA;
      //B_old=tempB;
      mpz_set (A_old, tempA);
      mpz_set (B_old, tempB);
    }
  mpz_set(A,A_old);
  mpz_set(B,B_old);
  
  mpz_clear (mm);
  mpz_clear (nn);
  mpz_clear (A_old);
  mpz_clear (A_new);
  mpz_clear (B_old);
  mpz_clear (B_new);
  mpz_clear (a);
  mpz_clear (temp);
  mpz_clear (tempA);
  mpz_clear (tempB);
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



