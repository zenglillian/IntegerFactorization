//Pollard Rho Method for integer factorization
//This algorithm finds a nontrivial factor of n


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;

void F(mpz_t n,mpz_t x,mpz_t a);
int main()
{
  mpz_t n,a,s,nsub2,g,V,U,UsubV, count;
    mpz_init_set_str (n,"50881147147252369",10); //input an composite number n
    mpz_init(a);
    mpz_init(s);
    mpz_init(nsub2);
    mpz_init(g);
    mpz_init(V);
    mpz_init(U);
    mpz_init(UsubV);
    mpz_init(count);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    
    if(mpz_divisible_ui_p (n,2)!=0){gmp_printf("The non-trivial factor of %Zd is 2 \n",n);}
    do{
        mpz_sub_ui(nsub2,n,2);
        mpz_urandomm (a,state,nsub2);
        if(mpz_cmp_ui(a,0)==0){mpz_add_ui(a,a,1);}
        mpz_urandomm (s,state,n);
        mpz_set(U,s);
        mpz_set(V,s);
        do{
            F(n,U,a);
            F(n,V,a);
            F(n,V,a);
            mpz_sub(UsubV, U,V);
            mpz_gcd (g, UsubV, n);
	    mpz_add_ui(count,count,1);
        } while(mpz_cmp_ui(g,1)==0);
    }while(mpz_cmp(g,n)==0);
    
    gmp_printf("The non-trivial factor of %Zd is %Zd \n",n,g);
    gmp_printf("Count: %Zd \n",count);
    mpz_clear(n);
    mpz_clear(a);
    mpz_clear(s);
    mpz_clear(nsub2);
    mpz_clear(g);
    mpz_clear(V);
    mpz_clear(U);
    mpz_clear(UsubV);
    return 0;
    
}


void F(mpz_t n,mpz_t x,mpz_t a)
{
    mpz_mul(x,x,x);
    mpz_add(x,x,a);
    mpz_mod(x, x,n);
}
