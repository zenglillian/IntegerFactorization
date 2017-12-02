//Fermat method for integer factorization
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <gmp.h>
#include <gmpxx.h>
using namespace std;

int main()
{
    mpz_t n,sqrtN,a,b,factor;
    mpz_t upbound;
    mpz_init(a);
    mpz_init(b);
    mpz_init(sqrtN);
    mpz_init(upbound);
    mpz_init(factor);
    
    //enter n
    mpz_init_set_str (n,"1591",10);
    //sqrtN: the floor of sqrt(n)
    mpz_sqrt(sqrtN,n);
    if(mpz_perfect_square_p(n)==0){mpz_add_ui(sqrtN,sqrtN,1);}
    //upperbound = (n+9)/6
    mpz_add_ui(upbound,n,9);
    mpz_div_ui(upbound,upbound,6);
    //lowerbound a is sqrt(n)
    mpz_set(a,sqrtN);
    
    while(mpz_cmp(a,upbound)<=0)
    {
        //b = a^2 - n
        mpz_pow_ui(b,a,2);
        mpz_sub(b,b,n);
        if(mpz_perfect_square_p(b)!=0)
        {
            //b = sqrt(a^2 - n)
            mpz_sqrt(b,b);
            //factor = a - b
            mpz_sub(factor,a,b);
            gmp_printf("The non-trivial factor of %Zd is %Zd \n",n,factor);
            
            mpz_clear(n);
            mpz_clear(sqrtN);
            mpz_clear(a);
            mpz_clear(b);
            mpz_clear(factor);
            mpz_clear(upbound);
            return 0;
        }
        //increment a by 1
        mpz_add_ui(a,a,1);
    }
    gmp_printf("%Zd is prime \n",n);
    
    mpz_clear(n);
    mpz_clear(sqrtN);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(factor);
    mpz_clear(upbound);
    return 0;
}

