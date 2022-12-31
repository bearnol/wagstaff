#include <iostream>
#include <gmpxx.h>
#include <gmp.h>

using namespace std;


int Rabin_Miller(mpz_class n)
/* given an integer n >= 3 returns 0 composite, 1 probably prime */
{
	return mpz_probab_prime_p(n.get_mpz_t(), 10);
}

mpz_class Wagstaff(mpz_class n)
/* calculate Wagstaff number from exponent */
{
	mpz_class i=1;
	mpz_class wagstaff=1;
	
	while (i <= n) {
		wagstaff*=2;
		i++;
	}
	wagstaff++;
	wagstaff/=3;

	return wagstaff;
}

mpz_class exp_mod(mpz_class x, mpz_class b, mpz_class n)
/* returns x ^ b mod n */
{
	mpz_t temp;
	mpz_init (temp);
	mpz_powm(temp, x.get_mpz_t(), b.get_mpz_t(), n.get_mpz_t());
	mpz_class temp_class(temp);
	mpz_clear (temp);
	return temp_class;
}

int main()
{
	mpz_class i=1;
	mpz_class wagstaff;
	
	while (i <= 10000) {
		i++;
		cout << i << '\r';
		if (!Rabin_Miller(i)) continue;
		wagstaff = Wagstaff(i);
		if (!(exp_mod(25, (wagstaff*3-2), (wagstaff*3))==25)) continue;
		cout << '\n';
	}
	return 0;
}
