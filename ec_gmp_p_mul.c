/*
	ec_gmp_p_mul.c
	2014 masterzorag@gmail.com

	entirely based on the gmplib implementation of elliptic curve scalar point 
	available at http://researchtrend.net/ijet32/6%20KULDEEP%20BHARDWAJ.pdf

	program sets curve domain parameters, an integer and a point to perform point 
	smultiplication by scalar, computing the derived point

	cryptographically speaking, it verifies private/public math correlation 

	Pass a first arg to verify known R point for first curve
	41da1a8f74ff8d3f1ce20ef3e9d8865c96014fe3
	73ca143c9badedf2d9d3c7573307115ccfe04f13
*/

// Point at Infinity is Denoted by (0,0)
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

struct Elliptic_Curve {
	mpz_t a;
	mpz_t b;
	mpz_t p;
};

struct Point {
	mpz_t x;
	mpz_t y;
};

struct Elliptic_Curve EC;

void Point_Doubling(struct Point P, struct Point *R)
{
	mpz_t slope, temp;
	mpz_init(temp);
	mpz_init(slope);
	
	if(mpz_cmp_ui(P.y, 0) != 0) {
		mpz_mul_ui(temp, P.y, 2);
		mpz_invert(temp, temp, EC.p);
		mpz_mul(slope, P.x, P.x);
		mpz_mul_ui(slope, slope, 3);
		mpz_add(slope, slope, EC.a);
		mpz_mul(slope, slope, temp);
		mpz_mod(slope, slope, EC.p);
		mpz_mul(R->x, slope, slope);
		mpz_sub(R->x, R->x, P.x);
		mpz_sub(R->x, R->x, P.x);
		mpz_mod(R->x, R->x, EC.p);
		mpz_sub(temp, P.x, R->x);
		mpz_mul(R->y, slope, temp);
		mpz_sub(R->y, R->y, P.y);
		mpz_mod(R->y, R->y, EC.p);
	} else {
		mpz_set_ui(R->x, 0);
		mpz_set_ui(R->y, 0);
	}
	mpz_clear(temp);
	mpz_clear(slope);
}

void Point_Addition(struct Point P, struct Point Q, struct Point *R)
{
	mpz_mod(P.x, P.x, EC.p);
	mpz_mod(P.y, P.y, EC.p);
	mpz_mod(Q.x, Q.x, EC.p);
	mpz_mod(Q.y, Q.y, EC.p);

	if(mpz_cmp_ui(P.x, 0) == 0 && mpz_cmp_ui(P.y, 0) == 0) {
		mpz_set(R->x, Q.x);
		mpz_set(R->y, Q.y);
		return;
	}

	if(mpz_cmp_ui(Q.x, 0) == 0 && mpz_cmp_ui(Q.y, 0) == 0) {
		mpz_set(R->x, P.x);
		mpz_set(R->y, P.y);
		return;
	}

	mpz_t temp;
	mpz_init(temp);

	if(mpz_cmp_ui(Q.y, 0) != 0) { 
		mpz_sub(temp, EC.p, Q.y);
		mpz_mod(temp, temp, EC.p);
	} else
		mpz_set_ui(temp, 0);

	//gmp_printf("\n temp=%Zd\n", temp);

	if(mpz_cmp(P.y, temp) == 0 && mpz_cmp(P.x, Q.x) == 0) {
		mpz_set_ui(R->x, 0);
		mpz_set_ui(R->y, 0);
		mpz_clear(temp);
		return;
	}
	
	if(mpz_cmp(P.x, Q.x) == 0 && mpz_cmp(P.y, Q.y) == 0)	{
		Point_Doubling(P, R);
		
		mpz_clear(temp);
		return;		
	} else {
		mpz_t slope;
		mpz_init_set_ui(slope, 0);

		mpz_sub(temp, P.x, Q.x);
		mpz_mod(temp, temp, EC.p);
		mpz_invert(temp, temp, EC.p);
		mpz_sub(slope, P.y, Q.y);
		mpz_mul(slope, slope, temp);
		mpz_mod(slope, slope, EC.p);
		mpz_mul(R->x, slope, slope);
		mpz_sub(R->x, R->x, P.x);
		mpz_sub(R->x, R->x, Q.x);
		mpz_mod(R->x, R->x, EC.p);
		mpz_sub(temp, P.x, R->x);
		mpz_mul(R->y, slope, temp);
		mpz_sub(R->y, R->y, P.y);
		mpz_mod(R->y, R->y, EC.p);
		
		mpz_clear(temp);
		mpz_clear(slope);
		return;
	}
}

void Scalar_Multiplication(struct Point P, struct Point *R, mpz_t m)
{
	struct Point Q, T;
	mpz_init(Q.x); mpz_init(Q.y);
	mpz_init(T.x); mpz_init(T.y);
	long no_of_bits, loop;
	
	no_of_bits = mpz_sizeinbase(m, 2);
	mpz_set_ui(R->x, 0);
	mpz_set_ui(R->y, 0);
	if(mpz_cmp_ui(m, 0) == 0)
		return;
		
	mpz_set(Q.x, P.x);
	mpz_set(Q.y, P.y);
	if(mpz_tstbit(m, 0) == 1){
		mpz_set(R->x, P.x);
		mpz_set(R->y, P.y);
	}

	for(loop = 1; loop < no_of_bits; loop++) {
		mpz_set_ui(T.x, 0);
		mpz_set_ui(T.y, 0);
		Point_Doubling(Q, &T);

		//gmp_printf("\n %Zd %Zd %Zd %Zd ", Q.x, Q.y, T.x, T.y);
		mpz_set(Q.x, T.x);
		mpz_set(Q.y, T.y);
		mpz_set(T.x, R->x);
		mpz_set(T.y, R->y);
		if(mpz_tstbit(m, loop))
			Point_Addition(T, Q, R);
	}
	
	mpz_clear(Q.x); mpz_clear(Q.y);
	mpz_clear(T.x); mpz_clear(T.y);
}		

int main(int argc, char *argv[])
{
	mpz_init(EC.a); 
	mpz_init(EC.b); 
	mpz_init(EC.p);
	
	struct Point P, R;
	mpz_init_set_ui(R.x, 0);
	mpz_init_set_ui(R.y, 0);
	mpz_init(P.x);
	mpz_init(P.y);
	
	mpz_t m;
	mpz_init(m);

	if(argv[1]){
	/*
		Valid test case
		----------------
		
		Curve domain parameters:
		p:	c1c627e1638fdc8e24299bb041e4e23af4bb5427
		a:	c1c627e1638fdc8e24299bb041e4e23af4bb5424
		b:	877a6d84155a1de374b72d9f9d93b36bb563b2ab		
	*/
		mpz_set_str(EC.p, "0xc1c627e1638fdc8e24299bb041e4e23af4bb5427", 0);
		mpz_set_str(EC.a, "0xc1c627e1638fdc8e24299bb041e4e23af4bb5424", 0);
		mpz_set_str(EC.b, "0x877a6d84155a1de374b72d9f9d93b36bb563b2ab", 0);
	/*	
		Base point:
		Gx: 010aff82b3ac72569ae645af3b527be133442131
		Gy: 46b8ec1e6d71e5ecb549614887d57a287df573cc
	*/
		mpz_set_str(P.x, "0x010aff82b3ac72569ae645af3b527be133442131", 0);
		mpz_set_str(P.y, "0x46b8ec1e6d71e5ecb549614887d57a287df573cc", 0);
	/*
		known verified R point for first curve
		R.x	41da1a8f74ff8d3f1ce20ef3e9d8865c96014fe3
		R.y	73ca143c9badedf2d9d3c7573307115ccfe04f13
		using this as 
		k: 00542d46e7b3daac8aeb81e533873aabd6d74bb710
	*/		
		mpz_set_str(m, "0x00542d46e7b3daac8aeb81e533873aabd6d74bb710", 0);
		
	} else {
	/*
		Curve domain parameters:
		p:	dfd7e09d5092e7a5d24fd2fec423f7012430ae9d
		a:	dfd7e09d5092e7a5d24fd2fec423f7012430ae9a
		b:	01914dc5f39d6da3b1fa841fdc891674fa439bd4
		N:	00dfd7e09d5092e7a5d25167ecfcfde992ebf8ecad
	*/
		mpz_set_str(EC.p, "0xdfd7e09d5092e7a5d24fd2fec423f7012430ae9d", 0);
		mpz_set_str(EC.a, "0xdfd7e09d5092e7a5d24fd2fec423f7012430ae9a", 0);
		mpz_set_str(EC.b, "0x01914dc5f39d6da3b1fa841fdc891674fa439bd4", 0);
	/*
		Base point:
		Gx:	70ee7b94f7d52ed6b1a1d3201e2d85d3b82a9810
		Gy:	0b23823cd6dc3df20979373e5662f7083f6aa56f
	*/
		mpz_set_str(P.x, "0x70ee7b94f7d52ed6b1a1d3201e2d85d3b82a9810", 0);
		mpz_set_str(P.y, "0x0b23823cd6dc3df20979373e5662f7083f6aa56f", 0);
	/*
		known verified R point for second curve
		R.x	5432bddd1f97418147aff016eaa6100834f2caa8
		R.y	c498b88965689ee44df349b066cd43cbf4f2c5d0
		problem is found discrete logaritm, unknown k
		so just set the same...
	*/
		mpz_set_str(m, "0x00542d46e7b3daac8aeb81e533873aabd6d74bb710", 0);
	}
	
	/*	p = k x G == R = m x P	*/
	Scalar_Multiplication(P, &R, m);
	
	mpz_out_str(stdout, 16, R.x); puts("");
	mpz_out_str(stdout, 16, R.y); puts("");

	// Free variables
	mpz_clear(EC.a); mpz_clear(EC.b); mpz_clear(EC.p);
	mpz_clear(R.x); mpz_clear(R.y);
	mpz_clear(P.x); mpz_clear(P.y);
	mpz_clear(m);
}