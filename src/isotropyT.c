#include "isotropyT.h"
#include "common.h"
#include "gfft.h"
#include "snoopy.h"

void isotropyT(struct Field fld, double *returnT){
	int i,j,k;	// Dummy counter variables

	// A_ijkl = <u_i u_j u_k u_l>
	double A1111=0.0, A2222=0.0, A3333=0.0;
	double A1112=0.0, A1122=0.0, A1222=0.0;
	double A1113=0.0, A1133=0.0, A1333=0.0;
	double A2223=0.0, A2233=0.0, A2333=0.0;
	double A1123=0.0, A1223=0.0, A1233=0.0;

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fld.vx[i];
		w2[i] =  fld.vy[i];
		w3[i] =  fld.vz[i];
	}

	projector(w1,w2,w3);
	
	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);

	// Compute v_i v_j components

#ifdef _OPENMP
	#pragma omp parallel for private(i) schedule(static)	
#endif
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
#ifndef WITH_2D
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
#endif
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
#ifndef WITH_2D
	gfft_r2c_t(wr6);
#endif
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	/* Done with Convolution for v_i v_j computation */
	
	// Now actually compute A_ijkl components
	for(i=0; i < NX_COMPLEX/NPROC; i++){
		for(j=0; j < NY_COMPLEX; j++){
			for(k=0; k < NZ_COMPLEX; k++){
#ifdef WITH_2D
					if(j==0) {
#else
					if(k==0) {
#endif
						// k=0, we have all of the modes
						A1111 = A1111 +     creal(conj(w4[IDX3D]) * w4[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2222 = A2222 +     creal(conj(w5[IDX3D]) * w5[IDX3D] / ((double) NTOTAL*NTOTAL));
						A3333 = A3333 +     creal(conj(w6[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1112 = A1112 +     creal(conj(w4[IDX3D]) * w7[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1122 = A1122 +     creal(conj(w4[IDX3D]) * w5[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1222 = A1222 +     creal(conj(w7[IDX3D]) * w5[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1113 = A1113 +     creal(conj(w4[IDX3D]) * w8[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1133 = A1133 +     creal(conj(w4[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1333 = A1333 +     creal(conj(w8[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2223 = A2223 +     creal(conj(w5[IDX3D]) * w9[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2233 = A2233 +     creal(conj(w5[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2333 = A2333 +     creal(conj(w9[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1123 = A1123 +     creal(conj(w4[IDX3D]) * w9[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1223 = A1223 +     creal(conj(w7[IDX3D]) * w9[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1233 = A1233 +     creal(conj(w7[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
					} else {
						// k>0, one half of the complex plane is represented
						A1111 = A1111 + 2.0*creal(conj(w4[IDX3D]) * w4[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2222 = A2222 + 2.0*creal(conj(w5[IDX3D]) * w5[IDX3D] / ((double) NTOTAL*NTOTAL));
						A3333 = A3333 + 2.0*creal(conj(w6[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1112 = A1112 + 2.0*creal(conj(w4[IDX3D]) * w7[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1122 = A1122 + 2.0*creal(conj(w4[IDX3D]) * w5[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1222 = A1222 + 2.0*creal(conj(w7[IDX3D]) * w5[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1113 = A1113 + 2.0*creal(conj(w4[IDX3D]) * w8[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1133 = A1133 + 2.0*creal(conj(w4[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1333 = A1333 + 2.0*creal(conj(w8[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2223 = A2223 + 2.0*creal(conj(w5[IDX3D]) * w9[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2233 = A2233 + 2.0*creal(conj(w5[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A2333 = A2333 + 2.0*creal(conj(w9[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1123 = A1123 + 2.0*creal(conj(w4[IDX3D]) * w9[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1223 = A1223 + 2.0*creal(conj(w7[IDX3D]) * w9[IDX3D] / ((double) NTOTAL*NTOTAL));
						A1233 = A1233 + 2.0*creal(conj(w7[IDX3D]) * w6[IDX3D] / ((double) NTOTAL*NTOTAL));
				}
			}
		}
	}

	// Get contributions from all the processes
	reduce(&A1111, 1);
	reduce(&A1112, 1);
	reduce(&A1222, 1);
	reduce(&A1122, 1);
	reduce(&A1123, 1);
	reduce(&A2222, 1);
	reduce(&A3333, 1);
	reduce(&A1113, 1);
	reduce(&A1133, 1);
	reduce(&A1333, 1);
	reduce(&A2223, 1);
	reduce(&A2233, 1);
	reduce(&A2333, 1);
	reduce(&A1223, 1);
	reduce(&A1233, 1);

	// <u^4> = <u_j u_j u_k u_k>
	double u4 = (A1111 + A1122 + A1133) + (A1122 + A2222 + A2233) + (A1133 + A2233 + A3333);

	// Normalise Aijkl
	A1111 = A1111 / u4;
	A2222 = A2222 / u4;
	A3333 = A3333 / u4;
	A1112 = A1112 / u4;
	A1122 = A1122 / u4;
	A1222 = A1222 / u4;
	A1113 = A1113 / u4;
	A1133 = A1133 / u4;
	A1333 = A1333 / u4;
	A2223 = A2223 / u4;
	A2233 = A2233 / u4;
	A2333 = A2333 / u4;
	A1123 = A1123 / u4;
	A1223 = A1223 / u4;
	A1233 = A1233 / u4;

	// A_ij = - A_ijll / 7
	double A11, A22, A33;
	double A12, A13, A23;

	A11 = -(A1111 + A1122 + A1133)/7.0;
	A22 = -(A1122 + A2222 + A2233)/7.0;
	A33 = -(A1133 + A2233 + A3333)/7.0;
	A12 = -(A1112 + A1222 + A1233)/7.0;
	A13 = -(A1113 + A1223 + A1333)/7.0;
	A23 = -(A1123 + A2223 + A2333)/7.0;

	// A = - A_ii/5
	double A = -(A11 + A22 + A33)/5.0;

	// T_ijkl = A_ijkl + A_ij d_kl + A_ik d_jl + A_il d_jk + A_jk d_il
	//                 + A_jl d_ik + A_kl d_ij + A (d_ij d_kl + d_ik d_jl + d_il d_jk)
	// where d_ij is the Kronecker delta

	// Assumming cyclic symmetry reduces the number of independent components
	// to T1111, T1112, T2221, T1122, T1123
	double T1111 = A1111 + 6.0*A11 + 3.0*A;
	double T1112 = A1112 + 3.0*A12;
	double T1222 = A1222 + 3.0*A12;
	double T1122 = A1122 + A11 + A22 + A;
	double T1123 = A1123 + A23;
	// But I'll compute the remaining ones anyway out of a sense of masochism...
	double T2222 = A2222 + 6.0*A22 + 3.0*A;
	double T3333 = A3333 + 6.0*A33 + 3.0*A;
	double T1113 = A1113 + 3.0*A13;
	double T1133 = A1133 + A11 + A33 + A;
	double T1333 = A1333 + 3.0*A13;
	double T2223 = A2223 + 3.0*A23;
	double T2233 = A2233 + A22 + A33 + A;
	double T2333 = A2333 + 3.0*A23;
	double T1223 = A1223 + A13;
	double T1233 = A1233 + A12;

	// Compute T_ijkl * T_ijkl
	// Note the combinatorial factors
	double Tmag =  T1111*T1111 + T2222*T2222 + T3333*T3333
	             + 4.0*T1112*T1112 + 6.0*T1122*T1122 + 4.0*T1222*T1222
							 + 4.0*T1113*T1113 + 6.0*T1133*T1133 + 4.0*T1333*T1333
							 + 4.0*T2223*T2223 + 6.0*T2233*T2233 + 4.0*T2333*T2333
							 + 12.0*T1123*T1123 + 12.0*T1223*T1223 + 12.0*T1233*T1233;

	// On these components, the fact that T_ijkl is traceless by construction
	// implies T1111 = - 2 T1122
	// and T1112 + T2221 + T1123 = 0
	// Hence we shall only return three components
	returnT[0] = T1111;
	returnT[1] = T1112;
	returnT[2] = T1123;
	// Ok, and return the scalar too...
	returnT[3] = Tmag;
	returnT[4] = u4;
}
