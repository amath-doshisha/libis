#ifndef IS_MACROS_H
#define IS_MACROS_H

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

typedef unsigned long int ulong;

#define MAT(A,I,J,LDA) ((A)[(I)+(J)*(LDA)])
#define MAT3(A,I,J,K,LDA1,LDA2) ((A)[(I)+(J)*(LDA1)+(K)*(LDA1)*(LDA2)])
#define COL(A,J,LDA)   ((A)[(J)*(LDA)])
#define MAX2(X,Y) (((X)>(Y))?(X):(Y))
#define MIN2(X,Y) (((X)<(Y))?(X):(Y))
#define STR_EQ(X,Y) (strcmp((X),(Y))==0)
#define STR_EQ_N(X,Y,N) (strncmp((X),(Y),(N))==0)
#define SWAP2(X,Y,T) { T a; a=(X); (X)=(Y); (Y)=a; }
#define ERROR_EXIT(S,I) { printf("Error in the function %s() at the %dth line in the file '%s'.\n", __func__, __LINE__, __FILE__); printf((S),(I)); exit(-1); }
#define ERROR_AT { printf("Error in the function %s() at the %dth line in the file '%s'.\n", __func__, __LINE__, __FILE__); }
#define NULL_EXC1(X){ if((X)==NULL) { printf("Null Exception occurs in the function %s() at the %dth line in the file '%s'.\n", __func__, __LINE__, __FILE__); exit(-1); } }
#define NULL_EXC2(X,Y){ if((X)==NULL || (Y)==NULL) { printf("Null Exception occurs in the function %s() at the %dth line in the file '%s'.\n", __func__, __LINE__, __FILE__); exit(-1); } }
#define NULL_EXC3(X,Y,Z){ if((X)==NULL || (Y)==NULL || (Z)==NULL) { printf("Null Exception occurs in the function %s() at the %dth line in the file '%s'.\n", __func__, __LINE__, __FILE__); exit(-1); } }
#define NULL_EXC4(X,Y,Z,T){ if((X)==NULL || (Y)==NULL || (Z)==NULL || (T)==NULL) { printf("Null Exception occurs in the function %s() at the %dth line in the file '%s'.\n", __func__, __LINE__, __FILE__); exit(-1); } }


// macros of complex operators
#define Z_REAL(X)          ((X).r)
#define Z_R(X)             ((X).r)
#define Z_IMAG(X)          ((X).i)
#define Z_I(X)             ((X).i)
#define Z_SET(X,R,I)       { Z_R(X) =(R);            Z_I(X) =(I); }
#define Z_LET(Z,X)         { Z_R(Z) =Z_R(X);         Z_I(Z) =Z_I(X); }         // Z =X
#define Z_ADD(Z,X)         { Z_R(Z)+=Z_R(X);         Z_I(Z)+=Z_I(X); }         // Z+=X
#define Z_PLUS_R(X,Y)      (Z_R(X)+Z_R(Y))                                     // (X+Y).r
#define Z_PLUS_I(X,Y)      (Z_I(X)+Z_I(Y))                                     // (X+Y).i
#define Z_SET_PLUS(Z,X,Y)  { Z_R(Z) =Z_PLUS_R(X,Y);  Z_I(Z) =Z_PLUS_I(X,Y); }  // Z =X+Y
#define Z_SUB(Z,X)         { Z_R(Z)-=Z_R(X);         Z_I(Z)-=Z_I(X); }         // Z-=X
#define Z_MINUS_R(X,Y)     (Z_R(X)-Z_R(Y))                                     // (X-Y).r
#define Z_MINUS_I(X,Y)     (Z_I(X)-Z_I(Y))                                     // (X-Y).i
#define Z_SET_MINUS(Z,X,Y) { Z_R(Z) =Z_MINUS_R(X,Y); Z_I(Z) =Z_MINUS_I(X,Y); } // Z =X-Y
#define Z_TIMES_R(X,Y)     (Z_R(X)*Z_R(Y)-Z_I(X)*Z_I(Y))                       // (X*Y).r
#define Z_TIMES_I(X,Y)     (Z_R(X)*Z_I(Y)+Z_I(X)*Z_R(Y))                       // (X*Y).i
#define Z_SET_TIMES(Z,X,Y) { Z_R(Z) =Z_TIMES_R(X,Y); Z_I(Z) =Z_TIMES_I(X,Y); } // Z =X*Y
#define Z_ADD_TIMES(Z,X,Y) { Z_R(Z)+=Z_TIMES_R(X,Y); Z_I(Z)+=Z_TIMES_I(X,Y); } // Z+=X*Y
#define Z_SUB_TIMES(Z,X,Y) { Z_R(Z)-=Z_TIMES_R(X,Y); Z_I(Z)-=Z_TIMES_I(X,Y); } // Z-=X*Y
#define Z_DOT_R(X,Y)       (Z_R(X)*Z_R(Y)+Z_I(X)*Z_I(Y))                       // (conj(X)*Y).r
#define Z_DOT_I(X,Y)       (Z_R(X)*Z_I(Y)-Z_I(X)*Z_R(Y))                       // (conj(X)*Y).i
#define Z_SET_DOT(Z,X,Y)   { Z_R(Z) =Z_DOT_R(X,Y);   Z_I(Z) =Z_DOT_I(X,Y); }   // Z =conj(X)*Y
#define Z_ADD_DOT(Z,X,Y)   { Z_R(Z)+=Z_DOT_R(X,Y);   Z_I(Z)+=Z_DOT_I(X,Y); }   // Z+=conj(X)*Y
#define Z_SUB_DOT(Z,X,Y)   { Z_R(Z)-=Z_DOT_R(X,Y);   Z_I(Z)-=Z_DOT_I(X,Y); }   // Z-=conj(X)*Y
#define Z_DIV_NUM_R(X,Y)   (Z_DOT_R(Y,X))                                      // num(X/Y).r
#define Z_DIV_NUM_I(X,Y)   (Z_DOT_I(Y,X))                                      // num(X/Y).i
#define Z_DIV_DEN(Y)       (Z_ABS2(Y))                                         // den(X/Y)
#define Z_DIV_R(X,Y)       ((Z_DIV_NUM_R(X,Y))/(Z_DIV_DEN(Y)))                 // (X/Y).r
#define Z_DIV_I(X,Y)       ((Z_DIV_NUM_I(X,Y))/(Z_DIV_DEN(Y)))                 // (X/Y).i
#define Z_SET_DIV(Z,X,Y)   { Z_R(Z) =Z_DIV_R(X,Y);   Z_I(Z) =Z_DIV_I(X,Y); }   // Z =X/Y
#define Z_INV_NUM_R(X)     (Z_CONJ_R(X))                                       // num((1/X)).r
#define Z_INV_NUM_I(X)     (Z_CONJ_I(X))                                       // num((1/X)).i
#define Z_INV_DEN(X)       (Z_ABS2(X))                                         // den((1/X))
#define Z_INV_R(X)         ((Z_INV_NUM_R(X))/(Z_INV_DEN(X)))                   // (1/X).r
#define Z_INV_I(X)         ((Z_INV_NUM_I(X))/(Z_INV_DEN(X)))                   // (1/X).r
#define Z_SET_INV(Z,X)     { Z_R(Z) =Z_INV_R(X);     Z_I(Z) =Z_INV_I(X); }     // Z=1/X
#define Z_SCALE(Z,A)       { Z_R(Z)*=(A);            Z_I(Z)*=(A); }            // Z*=A
#define Z_NEG_R(X)         (-(Z_R(X)))
#define Z_NEG_I(X)         (-(Z_I(X)))
#define Z_NEG(Z)           { Z_R(Z) =Z_NEG_R(Z);     Z_I(Z) =Z_NEG_I(Z); }     // Z=-Z
#define Z_SET_NEG(Z,X)     { Z_R(Z) =Z_NEG_R(X);     Z_I(Z) =Z_NEG_I(X); }     // Z=-X
#define Z_CONJ_R(X)        (Z_R(X))
#define Z_CONJ_I(X)        (-(Z_I(X)))
#define Z_CONJ(Z)          { Z_R(Z) =Z_CONJ_R(Z);    Z_I(Z) =Z_CONJ_I(Z); }    // Z=conj(Z)
#define Z_SET_CONJ(Z,X)    { Z_R(Z) =Z_CONJ_R(X);    Z_I(Z) =Z_CONJ_I(X); }    // Z=conj(X)
#define Z_ABS2(X)          (Z_R(X)*Z_R(X)+Z_I(X)*Z_I(X))                       // abs(X)^2
#define Z_ABS(X)           (sqrt(Z_ABS2(X)))                                   // abs(X)
#define Z_ARG(X)           (atan2((Z_I(X)),(Z_R(X))))                          // angle(X)
#define Z_DIST(X,Y)        (sqrt((Z_R(X)-Z_R(Y))*(Z_R(X)-Z_R(Y))+(Z_I(X)-Z_I(Y))*(Z_I(X)-Z_I(Y))))
#define Z_DIST2(X,Y)       ((Z_R(X)-Z_R(Y))*(Z_R(X)-Z_R(Y))+(Z_I(X)-Z_I(Y))*(Z_I(X)-Z_I(Y)))

#endif
//EOF
