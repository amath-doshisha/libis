#include<isys.h>

int main()
{
  bigint *a=NULL,*b=NULL,*c=NULL;

  a=bigint_allocate_int(30,60);
  b=bigint_allocate_str("50","-30");
  c=bigint_allocate();

  printf("a="); bigint_print(a); printf("\n");
  printf("b="); bigint_print(b); printf("\n");
  bigint_mul(c,a,b); printf("c=a*b="); bigint_print(c); printf("\n");
  bigint_div(c,a,b); printf("c=a/b="); bigint_print(c); printf("\n");
  bigint_inv(c,a); printf("c=1/a="); bigint_print(c); printf("\n");
  bigint_inv(c,b); printf("c=1/b="); bigint_print(c); printf("\n");
  bigint_add(c,a,b); printf("c=a+b="); bigint_print(c); printf("\n");
  bigint_sub(c,a,b); printf("c=a-b="); bigint_print(c); printf("\n");
  bigint_sub(c,b,a); printf("c=b-a="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,0); printf("c=b^0="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,1); printf("c=b^0="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,2); printf("c=b^2="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,3); printf("c=b^3="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,-1); printf("c=b^(-1)="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,-2); printf("c=b^(-2)="); bigint_print(c); printf("\n");
  bigint_pow_n(c,b,-3); printf("c=b^(-3)="); bigint_print(c); printf("\n");

  bigint_mul(a,a,b); printf("a=a*b="); bigint_print(a); printf("\n");
  bigint_add(a,a,b); printf("a=a+b="); bigint_print(a); printf("\n");
  bigint_sub(a,a,b); printf("a=a-b="); bigint_print(a); printf("\n");
  bigint_pow_n(a,a,2); printf("a=a^2="); bigint_print(a); printf("\n");
  bigint_inv(a,a); printf("a=1/a="); bigint_print(a); printf("\n");
  bigint_div(a,a,b); printf("a=a/b="); bigint_print(a); printf("\n");

  a=bigint_free(a);
  b=bigint_free(b);
  c=bigint_free(c);

  return 0;
}
