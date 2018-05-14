#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdarg.h>
#include"is_macros.h"
#include"is_rmulti.h"
#include"is_irmulti.h"
#include"is_irmat3.h"
#include"is_irmat.h"

/**
 @file  irmat3.c
 @brief 多倍長精度実数型irmultiの3次元配列に関する関数の定義
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */

/** @name irmulti型の3次元配列値の設定に関する関数 */
/** @{ */

/**
 @brief rmulti型の3次元配列の値のコピー [B0,B1]=[A0,A1]
 */
void irmat3_copy(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irset_r(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief rmulti型の3次元配列の値を倍精度実数から設定 B=rmulti(A).
 */
void irmat3_set_d(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, double *A, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
        irset_d(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief rmulti型の3次元配列の値を倍精度実数から設定 C=rmulti(A,B).
 */
void irmat3_set_dd(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, double *A, int LDA1, int LDA2, double *B, int LDB1, int LDB2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
        irset_d(MAT3(C0,i,j,k,LDC1,LDC2),MAT3(C1,i,j,k,LDC1,LDC2),MAT3(A,i,j,k,LDA1,LDA2),MAT3(B,i,j,k,LDB1,LDB2));
      }
    }
  }
}

/** @} */
/** @name irmulti型の3次元配列に関する関数 */
/** @{ */

/**
 @brief 区間の中心 [m-r,m+r]=[A0,A1]
 */
void irmat3_mid(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irmid(MAT3(mid,i,j,k,LD1mid,LD2mid),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief 区間の半径 [m-r,m+r]=[A0,A1]
 */
void irmat3_rad(int m, int n, int l, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irrad(MAT3(rad,i,j,k,LD1rad,LD2rad),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[A0,A1]
 */
void irmat3_mr(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irmr(MAT3(mid,i,j,k,LD1mid,LD2mid),MAT3(rad,i,j,k,LD1rad,LD2rad),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief irmulti型の3次元配列の符号反転 B=-A.
 */
void irmat3_neg(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irneg_r(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief 足し算 [Z0,Z1]=[X0,X1]+[Y0,Y1]
 */
void irmat3_add(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k; 
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	iradd_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }  
}

/**
 @brief 足し算 [Z0,Z1]=[X0,X1]+b
 */
void irmat3_add_r(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	iradd_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
}

/**
 @brief 引き算 [Z0,Z1]=[X0,X1]-[Y0,Y1]
 */
void irmat3_sub(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irsub_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
}

/**
 @brief 引き算 [Z0,Z1]=[x0,x1]-[Y0,Y1]
 */
void irmat3_sub_r1(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irsub_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),X0,X1,MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
}

/**
 @brief 引き算 [Z0,Z1]=[X0,X1]-[y0,y1]
 */
void irmat3_sub_r2(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irsub_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
}

/**
 @brief 掛け算 [Z0,Z1]=[X0,X1]*[Y0,Y1]
 */
void irmat3_mul(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irmul_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
}

/**
 @brief 掛け算 [Z0,Z1]=[X0,X1]*[y0,y1]
 */
void irmat3_mul_r(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irmul_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
}

/**
 @brief 割り算 [Z0,Z1]=[X0,X1]./[Y0,Y1]
 */
void irmat3_div(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irdiv_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
}

/**
 @brief 割り算 [Z0,Z1]=[x0,y1]./[Y0,Y1]
 */
void irmat3_div_r1(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irdiv_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),X0,X1,MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
}


/**
 @brief 割り算 [Z0,Z1]=[X0,Y1]./[x0,y1]
 */
void irmat3_div_r2(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irdiv_rr(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
}

/**
 @brief irmulti型の3次元配列の絶対値 C=abs(A).
 */
void irmat3_abs(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irabs_r(MAT3(C0,i,j,k,LDC1,LDC2),MAT3(C1,i,j,k,LDC1,LDC2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}


/**
 @brief 3次元配列の平方根 C=sqrt(A).
 */
void irmat3_sqrt(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	irsqrt_r(MAT3(C0,i,j,k,LDC1,LDC2),MAT3(C1,i,j,k,LDC1,LDC2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
}

/**
 @brief irmulti型の3次元配列の和 B=sum(A).
 */
void irmat3_sum(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti**A1, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    irvec_sum_irmat(m,n,&MAT3(B0,0,0,k,LDB1,LDB2),&MAT3(B1,0,0,k,LDB1,LDB2),&MAT3(A0,0,0,k,LDA1,LDA2),&MAT3(A1,0,0,k,LDA1,LDA2),LDA1);
  }
}

/**
 @brief irmulti型の3次元配列の最大値 [B0,B1]=[max(A0),max(A1)].
 */
void irmat3_max(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    irvec_max_irmat(m,n,&MAT3(B1,0,0,k,LDB1,LDB2),&MAT3(B0,0,0,k,LDB1,LDB2),&MAT3(A0,0,0,k,LDA1,LDA2),&MAT3(A1,0,0,k,LDA1,LDA2),LDA1);
  }
}

/**
 @brief irmulti型の3次元配列の最大値 [B0,B1]=[A0,max(A1)].
 */
void irmat3_umax(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    irvec_umax_irmat(m,n,&MAT3(B1,0,0,k,LDB1,LDB2),&MAT3(B0,0,0,k,LDB1,LDB2),&MAT3(A0,0,0,k,LDA1,LDA2),&MAT3(A1,0,0,k,LDA1,LDA2),LDA1);
  }
}

/**
 @brief irmulti型の3次元配列の最小値 [B0,B1]=[min(A0),min(A1)].
 */
void irmat3_min(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int k;
  for(k=0; k<l; k++){
    irvec_min_irmat(m,n,&MAT3(B1,0,0,k,LDB1,LDB2),&MAT3(B0,0,0,k,LDB1,LDB2),&MAT3(A0,0,0,k,LDA1,LDA2),&MAT3(A1,0,0,k,LDA1,LDA2),LDA1);
  }
}

/** @} */
