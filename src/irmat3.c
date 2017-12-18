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
 @file  irmat.c
 @brief 多倍長精度実数型rmultiの3次元配列に関する関数の定義
 @details rmulti型の関数に関する定義は@link rmulti.c@endlinkを参照のこと.
          rmulti型のベクトルの関数に関する定義は@link rvec.c@endlinkを参照のこと.
          rmulti型の行列の関数に関する定義は@link rmat.c@endlinkを参照のこと.
 */

/** @name rmulti型の3次元配列の初期化に関する関数 */
/** @{ */


/** @} */


/////////////////////////////////////////////////////

/** @name rmulti型の3次元配列のメンバ変数に関する関数 */
/** @{ */

/** @} */


////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の入出力に関する関数 */
/** @{ */

/** @} */

/////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の値の設定に関する関数 */

/**
 @brief rmulti型の3次元配列の値のコピー [B0,B1]=[A0,A1]
*/
int irmat3_copy(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=ircopy(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief rmulti型の3次元配列の値を倍精度実数から設定 B=rmulti(A).
*/
int irmat3_set_d(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, double *A, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irset_d(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/** @} */


////////////////////////////////////////////////////////////////////////


/** @name rmulti型の3次元配列の型変換に関する関数 */
/** @{ */

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name その他 */
/** @{ */


/**
 @brief 区間の中心 [m-r,m+r]=[A0,A1]
*/
int irmat3_mid(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irmid(MAT3(mid,i,j,k,LD1mid,LD2mid),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief 区間の半径 [m-r,m+r]=[A0,A1]
*/
int irmat3_rad(int m, int n, int l, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irrad(MAT3(rad,i,j,k,LD1rad,LD2rad),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief 区間の中心と半径 [m-r,m+r]=[A0,A1]
*/
int irmat3_mr(int m, int n, int l, rmulti **mid, int LD1mid, int LD2mid, rmulti **rad, int LD1rad, int LD2rad, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irmr(MAT3(mid,i,j,k,LD1mid,LD2mid),MAT3(rad,i,j,k,LD1rad,LD2rad),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

//編集済み
/**
 @brief irmulti型の3次元配列の符号反転 B=-A.
*/
int irmat3_neg(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irneg(MAT3(B0,i,j,k,LDB1,LDB2),MAT3(B1,i,j,k,LDB1,LDB2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}
//ここまで
/**
 @brief 足し算 [Z0,Z1]=[X0,X1]+[Y0,Y1]
*/
int irmat3_add(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=iradd(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
  return e;
}
//編集済み
/**
@brief 足し算 [Z0,Z1]=[X0,X1]+b
*/
int irmat3_add_r(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=iradd(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
  return e;
}

/**
 @brief 引き算 [Z0,Z1]=[X0,X1]-[Y0,Y1]
*/
int irmat3_sub(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irsub(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
  return e;
}

/**
@brief 引き算 [Z0,Z1]=[x0,x1]-[Y0,Y1]
*/
int irmat3_sub_r1(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irsub(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),X0,X1,MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
  return e;
}

/**
@brief 引き算 [Z0,Z1]=[X0,X1]-[y0,y1]
*/
int irmat3_sub_r2(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irsub(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
  return e;
}

/**
@brief 掛け算 [Z0,Z1]=[X0,X1]*[Y0,Y1]
*/
int irmat3_mul(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irmul(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
  return e;
}

/**
@brief 掛け算 [Z0,Z1]=[X0,X1]*[y0,y1]
*/
int irmat3_mul_r(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irmul(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
  return e;
}

/**
 @brief 割り算 [Z0,Z1]=[X0,X1]./[Y0,Y1]
*/
int irmat3_div(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irdiv(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
  return e;
}

/**
 @brief 割り算 [Z0,Z1]=[x0,y1]./[Y0,Y1]
*/
int irmat3_div_r1(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti *X0, rmulti *X1, rmulti **Y0, rmulti **Y1, int LDY1, int LDY2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irdiv(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),X0,X1,MAT3(Y0,i,j,k,LDY1,LDY2),MAT3(Y1,i,j,k,LDY1,LDY2));
      }
    }
  }
  return e;
}


/**
 @brief 割り算 [Z0,Z1]=[X0,Y1]./[x0,y1]
*/
int irmat3_div_r2(int m, int n, int l, rmulti **Z0, rmulti **Z1, int LDZ1, int LDZ2, rmulti **X0, rmulti **X1, int LDX1, int LDX2, rmulti *Y0, rmulti *Y1)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irdiv(MAT3(Z0,i,j,k,LDZ1,LDZ2),MAT3(Z1,i,j,k,LDZ1,LDZ2),MAT3(X0,i,j,k,LDX1,LDX2),MAT3(X1,i,j,k,LDX1,LDX2),Y0,Y1);
      }
    }
  }
  return e;
}

/** @} */

////////////////////////////////////////////////////////////////////////

/** @name rmulti型の3次元配列の数学関数に関する関数 */
/** @{ */

/**
 @brief 3次元配列の平方根 C=sqrt(A).
*/
int irmat3_sqrt(int m, int n, int l, rmulti **C0, rmulti **C1, int LDC1, int LDC2, rmulti **A0, rmulti **A1, int LDA1, int LDA2)
{
  int i,j,k,e=0;
  for(k=0; k<l; k++){
    for(j=0; j<n; j++){
      for(i=0; i<m; i++){
	e+=irsqrt(MAT3(C0,i,j,k,LDC1,LDC2),MAT3(C1,i,j,k,LDC1,LDC2),MAT3(A0,i,j,k,LDA1,LDA2),MAT3(A1,i,j,k,LDA1,LDA2));
      }
    }
  }
  return e;
}

/**
 @brief irmulti型の3次元配列の和 B=sum(A).
*/
int irmat3_sum(int m, int n, int l, rmulti **B0, rmulti **B1, int LDB1, int LDB2, rmulti **A0, rmulti**A1, int LDA1, int LDA2)
{
  int k,e=0;
  for(k=0; k<l; k++){
    e+=irvec_sum_irmat(m,n,&MAT3(B0,0,0,k,LDB1,LDB2),&MAT3(B1,0,0,k,LDB1,LDB2),&MAT3(A0,0,0,k,LDA1,LDA2),&MAT3(A1,0,0,k,LDA1,LDA2),LDA1);
  }
  return e;
}

//ここまで
/** @} */
