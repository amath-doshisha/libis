// 移行後に削除
/**
 * @brief 多倍長型のメモリの割当
 */
array *multi_allocate(char type, int m, int n, int l)
{
  int dim[3];
  dim[0]=m;
  dim[1]=n;
  dim[2]=l;
  return array_allocate(type,3,dim);
}

//EOF
