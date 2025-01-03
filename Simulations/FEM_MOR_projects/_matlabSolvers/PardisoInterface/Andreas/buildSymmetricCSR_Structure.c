

void buildSymmetricCSR_Structure(double* A_values_sym, int* A_ja_sym, int* A_ia_sym, int n_col, 
                                 double* A_valuesR, mwIndex* A_ja, mwIndex* A_ia)
{
  size_t i, j, k, kk;
  bool diagExist;
  /* build symmetric CSR structure */  
  k = 0;
  A_ia_sym[0] = 1;
  for (i = 0 ; i < n_col; i++) {
    diagExist = false;
    kk = 0;
    for (j = A_ia[i]; j < A_ia[i + 1]; j++) {
    /* upper part */
      if (A_ja[j] == i)
        /* diagonal element */
        diagExist = true;
      if (A_ja[j] >= i) {
        if (diagExist == false) {
          A_ja_sym[k] = (int)(i + 1);
          A_values_sym[k] = 0.0;          
          k += 1;
          kk += 1;
          diagExist = true;
        }
        A_ja_sym[k] = (int)(A_ja[j] + 1);
        A_values_sym[k] = A_valuesR[j];
        k += 1;
        kk += 1;
      }
    }
    if (diagExist == false) {
      A_ja_sym[k] = (int)(i + 1);
      A_values_sym[k] = 0.0;
      k += 1;
      kk += 1;
    }
    A_ia_sym[i + 1] =  (int)(A_ia_sym[i] + kk);
  }
}


void buildSymmetricCSR_StructureComplex(MKL_Complex16* A_values_sym, int* A_ja_sym, int* A_ia_sym, int n_col, 
                                        double* A_valuesR, double* A_valuesI, mwIndex* A_ja, mwIndex* A_ia)
{
  size_t i, j, k, kk;
  bool diagExist;
  /* build symmetric CSR structure */  
  k = 0;
  A_ia_sym[0] = 1;
  for (i = 0; i < n_col; i++) {
    diagExist = false;
    kk = 0;
    for (j = A_ia[i]; j < A_ia[i + 1]; j++) {
    /* upper part */
      if (A_ja[j] == i)
        /* diagonal element */
        diagExist = true;
      if (A_ja[j] >= i) {
        if (diagExist == false) {
          A_ja_sym[k] = (int)(i + 1);
          A_values_sym[k].real = 0.0;          
          A_values_sym[k].imag = 0.0;          
          k += 1;
          kk += 1;
          diagExist = true;
        }
        A_ja_sym[k] = (int)(A_ja[j] + 1);
        A_values_sym[k].real = A_valuesR[j];
        A_values_sym[k].imag = A_valuesI[j];
        k += 1;
        kk += 1;
      }
    }
    if (diagExist == false) {
      A_ja_sym[k] = (int)(i + 1);
      A_values_sym[k].real = 0.0;
      A_values_sym[k].imag = 0.0;
      k += 1;
      kk += 1;
    }
    A_ia_sym[i + 1] =  (int)(A_ia_sym[i] + kk);
  }
}


MKL_INT numbersOfZerosOnDiag(mwIndex* A_ja, mwIndex* A_ia, MKL_INT nnz, MKL_INT n_row)
{
  MKL_INT rowCnt, colCnt, nnzCnt, numZerosOnDiag = 0, nnzOnDiag = 0;
  #if DEBUG
  for (nnzCnt = 0; nnzCnt < nnz; nnzCnt++)
    printf("A_ja[%d] = %d\n", nnzCnt, A_ja[nnzCnt]);
  for (rowCnt = 0; rowCnt <= n_row; rowCnt++)
    printf("A_ia[%d] = %d\n", rowCnt, A_ia[rowCnt]);
  #endif
  for (rowCnt = 0; rowCnt < n_row; rowCnt++)
    for (colCnt = A_ia[rowCnt]; colCnt < A_ia[rowCnt + 1]; colCnt++)
      if (A_ja[colCnt] == rowCnt)
        nnzOnDiag++;
  numZerosOnDiag = n_row - nnzOnDiag;
  #if DEBUG
  printf("numZerosOnDiag = %d\n", numZerosOnDiag);
  #endif
  return numZerosOnDiag;
}

