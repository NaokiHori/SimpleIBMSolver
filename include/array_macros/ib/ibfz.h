#if !defined(INCLUDE_ARRAY_MACROS_IB_IBFZ_H)
#define INCLUDE_ARRAY_MACROS_IB_IBFZ_H

// This file is generated by tools/define_arrays.py

#if NDIMS == 3
// [1 : isize+0], [0 : jsize+1], [0 : ksize+1]
#define IBFZ(I, J, K) (ibfz[(I-1) + (isize+0) * ((J  ) + (jsize+2) * (K  ))])
#define IBFZ_NADDS (int [NDIMS][2]){ {0, 0}, {1, 1}, {1, 1}, }
#endif

#endif // INCLUDE_ARRAY_MACROS_IB_IBFZ_H