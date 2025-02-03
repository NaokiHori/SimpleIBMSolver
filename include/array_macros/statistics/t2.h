#if !defined(INCLUDE_ARRAY_MACROS_STATISTICS_T2_H)
#define INCLUDE_ARRAY_MACROS_STATISTICS_T2_H

// This file is generated by tools/define_arrays.py

// [0 : isize+1], [1 : jsize+0], [1 : ksize+0]
#define T2(I, J, K) (t2[(I  ) + (isize+2) * ((J-1) + (jsize+0) * (K-1))])
#define T2_NADDS (int [NDIMS][2]){ {1, 1}, {0, 0}, {0, 0}, }

#endif // INCLUDE_ARRAY_MACROS_STATISTICS_T2_H
