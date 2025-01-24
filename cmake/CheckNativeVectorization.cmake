# The C source code in this file is subject to the following license:
# MIT License

# Copyright (c) Microsoft Corporation. All rights reserved.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

macro(check_native_vectorization)
  # Check for native vectorization
  include(CheckCSourceRuns)
  check_c_source_runs("
  #include <immintrin.h>
  int main ()
  {
    volatile __m128 a, b, c;
    const float src[4] = { 1.0f, 2.0f, 3.0f, 4.0f };
    float dst[4];
    a = _mm_loadu_ps( src );
    b = _mm_loadu_ps( src );
    c = _mm_hadd_ps( a, b );
    _mm_storeu_ps( dst, c );
    if( dst[0] != 3.0f || dst[1] != 7.0f || dst[2] != 3.0f || dst[3] != 7.0f ){
      return -1;
    }
    return 0;
  }"
  SSE3_WORKS)

    check_c_source_runs("
    #include <immintrin.h>
    int main()
    {
        volatile __m256 a, b, c;
        const float src[8] = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
        float dst[8];
        a = _mm256_loadu_ps( src );
        b = _mm256_loadu_ps( src );
        c = _mm256_add_ps( a, b );
        _mm256_storeu_ps( dst, c );
        for( int i = 0; i < 8; i++ ){
            if( ( src[i] + src[i] ) != dst[i] ){
                return -1;
            }
        }
        return 0;
    }"
    AVX_WORKS)

    check_c_source_runs("
    #include <immintrin.h>
    int main()
    {
      volatile __m256i a, b, c;
      const int src[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
      int dst[8];
      a =  _mm256_loadu_si256( (__m256i*)src );
      b =  _mm256_loadu_si256( (__m256i*)src );
      c = _mm256_add_epi32( a, b );
      _mm256_storeu_si256( (__m256i*)dst, c );
      for( int i = 0; i < 8; i++ ){
        if( ( src[i] + src[i] ) != dst[i] ){
          return -1;
        }
      }
      return 0;
    }"
    AVX2_WORKS)

    check_c_source_runs("
        #include <immintrin.h>
        int main()
        {
            volatile __m512d a, b, c;
            const double src[8] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
            double dst[8];
            a = _mm512_loadu_pd( src );
            b = _mm512_loadu_pd( src );
            c = _mm512_add_pd( a, b );
            _mm512_storeu_pd( (__m512d*) dst, c );
            for( int i = 0; i < 8; i++ ){
                if( ( src[i] + src[i] ) != dst[i] ){
                    return -1;
                }
            }
            return 0;
        }"
    AVX512F_WORKS)
endmacro()
