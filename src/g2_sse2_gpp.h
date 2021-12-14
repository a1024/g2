#ifndef     G2_SSE2_GPP_H
#define     G2_SSE2_GPP_H
#include    <smmintrin.h>//SSE4.1
#ifdef      __GNUC__
inline int& m128i_i32(__m128i &v){return *(int*)&v;}
inline int m128i_i32(__m128i const &v){return *(int*)&v;}
inline int& m128i_i32(__m128i &v, int component){return ((int*)&v)[component];}
inline int m128i_i32(__m128i const &v, int component){return ((int*)&v)[component];}

inline float& m128_f32(__m128 &v){return *(float*)&v;}
inline float m128_f32(__m128 const &v){return *(float*)&v;}
inline float& m128_f32(__m128 &v, int component){return ((float*)&v)[component];}
inline float m128_f32(__m128 const &v, int component){return ((float*)&v)[component];}

inline double& m128d_v0(__m128d &v){return *(double*)&v;}
inline double m128d_v0(__m128d const &v){return *(double*)&v;}
inline double& m128d_vn(__m128d &v, int component){return ((double*)&v)[component];}
inline double m128d_vn(__m128d const &v, int component){return ((double*)&v)[component];}
#else
inline int& m128i_i32(__m128i &v){return v.m128i_i32[0];}
inline int m128i_i32(__m128i const &v){return v.m128i_i32[0];}
inline int& m128i_i32(__m128i &v, int component){return v.m128i_i32[component];}
inline int m128i_i32(__m128i const &v, int component){return v.m128i_i32[component];}

inline float& m128_f32(__m128 &v){return v.m128_f32[0];}
inline float m128_f32(__m128 const &v){return v.m128_f32[0];}
inline float& m128_f32(__m128 &v, int component){return v.m128_f32[component];}
inline float m128_f32(__m128 const &v, int component){return v.m128_f32[component];}

inline double& m128d_v0(__m128d &v){return v.m128d_f64[0];}
inline double m128d_v0(__m128d const &v){return v.m128d_f64[0];}
inline double& m128d_vn(__m128d &v, int component){return v.m128d_f64[component];}
inline double m128d_vn(__m128d const &v, int component){return v.m128d_f64[component];}
#endif
#endif