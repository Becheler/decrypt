Notes by Ziheng Yang
18 March 2017

(1) To compile, try one of the following

   Linux/UNIX gcc compiler:
      gcc -o bpp -O3 bpp.c tools.c -lm
      gcc -o MCcoal -DSIMULATION bpp.c tools.c -lm
      gcc -o bpp_sse -O3 -DUSE_SSE -msse3 bpp.c tools.c -lm
      gcc -o bpp_avx -O3 -DUSE_AVX -mavx bpp.c tools.c -lm

      gcc -o bpp -O3 bpp.c tools.c -lm
      gcc -o MCcoal -DSIMULATION bpp.c tools.c -lm

   INTEL icc compiler:
      icc -o bpp -fast bpp.c tools.c -lm
      icc -o MCcoal -DSIMULATION -fast bpp.c tools.c -lm

   MAC OSX intel:
      cc -o bpp -O3 bpp.c tools.c -lm
      cc -o MCcoal -DSIMULATION -O3 bpp.c tools.c -lm

   Windows Microsoft Visual C++:
      cl -Febpp_NoSSE.exe -O2 -W3 -D_CRT_SECURE_NO_WARNINGS bpp.c tools.c /F10000000
      cl -O2 -FeMCcoal.exe -DSIMULATION bpp.c tools.c /F10000000
      cl -Febpp_sse.exe -DUSE_SSE -arch:SSE2 -Ox -W3 -D_CRT_SECURE_NO_WARNINGS bpp.c tools.c /F10000000
      cl -Febpp.exe -DUSE_AVX -arch:AVX -Ox -W3 -D_CRT_SECURE_NO_WARNINGS bpp.c tools.c /F10000000

(2) To run an example analysis, try 

   cd examples

   ../bpp yu2001.bpp.ctl

   ../bpp ChenLi2001.bpp.ctl

   ../bpp lizard.bpp.ctl

(3) Copyright notice and disclaimer

The software package is provided "as is" without warranty of any
kind. In no event shall the author or his employer be held responsible
for any damage resulting from the use of this software, including but
not limited to the frustration that you may experience in using the
package.  The program package, including source codes, example data
sets, executables, and this documentation, is maintained by Ziheng
Yang and distributed under the GNU GPL v3.
