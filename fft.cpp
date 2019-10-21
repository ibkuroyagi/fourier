#include <stdio.h>
#include <math.h>
const double PI = 3.141592653589793;

/*
高速フーリエ変換
int n ：データ数（２のべき乗）
int flg ：順変換:-1, 逆変換:1
double *ar ：データ配列（実部；要素数 n)
double *ai ：データ配列（虚部；要素数 n）
*/
int fft(int n, int flg, double *ar, double *ai)
{
  long m, mh, i, j, k, irev;
  double wr, wi, xr, xi;
  double theta;

  theta = flg * 2 * PI / n;

  i = 0;
  for (j = 1; j < n - 1; j++)
  {
    for (k = n >> 1; k > (i ^= k); k >>= 1)
      ;
    if (j < i)
    {
      xr = ar[j];
      xi = ai[j];
      ar[j] = ar[i];
      ai[j] = ai[i];
      ar[i] = xr;
      ai[i] = xi;
    }
  }
  for (mh = 1; (m = mh << 1) <= n; mh = m)
  {
    irev = 0;
    for (i = 0; i < n; i += m)
    {
      wr = cos(theta * irev);
      wi = sin(theta * irev);
      for (k = n >> 2; k > (irev ^= k); k >>= 1)
        ;
      for (j = i; j < mh + i; j++)
      {
        k = j + mh;
        xr = ar[j] - ar[k];
        xi = ai[j] - ai[k];
        ar[j] += ar[k];
        ai[j] += ai[k];
        ar[k] = wr * xr - wi * xi;
        ai[k] = wr * xi + wi * xr;
      }
    }
  }

  if (flg == -1)
  {
    for (i = 0; i<n; i++)
    {
      ar[i] /= n;
      ai[i] /= n;
    }
  }
  return 0;
}
