package spline

import tridiagonal "github.com/AntonCkya/numeric/Tridiagonal"

func Spline(x, f []float64, n int) ([]float64, []float64, []float64, []float64) {
	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	d := make([]float64, n)

	h := make([]float64, n)
	for i := 1; i < n; i++ {
		h[i] = x[i] - x[i-1]
	}

	ta := make([]float64, n-2)
	tb := make([]float64, n-2)
	tc := make([]float64, n-2)
	td := make([]float64, n-2)
	for i := 0; i < n-2; i++ {
		ta[i] = h[i+1]
		tb[i] = 2 * (h[i+1] + h[i+2])
		tc[i] = h[i+2]
		td[i] = 3 * ((f[i+2]-f[i+1])/h[i+2] - (f[i+1]-f[i])/h[i+1])
	}
	ta[0] = 0
	tc[n-3] = 0
	solvedC := tridiagonal.SolveTridiagonal(ta, tb, tc, td, n-2)

	for i := 1; i < n; i++ {
		a[i] = f[i-1]
		if i > 1 {
			c[i] = solvedC[i-2]
		}
	}
	c[1] = 0.0
	for i := 1; i < n-1; i++ {
		b[i] = (f[i]-f[i-1])/h[i] - h[i]*(c[i+1]+2*c[i])/3.0
		d[i] = (c[i+1] - c[i]) / (3 * h[i])
	}
	b[n-1] = ((f[n-1] - f[n-2]) / h[n-1]) - ((2.0 / 3.0) * h[n-1] * c[n-1])
	d[n-1] = -c[n-1] / (3 * h[n-1])
	return a, b, c, d
}
