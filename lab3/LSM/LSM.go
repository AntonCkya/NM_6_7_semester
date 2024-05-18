package LSM

import (
	"math"

	"github.com/AntonCkya/numeric/LU"
)

func MatrixTranspose(m [][]float64) [][]float64 {
	var c [][]float64
	for i := 0; i < len(m[0]); i++ {
		c = append(c, make([]float64, len(m)))
	}

	for i := 0; i < len(c); i++ {
		for j := 0; j < len(c[0]); j++ {
			c[i][j] = m[j][i]
		}
	}
	return c
}

func MatrixMult(a [][]float64, b [][]float64) [][]float64 {
	var c [][]float64
	for i := 0; i < len(a); i++ {
		c = append(c, make([]float64, len(b[0])))
	}

	for i := 0; i < len(a); i++ {
		for j := 0; j < len(b[0]); j++ {
			for k := 0; k < len(b); k++ {
				c[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return c
}

func MatrixVectorMult(a [][]float64, b []float64) []float64 {
	c := make([]float64, len(a))
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a[0]); j++ {
			c[i] += b[j] * a[i][j]
		}
	}
	return c
}

func A(x, y []float64, p, n int) []float64 {
	var a [][]float64
	for i := 0; i < n; i++ {
		a = append(a, make([]float64, p+1))
	}
	for i := 0; i < n; i++ {
		for j := 0; j < p+1; j++ {
			a[i][j] = math.Pow(x[i], float64(j))
		}
	}
	at := MatrixTranspose(a)
	r := MatrixVectorMult(at, y)
	ata := MatrixMult(at, a)
	return LU.SolveLU(ata, r, len(r))
}

func Error(x, y, a []float64) float64 {
	POLY := func(x float64, p []float64) float64 {
		res := 0.0
		for i := 0; i < len(p); i++ {
			res += math.Pow(x, float64(i)) * p[i]
		}
		return res
	}
	delta := 0.0
	for i := 0; i < len(x); i++ {
		delta += math.Pow(POLY(x[i], a)-y[i], 2.0)
	}
	return delta
}
