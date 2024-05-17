package interpolation

import "math"

func PolyMult(x1, x2 []float64) []float64 {
	res := make([]float64, len(x1)+len(x2))
	for i := 0; i < len(x1); i++ {
		for j := 0; j < len(x2); j++ {
			res[i+j] += x1[i] * x2[j]
		}
	}
	return res
}

func PolyPlus(x1, x2 []float64, n int) []float64 {
	res := make([]float64, n)
	for i := 0; i < n; i++ {
		if i >= len(x2) {
			res[i] = x1[i]
		} else if i >= len(x1) {
			res[i] = x2[i]
		} else {
			res[i] = x1[i] + x2[i]
		}
	}
	return res
}

func PolyScalarMult(x1 []float64, k float64, n int) []float64 {
	res := make([]float64, n)
	for i := 0; i < n; i++ {
		res[i] = x1[i] * k
	}
	return res
}

func PolynomLagrange(x, y []float64, n int) []float64 {
	var res []float64
	res = append(res, 0)
	for i := 0; i < n; i++ {
		li := make([]float64, 1)
		li[0] = 1
		for j := 0; j < n; j++ {
			if i == j {
				continue
			}
			xx := make([]float64, 2)
			xx[0] = -x[j]
			xx[1] = 1
			li = PolyMult(li, xx)
			li = PolyScalarMult(li, 1/(x[i]-x[j]), len(li))
		}
		res = PolyPlus(res, PolyScalarMult(li, y[i], len(li)), int(math.Max(float64(len(li)), float64(len(res)))))
	}
	return res
}

func diff(x, y []float64, l, r int) float64 {
	if l+1 == r {
		return (y[l] - y[r]) / (x[l] - x[r])
	} else {
		return (diff(x, y, l, r-1) - diff(x, y, l+1, r)) / (x[l] - x[r])
	}
}

func PolynomNewton(x, y []float64, n int) []float64 {
	j := 1
	var res []float64
	res = append(res, y[0])
	var li []float64
	li = append(li, -x[0], 1)
	for i := 1; i < n; i++ {
		res = PolyPlus(res, PolyScalarMult(li, diff(x, y, 0, j), len(li)), int(math.Max(float64(len(li)), float64(len(res)))))
		xx := make([]float64, 2)
		xx[0] = -x[i]
		xx[1] = 1
		li = PolyMult(li, xx)
		j++
	}
	return res
}
