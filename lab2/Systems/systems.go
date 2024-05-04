package systems

import (
	"math"

	"github.com/AntonCkya/numeric_methods/LU"
)

func GetQ(dphi [][]func([]float64) float64, a1, b1, a2, b2 float64) float64 {
	x1 := (a1+b1)/2 + math.Abs(b1-a1)
	x2 := (a2+b2)/2 + math.Abs(b2-a2)
	var x []float64
	x = append(x, x1, x2)
	return math.Max(math.Abs(dphi[0][0](x))+math.Abs(dphi[0][1](x)), math.Abs(dphi[1][0](x))+math.Abs(dphi[1][1](x)))
}

func VectorNormal(v []float64) float64 {
	res := math.Inf(-1)
	for i := 0; i < len(v); i++ {
		res = math.Max(res, math.Abs(v[i]))
	}
	return res
}

func VectorSum(a []float64, b []float64) []float64 {
	vSize := len(b)
	c := make([]float64, vSize)
	for i := 0; i < vSize; i++ {
		c[i] = a[i] + b[i]
	}
	return c
}

func VectorSubstr(a []float64, b []float64) []float64 {
	vSize := len(b)
	for i := 0; i < vSize; i++ {
		b[i] = -b[i]
	}
	return VectorSum(a, b)
}

func SimpleIterationsMethod(phi []func([]float64) float64, dphi [][]func([]float64) float64, a1, b1, a2, b2, eps float64) ([]float64, int) {
	q := GetQ(dphi, a1, b1, a2, b2)
	x := append(make([]float64, 0), a1, b1)
	dx := math.Inf(+1)
	count := 0
	for eps < dx*(q/(1-q)) {
		xNew := append(make([]float64, 0), phi[0](x), phi[1](x))
		dx = VectorNormal(VectorSubstr(xNew, x))
		x = xNew
		count++
	}
	return x, count
}

func GetMatrixJ(df [][]func([]float64) float64, x []float64) [][]float64 {
	var J [][]float64
	for i := 0; i < 2; i++ {
		J = append(J, append(make([]float64, 0), df[i][0](x), df[i][1](x)))
	}
	return J
}

func GetVectorB(f []func([]float64) float64, x []float64) []float64 {
	var b []float64
	for i := 0; i < 2; i++ {
		b = append(b, f[i](x))
	}
	return b
}

func NewtonMethod(f []func([]float64) float64, df [][]func([]float64) float64, a1, b1, eps float64) ([]float64, int) {
	dx := math.Inf(+1)
	x := append(make([]float64, 0), a1, b1)
	count := 0
	for eps < dx {
		J := GetMatrixJ(df, x)
		b := GetVectorB(f, x)
		delta := LU.SolveLU(J, b, len(b))
		xNew := VectorSubstr(x, delta)
		dx = VectorNormal(VectorSubstr(xNew, x))
		x = xNew
		count++
	}
	return x, count
}
