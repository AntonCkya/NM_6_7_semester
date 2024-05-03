package iterations

import (
	"math"
)

func MatrixNormal(m [][]float64) float64 {
	res := math.Inf(-1)
	for i := 0; i < len(m); i++ {
		s := 0.0
		for j := 0; j < len(m); j++ {
			s += math.Abs(m[i][j])
		}
		res = math.Max(res, s)
	}
	return res
}

func VectorNormal(v []float64) float64 {
	res := math.Inf(-1)
	for i := 0; i < len(v); i++ {
		res = math.Max(res, math.Abs(v[i]))
	}
	return res
}

func Jakobi(a [][]float64, b []float64, mSize int) ([][]float64, []float64) {
	var alpha [][]float64
	for i := 0; i < mSize; i++ {
		alpha = append(alpha, make([]float64, mSize))
	}
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			if i != j {
				alpha[i][j] = -a[i][j] / a[i][i]
			}
		}
	}

	beta := make([]float64, mSize)
	for i := 0; i < mSize; i++ {
		beta[i] = b[i] / a[i][i]
	}
	return alpha, beta
}

func MatrixVectorMult(a [][]float64, b []float64) []float64 {
	//pls matching matrix & vector
	mSize := len(b)
	c := make([]float64, mSize)
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			c[i] += b[j] * a[i][j]
		}
	}
	return c
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

func SolveSimpleIt(m [][]float64, b []float64, eps float64, mSize int) ([]float64, int) {
	alpha, beta := Jakobi(m, b, mSize)
	currEps := math.Inf(1)
	aNorm := MatrixNormal(alpha)
	count := 0
	x := VectorSum(beta, make([]float64, mSize))
	for eps < currEps {
		c := MatrixVectorMult(alpha, x)
		newX := VectorSum(beta, c)
		currEps = (aNorm / (1 - aNorm)) * VectorNormal(VectorSubstr(newX, x))
		x = newX
		count++
	}
	return x, count
}

func SolveZeidel(m [][]float64, b []float64, eps float64, mSize int) ([]float64, int) {
	alpha, beta := Jakobi(m, b, mSize)

	var c [][]float64
	for i := 0; i < mSize; i++ {
		c = append(c, make([]float64, mSize))
	}
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			c[i][j] = alpha[i][j]
		}
	}

	currEps := math.Inf(1)
	aNorm := MatrixNormal(alpha)
	cNorm := MatrixNormal(c)
	count := 0
	x := VectorSum(beta, make([]float64, mSize))
	for eps < currEps {
		newX := VectorSum(beta, make([]float64, mSize))
		for i := 0; i < mSize; i++ {
			for j := 0; j < i; j++ {
				newX[i] += newX[j] * alpha[i][j]
			}
			for j := i; j < mSize; j++ {
				newX[i] += x[j] * alpha[i][j]
			}
		}
		currEps = (cNorm / (1 - aNorm)) * VectorNormal(VectorSubstr(newX, x))
		x = newX
		count++
	}
	return x, count
}
