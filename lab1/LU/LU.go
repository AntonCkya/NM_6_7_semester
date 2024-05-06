package LU

import (
	"math"

	iterations "github.com/AntonCkya/numeric_methods/Iterations"
)

func GetLU(m [][]float64, mSize int) ([][]float64, [][]float64, [][]float64) {
	var L [][]float64
	var U [][]float64
	var P [][]float64

	for i := 0; i < mSize; i++ {
		L = append(L, make([]float64, mSize))
		U = append(U, make([]float64, mSize))
		P = append(P, make([]float64, mSize))
	}

	for i := 0; i < mSize; i++ {
		L[i][i] = 1
		P[i][i] = 1
	}

	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			U[i][j] = m[i][j]
		}
	}

	for i := 0; i < mSize; i++ {
		kMax := i
		maxx := math.Inf(-1)
		for k := i; k < mSize; k++ {
			if U[i][k] > maxx {
				maxx = U[i][k]
				kMax = k
			}
		}
		for k := 0; k < mSize; k++ {
			temp := U[i][k]
			U[i][k] = U[kMax][k]
			U[kMax][k] = temp
		}
		for k := 0; k < mSize; k++ {
			temp := L[i][k]
			L[i][k] = L[kMax][k]
			L[kMax][k] = temp
		}
		for k := 0; k < mSize; k++ {
			temp := L[k][i]
			L[k][i] = L[k][kMax]
			L[k][kMax] = temp
		}
		for k := 0; k < mSize; k++ {
			temp := P[k][i]
			P[k][i] = P[k][kMax]
			P[k][kMax] = temp
		}
		for j := i + 1; j < mSize; j++ {
			l := U[j][i] / U[i][i]
			L[j][i] = l
			for k := 0; k < mSize; k++ {
				U[j][k] = U[j][k] - (U[i][k] * L[j][i])
			}
		}
	}

	return L, U, P
}

func DetLU(m [][]float64, mSize int) float64 {
	_, U, P := GetLU(m, mSize)
	sign := 1
	for i := 0; i < mSize; i++ {
		if P[i][i] != 1 {
			sign *= -1
		}
	}
	res := 1.0
	for i := 0; i < mSize; i++ {
		res *= U[i][i]
	}
	return res * float64(sign)
}

func SolveLU(m [][]float64, b []float64, mSize int) []float64 {
	L, U, P := GetLU(m, mSize)
	return SolveWithLU(L, U, P, b, mSize)
}

func SolveWithLU(L [][]float64, U [][]float64, P [][]float64, b []float64, mSize int) []float64 {
	y := make([]float64, mSize)
	b = iterations.MatrixVectorMult(P, b)
	for i := 0; i < mSize; i++ {
		temp := 0.0
		for j := 0; j < i; j++ {
			temp += L[i][j] * y[j]
		}
		y[i] = (b[i] - temp) / L[i][i]
	}
	x := make([]float64, mSize)
	for i := mSize - 1; i >= 0; i-- {
		temp := 0.0
		for j := mSize - 1; j >= i; j-- {
			temp += U[i][j] * x[j]
		}
		x[i] = (y[i] - temp) / U[i][i]
	}
	return x
}

func InvertLU(m [][]float64, mSize int) [][]float64 {
	e := make([]float64, mSize)
	e[0] = 1
	var res [][]float64
	for i := 0; i < mSize; i++ {
		res = append(res, make([]float64, mSize))
	}
	L, U, P := GetLU(m, mSize)
	for i := 1; i <= mSize; i++ {
		st := SolveWithLU(L, U, P, e, mSize)
		for j := 0; j < mSize; j++ {
			res[j][i-1] = st[j]
		}
		e[i-1] = 0
		if i != mSize {
			e[i] = 1
		}
	}

	return res
}
