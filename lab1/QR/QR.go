package QR

import (
	"math"
)

func VectorVectorMatrix(a []float64, b []float64) [][]float64 {
	var c [][]float64
	for i := 0; i < len(a); i++ {
		c = append(c, make([]float64, len(a)))
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a); j++ {
			c[i][j] = a[i] * b[j]
		}
	}
	return c
}

func MatrixMult(a [][]float64, b [][]float64, mSize int) [][]float64 {
	var c [][]float64
	for i := 0; i < mSize; i++ {
		c = append(c, make([]float64, mSize))
	}

	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			for k := 0; k < mSize; k++ {
				c[i][j] += a[i][k] * b[k][j]
			}
		}
	}
	return c
}

func MatrixTranspose(m [][]float64, mSize int) [][]float64 {
	var c [][]float64
	for i := 0; i < mSize; i++ {
		c = append(c, make([]float64, mSize))
	}

	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			c[i][j] = m[j][i]
		}
	}
	return c
}

func VectorVectorNumber(a []float64, b []float64) float64 {
	c := 0.0
	for i := 0; i < len(a); i++ {
		c += a[i] * b[i]
	}
	return c
}

func MatrixNumberMult(a [][]float64, n float64) [][]float64 {
	var c [][]float64
	for i := 0; i < len(a); i++ {
		c = append(c, make([]float64, len(a)))
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a); j++ {
			c[i][j] = a[i][j] * n
		}
	}
	return c
}

func MatrixSum(a [][]float64, b [][]float64) [][]float64 {
	var c [][]float64
	for i := 0; i < len(a); i++ {
		c = append(c, make([]float64, len(a)))
	}
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a); j++ {
			c[i][j] = a[i][j] + b[i][j]
		}
	}
	return c
}

func MatrixSubstr(a [][]float64, b [][]float64) [][]float64 {
	for i := 0; i < len(a); i++ {
		for j := 0; j < len(a); j++ {
			b[i][j] = -b[i][j]
		}
	}
	return MatrixSum(a, b)
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

func VectorSquares(v []float64) float64 {
	res := 0.0
	for i := 0; i < len(v); i++ {
		res += v[i] * v[i]
	}
	return math.Sqrt(res)
}

func Sign(x float64) int {
	switch {
	case x > 0:
		return 1
	case x == 0:
		return 0
	default:
		return -1
	}
}

func GetHouseholderMatrix(mSize int, v []float64, i int) [][]float64 {
	vNew := VectorSum(v, make([]float64, mSize))
	vNew[i] += float64(Sign(v[i])) * VectorSquares(v)
	var e [][]float64
	for k := 0; k < mSize; k++ {
		e = append(e, make([]float64, mSize))
	}
	for k := 0; k < mSize; k++ {
		e[k][k] = 1
	}
	return MatrixSubstr(e, MatrixNumberMult(VectorVectorMatrix(vNew, vNew), 2/(VectorVectorNumber(vNew, vNew))))
}

func GetQR(m [][]float64, mSize int) ([][]float64, [][]float64) {
	var Q [][]float64
	for i := 0; i < mSize; i++ {
		Q = append(Q, make([]float64, mSize))
	}
	for i := 0; i < mSize; i++ {
		Q[i][i] = 1.0
	}
	var R [][]float64
	for i := 0; i < mSize; i++ {
		R = append(R, m[i])
	}
	for i := 0; i < mSize-1; i++ {
		b := make([]float64, mSize)
		for j := i; j < mSize; j++ {
			b[j] = R[j][i]
		}
		H := GetHouseholderMatrix(mSize, b, i)
		Q = MatrixMult(Q, H, mSize)
		R = MatrixMult(H, R, mSize)
	}
	return Q, R
}

func ComplexSolve(a11, a12, a21, a22, eps float64) (complex128, complex128, bool) {
	a := 1.0
	b := -a11 - a22
	c := a11*a22 - a12*a21
	d := b*b - 4*a*c
	if d > eps {
		return complex(0, 0), complex(0, 0), false
	} else {
		return (complex(-b, 0) + complex(0, math.Sqrt(-d))) / complex(2*a, 0),
			(complex(-b, 0) - complex(0, math.Sqrt(-d))) / complex(2*a, 0), true
	}

}

func QRAlgo(m [][]float64, mSize int, eps float64) ([]complex128, int) {
	var mNew [][]float64
	for i := 0; i < mSize; i++ {
		mNew = append(mNew, m[i])
	}
	count := 0
	res := make([]complex128, mSize)
	for true {
		count++
		Q, R := GetQR(mNew, mSize)
		mNew = MatrixMult(R, Q, mSize)
		br := true
		i := 0
		for i < mSize {
			if i < mSize-1 && math.Abs(mNew[i+1][i]) > eps {
				e1, e2, f := ComplexSolve(mNew[i][i], mNew[i][i+1], mNew[i+1][i], mNew[i+1][i+1], eps)
				if f {
					res[i] = e1
					res[i+1] = e2
					i++
				} else {
					br = false
				}
			} else {
				res[i] = complex(mNew[i][i], 0)
			}
			i++
		}
		if br {
			break
		}
	}
	return res, count
}
