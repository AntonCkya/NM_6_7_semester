package rotations

import (
	"math"
)

func MatrixSquares(m [][]float64) float64 {
	res := 0.0
	for i := 0; i < len(m); i++ {
		for j := i + 1; j < len(m); j++ {
			res += m[i][j] * m[i][j]
		}
	}
	return math.Sqrt(res)
}

func GetMaxNoDiag(m [][]float64, mSize int) (int, int) {
	iMax := 0
	jMax := 0
	maxElem := math.Inf(-1)
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			if i != j {
				if math.Abs(m[i][j]) > maxElem {
					iMax = i
					jMax = j
					maxElem = math.Abs(m[i][j])
				}
			}
		}
	}
	return iMax, jMax
}

func GetRotationMatrix(mSize int, phi float64, i int, j int) [][]float64 {
	var u [][]float64
	for k := 0; k < mSize; k++ {
		u = append(u, make([]float64, mSize))
	}
	for k := 0; k < mSize; k++ {
		u[k][k] = 1
	}
	u[i][i] = math.Cos(phi)
	u[j][j] = math.Cos(phi)
	u[i][j] = -math.Sin(phi)
	u[j][i] = math.Sin(phi)
	return u
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

func Rotation(m [][]float64, mSize int, eps float64) ([][]float64, [][]float64, int) {
	var mNew [][]float64
	for i := 0; i < mSize; i++ {
		mNew = append(mNew, m[i])
	}
	var u [][]float64
	for i := 0; i < mSize; i++ {
		u = append(u, make([]float64, mSize))
	}
	for i := 0; i < mSize; i++ {
		u[i][i] = 1
	}
	count := 0
	var phi float64
	for MatrixSquares(mNew) > eps {
		i, j := GetMaxNoDiag(mNew, mSize)
		if mNew[i][i] == mNew[j][j] {
			phi = math.Pi / 4
		} else {
			phi = 0.5 * math.Atan2(2*mNew[i][j], (mNew[i][i]-mNew[j][j]))
		}
		uNew := GetRotationMatrix(mSize, phi, i, j)
		u = MatrixMult(u, uNew, mSize)
		mNew = MatrixMult(MatrixMult(MatrixTranspose(uNew, mSize), mNew, mSize), uNew, mSize)
		count++
	}
	return mNew, u, count
}
