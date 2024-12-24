package rotations

import (
	"fmt"
	"math"
	"sync"
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

func GetMaxNoDiagOneThread(m [][]float64, mSize int) (int, int) {
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

func GetMaxConcurrent(m []float64, mSize, i0, i1 int) int {
	//fmt.Println("max", i0, i1)
	if i0-i1 == 0 {
		return i0
	} else if i0-i1 == -1 {
		if math.Abs(m[i0]) >= math.Abs(m[i1]) {
			return i0
		} else {
			return i1
		}
	} else {
		i01 := (i0 + i1) / 2

		var m1, m0 int
		var wg sync.WaitGroup
		wg.Add(2)
		go func() {
			defer wg.Done()
			m0 = GetMaxConcurrent(m, mSize, i0, i01)
		}()
		go func() {
			defer wg.Done()
			m1 = GetMaxConcurrent(m, mSize, i01+1, i1)
		}()
		wg.Wait()
		if math.Abs(m[m0]) >= math.Abs(m[m1]) {
			return m0
		} else {
			return m1
		}
	}
}

func GetMaxNoDiag(m [][]float64, mSize int) (int, int) {
	mm := make([]float64, mSize*mSize)
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			if i == j {
				mm[i*mSize+j] = 0
			} else {
				mm[i*mSize+j] = m[i][j]
			}
		}
	}
	maxmm := GetMaxConcurrent(mm, mSize, 0, mSize*mSize-1)
	jj := maxmm % mSize
	ii := maxmm / mSize
	return ii, jj
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

	var wg sync.WaitGroup
	wg.Add(mSize * mSize)
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			go func() {
				//fmt.Println("MX Mult", i, j)
				defer wg.Done()
				for k := 0; k < mSize; k++ {
					// вынесение умножений в потоки не даст времени
					c[i][j] += a[i][k] * b[k][j]
				}
			}()
		}
	}
	wg.Wait()
	return c
}

func MatrixVectorMult(a [][]float64, b []float64) []float64 {
	//for error counting
	mSize := len(b)
	c := make([]float64, mSize)

	var wg sync.WaitGroup
	wg.Add(mSize)
	for i := 0; i < mSize; i++ {
		go func() {
			//fmt.Println("MV Mult", i)
			defer wg.Done()
			for j := 0; j < mSize; j++ {
				c[i] += b[j] * a[i][j]
			}
		}()
	}
	wg.Wait()
	return c
}

func MatrixTranspose(m [][]float64, mSize int) [][]float64 {
	var c [][]float64
	for i := 0; i < mSize; i++ {
		c = append(c, make([]float64, mSize))
	}

	var wg sync.WaitGroup
	wg.Add(mSize * mSize)
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			go func() {
				//fmt.Println("Transpose", i, j)
				defer wg.Done()
				c[i][j] = m[j][i]
			}()
		}
	}
	wg.Wait()
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

		xxx := MatrixMult(mNew, uNew, mSize)
		yyy := MatrixMult(MatrixTranspose(uNew, mSize), mNew, mSize)

		mNew = MatrixMult(MatrixMult(MatrixTranspose(uNew, mSize), mNew, mSize), uNew, mSize)

		fmt.Println("U:")
		for i := 0; i < mSize; i++ {
			for j := 0; j < mSize; j++ {
				fmt.Printf("%.4f ", uNew[j][i])
			}
			fmt.Println()
		}
		fmt.Println("Mnew:")
		for i := 0; i < mSize; i++ {
			for j := 0; j < mSize; j++ {
				fmt.Printf("%.4f ", mNew[j][i])
			}
			fmt.Println()
		}
		fmt.Println("MU:")
		for i := 0; i < mSize; i++ {
			for j := 0; j < mSize; j++ {
				fmt.Printf("%.4f ", xxx[j][i])
			}
			fmt.Println()
		}
		fmt.Println("UtM:")
		for i := 0; i < mSize; i++ {
			for j := 0; j < mSize; j++ {
				fmt.Printf("%.4f ", yyy[j][i])
			}
			fmt.Println()
		}

		count++
	}
	return mNew, u, count
}
