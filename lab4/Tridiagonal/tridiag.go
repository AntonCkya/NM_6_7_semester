package tridiagonal

func MatrixToTridiagonal(m [][]float64, mSize int) ([]float64, []float64, []float64) {
	A := make([]float64, mSize)
	B := make([]float64, mSize)
	C := make([]float64, mSize)
	for i := 0; i < mSize; i++ {
		if i > 0 {
			A[i] = m[i][i-1]
		}
		B[i] = m[i][i]
		if i < mSize-1 {
			C[i] = m[i][i+1]
		}
	}
	return A, B, C
}

func SolveTridiagonal(a, b, c, d []float64, mSize int) []float64 {
	p := make([]float64, mSize)
	q := make([]float64, mSize)

	p[0] = -c[0] / b[0]
	q[0] = d[0] / b[0]

	for i := 1; i < mSize; i++ {
		p[i] = -c[i] / (b[i] + a[i]*p[i-1])
		q[i] = (d[i] - a[i]*q[i-1]) / (b[i] + a[i]*p[i-1])
	}

	x := make([]float64, mSize)
	x[mSize-1] = q[mSize-1]

	for i := mSize - 2; i >= 0; i-- {
		x[i] = p[i]*x[i+1] + q[i]
	}

	return x
}
