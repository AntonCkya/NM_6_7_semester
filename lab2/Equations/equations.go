package equations

import "math"

func SimpleIterationsMethod(phi, dphi func(float64) float64, a, b, eps float64) (float64, int) {
	q := math.Min(math.Abs(dphi(a)), math.Abs(dphi(b)))
	x := b
	dx := math.Inf(+1)
	count := 0
	for eps < dx*(q/(1-q)) {
		xNew := phi(x)
		dx = math.Abs(xNew - x)
		x = xNew
		count++
	}
	return x, count
}

func NewtonMethod(f, df, d2f func(float64) float64, eps, x0 float64) (float64, int) {
	if f(x0)*d2f(x0) <= 0 {
		return 0.0, -1
	}
	x := x0
	dx := math.Inf(+1)
	count := 0
	for eps < dx {
		xNew := x - f(x)/df(x)
		dx = math.Abs(xNew - x)
		x = xNew
		count++
	}
	return x, count
}
