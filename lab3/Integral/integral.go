package integral

import "math"

func IntegralRect(f func(float64) float64, a, b, h float64) float64 {
	x1 := a
	x2 := a + h
	res := 0.0
	for x1 < b {
		res += f((x1 + x2) / 2)
		x1 = x2
		x2 += h
	}
	return res * h
}

func IntegralTrapezoid(f func(float64) float64, a, b, h float64) float64 {
	x1 := a
	x2 := a + h
	res := 0.0
	for x1 < b {
		res += (f(x1) + f(x2))
		x1 = x2
		x2 += h
	}
	return res * (h / 2)
}

func IntegralSimpson(f func(float64) float64, a, b, h float64) float64 {
	x := a + h
	res := 0.0
	for x < b {
		res += f(x-h) + 4*f(x) + f(x+h)
		x += h * 2
	}
	return res * (h / 3)
}

func RRRMethod(Fh, Fkh, k, p float64) float64 {
	return (Fh - Fkh) / (math.Pow(k, p) - 1.0)
}
