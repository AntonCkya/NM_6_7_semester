package boundary

import (
	"math"

	cauchy "github.com/AntonCkya/numeric_methods/Cauchy"
	tridiagonal "github.com/AntonCkya/numeric_methods/Tridiagonal"
)

func DiffShooting(system []func(x, y, z float64) float64, xl, xr, al, bl, yl, ar, br, yr, h, eps float64) ([]float64, []float64) {
	n0, n1, n2 := 1.0, 0.8, 2.0
	var x1, y1, z1 []float64
	for math.Abs(n2) > eps {
		_, y0, z0 := cauchy.DiffAdams(
			system, append(make([]float64, 0, 2), xl, n0),
			append(make([]float64, 0, 2), xl, (yl-(al*n0))/bl), xl, xr, h)
		x1, y1, z1 = cauchy.DiffAdams(
			system, append(make([]float64, 0, 2), xl, n1),
			append(make([]float64, 0, 2), xl, (yl-(al*n1))/bl), xl, xr, h)
		f0 := ar*y0[len(y0)-1] + br*z0[len(z0)-1] - yr
		f1 := ar*y1[len(y1)-1] + br*z1[len(z1)-1] - yr
		n2 = n1 - (n1-n0)/(f1-f0)*f1
		n0 = n1
		n1 = n2
	}
	return x1, y1
}

func DiffFinite(system []func(x, y, z float64) float64, xl, xr, al, bl, yl, ar, br, yr, h float64) ([]float64, []float64) {
	p := func(x float64) float64 {
		return -1 * system[1](x, 0, 1)
	}
	q := func(x float64) float64 {
		return -1 * system[1](x, 1, 0)
	}
	f := func(x float64) float64 {
		return system[0](x, 0, 0)
	}

	var a, b, c, d, x []float64
	x = append(x, xl)
	count := 0

	a = append(a, 0)
	b = append(b, al-bl/h)
	c = append(c, bl/h)
	d = append(d, yl)

	for x[count]+0.0000001 < xr-h {
		x = append(x, x[count]+h)
		count++
		a = append(a, 1-p(x[count])*h*0.5)
		b = append(b, -2+h*h*q(x[count]))
		c = append(c, 1+p(x[count])*h*0.5)
		d = append(d, h*h*f(x[count]))
	}

	x = append(x, xr)
	a = append(a, -br/h)
	b = append(b, ar+br/h)
	c = append(c, 0)
	d = append(d, yr)

	y := tridiagonal.SolveTridiagonal(a, b, c, d, len(b))
	return x, y
}
