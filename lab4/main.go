package main

import (
	"fmt"
	"math"

	cauchy "github.com/AntonCkya/numeric_methods/Cauchy"
	plotter "github.com/AntonCkya/numeric_methods/Plotter"
)

func CauchyRunner() {
	SYS := []func(x, y, z float64) float64{
		func(x, y, z float64) float64 {
			return z
		},
		func(x, y, z float64) float64 {
			return ((x+1)/x)*z - (y / x)
		},
	}
	SOLVE := func(x float64) float64 {
		return x + 1 + math.Pow(math.E, x)
	}
	yKey := append(make([]float64, 0, 2), 1.0, 2.0+math.E)
	zKey := append(make([]float64, 0, 2), 1.0, 1.0+math.E)
	xl := 1.0
	xr := 2.0
	h := 0.1

	xEuler, yEuler := cauchy.DiffEuler(SYS, yKey, zKey, xl, xr, h)
	_, y2Euler := cauchy.DiffEuler(SYS, yKey, zKey, xl, xr, h/2)
	var yTrue []float64
	for i := 0; i < len(xEuler); i++ {
		yTrue = append(yTrue, SOLVE(xEuler[i]))
	}

	plotter.Plot2(xEuler, yEuler, xEuler, yTrue, "euler", []string{"euler", "function"})

	xRunge, yRunge, _ := cauchy.DiffRungeKutt(SYS, yKey, zKey, xl, xr, h)
	_, y2Runge, _ := cauchy.DiffRungeKutt(SYS, yKey, zKey, xl, xr, h/2)
	plotter.Plot2(xRunge, yRunge, xRunge, yTrue, "runge-kutt", []string{"runge-kutt", "function"})

	xAdams, yAdams := cauchy.DiffAdams(SYS, yKey, zKey, xl, xr, h)
	_, y2Adams := cauchy.DiffAdams(SYS, yKey, zKey, xl, xr, h/2)
	plotter.Plot2(xAdams, yAdams, xAdams, yTrue, "adams", []string{"adams", "function"})

	fmt.Println("Euler error (RRR): ", cauchy.RRRmethod(yEuler, y2Euler, 1))
	fmt.Println("Runge-Kutt error (RRR): ", cauchy.RRRmethod(yRunge, y2Runge, 4))
	fmt.Println("Adams error(RRR): ", cauchy.RRRmethod(yAdams, y2Adams, 4))
	fmt.Println()
	fmt.Println("Euler error: ", cauchy.AbsoluteError(yEuler, yTrue))
	fmt.Println("Runge-Kutt error: ", cauchy.AbsoluteError(yRunge, yTrue))
	fmt.Println("Adams error: ", cauchy.AbsoluteError(yAdams, yTrue))
}

func BoundaryRunner() {

}

func main() {
	fmt.Println("Select lab 4.X:")
	fmt.Println("1: Cauchy")
	fmt.Println("2: boundary")

	var point int
	fmt.Scan(&point)
	switch point {
	case 1:
		CauchyRunner()
	case 2:
		BoundaryRunner()
	default:
		fmt.Println("I don't know u wrong")
	}
}
