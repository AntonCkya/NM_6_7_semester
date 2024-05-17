package main

import (
	"fmt"
	"math"

	derivative "github.com/AntonCkya/numeric/Derivative"
	integral "github.com/AntonCkya/numeric/Integral"
	interpolation "github.com/AntonCkya/numeric/Interpolation"
	plotter "github.com/AntonCkya/numeric/Plotter"
)

func InterpolationRunner() {
	F := func(x float64) float64 {
		return math.Pow(math.E, x) + x
	}

	POLY := func(x float64, p []float64) float64 {
		res := 0.0
		for i := 0; i < len(p); i++ {
			res += math.Pow(x, float64(i)) * p[i]
		}
		return res
	}

	fmt.Println("Size:")
	var n int
	fmt.Scan(&n)

	var x []float64
	var y []float64
	fmt.Println("X:")
	for i := 0; i < n; i++ {
		var xx float64
		fmt.Scan(&xx)
		x = append(x, xx)
		y = append(y, F(xx))
	}

	fmt.Println("X*:")
	var xx float64
	fmt.Scan(&xx)

	lagrange := interpolation.PolynomLagrange(x, y, n)
	fmt.Println("\nLagrange:")
	for i := 0; i < n; i++ {
		if i == 0 {
			fmt.Print(lagrange[i])
		} else {
			fmt.Print(lagrange[i], "*x^", i)
		}
		if i != n-1 {
			fmt.Print(" + ")
		}
	}
	fmt.Println()
	fmt.Println("f(x*) = ", POLY(xx, lagrange))
	fmt.Println("error = ", math.Abs(POLY(xx, lagrange)-F(xx)))

	var lx1, lx2, ly1, ly2 []float64
	for i := x[0] * 2; i < x[n-1]*2; i += 0.01 {
		lx1 = append(lx1, i)
		lx2 = append(lx2, i)
		ly1 = append(ly1, F(i))
		ly2 = append(ly2, POLY(i, lagrange))
	}
	plotter.Plot2(lx1, ly1, lx2, ly2, "lagrange", []string{"function", "lagrange"})

	newton := interpolation.PolynomNewton(x, y, n)
	fmt.Println("\nNewton:")
	for i := 0; i < n; i++ {
		if i == 0 {
			fmt.Print(newton[i])
		} else {
			fmt.Print(newton[i], "*x^", i)
		}
		if i != n-1 {
			fmt.Print(" + ")
		}
	}
	fmt.Println()
	fmt.Println("f(x*) = ", POLY(xx, newton))
	fmt.Println("error = ", math.Abs(POLY(xx, newton)-F(xx)))

	var nx1, nx2, ny1, ny2 []float64
	for i := x[0] * 2; i < x[n-1]*2; i += 0.01 {
		nx1 = append(nx1, i)
		nx2 = append(nx2, i)
		ny1 = append(ny1, F(i))
		ny2 = append(ny2, POLY(i, lagrange))
	}
	plotter.Plot2(nx1, ny1, nx2, ny2, "newton", []string{"function", "newton"})
}

func SplineRunner() {

}

func LSMRunner() {

}

func DerivativeRunner() {
	fmt.Println("3.4-----Derivative-----")

	fmt.Println("Size:")
	var tableSize int
	fmt.Scan(&tableSize)

	var table [][]float64
	table = append(table, make([]float64, tableSize))
	table = append(table, make([]float64, tableSize))

	fmt.Println("X:")
	for i := 0; i < tableSize; i++ {
		var x float64
		fmt.Scan(&x)
		table[0][i] = x
	}

	fmt.Println("Y:")
	for i := 0; i < tableSize; i++ {
		var y float64
		fmt.Scan(&y)
		table[1][i] = y
	}

	fmt.Println("X*:")
	var xx float64
	fmt.Scan(&xx)

	xxIndex := -1
	for i := 0; i < tableSize-1; i++ {
		if table[0][i] <= xx && table[0][i+1] >= xx {
			xxIndex = i
			break
		}
	}

	dx := derivative.DX(table, xx, xxIndex)
	d2x := derivative.D2X(table, xx, xxIndex)

	fmt.Println("D1: ", dx)
	fmt.Println("D2: ", d2x)
}

func IntegralRunner() {
	fmt.Println("3.5-----Integral-----")

	F := func(x float64) float64 {
		return 1 / (256 - math.Pow(x, 4))
	}
	a := -2.0
	b := 2.0
	h1 := 1.0
	h2 := 0.5

	rect1 := integral.IntegralRect(F, a, b, h1)
	rect2 := integral.IntegralRect(F, a, b, h2)

	trap1 := integral.IntegralTrapezoid(F, a, b, h1)
	trap2 := integral.IntegralTrapezoid(F, a, b, h2)

	simp1 := integral.IntegralSimpson(F, a, b, h1)
	simp2 := integral.IntegralSimpson(F, a, b, h2)

	fmt.Println("Rect method (h1): ", rect1)
	fmt.Println("Rect method (h2): ", rect2)
	fmt.Println()
	fmt.Println("Trapezoid method (h1): ", trap1)
	fmt.Println("Trapezoid method (h2): ", trap2)
	fmt.Println()
	fmt.Println("Simpson method (h1): ", simp1)
	fmt.Println("Simpson method (h2): ", simp2)
	fmt.Println()
	k := h2 / h1
	p := 2.0
	fmt.Println("Rect err: ", integral.RRRMethod(rect1, rect2, k, p))
	fmt.Println("Trapezoid err: ", integral.RRRMethod(trap1, trap2, k, p))
	fmt.Println("Simpson err: ", integral.RRRMethod(simp1, simp2, k, p))
}

func main() {
	fmt.Println("Select lab 3.X:")
	fmt.Println("1: interpolation")
	fmt.Println("2: spline")
	fmt.Println("3: LSM")
	fmt.Println("4: derivative")
	fmt.Println("5: integral")

	var point int
	fmt.Scan(&point)
	switch point {
	case 1:
		InterpolationRunner()
	case 2:
		SplineRunner()
	case 3:
		LSMRunner()
	case 4:
		DerivativeRunner()
	case 5:
		IntegralRunner()
	default:
		fmt.Println("I don't know u wrong")
	}
}
