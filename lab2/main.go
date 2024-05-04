package main

import (
	"fmt"
	"math"

	equations "github.com/AntonCkya/numeric_methods/Equations"
	systems "github.com/AntonCkya/numeric_methods/Systems"
)

func EquationRunner() {
	PHI := func(x float64) float64 {
		return 0.5 * math.Log2(5*x+2)
		// return math.Log(4-3*x) / 2
	}
	DPHI := func(x float64) float64 {
		return 5 / (2 * math.Log(2) * (5*x + 2))
		// return -3 / (2 * (4 - 3*x))
	}
	F := func(x float64) float64 {
		return math.Pow(4, x) - 5*x - 2
		// return math.Pow(math.E, 2*x) + 3*x - 4
	}
	DF := func(x float64) float64 {
		return math.Pow(4, x)*math.Log(4) - 5
		// return 2*math.Pow(math.E, 2*x) + 3
	}
	D2F := func(x float64) float64 {
		return math.Pow(4, x) * math.Log(4) * math.Log(4)
		// return 4 * math.Pow(math.E, 2*x)
	}

	fmt.Println("2.1-----Equations-----")
	fmt.Println("Simple Iterations:")

	res, count := equations.SimpleIterationsMethod(PHI, DPHI, 1.0, 2.0, 0.001)
	fmt.Println("answer: ", res)
	fmt.Println("count: ", count)

	fmt.Println("\nNewton:")

	res2, count2 := equations.NewtonMethod(F, DF, D2F, 0.001, 2)
	fmt.Println("answer: ", res2)
	fmt.Println("count: ", count2)
}

func SystemRunner() {
	fmt.Println("2.2-----Systems-----")

	PHI := []func([]float64) float64{
		func(x []float64) float64 {
			return math.Cos(x[1]) / 3
			//return 0.3 - 0.1*x[0]*x[0] - 0.2*x[1]*x[1]
		},
		func(x []float64) float64 {
			return math.Pow(math.E, x[0]) / 3
			//return 0.7 - 0.2*x[0]*x[0] + 0.1*x[0]*x[1]
		},
	}
	DPHI := [][]func([]float64) float64{
		{
			func(x []float64) float64 {
				return 0
				//return -0.2 * x[0]
			},
			func(x []float64) float64 {
				return -math.Sin(x[1]) / 3
				//return -0.4 * x[1]
			},
		},
		{
			func(x []float64) float64 {
				return math.Pow(math.E, x[0]) / 3
				//return -0.4*x[0] + 0.1*x[1]
			},
			func(x []float64) float64 {
				return 0
				//return 0.1 * x[0]
			},
		},
	}
	F := []func([]float64) float64{
		func(x []float64) float64 {
			return 3*x[0] - math.Cos(x[1])
			//return 0.1*x[0]*x[0] + x[0] + 0.2*x[1]*x[1] - 0.3
		},
		func(x []float64) float64 {
			return 3*x[1] - math.Pow(math.E, x[0])
			//return 0.2*x[0]*x[0] + x[1] - 0.1*x[0]*x[1] - 0.7
		},
	}
	DF := [][]func([]float64) float64{
		{
			func(x []float64) float64 {
				return 3
				//return 0.2*x[0] + 1
			},
			func(x []float64) float64 {
				return math.Sin(x[1])
				//return 0.4 * x[1]
			},
		},
		{
			func(x []float64) float64 {
				return -math.Pow(math.E, x[0])
				//return 0.4*x[0] - 0.1*x[1]
			},
			func(x []float64) float64 {
				return 3
				//return 1 - 0.1*x[0]
			},
		},
	}

	fmt.Println("Simple Iterations:")
	res, count := systems.SimpleIterationsMethod(PHI, DPHI, 0.0, 0.5, 0.5, 1, 0.0001)
	for i := 0; i < 2; i++ {
		fmt.Printf("x%d = %f\n", i+1, res[i])
	}
	fmt.Println("count: ", count)

	fmt.Println("\nNewton:")
	res2, count2 := systems.NewtonMethod(F, DF, 0.0, 0.0, 0.0001)
	for i := 0; i < 2; i++ {
		fmt.Printf("x%d = %f\n", i+1, res2[i])
	}
	fmt.Println("count: ", count2)

}

func main() {
	fmt.Println("Select lab 2.X:")
	fmt.Println("1: equation")
	fmt.Println("2: system")

	var point int
	fmt.Scan(&point)
	switch point {
	case 1:
		EquationRunner()
	case 2:
		SystemRunner()
	default:
		fmt.Println("I don't know u wrong")
	}
}
