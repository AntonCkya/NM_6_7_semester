package main

import (
	"fmt"

	iterations "github.com/AntonCkya/numeric_methods/Iterations"
	"github.com/AntonCkya/numeric_methods/LU"
	"github.com/AntonCkya/numeric_methods/QR"
	rotations "github.com/AntonCkya/numeric_methods/Rotations"
	tridiagonal "github.com/AntonCkya/numeric_methods/Tridiagonal"
)

func LURunner() {
	fmt.Println("1.1-----LU-----")
	fmt.Println("Size:")
	var M [][]float64

	var mSize int
	fmt.Scan(&mSize)

	for i := 0; i < mSize; i++ {
		M = append(M, make([]float64, mSize))
	}

	fmt.Println("Matrix:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			var ij float64
			fmt.Scan(&ij)
			M[i][j] = ij
		}
	}

	fmt.Println("b:")
	b := make([]float64, mSize)
	for i := 0; i < mSize; i++ {
		var bb float64
		fmt.Scan(&bb)
		b[i] = bb
	}

	L, U, P := LU.GetLU(M, mSize)

	fmt.Println("L:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.2f ", L[i][j])
		}
		fmt.Print("\n")
	}

	fmt.Println("U:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.2f ", U[i][j])
		}
		fmt.Print("\n")
	}

	fmt.Println("P:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.0f ", P[i][j])
		}
		fmt.Print("\n")
	}

	det := LU.DetLU(M, mSize)
	fmt.Println("det: ", det)

	INV := LU.InvertLU(M, mSize)
	fmt.Println("inv:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.2f ", INV[i][j])
		}
		fmt.Print("\n")
	}

	solve := LU.SolveLU(M, b, mSize)
	fmt.Println("solve:")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.6f\n", solve[i])
	}
}

func TridiagRunner() {
	fmt.Println("1.2-----TRIDIAGONAL-----")
	fmt.Println("Size:")
	var M [][]float64

	var mSize int
	fmt.Scan(&mSize)

	for i := 0; i < mSize; i++ {
		M = append(M, make([]float64, mSize))
	}

	fmt.Println("Matrix:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			var ij float64
			fmt.Scan(&ij)
			M[i][j] = ij
		}
	}

	fmt.Println("b:")
	b := make([]float64, mSize)
	for i := 0; i < mSize; i++ {
		var bb float64
		fmt.Scan(&bb)
		b[i] = bb
	}

	A, B, C := tridiagonal.MatrixToTridiagonal(M, mSize)

	fmt.Print("A: ")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.2f ", A[i])
	}
	fmt.Print("\nB: ")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.2f ", B[i])
	}
	fmt.Print("\nC: ")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.2f ", C[i])
	}
	fmt.Println()

	solve2 := tridiagonal.SolveTridiagonal(M, b, mSize)
	fmt.Println("solve:")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.6f\n", solve2[i])
	}
}

func IterationsRunner() {
	fmt.Println("1.3-----ITERATIONS-----")
	fmt.Println("Size:")
	var M [][]float64

	var mSize int
	fmt.Scan(&mSize)

	for i := 0; i < mSize; i++ {
		M = append(M, make([]float64, mSize))
	}

	fmt.Println("Matrix:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			var ij float64
			fmt.Scan(&ij)
			M[i][j] = ij
		}
	}

	fmt.Println("b:")
	b := make([]float64, mSize)
	for i := 0; i < mSize; i++ {
		var bb float64
		fmt.Scan(&bb)
		b[i] = bb
	}

	fmt.Println("eps:")
	var eps float64
	fmt.Scan(&eps)

	solve, count := iterations.SolveSimpleIt(M, b, eps, mSize)

	fmt.Println("iterations solve:")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.6f\n", solve[i])
	}

	fmt.Println("count:")
	fmt.Printf("%d\n", count)

	solve2, count2 := iterations.SolveZeidel(M, b, eps, mSize)
	fmt.Println("\nZeidel solve:")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.6f\n", solve2[i])
	}

	fmt.Println("count:")
	fmt.Printf("%d\n", count2)
}

func RotationRunner() {
	fmt.Println("1.4-----ROTATION-----")
	fmt.Println("Size:")
	var M [][]float64

	var mSize int
	fmt.Scan(&mSize)

	for i := 0; i < mSize; i++ {
		M = append(M, make([]float64, mSize))
	}

	fmt.Println("Matrix:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			var ij float64
			fmt.Scan(&ij)
			M[i][j] = ij
		}
	}

	fmt.Println("eps:")
	var eps float64
	fmt.Scan(&eps)

	eValue, eVector, count := rotations.Rotation(M, mSize, eps)

	fmt.Println("eigenvectors:")
	for i := 0; i < mSize; i++ {
		fmt.Print("(")
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.4f ", eVector[j][i])
		}
		fmt.Print(")\n")
	}

	fmt.Println("\neigenvalues:")
	for i := 0; i < mSize; i++ {
		fmt.Printf("%.6f\n", eValue[i][i])
	}

	fmt.Println()
	for i := 0; i < mSize; i++ {
		var eVecVal []float64
		var eCurrVec []float64
		for j := 0; j < mSize; j++ {
			eVecVal = append(eVecVal, eValue[i][i]*eVector[j][i])
			eCurrVec = append(eCurrVec, eVector[j][i])
		}
		fmt.Printf("Check %d:\n", i+1)
		fmt.Print("(")
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.4f ", eVecVal[j])
		}
		fmt.Print(")\n(")
		MxEV := iterations.MatrixVectorMult(M, eCurrVec)
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.4f ", MxEV[j])
		}
		fmt.Print(")\n")
	}

	fmt.Println("\ncount:")
	fmt.Printf("%d\n", count)
}

func QRRunner() {
	fmt.Println("1.5-----QR-----")
	fmt.Println("Size:")
	var M [][]float64

	var mSize int
	fmt.Scan(&mSize)

	for i := 0; i < mSize; i++ {
		M = append(M, make([]float64, mSize))
	}

	fmt.Println("Matrix:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			var ij float64
			fmt.Scan(&ij)
			M[i][j] = ij
		}
	}

	fmt.Println("eps:")
	var eps float64
	fmt.Scan(&eps)

	Q, R := QR.GetQR(M, mSize)
	fmt.Println("Q:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.4f ", Q[i][j])
		}
		fmt.Print("\n")
	}
	fmt.Println("\nR:")
	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.4f ", R[i][j])
		}
		fmt.Print("\n")
	}

	solve, count := QR.QRAlgo(M, mSize, eps)
	fmt.Println("\neigenvalues:")
	for j := 0; j < mSize; j++ {
		if imag(solve[j]) == 0.0 {
			fmt.Printf("%.4f ", real(solve[j]))
		} else {
			if imag(solve[j]) > 0.0 {
				fmt.Printf("%.4f+%.4fj ", real(solve[j]), imag(solve[j]))
			} else {
				fmt.Printf("%.4f%.4fj ", real(solve[j]), imag(solve[j]))
			}
		}
	}
	fmt.Println("\ncount:")
	fmt.Printf("%d\n", count)
}

func main() {
	fmt.Println("Select lab 1.X:")
	fmt.Println("1: LU-decomposition")
	fmt.Println("2: Tridiagonal matrix algo")
	fmt.Println("3: Iterations/Zeidel algo")
	fmt.Println("4: Rotations method")
	fmt.Println("5: QR-decomposition")

	var point int
	fmt.Scan(&point)
	switch point {
	case 1:
		LURunner()
	case 2:
		TridiagRunner()
	case 3:
		IterationsRunner()
	case 4:
		RotationRunner()
	case 5:
		QRRunner()
	default:
		fmt.Println("I don't know u wrong")
	}
}
