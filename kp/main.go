package main

import (
	"fmt"

	rotations "github.com/AntonCkya/numeric_methods/rotations"
)

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
		MxEV := rotations.MatrixVectorMult(M, eCurrVec)
		for j := 0; j < mSize; j++ {
			fmt.Printf("%.4f ", MxEV[j])
		}
		fmt.Print(")\n")
	}

	fmt.Println("\ncount:")
	fmt.Printf("%d\n", count)
}

func main() {
	RotationRunner()
}
