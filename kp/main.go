package main

import (
	"fmt"
	"runtime"
	"time"

	rotations "github.com/AntonCkya/numeric_methods/rotations"
)

func RotationRunner(M [][]float64, mSize int, eps float64) {
	rotations.Rotation(M, mSize, eps)
}

func main() {
	var M [][]float64

	var mSize int
	fmt.Scan(&mSize)

	for i := 0; i < mSize; i++ {
		M = append(M, make([]float64, mSize))
	}

	for i := 0; i < mSize; i++ {
		for j := 0; j < mSize; j++ {
			var ij float64
			fmt.Scan(&ij)
			M[i][j] = ij
		}
	}

	var eps float64
	fmt.Scan(&eps)
	fmt.Println(mSize)
	var count int
	fmt.Scan(&count)
	for i := 0; i < count; i++ {
		var procs int
		fmt.Scan(&procs)
		runtime.GOMAXPROCS(procs)
		start := time.Now()
		RotationRunner(M, mSize, eps)
		elapsed := time.Since(start)
		fmt.Println(procs, ":", elapsed)
	}
}

