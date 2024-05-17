package derivative

func DX(table [][]float64, x float64, i int) float64 {
	return (table[1][i+1]-table[1][i])/(table[0][i+1]-table[0][i]) + ((table[1][i+2]-table[1][i+1])/(table[0][i+2]-table[0][i+1])-(table[1][i+1]-table[1][i])/(table[0][i+1]-table[0][i]))/(table[0][i+2]-table[0][i])*(2*x-table[0][i]-table[0][i+1])
}

func D2X(table [][]float64, x float64, i int) float64 {
	return 2 * ((table[1][i+2]-table[1][i+1])/(table[0][i+2]-table[0][i+1]) - (table[1][i+1]-table[1][i])/(table[0][i+1]-table[0][i])) / (table[0][i+2] - table[0][i])
}
