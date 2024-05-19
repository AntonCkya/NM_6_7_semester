package cauchy

import (
	"math"
)

func DiffEuler(system []func(x, y, z float64) float64, yKey, zKey []float64, xl, xr, h float64) ([]float64, []float64, []float64) {
	var x, y, z []float64
	x = append(x, xl)
	y = append(y, yKey[1])
	z = append(z, zKey[1])
	count := 0
	for x[count]+0.0000001 < xr {
		x = append(x, x[count]+h)
		y = append(y, y[count]+h*system[0](x[count], y[count], z[count]))
		z = append(z, z[count]+h*system[1](x[count], y[count], z[count]))
		count++
	}
	return x, y, z
}

func DiffRungeKutt(system []func(x, y, z float64) float64, yKey, zKey []float64, xl, xr, h float64) ([]float64, []float64, []float64) {
	var x, y, z []float64
	x = append(x, xl)
	y = append(y, yKey[1])
	z = append(z, zKey[1])
	count := 0
	for x[count]+0.0000001 < xr {
		K1 := h * system[0](x[count], y[count], z[count])
		L1 := h * system[1](x[count], y[count], z[count])
		K2 := h * system[0](x[count]+h/2, y[count]+K1/2, z[count]+L1/2)
		L2 := h * system[1](x[count]+h/2, y[count]+K1/2, z[count]+L1/2)
		K3 := h * system[0](x[count]+h/2, y[count]+K2/2, z[count]+L2/2)
		L3 := h * system[1](x[count]+h/2, y[count]+K2/2, z[count]+L2/2)
		K4 := h * system[0](x[count]+h, y[count]+K3, z[count]+L3)
		L4 := h * system[1](x[count]+h, y[count]+K3, z[count]+L3)
		deltaY := (K1 + 2*K2 + 2*K3 + K4) / 6
		deltaZ := (L1 + 2*L2 + 2*L3 + L4) / 6
		x = append(x, x[count]+h)
		y = append(y, y[count]+deltaY)
		z = append(z, z[count]+deltaZ)
		count++
	}
	return x, y, z
}

func DiffAdams(system []func(x, y, z float64) float64, yKey, zKey []float64, xl, xr, h float64) ([]float64, []float64, []float64) {
	x, y, z := DiffRungeKutt(system, yKey, zKey, xl, xl+4*h, h)
	count := len(x) - 1
	for x[count]+0.0000001 < xr {
		x = append(x, x[count]+h)
		y = append(y, y[count]+(h/24)*(55*system[0](x[count], y[count], z[count])-59*system[0](x[count-1], y[count-1], z[count-1])+37*system[0](x[count-2], y[count-2], z[count-2])-9*system[0](x[count-3], y[count-3], z[count-3])))
		z = append(z, z[count]+(h/24)*(55*system[1](x[count], y[count], z[count])-59*system[1](x[count-1], y[count-1], z[count-1])+37*system[1](x[count-2], y[count-2], z[count-2])-9*system[1](x[count-3], y[count-3], z[count-3])))
		count++
	}
	return x, y, z
}

func RRRmethod(y1, y2 []float64, m float64) float64 {
	res := 0.0
	for i := 0; i < len(y1); i++ {
		res = math.Max(res, math.Abs(y1[i]-y2[i*2])/math.Pow(2.0, m))
	}
	return res
}

func AbsoluteError(y1, y2 []float64) float64 {
	res := 0.0
	for i := 0; i < len(y1); i++ {
		res = math.Max(res, math.Abs(y1[i]-y2[i]))
	}
	return res
}
