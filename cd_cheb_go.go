package cd_cheb_go

//Chebyshev polynomials for calculating position and velocity of celestial objects
// using NASA JPL ephemerides

func Chebyshev(order int, x float64, data []float64) float64 {

	// Evaluate a Chebyshev polynomial
	var bk float64
	two_x := 2 * x
	bkp2 := data[order]
	bkp1 := two_x*bkp2 + data[order-1]

	for n := order - 2; n > 0; n-- {
		bk = data[n] + two_x*bkp1 - bkp2
		bkp2 = bkp1
		bkp1 = bk
	}
	return data[0] + x*bkp1 - bkp2
}

func DerChebyshev(order int, x float64, data []float64) float64 {
	//Evaluate the derivative of a Chebyshev polynomial

	var bk float64
	two_x := 2 * x
	bkp2 := float64(order) * data[order]
	bkp1 := two_x*bkp2 + float64(order-1)*data[order-1]

	for n := order - 2; n > 1; n-- {
		bk = float64(n)*data[n] + two_x*bkp1 - bkp2
		bkp2 = bkp1
		bkp1 = bk
	}

	return data[1] + two_x*bkp1 - bkp2
}
