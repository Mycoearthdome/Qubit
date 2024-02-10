package main

import (
	"crypto/rand"
	"encoding/binary"
	"fmt"
	"math"
	"math/cmplx"

	"github.com/shopspring/decimal"
)

type Qubit struct { //
	//gives a qubit structures as :                           |	A O\		/A O		/|
	Alpha complex128 // probability amplitude for |0⟩    0	  |			+ 			=  / |
	Beta  complex128 // probability amplitude for |1⟩    1    | B T/		\B T		 |
	Omega complex128 // probability amplitude for 0|⟩      0  |							 |
	Theta complex128 // probability amplitude for 1|)      1  |                        -----

}

func AddComplex(s0 complex128, s1 complex128) complex128 {

	var retAddedReal float64
	var retAddedImag float64

	realPart1 := real(s0)
	imagPart1 := imag(s0)
	realPart2 := real(s1)
	imagPart2 := imag(s1)

	if !math.IsNaN(realPart1) && !math.IsNaN(realPart2) && !math.IsNaN(imagPart1) && !math.IsNaN(imagPart2) {
		if !math.IsInf(realPart1, 0) && !math.IsInf(realPart2, 0) && !math.IsInf(imagPart1, 0) && !math.IsInf(imagPart2, 0) {
			Real1 := decimal.NewFromFloat(realPart1)
			Imag1 := decimal.NewFromFloat(imagPart1)
			Real2 := decimal.NewFromFloat(realPart2)
			Imag2 := decimal.NewFromFloat(imagPart2)

			AddedReal := Real1.Add(Real2).String()
			AddedImag := Imag1.Add(Imag2).String()

			if len(AddedReal) == 1 && int(AddedReal[0]) == 49 { // 1
				retAddedReal = 1.0
			} else {
				retAddedReal = float64(0)
			}
			if len(AddedImag) == 1 && int(AddedReal[0]) == 49 { // 1
				retAddedImag = 1.0
			} else {
				retAddedImag = float64(0)
			}
		}
	}
	return complex(retAddedReal, retAddedImag)
}

func CollapsingQubits(q1 Qubit, q2 Qubit) bool {
	var magnitude bool = false
	var phase bool = false
	CollapseThreshold := complex(1, 1)

	Q1A := cmplx.Pow(q1.Alpha, 2)
	Q2A := cmplx.Pow(q2.Alpha, 2)
	Q1O := cmplx.Pow(q1.Omega, 2)
	Q2O := cmplx.Pow(q2.Omega, 2)
	Q1B := cmplx.Pow(q1.Beta, 2)
	Q2B := cmplx.Pow(q2.Beta, 2)
	Q1T := cmplx.Pow(q1.Theta, 2)
	Q2T := cmplx.Pow(q2.Theta, 2)

	if cmplx.Abs(AddComplex(Q1A, Q2A)) == cmplx.Abs(CollapseThreshold) {
		if cmplx.Abs(AddComplex(Q1O, Q2O)) == cmplx.Abs(CollapseThreshold) {
			if cmplx.Abs(AddComplex(Q1B, Q2B)) == cmplx.Abs(CollapseThreshold) {
				if cmplx.Abs(AddComplex(Q1T, Q2T)) == cmplx.Abs(CollapseThreshold) {
					fmt.Println("Magnitude match")
					magnitude = true
				}
			}
		}
	}

	if cmplx.Phase(AddComplex(Q1A, Q2A)) == cmplx.Phase(CollapseThreshold) {
		if cmplx.Phase(AddComplex(Q1O, Q2O)) == cmplx.Phase(CollapseThreshold) {
			if cmplx.Phase(AddComplex(Q1B, Q2B)) == cmplx.Phase(CollapseThreshold) {
				if cmplx.Phase(AddComplex(Q1T, Q2T)) == cmplx.Phase(CollapseThreshold) {
					fmt.Println("Phase match")
					phase = true
				}
			}
		}
	}

	if magnitude && phase {
		return true
	} else {
		return false
	}
}

func secureRandomFloat64() float64 {

	var randomBytes [8]byte

	// Read 8 random bytes from crypto/rand
	_, err := rand.Read(randomBytes[:])
	if err != nil {
		panic(err)
	}

	// Convert the random bytes to a float64 in the range [0.0, 1.0)
	randomFloat := math.Float64frombits(binary.BigEndian.Uint64(randomBytes[:]))

	return randomFloat
}

func main() {

	for {
		// QBit 1
		var A1 float64 = secureRandomFloat64()
		var A1i float64 = secureRandomFloat64()
		var O1 float64 = secureRandomFloat64()
		var O1i float64 = secureRandomFloat64()
		var B1 float64 = secureRandomFloat64()
		var B1i float64 = secureRandomFloat64()
		var T1 float64 = secureRandomFloat64()
		var T1i float64 = secureRandomFloat64()

		//Qbit 2
		var A2 float64 = secureRandomFloat64()
		var A2i float64 = secureRandomFloat64()
		var O2 float64 = secureRandomFloat64()
		var O2i float64 = secureRandomFloat64()
		var B2 float64 = secureRandomFloat64()
		var B2i float64 = secureRandomFloat64()
		var T2 float64 = secureRandomFloat64()
		var T2i float64 = secureRandomFloat64()

		// Initialize a qubit with some probability amplitudes
		q1 := Qubit{Alpha: cmplx.Sqrt(complex(A1, A1i)), Omega: cmplx.Sqrt(complex(O1, O1i)), Beta: cmplx.Sqrt(complex(B1, B1i)), Theta: cmplx.Sqrt(complex(T1, T1i))}
		q2 := Qubit{Alpha: cmplx.Sqrt(complex(A2, A2i)), Omega: cmplx.Sqrt(complex(O2, O2i)), Beta: cmplx.Sqrt(complex(B2, B2i)), Theta: cmplx.Sqrt(complex(T2, T2i))}

		if CollapsingQubits(q1, q2) {
			fmt.Println("They Collapsed to their respective states.")
			fmt.Println("q1=", q1)
			fmt.Println("q2=", q2)
			break
		}

	}
}
