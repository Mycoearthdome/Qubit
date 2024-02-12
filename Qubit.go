package main

import (
	"crypto/rand"
	"encoding/binary"
	"encoding/json"
	"flag"
	"fmt"
	"math"
	"math/cmplx"
	"os"
	"runtime"
	"sync"
	"time"

	"github.com/shopspring/decimal"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

type Qubit struct {
	Alpha complex128 // probability amplitude for |0⟩
	Beta  complex128 // probability amplitude for |1⟩
	Omega complex128 // probability amplitude for 0|⟩
	Theta complex128 // probability amplitude for 1|⟩
}

type QubitRI struct {
	AlphaReal float64
	AlphaImag float64
	BetaReal  float64
	BetaImag  float64
	OmegaReal float64
	OmegaImag float64
	ThetaReal float64
	ThetaImag float64
	Tick      int
}

func normalizeComplex(s0 complex128) complex128 {
	magnitude := cmplx.Abs(s0)
	return s0 / complex(magnitude, 0)
}

func AddComplex(s0 complex128, s1 complex128) complex128 {
	realPart1, imagPart1 := real(s0), imag(s0)
	realPart2, imagPart2 := real(s1), imag(s1)

	if !math.IsNaN(realPart1) && !math.IsNaN(realPart2) && !math.IsNaN(imagPart1) && !math.IsNaN(imagPart2) {
		if !math.IsInf(realPart1, 0) && !math.IsInf(realPart2, 0) && !math.IsInf(imagPart1, 0) && !math.IsInf(imagPart2, 0) {

			Real1 := decimal.NewFromFloat(realPart1)
			Imag1 := decimal.NewFromFloat(imagPart1)
			Real2 := decimal.NewFromFloat(realPart2)
			Imag2 := decimal.NewFromFloat(imagPart2)

			AddedReal := Real1.Add(Real2)
			AddedImag := Imag1.Add(Imag2)

			AddedRealF, _ := AddedReal.Float64()
			AddedImagF, _ := AddedImag.Float64()

			return complex(AddedRealF, AddedImagF)
		}
	}
	return s1
}

func Tick(wg *sync.WaitGroup, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T complex128, Moment int, magnitude, phase chan bool) {
	defer wg.Done()

	var mag bool = false
	var pha bool = false

	CollapseThreshold := complex(1, 1)

	switch Moment {
	case 1:
		if cmplx.Abs(AddComplex(Q1A, Q2A)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1O, Q2O)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1B, Q2B)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1T, Q2T)) == cmplx.Abs(CollapseThreshold) {
			//fmt.Println("Magnitude match")
			mag = true
		}

		if cmplx.Phase(AddComplex(Q1A, Q2A)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1O, Q2O)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1B, Q2B)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1T, Q2T)) == cmplx.Phase(CollapseThreshold) {
			//fmt.Println("Phase match")
			pha = true
		}
	case 2:
		if cmplx.Abs(AddComplex(Q1T, Q2A)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1A, Q2O)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1O, Q2B)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1B, Q2T)) == cmplx.Abs(CollapseThreshold) {
			//fmt.Println("Magnitude match")
			mag = true
		}

		if cmplx.Phase(AddComplex(Q1T, Q2A)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1A, Q2O)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1O, Q2B)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1B, Q2T)) == cmplx.Phase(CollapseThreshold) {
			//fmt.Println("Phase match")
			pha = true
		}
	case 3:
		if cmplx.Abs(AddComplex(Q1B, Q2A)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1T, Q2O)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1A, Q2B)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1O, Q2T)) == cmplx.Abs(CollapseThreshold) {
			//fmt.Println("Magnitude match")
			mag = true
		}

		if cmplx.Phase(AddComplex(Q1B, Q2A)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1T, Q2O)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1A, Q2B)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1O, Q2T)) == cmplx.Phase(CollapseThreshold) {
			//fmt.Println("Phase match")
			pha = true
		}
	case 4:
		if cmplx.Abs(AddComplex(Q1O, Q2A)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1B, Q2O)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1T, Q2B)) == cmplx.Abs(CollapseThreshold) &&
			cmplx.Abs(AddComplex(Q1A, Q2T)) == cmplx.Abs(CollapseThreshold) {
			//fmt.Println("Magnitude match")
			mag = true
		}

		if cmplx.Phase(AddComplex(Q1O, Q2A)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1B, Q2O)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1T, Q2B)) == cmplx.Phase(CollapseThreshold) &&
			cmplx.Phase(AddComplex(Q1A, Q2T)) == cmplx.Phase(CollapseThreshold) {
			//fmt.Println("Phase match")
			pha = true
		}
	}

	if !mag && !pha {
		magnitude <- false
		phase <- false
	} else {
		if !mag && pha {
			magnitude <- false
			phase <- false
		}
		if mag && !pha {
			magnitude <- false
			phase <- false
		}
		if mag && pha {
			magnitude <- true
			phase <- true
		}
	}
}

func CollapsingQubits(q1 Qubit, q2 Qubit) (bool, int) {
	var wg sync.WaitGroup

	magnitude1 := make(chan bool, 1)
	magnitude2 := make(chan bool, 1)
	magnitude3 := make(chan bool, 1)
	magnitude4 := make(chan bool, 1)
	phase1 := make(chan bool, 1)
	phase2 := make(chan bool, 1)
	phase3 := make(chan bool, 1)
	phase4 := make(chan bool, 1)

	Q1A, Q2A := cmplx.Pow(q1.Alpha, 2), cmplx.Pow(q2.Alpha, 2)
	Q1O, Q2O := cmplx.Pow(q1.Omega, 2), cmplx.Pow(q2.Omega, 2)
	Q1B, Q2B := cmplx.Pow(q1.Beta, 2), cmplx.Pow(q2.Beta, 2)
	Q1T, Q2T := cmplx.Pow(q1.Theta, 2), cmplx.Pow(q2.Theta, 2)

	Q1A, Q2A = normalizeComplex(Q1A), normalizeComplex(Q2A)
	Q1O, Q2O = normalizeComplex(Q1O), normalizeComplex(Q2O)
	Q1B, Q2B = normalizeComplex(Q1B), normalizeComplex(Q2B)
	Q1T, Q2T = normalizeComplex(Q1T), normalizeComplex(Q2T)

	wg.Add(4)

	go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 1, magnitude1, phase1)
	go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 2, magnitude2, phase2)
	go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 3, magnitude3, phase3)
	go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 4, magnitude4, phase4)

	wg.Wait()

	if <-magnitude1 && <-phase1 {
		return true, 1
	}
	if <-magnitude2 && <-phase2 {
		return true, 2
	}
	if <-magnitude3 && <-phase3 {
		return true, 3
	}
	if <-magnitude4 && <-phase4 {
		return true, 4
	}
	return false, 0
}

func secureRandomFloat64() float64 {
	var randomBytes [8]byte

	_, err := rand.Read(randomBytes[:])
	if err != nil {
		panic(err)
	}

	randomFloat := math.Float64frombits(binary.BigEndian.Uint64(randomBytes[:]))

	return randomFloat
}

func RebuildFromTick(tick int, q1 Qubit) Qubit {

	var q3 Qubit

	switch tick {
	case 1:
		return q1
	case 2:
		q3.Alpha = q1.Theta
		q3.Omega = q1.Alpha
		q3.Beta = q1.Omega
		q3.Theta = q1.Beta
		return q3
	case 3:
		q3.Alpha = q1.Beta
		q3.Omega = q1.Theta
		q3.Beta = q1.Alpha
		q3.Theta = q1.Omega
		return q3
	case 4:
		q3.Alpha = q1.Omega
		q3.Omega = q1.Beta
		q3.Beta = q1.Theta
		q3.Theta = q1.Alpha
		return q3
	default:
		return q1
	}
}

func WriteQubits(JSONfilename string, WriteBuffer int64) map[int64][]QubitRI {
	var CollectedQubits int64
	var TotalCollectedQubits int64
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	dictQubit := make(map[int64][]QubitRI)
	var Duration int64
	for {
		if m.Sys-m.Alloc < 102400 { // 100MB
			break
		}
		for i := 0; CollectedQubits < WriteBuffer; i++ {
			var listQubit []QubitRI
			var Then int64 = time.Now().UnixNano()
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

			q1 := Qubit{Alpha: cmplx.Sqrt(normalizeComplex(complex(A1, A1i))),
				Omega: cmplx.Sqrt(normalizeComplex(complex(O1, O1i))),
				Beta:  cmplx.Sqrt(normalizeComplex(complex(B1, B1i))),
				Theta: cmplx.Sqrt(normalizeComplex(complex(T1, T1i)))}

			q2 := Qubit{Alpha: cmplx.Sqrt(normalizeComplex(complex(A2, A2i))),
				Omega: cmplx.Sqrt(normalizeComplex(complex(O2, O2i))),
				Beta:  cmplx.Sqrt(normalizeComplex(complex(B2, B2i))),
				Theta: cmplx.Sqrt(normalizeComplex(complex(T2, T2i)))}

			Result, tick := CollapsingQubits(q1, q2)
			if Result {

				q3 := RebuildFromTick(tick, q1)

				q3RI := QubitRI{
					AlphaReal: real(q3.Alpha),
					AlphaImag: imag(q3.Alpha),
					OmegaReal: real(q3.Omega),
					OmegaImag: imag(q3.Omega),
					BetaReal:  real(q3.Beta),
					BetaImag:  imag(q3.Beta),
					ThetaReal: real(q3.Theta),
					ThetaImag: imag(q3.Theta),
					Tick:      tick,
				}
				q2RI := QubitRI{
					AlphaReal: real(q2.Alpha),
					AlphaImag: imag(q2.Alpha),
					OmegaReal: real(q2.Omega),
					OmegaImag: imag(q2.Omega),
					BetaReal:  real(q2.Beta),
					BetaImag:  imag(q2.Beta),
					ThetaReal: real(q2.Theta),
					ThetaImag: imag(q2.Theta),
					Tick:      tick,
				}
				Duration = time.Now().UnixNano() - Then
				listQubit = append(listQubit, q3RI, q2RI)
				dictQubit[Duration] = listQubit
				fmt.Printf("\rFound Qubits = %d", TotalCollectedQubits)
				CollectedQubits++
				TotalCollectedQubits++
			}
		}
		CollectedQubits = 0
		fmt.Printf("...Saving ... %d Qubits to file %s\n", TotalCollectedQubits, JSONfilename)
		jsonData, err := json.Marshal(dictQubit)
		if err != nil {
			panic(err)
		}

		file, err := os.Create(JSONfilename)
		if err != nil {
			panic(err)
		}

		_, err = file.Write(jsonData)
		if err != nil {
			panic(err)
		}

		file.Close()
	}
	return dictQubit
}

func readFileIntoBytes(filePath string) ([]byte, error) {
	// Open the file
	file, err := os.Open(filePath)
	if err != nil {
		return nil, fmt.Errorf("ERROR OPENNING FILE: %v", err)
	}
	defer file.Close()

	// Get the file size
	fileInfo, err := file.Stat()
	if err != nil {
		return nil, fmt.Errorf("ERROR GETTING FILE INFO: %v", err)
	}
	fileSize := fileInfo.Size()

	// Read the entire file into a byte slice
	fileContent := make([]byte, fileSize)
	_, err = file.Read(fileContent)
	if err != nil {
		return nil, fmt.Errorf("ERROR READING FILE: %v", err)
	}

	return fileContent, nil
}

func ReadJSONtoDict(JSONfilename string, dictQubit map[int64][]QubitRI) {
	data, err := readFileIntoBytes(JSONfilename)
	if err != nil {
		panic(err)
	}

	err = json.Unmarshal(data, &dictQubit)
	if err != nil {
		panic(err)
	}
}

func cumulativeSum(numbers []int64) []int64 {
	var cumSum []int64
	var sum int64 = 0

	for _, num := range numbers {
		sum += num
		cumSum = append(cumSum, sum)
	}

	return cumSum
}

func Graph(lstTickOrQubits []int64, tick bool) {

	// Calculate the cumulative sum
	cumSum := cumulativeSum(lstTickOrQubits)

	// Create a new plot
	p := plot.New()

	// Create scatter plot for original data
	pts := make(plotter.XYs, len(lstTickOrQubits))
	for i, num := range lstTickOrQubits {
		pts[i].X = float64(i + 1)
		pts[i].Y = float64(num)
	}
	scatter, err := plotter.NewScatter(pts)
	if err != nil {
		panic(err)
	}

	// Create line plot for cumulative sum
	line, err := plotter.NewLine(plotter.XYs{
		{X: 1, Y: float64(cumSum[0])},
		{X: float64(len(lstTickOrQubits)), Y: float64(cumSum[len(lstTickOrQubits)-1])},
	})
	if err != nil {
		panic(err)
	}

	// Add the scatter and line plots to the plot
	p.Add(scatter, line)

	// Set axis labels
	p.X.Label.Text = "Index"
	p.Y.Label.Text = "Values"

	// Save the plot to a PNG file

	if !tick {
		if err := p.Save(12*vg.Inch, 8*vg.Inch, "Qubits_cumulative_sum_plot.png"); err != nil {
			panic(err)
		}
		printf("Saved Qubits graph... Qubits_cumulative_sum_plot.png")
	} else {
		if err := p.Save(12*vg.Inch, 8*vg.Inch, "Ticks_cumulative_sum_plot.png"); err != nil {
			panic(err)
		}
		printf("Saved Ticks graph... Ticks_cumulative_sum_plot.png")
	}

}

func main() {

	var WriteBufferThreshold int64 = 1000
	dictQubit := make(map[int64][]QubitRI)
	var listNanoseconds []int64
	var listTick []int64

	var w string
	var r string

	flag.StringVar(&w, "w", "", "Write Qubits")
	flag.StringVar(&r, "r", "", "Read Qubits")

	flag.Parse()

	if len(os.Args) < 3 {
		fmt.Println("Usage: Qubit -w/-r <Qubits file>")
		os.Exit(1)
	}

	if w != "" {
		JSON_OUT_FILENAME := w
		dictQubit = WriteQubits(JSON_OUT_FILENAME, WriteBufferThreshold)
	} else {
		if r != "" {
			JSON_OUT_FILENAME := r
			ReadJSONtoDict(JSON_OUT_FILENAME, dictQubit)
		} else {
			fmt.Println("Usage: Qubit -w/-r <Qubits file>")
			os.Exit(1)
		}
	}

	for timestamp, listQubit := range dictQubit {
		/*
			q1 := Qubit{Alpha: complex(listQubit[0].AlphaReal, listQubit[0].AlphaImag),
				Beta:  complex(listQubit[0].BetaReal, listQubit[0].BetaImag),
				Omega: complex(listQubit[0].OmegaReal, listQubit[0].OmegaImag),
				Theta: complex(listQubit[0].ThetaReal, listQubit[0].ThetaImag),
			}

			q2 := Qubit{Alpha: complex(listQubit[1].AlphaReal, listQubit[1].AlphaImag),
				Beta:  complex(listQubit[1].BetaReal, listQubit[1].BetaImag),
				Omega: complex(listQubit[1].OmegaReal, listQubit[1].OmegaImag),
				Theta: complex(listQubit[1].ThetaReal, listQubit[1].ThetaImag),
			}
		*/
		listTick = append(listTick, int64(listQubit[0].Tick))
		listNanoseconds = append(listNanoseconds, timestamp)
	}

	Graph(listNanoseconds, false)
	Graph(listTick, true)
}
