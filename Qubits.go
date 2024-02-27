package main

import (
	"bytes"
	"encoding/binary"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"math/cmplx"
	"math/rand"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/biogo/biogo/alphabet"
)

const apiToken = "your_api_token"
const apiEndpoint = "https://quantum-computing.ibm.com/api/YOUR_API_VERSION"

type QubitRI struct {
	AlphaReal float64
	AlphaImag float64
	BetaReal  float64
	BetaImag  float64
	OmegaReal float64
	OmegaImag float64
	ThetaReal float64
	ThetaImag float64
}

///// CREDITS ----> github.com/DylanMeeus/GoAudio/

func generate(datafile string, reqFrames int64, wfmt WaveFmt) ([]float64, []float64) {
	MinFloat64 := 2.2 * math.Pow10(-308)
	frames := make([]float64, reqFrames)
	var OutOfBoundsFrames []float64
	//osc, err := NewOscillator(wfmt.SampleRate, shape)
	//if err != nil {
	//	panic(err)
	//}
	//rand.Seed(time.Now().Unix())
	//var data2 [][]byte
	data, err := ReadDataFile(datafile)
	if err != nil {
		panic(err)
	}

	size := wfmt.BitsPerSample / 8 //32 = 4
	var FrameSamples [][]byte
	for i := 0; i < len(data); i = i + size {
		FrameSamples = append(FrameSamples, data[i:i+size])
	}
	j := 0
	for i := range frames {
		if i >= len(FrameSamples) {
			frames[i] = MinFloat64 //hexadecimal representation = NULL BYTES PADDING
		} else {
			frame := bitsToFloat(FrameSamples[j]) * math.Pow10(307) * math.Pow10(1)
			if frame >= MinFloat64 && frame <= math.MaxFloat64 {
				frames[i] = frame
			} else {
				frames[i] = 0
				OutOfBoundsFrames = append(OutOfBoundsFrames, frame)
			}
			//fmt.Println(frames[i])
			//data2 = append(data2, floatToBytes(float64(frames[i]), size))
			j++
		}
	}
	return frames, OutOfBoundsFrames
}

// Breakpoints is a collection of Breakpoint (time:value) pairs
type Breakpoints []Breakpoint

// Breakpoint represents a time:value pair
type Breakpoint struct {
	Time, Value float64
}

// BreakpointStream can be used to to treat breakpoints as a stream of data
// Each 'tick' can manipulate the state of the breakpoint stream
type BreakpointStream struct {
	Breakpoints     Breakpoints
	Left            Breakpoint
	Right           Breakpoint
	IndexLeft       int
	IndexRight      int
	CurrentPosition float64 // current position in timeframes
	Increment       float64
	Width           float64
	Height          float64
	HasMore         bool
}

// Tick returns the next value in the breakpoint stream
func (b *BreakpointStream) Tick() (out float64) {
	//if !b.HasMore {
	// permanently the last value
	//	return b.Right.Value
	//}
	if b.Width == 0.0 {
		out = b.Right.Value
	} else {
		// figure out value from linear interpolation
		frac := (float64(b.CurrentPosition) - b.Left.Time) / b.Width
		out = b.Left.Value + (b.Height * frac)
	}

	// prepare for next frame
	b.CurrentPosition += b.Increment
	if b.CurrentPosition > b.Right.Time {
		// move to next span
		b.IndexLeft++
		b.IndexRight++
		if b.IndexRight < len(b.Breakpoints) {
			b.Left = b.Breakpoints[b.IndexLeft]
			b.Right = b.Breakpoints[b.IndexRight]
			b.Width = b.Right.Time - b.Left.Time
			b.Height = b.Right.Value - b.Left.Value
		} else {
			// no more points
			b.HasMore = false
		}
	}
	return out
}

// NewBreakpointStream represents a slice of breakpoints streamed at a given sample rate
func NewBreakpointStream(bs []Breakpoint, sr int) (*BreakpointStream, error) {
	if len(bs) == 0 {
		return nil, errors.New("Need at least two points to create a stream")
	}
	left, right := bs[0], bs[1]
	return &BreakpointStream{
		Breakpoints:     Breakpoints(bs),
		Increment:       1.0 / float64(sr),
		IndexLeft:       0,
		IndexRight:      1,
		CurrentPosition: 0,
		Left:            left,
		Right:           right,
		Width:           right.Time - left.Time,   // first span
		Height:          right.Value - left.Value, // diff of first span
		HasMore:         len(bs) > 0,
	}, nil
}

// ParseBreakpoints reads the breakpoints from an io.Reader
// and turns them into a slice.
// A file is expected to be [time: value] formatted
// Will panic if file format is wrong
// TODO: don't panic
func ParseBreakpoints(in io.Reader) ([]Breakpoint, error) {
	data, err := io.ReadAll(in)
	if err != nil {
		return nil, err
	}

	lines := strings.Split(string(data), "\n")

	brkpnts := []Breakpoint{}
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		parts := strings.Split(line, ":")
		if len(parts) != 2 {
			return brkpnts, err
		}
		time := parts[0]
		value := parts[1]

		tf, err := strconv.ParseFloat(time, 64)
		if err != nil {
			return brkpnts, err
		}
		vf, err := strconv.ParseFloat(value, 64)
		if err != nil {
			return brkpnts, err
		}

		brkpnts = append(brkpnts, Breakpoint{
			Time:  tf,
			Value: vf,
		})

	}
	return brkpnts, nil
}

// ValueAt returns the expected value at a given time (expressed as float64) by linear interpolation
// Returns the index at which we found our value as well as the value itself.
func ValueAt(bs []Breakpoint, time float64, startIndex int) (index int, value float64) {
	if len(bs) == 0 {
		return 0, 0
	}
	npoints := len(bs)

	// first we need to find a span containing our timeslot
	startSpan := startIndex // start of span
	for _, b := range bs[startSpan:] {
		if b.Time > time {
			break
		}
		startSpan++
	}

	// Our span is never-ending (the last point in our breakpoint file was hit)
	if startSpan == npoints {
		return startSpan, bs[startSpan-1].Value
	}

	left := bs[startSpan-1]
	right := bs[startSpan]

	// check for instant jump
	// 2 points having the same time...
	width := right.Time - left.Time

	if width == 0 {
		return startSpan, right.Value
	}

	frac := (time - left.Time) / width

	val := left.Value + ((right.Value - left.Value) * frac)

	return startSpan, val
}

// MinMaxValue returns the smallest and largest value found in the breakpoint file
func MinMaxValue(bs []Breakpoint) (smallest float64, largest float64) {
	// TODO: implement as SORT and return the first and last element
	if len(bs) == 0 {
		return
	}
	smallest = bs[0].Value
	largest = bs[0].Value
	for _, b := range bs[1:] {
		if b.Value < smallest {
			smallest = b.Value
		} else if b.Value > largest {
			largest = b.Value
		} else {
			// no op
		}
	}
	return
}

// Any returns true if any breakpoint matches the filter.
func (bs Breakpoints) Any(f func(Breakpoint) bool) bool {
	for _, b := range bs {
		if f(b) {
			return true
		}
	}
	return false
}

// DFT is a discrete fourier transformation on the input frames
// DEPRECATED
// Please use FFT unless you are sure  you want this one..
func DFT(input []Frame) []complex128 {
	N := len(input)

	output := make([]complex128, len(input))

	reals := make([]float64, len(input))
	imgs := make([]float64, len(input))
	for i, frame := range input {
		for n := 0; n < N; n++ {
			reals[i] += float64(frame) * math.Cos(float64(i*n)*tau/float64(N))
			imgs[i] += float64(frame) * math.Sin(float64(i*n)*tau/float64(N))
		}

		reals[i] /= float64(N)
		imgs[i] /= float64(N)
	}

	for i := 0; i < len(reals); i++ {
		output[i] = complex(reals[i], imgs[i])
	}

	return output
}

// HFFT mutates freqs!
func HFFT(input []Frame, freqs []complex128, n, step int) {
	if n == 1 {
		freqs[0] = complex(input[0], 0)
		return
	}

	h := n / 2

	HFFT(input, freqs, h, 2*step)
	HFFT(input[step:], freqs[h:], h, 2*step)

	for k := 0; k < h; k++ {
		a := -2 * math.Pi * float64(k) * float64(n)
		e := cmplx.Rect(1, a) * freqs[k+h]
		freqs[k], freqs[k+h] = freqs[k]+e, freqs[k]-e
	}
}

// FFT (Fast Fourier Transform) implementation
func FFT(input []Frame) []complex128 {
	freqs := make([]complex128, len(input))
	HFFT(input, freqs, len(input), 1)
	return freqs
}

// EqualTemperedNote
type EqualTemperedNote float64

// These constants represent the frequncies for the equal-tempered scale tuned to A4 = 440Hz
const (
	C0  EqualTemperedNote = 16.35
	C0S EqualTemperedNote = 17.32
	D0  EqualTemperedNote = 18.35
	D0S EqualTemperedNote = 19.45
	E0  EqualTemperedNote = 20.60
	F0  EqualTemperedNote = 21.83
	F0S EqualTemperedNote = 23.12
	G0  EqualTemperedNote = 24.50
	G0S EqualTemperedNote = 25.96
	A0  EqualTemperedNote = 27.50
	A0S EqualTemperedNote = 29.14
	B0  EqualTemperedNote = 30.87
	C1  EqualTemperedNote = 32.70
	C1S EqualTemperedNote = 34.65
	D1  EqualTemperedNote = 36.71
	D1S EqualTemperedNote = 38.89
	E1  EqualTemperedNote = 41.20
	F1  EqualTemperedNote = 43.65
	F1S EqualTemperedNote = 46.25
	G1  EqualTemperedNote = 49.00
	G1S EqualTemperedNote = 51.91
	A1  EqualTemperedNote = 55.00
	A1S EqualTemperedNote = 58.27
	B1  EqualTemperedNote = 61.74
	C2  EqualTemperedNote = 65.41
	C2S EqualTemperedNote = 69.30
	D2  EqualTemperedNote = 73.42
	D2S EqualTemperedNote = 77.78
	E2  EqualTemperedNote = 82.41
	F2  EqualTemperedNote = 87.31
	F2S EqualTemperedNote = 92.50
	G2  EqualTemperedNote = 98.00
	G2S EqualTemperedNote = 103.83
	A2  EqualTemperedNote = 110.00
	A2S EqualTemperedNote = 116.54
	B2  EqualTemperedNote = 123.47
	C3  EqualTemperedNote = 130.81
	C3S EqualTemperedNote = 138.59
	D3  EqualTemperedNote = 146.83
	D3S EqualTemperedNote = 155.56
	E3  EqualTemperedNote = 164.81
	F3  EqualTemperedNote = 174.61
	F3S EqualTemperedNote = 185.00
	G3  EqualTemperedNote = 196.00
	G3S EqualTemperedNote = 207.65
	A3  EqualTemperedNote = 220.00
	A3S EqualTemperedNote = 233.08
	B3  EqualTemperedNote = 246.94
	C4  EqualTemperedNote = 261.63
	C4S EqualTemperedNote = 277.18
	D4  EqualTemperedNote = 293.66
	D4S EqualTemperedNote = 311.13
	E4  EqualTemperedNote = 329.63
	F4  EqualTemperedNote = 349.23
	F4S EqualTemperedNote = 369.99
	G4  EqualTemperedNote = 392.00
	G4S EqualTemperedNote = 415.30
	A4  EqualTemperedNote = 440.00
	A4S EqualTemperedNote = 466.16
	B4  EqualTemperedNote = 493.88
	C5  EqualTemperedNote = 523.25
	C5S EqualTemperedNote = 554.37
	D5  EqualTemperedNote = 587.33
	D5S EqualTemperedNote = 622.25
	E5  EqualTemperedNote = 659.25
	F5  EqualTemperedNote = 698.46
	F5S EqualTemperedNote = 739.99
	G5  EqualTemperedNote = 783.99
	G5S EqualTemperedNote = 830.61
	A5  EqualTemperedNote = 880.00
	A5S EqualTemperedNote = 932.33
	B5  EqualTemperedNote = 987.77
	C6  EqualTemperedNote = 1046.50
	C6S EqualTemperedNote = 1108.73
	D6  EqualTemperedNote = 1174.66
	D6S EqualTemperedNote = 1244.51
	E6  EqualTemperedNote = 1318.51
	F6  EqualTemperedNote = 1396.91
	F6S EqualTemperedNote = 1479.98
	G6  EqualTemperedNote = 1567.98
	G6S EqualTemperedNote = 1661.22
	A6  EqualTemperedNote = 1760.00
	A6S EqualTemperedNote = 1864.66
	B6  EqualTemperedNote = 1975.53
	C7  EqualTemperedNote = 2093.00
	C7S EqualTemperedNote = 2217.46
	D7  EqualTemperedNote = 2349.32
	D7S EqualTemperedNote = 2489.02
	E7  EqualTemperedNote = 2636.02
	F7  EqualTemperedNote = 2793.83
	F7S EqualTemperedNote = 2959.96
	G7  EqualTemperedNote = 3135.96
	G7S EqualTemperedNote = 3322.44
	A7  EqualTemperedNote = 3520.00
	A7S EqualTemperedNote = 3729.31
	B7  EqualTemperedNote = 3951.07
	C8  EqualTemperedNote = 4186.01
	C8S EqualTemperedNote = 4434.92
	D8  EqualTemperedNote = 4698.63
	D8S EqualTemperedNote = 4978.03
	E8  EqualTemperedNote = 5274.04
	F8  EqualTemperedNote = 5587.65
	F8S EqualTemperedNote = 5919.91
	G8  EqualTemperedNote = 6271.93
	G8S EqualTemperedNote = 6644.88
	A8  EqualTemperedNote = 7040.00
	A8S EqualTemperedNote = 7458.62
	B8  EqualTemperedNote = 7902.13
)

// Lowpass applies a low-pass filter to the frames
// Does not modify the input signal
func Lowpass(fs []float64, freq, delay, sr float64) []float64 {
	output := make([]float64, len(fs))
	copy(output, fs)

	costh := 2. - math.Cos((tau*freq)/sr)
	coef := math.Sqrt(costh*costh-1.) - costh

	for i, a := range output {
		output[i] = a*(1+coef) - delay*coef
		delay = output[i]
	}

	return output
}

// Highpass applies a high-pass filter to the frames.
// Does not modify the input signal
func Highpass(fs []float64, freq, delay, sr float64) []float64 {
	output := make([]float64, len(fs))
	copy(output, fs)

	b := 2. - math.Cos(tau*freq/sr)
	coef := b - math.Sqrt(b*b-1.)

	for i, a := range output {
		output[i] = a*(1.-coef) - delay*coef
		delay = output[i]
	}

	return output
}

// Balance a signal (rescale output signal)
func Balance(signal, comparator, delay []float64, frequency, samplerate float64) []float64 {
	c := make([]float64, len(signal))
	copy(signal, c)

	costh := 2. - math.Cos(tau*frequency/samplerate)
	coef := math.Sqrt(costh*costh-1.) - costh

	for i, s := range signal {
		ss := signal[i]
		if signal[i] < 0 {
			ss = -s
		}
		delay[0] = ss*(1+coef) - (delay[0] * coef)

		if comparator[i] < 0 {
			comparator[i] = -comparator[i]
		}
		delay[1] = comparator[i]*(1+coef) - (delay[1] * coef)
		if delay[0] != 0 {
			c[i] = s * (delay[0] / delay[1])
		} else {
			c[i] = s * delay[1]
		}
	}
	return c
}

// Gtable is a Guard-table for oscillator lookup
type Gtable struct {
	data []float64
}

// Len returns the length of the data segment without the guard point
func Len(g *Gtable) int {
	return len(g.data) - 1
}

func NewGtable(data []float64) *Gtable {
	return &Gtable{data}
}

// NewSineTable returns a lookup table populated for sine-wave generation.
func NewSineTable(length int) *Gtable {
	g := &Gtable{}
	if length == 0 {
		return g
	}
	g.data = make([]float64, length+1) // one extra for the guard point.
	step := tau / float64(Len(g))
	for i := 0; i < Len(g); i++ {
		g.data[i] = math.Sin(step * float64(i))
	}
	// store a guard point
	g.data[len(g.data)-1] = g.data[0]
	return g
}

// NewTriangleTable generates a lookup table for a triangle wave
// of the specified length and with the requested number of harmonics.
func NewTriangleTable(length int, nharmonics int) (*Gtable, error) {
	if length == 0 || nharmonics == 0 || nharmonics >= length/2 {
		return nil, errors.New("invalid arguments for creation of Triangle Table")
	}

	g := &Gtable{}
	g.data = make([]float64, length+1)

	step := tau / float64(length)

	// generate triangle waveform
	harmonic := 1.0
	for i := 0; i < nharmonics; i++ {
		amp := 1.0 / (harmonic * harmonic)
		for j := 0; j < length; j++ {
			g.data[j] += amp * math.Cos(step*harmonic*float64(j))
		}
		harmonic += 2 // triangle wave has only odd harmonics
	}
	// normalize the values to be in the [-1;1] range
	g.data = normalize(g.data)
	return g, nil
}

// normalize the functions to the range -1, 1
func normalize(xs []float64) []float64 {
	length := len(xs)
	maxamp := 0.0
	for i := 0; i < length; i++ {
		amp := math.Abs(xs[i])
		if amp > maxamp {
			maxamp = amp
		}
	}

	maxamp = 1.0 / maxamp
	for i := 0; i < length; i++ {
		xs[i] *= maxamp
	}
	xs[len(xs)-1] = xs[0]
	return xs
}

// LookupOscillator is an oscillator that's more gentle on your CPU
// By performing a table lookup to generate the required waveform..
type LookupOscillator struct {
	Oscillator
	Table      *Gtable
	SizeOverSr float64 // convenience variable for calculations
}

// NewLookupOscillator creates a new oscillator which
// performs a table-lookup to generate the required waveform
func NewLookupOscillator(sr int, t *Gtable, phase float64) (*LookupOscillator, error) {
	if t == nil || len(t.data) == 0 {
		return nil, errors.New("invalid table provided for lookup oscillator")
	}

	return &LookupOscillator{
		Oscillator: Oscillator{
			curfreq:  0.0,
			curphase: float64(Len(t)) * phase,
			incr:     0.0,
		},
		Table:      t,
		SizeOverSr: float64(Len(t)) / float64(sr),
	}, nil

}

// TruncateTick performs a lookup and truncates the value
// index down (if the index for lookup = 10.5, return index 10)
func (l *LookupOscillator) TruncateTick(freq float64) float64 {
	return l.BatchTruncateTick(freq, 1)[0]
}

// BatchTruncateTick returns a slice of samples from the oscillator of the requested length
func (l *LookupOscillator) BatchTruncateTick(freq float64, nframes int) []float64 {
	out := make([]float64, nframes)
	for i := 0; i < nframes; i++ {
		index := l.curphase
		if l.curfreq != freq {
			l.curfreq = freq
			l.incr = l.SizeOverSr * l.curfreq
		}
		curphase := l.curphase
		curphase += l.incr
		for curphase > float64(Len(l.Table)) {
			curphase -= float64(Len(l.Table))
		}
		for curphase < 0.0 {
			curphase += float64(Len(l.Table))
		}
		l.curphase = curphase
		out[i] = l.Table.data[int(index)]
	}
	return out
}

// InterpolateTick performs a lookup but interpolates the value if the
// requested index does not appear in the table.
func (l *LookupOscillator) InterpolateTick(freq float64) float64 {
	return l.BatchInterpolateTick(freq, 1)[0]
}

// BatchInterpolateTick performs a lookup for N frames, and interpolates the value if the
// requested index does not appear in the table.
func (l *LookupOscillator) BatchInterpolateTick(freq float64, nframes int) []float64 {
	out := make([]float64, nframes)
	for i := 0; i < nframes; i++ {
		baseIndex := int(l.curphase)
		nextIndex := baseIndex + 1
		if l.curfreq != freq {
			l.curfreq = freq
			l.incr = l.SizeOverSr * l.curfreq
		}
		curphase := l.curphase
		frac := curphase - float64(baseIndex)
		val := l.Table.data[baseIndex]
		slope := l.Table.data[nextIndex] - val
		val += frac * slope
		curphase += l.incr

		for curphase > float64(Len(l.Table)) {
			curphase -= float64(Len(l.Table))
		}
		for curphase < 0.0 {
			curphase += float64(Len(l.Table))
		}

		l.curphase = curphase
		out[i] = val
	}
	return out
}

const tau = (2 * math.Pi)

// Shape for defining the different possible waveform shapes for use with the Oscillator
type Shape int

// Shapes for which we can generate waveforms
const (
	SINE Shape = iota
	SQUARE
	DOWNWARD_SAWTOOTH
	UPWARD_SAWTOOTH
	TRIANGLE
)

var (
	shapeCalcFunc = map[Shape]func(float64) float64{
		SINE:              sineCalc,
		SQUARE:            squareCalc,
		TRIANGLE:          triangleCalc,
		DOWNWARD_SAWTOOTH: downSawtoothCalc,
		UPWARD_SAWTOOTH:   upwSawtoothCalc,
	}
)

// Oscillator represents a wave-oscillator where each tick is calculated in the moment.
type Oscillator struct {
	curfreq  float64
	curphase float64
	incr     float64
	twopiosr float64 // (2*PI) / samplerate
	tickfunc func(float64) float64
}

// NewOscillator set to a given sample rate
func NewOscillator(sr int, shape Shape) (*Oscillator, error) {
	cf, ok := shapeCalcFunc[shape]
	if !ok {
		return nil, fmt.Errorf("Shape type %v not supported", shape)
	}
	return &Oscillator{
		twopiosr: tau / float64(sr),
		tickfunc: cf,
	}, nil
}

// NewPhaseOscillator creates a new oscillator where the initial phase is offset
// by a given phase
func NewPhaseOscillator(sr int, phase float64, shape Shape) (*Oscillator, error) {
	cf, ok := shapeCalcFunc[shape]
	if !ok {
		return nil, fmt.Errorf("Shape type %v not supported", shape)
	}
	return &Oscillator{
		twopiosr: tau / float64(sr),
		tickfunc: cf,
		curphase: tau * phase,
	}, nil
}

// Tick generates the next value of the oscillator waveform at a given frequency in Hz
func (o *Oscillator) Tick(freq float64) float64 {
	if o.curfreq != freq {
		o.curfreq = freq
		o.incr = o.twopiosr * freq
	}
	val := o.tickfunc(o.curphase)
	o.curphase += o.incr
	if o.curphase >= tau {
		o.curphase -= tau
	}
	if o.curphase < 0 {
		o.curphase = tau
	}
	return val

}

func triangleCalc(phase float64) float64 {
	val := 2.0*(phase*(1.0/tau)) - 1.0
	if val < 0.0 {
		val = -val
	}
	val = 2.0 * (val - 0.5)
	return val
}

func upwSawtoothCalc(phase float64) float64 {
	val := 2.0*(phase*(1.0/tau)) - 1.0
	return val
}

func downSawtoothCalc(phase float64) float64 {
	val := 1.0 - 2.0*(phase*(1.0/tau))
	return val
}

func squareCalc(phase float64) float64 {
	val := -1.0
	if phase <= math.Pi {
		val = 1.0
	}
	return val
}

func sineCalc(phase float64) float64 {
	return math.Sin(phase)
}

// type aliases for conversion functions
type (
	bytesToIntF   func([]byte) int
	bytesToFloatF func([]byte) float64
)

var (
	// figure out which 'to int' function to use..
	byteSizeToIntFunc = map[int]bytesToIntF{
		16: bits16ToInt,
		24: bits24ToInt,
		32: bits32ToInt,
	}

	byteSizeToFloatFunc = map[int]bytesToFloatF{
		16: bitsToFloat,
		32: bitsToFloat,
		64: bitsToFloat,
	}

	// max value depending on the bit size
	maxValues = map[int]int{
		8:  math.MaxInt8,
		16: math.MaxInt16,
		32: math.MaxInt32,
		64: math.MaxInt64,
	}
)

func ReadDataFile(f string) ([]byte, error) {
	// open as read-only file
	file, err := os.Open(f)
	if err != nil {
		return []byte{}, err
	}
	defer file.Close()

	return io.ReadAll(file)
}

// ReadWaveFile parses a .wave file into a Wave struct
func ReadWaveFile(f string) (Wave, error) {
	// open as read-only file
	file, err := os.Open(f)
	if err != nil {
		return Wave{}, err
	}
	defer file.Close()

	return ReadWaveFromReader(file)
}

// ReadWaveFromReader parses an io.Reader into a Wave struct
func ReadWaveFromReader(reader io.Reader) (Wave, error) {
	data, err := io.ReadAll(reader)
	if err != nil {
		return Wave{}, err
	}

	data = deleteJunk(data)

	hdr := readHeader(data)

	wfmt := readFmt(data)

	wavdata := readData(data, wfmt)

	frames := parseRawData(wfmt, wavdata.RawData)
	wavdata.Frames = frames

	return Wave{
		WaveHeader: hdr,
		WaveFmt:    wfmt,
		WaveData:   wavdata,
	}, nil
}

// for our wave format we expect double precision floats
func bitsToFloat(b []byte) float64 {
	var bits uint64
	switch len(b) {
	case 2:
		bits = uint64(binary.LittleEndian.Uint16(b))
	case 4:
		bits = uint64(binary.LittleEndian.Uint32(b))
	case 8:
		bits = binary.LittleEndian.Uint64(b)
	default:
		panic("Can't parse to float..")
	}
	float := math.Float64frombits(bits)
	return float
}

func bits16ToInt(b []byte) int {
	if len(b) != 2 {
		panic("Expected size 4!")
	}
	var payload int16
	buf := bytes.NewReader(b)
	err := binary.Read(buf, binary.LittleEndian, &payload)
	if err != nil {
		// TODO: make safe
		panic(err)
	}
	return int(payload) // easier to work with ints
}

func bits24ToInt(b []byte) int {
	if len(b) != 3 {
		panic("Expected size 3!")
	}
	// add some padding to turn a 24-bit integer into a 32-bit integer
	b = append([]byte{0x00}, b...)
	var payload int32
	buf := bytes.NewReader(b)
	err := binary.Read(buf, binary.LittleEndian, &payload)
	if err != nil {
		// TODO: make safe
		panic(err)
	}
	return int(payload) // easier to work with ints
}

// turn a 32-bit byte array into an int
func bits32ToInt(b []byte) int {
	if len(b) != 4 {
		panic("Expected size 4!")
	}
	var payload int32
	buf := bytes.NewReader(b)
	err := binary.Read(buf, binary.LittleEndian, &payload)
	if err != nil {
		// TODO: make safe
		panic(err)
	}
	return int(payload) // easier to work with ints
}

func readData(b []byte, wfmt WaveFmt) WaveData {
	wd := WaveData{}

	start := 36 + wfmt.ExtraParamSize
	subchunk2ID := b[start : start+4]
	wd.Subchunk2ID = subchunk2ID

	subsize := bits32ToInt(b[start+4 : start+8])
	wd.Subchunk2Size = subsize

	wd.RawData = b[start+8:]

	return wd
}

// Should we do n-channel separation at this point?
func parseRawData(wfmt WaveFmt, rawdata []byte) []Frame {
	bytesSampleSize := wfmt.BitsPerSample / 8
	// TODO: sanity-check that this is a power of 2? I think only those sample sizes are
	// possible

	frames := []Frame{}
	// read the chunks
	for i := 0; i < len(rawdata); i += bytesSampleSize {
		rawFrame := rawdata[i : i+bytesSampleSize]
		unscaledFrame := byteSizeToIntFunc[wfmt.BitsPerSample](rawFrame)
		scaled := scaleFrame(unscaledFrame, wfmt.BitsPerSample)
		frames = append(frames, scaled)
	}
	return frames
}

func scaleFrame(unscaled, bits int) Frame {
	maxV := maxValues[bits]
	return Frame(float64(unscaled) / float64(maxV))

}

// deleteJunk will remove the JUNK chunks if they are present
func deleteJunk(b []byte) []byte {
	var junkStart, junkEnd int

	for i := 0; i < len(b)-4; i++ {
		if strings.ToLower(string(b[i:i+4])) == "junk" {
			junkStart = i
		}

		if strings.ToLower(string(b[i:i+3])) == "fmt" {
			junkEnd = i
		}
	}

	if junkStart != 0 {
		cpy := make([]byte, len(b[0:junkStart]))
		copy(cpy, b[0:junkStart])
		cpy = append(cpy, b[junkEnd:]...)
		return cpy
	}

	return b
}

// readFmt parses the FMT portion of the WAVE file
// assumes the entire binary representation is passed!
func readFmt(b []byte) WaveFmt {
	wfmt := WaveFmt{}
	subchunk1ID := b[12:16]
	wfmt.Subchunk1ID = subchunk1ID

	subchunksize := bits32ToInt(b[16:20])
	wfmt.Subchunk1Size = subchunksize

	format := bits16ToInt(b[20:22])
	wfmt.AudioFormat = format

	numChannels := bits16ToInt(b[22:24])
	wfmt.NumChannels = numChannels

	sr := bits32ToInt(b[24:28])
	wfmt.SampleRate = sr

	br := bits32ToInt(b[28:32])
	wfmt.ByteRate = br

	ba := bits16ToInt(b[32:34])
	wfmt.BlockAlign = ba

	bps := bits16ToInt(b[34:36])
	wfmt.BitsPerSample = bps

	// parse extra (optional) elements..

	if subchunksize != 16 {
		// only for compressed files (non-PCM)
		extraSize := bits16ToInt(b[36:38])
		wfmt.ExtraParamSize = extraSize
		wfmt.ExtraParams = b[38 : 38+extraSize]
	}

	return wfmt
}

// TODO: make safe.
func readHeader(b []byte) WaveHeader {
	// the start of the bte slice..
	hdr := WaveHeader{}
	hdr.ChunkID = b[0:4]
	if string(hdr.ChunkID) != "RIFF" {
		panic("Invalid file")
	}

	chunkSize := b[4:8]
	var size uint32
	buf := bytes.NewReader(chunkSize)
	err := binary.Read(buf, binary.LittleEndian, &size)
	if err != nil {
		panic(err)
	}
	hdr.ChunkSize = int(size) // easier to work with ints

	format := b[8:12]
	if string(format) != "WAVE" {
		panic("Format should be WAVE")
	}
	hdr.Format = string(format)
	return hdr
}

var (
	// noteIndex for use in calculations where a user passes a note
	noteIndex = map[string]int{
		"a":  0,
		"a#": 1,
		"bb": 1,
		"b":  2,
		"c":  3,
		"c#": 4,
		"db": 4,
		"d":  5,
		"d#": 6,
		"eb": 6,
		"e":  7,
		"f":  8,
		"f#": 9,
		"gb": 9,
		"g":  10,
		"g#": 11,
		"ab": 11,
	}
)

var (
	s      = struct{}{}
	valid  = map[string]interface{}{"a": s, "b": s, "c": s, "d": s, "e": s, "f": s, "g": s, "#": s}
	digits = map[string]interface{}{"0": s, "1": s, "2": s, "3": s, "4": s, "5": s, "6": s, "7": s, "8": s, "9": s}
)

// ADSR creates an attack -> decay -> sustain -> release envelope
// time durations are passes as seconds.
// returns the value + the current time
func ADSR(maxamp, duration, attacktime, decaytime, sus, releasetime, controlrate float64, currentframe int) float64 {
	dur := duration * controlrate
	at := attacktime * controlrate
	dt := decaytime * controlrate
	rt := releasetime * controlrate
	cnt := float64(currentframe)

	amp := 0.0
	if cnt < dur {
		if cnt <= at {
			// attack
			amp = cnt * (maxamp / at)
		} else if cnt <= (at + dt) {
			// decay
			amp = ((sus-maxamp)/dt)*(cnt-at) + maxamp
		} else if cnt <= dur-rt {
			// sustain
			amp = sus
		} else if cnt > (dur - rt) {
			// release
			amp = -(sus/rt)*(cnt-(dur-rt)) + sus
		}
	}

	return amp
}

// NoteToFrequency turns a given note & octave into a frequency
// using Equal-Tempered tuning with reference pitch = A440
func NoteToFrequency(note string, octave int) float64 {
	// TODO: Allow for tuning systems other than Equal-Tempered A440?
	// clean the input
	note = strings.ToLower(strings.TrimSpace(note))
	ni := noteIndex[note]
	if ni >= 3 {
		// correct for octaves starting at C, not A.
		octave--
	}
	FR := 440.
	// we adjust the octave (-4) as the reference frequency is in the fourth octave
	// this effectively allows us to generate any octave above or below the reference octave
	return FR * math.Pow(2, float64(octave-4)+(float64(ni)/12.))
}

// parseNoteOctave returns the note + octave value
func parseNoteOctave(note string) (string, int, error) {
	note = strings.ToLower(note)
	notePart := strings.Map(func(r rune) rune {
		if _, ok := valid[string(r)]; !ok {
			return rune(-1)
		}
		return r
	}, note)

	digitPart := strings.Map(func(r rune) rune {
		if _, ok := digits[string(r)]; !ok {
			return rune(-1)
		}
		return r
	}, note[len(notePart):])

	octave, err := strconv.Atoi(digitPart)
	if err != nil {
		return "", 0, err
	}

	return notePart, octave, nil
}

// ParseNoteToFrequency tries to parse a string representation of a note+octave (e.g C#4)
// and will return a float64 frequency value using 'NoteToFrequency'
func ParseNoteToFrequency(note string) (float64, error) {
	nt, oct, err := parseNoteOctave(note)
	if err != nil {
		return -1, err
	}
	return NoteToFrequency(nt, oct), nil
}

// representation of the wave file, used by reader.go and writer.go

// Frame is a single float64 value of raw audio data
type Frame float64

/*

╔════════╤════════════════╤══════╤═══════════════════════════════════════════════════╗
║ Offset │ Field          │ Size │ -- start of header                                ║
╠════════╪════════════════╪══════╪═══════════════════════════════════════════════════╣
║ 0      │ ChunkID        │ 4    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 4      │ ChunkSize      │ 4    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 8      │ Format         │ 8    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ --     │ --             │ --   │ -- start of fmt                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 12     │ SubchunkID     │ 4    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 16     │ SubchunkSize   │ 4    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 20     │ AudioFormat    │ 2    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 22     │ NumChannels    │ 2    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 24     │ SampleRate     │ 4    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 28     │ ByteRate       │ 4    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 32     │ BlockAlign     │ 2    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 34     │ BitsPerSample  │ 2    │                                                   ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ * 36   │ ExtraParamSize │ 2    │ Optional! Only when not PCM                       ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ * 38   │ ExtraParams    │ *    │ Optional! Only when not PCM                       ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ --     │ --             │ --   │ -- start of data, assuming PCM                    ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 36     │ Subchunk2ID    │ 4    │ (offset by extra params of subchunk 1 if not PCM) ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 40     │ SubchunkSize   │ 4    │ (offset by extra params of subchunk 1 if not PCM) ║
╟────────┼────────────────┼──────┼───────────────────────────────────────────────────╢
║ 44     │ Data           │ *    │ (offset by extra params of subchunk 1 if not PCM) ║
╚════════╧════════════════╧══════╧═══════════════════════════════════════════════════╝


*/

// Wave represents an entire .wav audio file
type Wave struct {
	WaveHeader
	WaveFmt
	WaveData
}

// WaveHeader describes the header each WAVE file should start with
type WaveHeader struct {
	ChunkID   []byte // should be RIFF on little-endian or RIFX on big-endian systems..
	ChunkSize int
	Format    string // sanity-check, should be WAVE (//TODO: keep as []byte?)
}

// WaveFmt describes the format of the sound-information in the data subchunks
type WaveFmt struct {
	Subchunk1ID    []byte // should contain "fmt"
	Subchunk1Size  int    // 16 for PCM
	AudioFormat    int    // PCM = 1 (Linear Quantization), if not 1, compression was used.
	NumChannels    int    // Mono 1, Stereo = 2, ..
	SampleRate     int    // 44100 for CD-Quality, etc..
	ByteRate       int    // SampleRate * NumChannels * BitsPerSample / 8
	BlockAlign     int    // NumChannels * BitsPerSample / 8 (number of bytes per sample)
	BitsPerSample  int    // 8 bits = 8, 16 bits = 16, .. :-)
	ExtraParamSize int    // if not PCM, can contain extra params
	ExtraParams    []byte // the actual extra params.
}

// WaveData contains the raw sound data
type WaveData struct {
	Subchunk2ID   []byte // Identifier of subchunk
	Subchunk2Size int    // size of raw sound data
	RawData       []byte // raw sound data itself
	Frames        []Frame
}

// NewWaveFmt can be used to generate a complete WaveFmt by calculating the remaining props
func NewWaveFmt(format, channels, samplerate, bitspersample int, extraparams []byte) WaveFmt {
	return WaveFmt{
		Subchunk1ID:    Format,
		Subchunk1Size:  16, // assume PCM for now //16
		AudioFormat:    format,
		NumChannels:    channels,
		SampleRate:     samplerate,
		ByteRate:       samplerate * channels * (bitspersample / 8.0),
		BlockAlign:     channels * (bitspersample / 8),
		BitsPerSample:  bitspersample,
		ExtraParamSize: len(extraparams),
		ExtraParams:    extraparams,
	}
}

// SetChannels changes the FMT to adapt to a new amount of channels
func (wfmt *WaveFmt) SetChannels(n uint) {
	wfmt.NumChannels = int(n)
	wfmt.ByteRate = (wfmt.SampleRate * wfmt.NumChannels * wfmt.BitsPerSample) / 8
	wfmt.BlockAlign = (wfmt.NumChannels * wfmt.BitsPerSample) / 8
}

// wavetable implementation

// FourierTable constructs a lookup table based on fourier addition with 'nharmns' harmonics
// If amps is provided, scales the harmonics by the provided amp
func FourierTable(nharms int, amps []float64, length int, phase float64) []float64 {
	table := make([]float64, length+2)
	phase *= tau

	for i := 0; i < nharms; i++ {
		for n := 0; n < len(table); n++ {
			amp := 1.0
			if i < len(amps) {
				amp = amps[i]
			}
			angle := float64(i+1) * (float64(n) * tau / float64(length))
			table[n] += (amp * math.Cos(angle+phase))
		}
	}
	return normalize(table)
}

// SawTable creates a sawtooth wavetable using Fourier addition
func SawTable(nharms, length int) []float64 {
	amps := make([]float64, nharms)
	for i := 0; i < len(amps); i++ {
		amps[i] = 1.0 / float64(i+1)
	}
	return FourierTable(nharms, amps, length, -0.25)
}

// SquareTable uses fourier addition to create a square waveform
func SquareTable(nharms, length int) []float64 {
	amps := make([]float64, nharms)
	for i := 0; i < len(amps); i += 2 {
		amps[i] = 1.0 / float64(i+1)
	}
	return FourierTable(nharms, amps, length, -0.25)
}

// TriangleTable uses fourier addition to create a triangle waveform
func TriangleTable(nharms, length int) []float64 {
	amps := make([]float64, nharms)
	for i := 0; i < nharms; i += 2 {
		amps[i] = 1.0 / (float64(i+1) * float64(i+1))
	}
	return FourierTable(nharms, amps, length, 0)
}

// Consts that appear in the .WAVE file format
var (
	ChunkID          = []byte{0x52, 0x49, 0x46, 0x46} // RIFF
	BigEndianChunkID = []byte{0x52, 0x49, 0x46, 0x58} // RIFX
	WaveID           = []byte{0x57, 0x41, 0x56, 0x45} // WAVE
	Format           = []byte{0x66, 0x6d, 0x74, 0x20} // FMT
	Subchunk2ID      = []byte{0x64, 0x61, 0x74, 0x61} // DATA
)

type intsToBytesFunc func(i int) []byte

var (
	// intsToBytesFm to map X-bit int to byte functions
	intsToBytesFm = map[int]intsToBytesFunc{
		16: int16ToBytes,
		32: int32ToBytes,
	}
)

// WriteFrames writes the slice to disk as a .wav file
// the WaveFmt metadata needs to be correct
// WaveData and WaveHeader are inferred from the samples however..
func WriteFrames(samples []Frame, wfmt WaveFmt, file string) error {
	return WriteWaveFile(samples, wfmt, file)
}

func WriteWaveFile(samples []Frame, wfmt WaveFmt, file string) error {
	f, err := os.Create(file)
	if err != nil {
		return err
	}
	defer f.Close()

	return WriteWaveToWriter(samples, wfmt, f)
}

func WriteWaveToWriter(samples []Frame, wfmt WaveFmt, writer io.Writer) error {
	wfb := fmtToBytes(wfmt)
	data, databits := framesToData(samples, wfmt)
	hdr := createHeader(data)

	_, err := writer.Write(hdr)
	if err != nil {
		return err
	}
	_, err = writer.Write(wfb)
	if err != nil {
		return err
	}
	_, err = writer.Write(databits)
	if err != nil {
		return err
	}

	return nil
}

func int16ToBytes(i int) []byte {
	b := make([]byte, 2)
	in := uint16(i)
	binary.LittleEndian.PutUint16(b, in)
	return b
}

func int32ToBytes(i int) []byte {
	b := make([]byte, 4)
	in := uint32(i)
	binary.LittleEndian.PutUint32(b, in)
	return b
}

func framesToData(frames []Frame, wfmt WaveFmt) (WaveData, []byte) {
	b := []byte{}
	raw := samplesToRawData(frames, wfmt)

	// We receive frames but have to store the size of the samples
	// The size of the samples is frames / channels..
	subchunksize := (len(frames) * wfmt.NumChannels * wfmt.BitsPerSample) / 8
	subBytes := int32ToBytes(subchunksize)

	// construct the data part..
	b = append(b, Subchunk2ID...)
	b = append(b, subBytes...)
	b = append(b, raw...)

	wd := WaveData{
		Subchunk2ID:   Subchunk2ID,
		Subchunk2Size: subchunksize,
		RawData:       raw,
		Frames:        frames,
	}
	return wd, b
}

func floatToBytes(f float64, nBytes int) []byte {
	bits := math.Float64bits(f * math.Pow10(-307) * math.Pow10(-1))
	bs := make([]byte, 8)
	binary.LittleEndian.PutUint64(bs, bits)
	// trim padding
	switch nBytes {
	case 2:
		return bs[:2]
	case 4:
		return bs[:4]
	}
	return bs
}

// Turn the samples into raw data...
func samplesToRawData(samples []Frame, props WaveFmt) []byte {
	raw := []byte{}
	for _, s := range samples {
		// the samples are scaled - rescale them?
		rescaled := rescaleFrame(s, props.BitsPerSample)
		bits := intsToBytesFm[props.BitsPerSample](rescaled)
		raw = append(raw, bits...)
	}
	return raw
}

// rescale frames back to the original values..
func rescaleFrame(s Frame, bits int) int {
	rescaled := float64(s) * float64(maxValues[bits])
	return int(rescaled)
}

func fmtToBytes(wfmt WaveFmt) []byte {
	b := []byte{}

	subchunksize := int32ToBytes(wfmt.Subchunk1Size)
	audioformat := int16ToBytes(wfmt.AudioFormat)
	numchans := int16ToBytes(wfmt.NumChannels)
	sr := int32ToBytes(wfmt.SampleRate)
	br := int32ToBytes(wfmt.ByteRate)
	blockalign := int16ToBytes(wfmt.BlockAlign)
	bitsPerSample := int16ToBytes(wfmt.BitsPerSample)

	b = append(b, wfmt.Subchunk1ID...)
	b = append(b, subchunksize...)
	b = append(b, audioformat...)
	b = append(b, numchans...)
	b = append(b, sr...)
	b = append(b, br...)
	b = append(b, blockalign...)
	b = append(b, bitsPerSample...)

	return b
}

// turn the sample to a valid header
func createHeader(wd WaveData) []byte {
	// write chunkID
	bits := []byte{}

	chunksize := 36 + wd.Subchunk2Size
	cb := int32ToBytes(chunksize)

	bits = append(bits, ChunkID...) // in theory switch on endianness..
	bits = append(bits, cb...)
	bits = append(bits, WaveID...)

	return bits
}

//////END CREDITS ----https://github.com/DylanMeeus/GoAudio/

func Float642Frame(floats []float64) []Frame {
	Frames := []Frame{}
	for i := 0; i < len(floats); i++ {
		Frames = append(Frames, Frame(floats[i]*math.Pow10(13)))
	}
	return Frames
}

func int64ToFloat64Slice(input []int64) []float64 {
	result := make([]float64, len(input))
	for i, v := range input {
		result[i] = float64(v)
	}
	return result
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

func FetchQubitsIBM(Qubitfilename string, NumberOfNewGeneratedFloats64 int64) ([][]QubitRI, int) {
	var q1RI QubitRI
	var q2RI QubitRI
	var Combo []QubitRI
	var Echoes [][]QubitRI
	var QuantumFloats []float64
	data, err := readFileIntoBytes(Qubitfilename)
	if err != nil {
		panic(err)
	}

	QuantumFloatsString := strings.Split(string(data), "\n")
	for i := 0; i < len(QuantumFloatsString); i++ {
		value, _ := strconv.ParseFloat(QuantumFloatsString[i], 32)
		QuantumFloats = append(QuantumFloats, value)
	}

	OriginalQubitsFloatsSpan := len(QuantumFloats)

	for i := NumberOfNewGeneratedFloats64; i < int64(len(QuantumFloatsString)); i = i * 16 {
		NumberOfNewGeneratedFloats64 = i
	}

	QuantumFloats = append(QuantumFloats, GenerateMoreValues(NumberOfNewGeneratedFloats64-int64(len(QuantumFloatsString)), QuantumFloats)...)

	Qte := int(len(QuantumFloats)/16) - 1 //TODO:REVIEW if this is still necessary. (-1)

	for i := 0; i < Qte; i += 16 {

		One := QuantumFloats[i]
		Two := QuantumFloats[i+1]
		Three := QuantumFloats[i+2]
		Four := QuantumFloats[i+3]
		Five := QuantumFloats[i+4]
		Six := QuantumFloats[i+5]
		Seven := QuantumFloats[i+6]
		Eight := QuantumFloats[i+7]
		Nine := QuantumFloats[i+8]
		Ten := QuantumFloats[i+9]
		Eleven := QuantumFloats[i+10]
		Twelve := QuantumFloats[i+11]
		Thirteen := QuantumFloats[i+12]
		Fourthteen := QuantumFloats[i+13]
		Fifthteen := QuantumFloats[i+14]
		Sixteen := QuantumFloats[i+15]

		q1RI = QubitRI{
			AlphaReal: One,
			AlphaImag: Two,
			OmegaReal: Three,
			OmegaImag: Four,
			BetaReal:  Five,
			BetaImag:  Six,
			ThetaReal: Seven,
			ThetaImag: Eight,
		}

		q2RI = QubitRI{
			AlphaReal: Nine,
			AlphaImag: Ten,
			OmegaReal: Eleven,
			OmegaImag: Twelve,
			BetaReal:  Thirteen,
			BetaImag:  Fourthteen,
			ThetaReal: Fifthteen,
			ThetaImag: Sixteen,
		}
		Combo = append(Combo, q1RI, q2RI)
		Echoes = append(Echoes, Combo)
		Combo = nil
	}
	return Echoes, OriginalQubitsFloatsSpan
}

func ReadJSON(JSONfilename string) [][]QubitRI {
	data, err := readFileIntoBytes(JSONfilename)
	var lstQubit [][]QubitRI
	if err != nil {
		panic(err)
	}

	err = json.Unmarshal(data, &lstQubit)
	if err != nil {
		panic(err)
	}
	return lstQubit
}

func WriteQubits(JSONfilename string, QubitsFilename string, NbToBeGeneratedATopFromMean int64) ([][]QubitRI, int) {

	Echoes, OriginalQubitsFloatsSpan := FetchQubitsIBM(QubitsFilename, NbToBeGeneratedATopFromMean)

	jsonData, err := json.Marshal(Echoes)
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

	return Echoes, OriginalQubitsFloatsSpan
}

func Float64sFromBytes(bytes []byte) []float64 {
	numFloats := len(bytes) / 8
	floats := make([]float64, numFloats)

	for i := 0; i < numFloats; i++ {
		start := i * 8
		end := (i + 1) * 8
		bits := binary.LittleEndian.Uint64(bytes[start:end])
		float := math.Float64frombits(bits)
		floats[i] = float
	}

	return floats
}

func float64ToBytes(data float64) []byte {
	var buf bytes.Buffer
	err := binary.Write(&buf, binary.LittleEndian, data)
	if err != nil {
		fmt.Println("binary.Write failed:", err)
	}
	return buf.Bytes()
}

func floats64ToBytes(data []float64) []byte {
	var buf bytes.Buffer
	for i := 0; i < len(data); i++ {
		err := binary.Write(&buf, binary.LittleEndian, data[i])
		if err != nil {
			fmt.Println("binary.Write failed:", err)
		}
	}
	return buf.Bytes()
}

func main() {
	var lstQubitsRI [][]QubitRI
	var w string
	var r string
	var i string
	var Realq1 []float64
	var Realq2 []float64
	var Imagq1 []float64
	var Imagq2 []float64
	var GenerateFromMean int64 = 16 * 16 * 16 * 16 * 16 * 16 // 16777216 max Qubits Floats64 (run time)
	var OriginalQubitsFloatsSpan int

	flag.StringVar(&w, "w", "", "Write Qubits")
	flag.StringVar(&i, "i", "", "IBM Qubits file")
	flag.StringVar(&r, "r", "", "Read Qubits")

	flag.Parse()

	if len(os.Args) < 3 {
		fmt.Println("Usage: Qubits -i <IBM_Qubits_file> -w <JSON_Out_Filename> /-r <JSON_In_Filename>")
		os.Exit(1)
	}

	if w != "" {
		JSON_OUT_FILENAME := w
		IBM_Qubits_File := i
		lstQubitsRI, OriginalQubitsFloatsSpan = WriteQubits(JSON_OUT_FILENAME, IBM_Qubits_File, GenerateFromMean)
	} else {
		if r != "" {
			JSON_OUT_FILENAME := r
			lstQubitsRI = ReadJSON(JSON_OUT_FILENAME)
		} else {
			fmt.Println("Usage: Qubit -w/-r <Qubits file>")
			os.Exit(1)
		}
	}

	for _, Qubits := range lstQubitsRI {
		//fmt.Println(Qubits)
		Realq1 = append(Realq1, Qubits[0].AlphaReal, Qubits[0].BetaReal, Qubits[0].OmegaReal, Qubits[0].ThetaReal)
		Realq2 = append(Realq2, Qubits[1].AlphaReal, Qubits[1].BetaReal, Qubits[1].OmegaReal, Qubits[1].ThetaReal)
		Imagq1 = append(Imagq1, Qubits[0].AlphaImag, Qubits[0].BetaImag, Qubits[0].OmegaImag, Qubits[0].ThetaImag)
		Imagq2 = append(Imagq2, Qubits[1].AlphaImag, Qubits[1].BetaImag, Qubits[1].OmegaImag, Qubits[1].ThetaImag)
		//fmt.Printf("%d:%v+%v\n\n", i, q1, q2) // DEBUG
	}

	/*
		fmt.Println(len(Realq1))
		fmt.Println(len(Realq2))
		fmt.Println(len(Imagq1))
		fmt.Println(len(Imagq2))
	*/

	var data []byte
	for i := 0; i < 4; i++ {
		switch i {
		case 0:
			data = floats64ToBytes(Realq1)
		case 1:
			data = floats64ToBytes(Realq2)
		case 2:
			data = floats64ToBytes(Imagq1)
		case 3:
			data = floats64ToBytes(Imagq2)
		}

		var addFirst uint8 = 0
		var addSecond uint8 = 0
		var addThird uint8 = 0
		var addFourth uint8 = 0
		var addFifth uint8 = 0
		var addSixth uint8 = 0
		var addSeventh uint8 = 0
		var addEight uint8 = 0
		var addNinth uint8 = 0
		var addTenth uint8 = 0
		var addEleventh uint8 = 0
		var addTwelfth uint8 = 0
		var addThirteenth uint8 = 0
		var addFourtheenth uint8 = 0
		var addFifteenth uint8 = 0
		var addSixteenth uint8 = 0
		var addSeventeenth uint8 = 0
		var addEighteenth uint8 = 0
		var addNineteenth uint8 = 0
		var addTwentieth uint8 = 0
		var addTwentyFirst uint8 = 0
		var PreviousIndex int = 0
		var TailRemoved = false
		for j := 0; j < len(data); j++ {
			var first uint8 = 1
			var second uint8 = 1
			var third uint8 = 1
			var fourth uint8 = 1
			var fifth uint8 = 1
			var sixth uint8 = 1
			var seventh uint8 = 1
			var eight uint8 = 1
			var nine uint8 = 1
			var ten uint8 = 1
			var eleven uint8 = 1
			var twelve uint8 = 1
			var thirteen uint8 = 1
			var fourtheen uint8 = 1
			var fifthteen uint8 = 1
			var sixteen uint8 = 1
			var seventeen uint8 = 1
			var eighteen uint8 = 1
			var nineteen uint8 = 1
			var twenty uint8 = 1
			var twentyOne uint8 = 1
			var AdjustementsDone = false
			if !TailRemoved {
				for i := 0; i < len(data); i++ {
					if i >= 7 {
						first = uint8(data[i-7])
						second = uint8(data[i-6])
						third = uint8(data[i-5])
						fourth = uint8(data[i-4])
						fifth = uint8(data[i-3])
						sixth = uint8(data[i-2])
						seventh = uint8(data[i-1])
						eight = uint8(data[i])

						if (first+second+third == 0) || (second+third+fourth == 0) || (third+fourth+fifth == 0) || (fourth+fifth+sixth == 0) || (fifth+sixth+seventh == 0) || (sixth+seventh+eight == 0) {
							i = i + 24
							first = uint8(data[i-3])
							second = uint8(data[i-2])
							third = uint8(data[i-1])
							if first+second+third == 0 {
								first = 1
								second = 1
								third = 1
								TailRemoved = true

								data = data[i:]                   // remove the trailing qubits information fron the previous DNA segment.
								f, err := os.Create("Qubits.DAT") // to build the reference map with bin2png
								if err != nil {
									panic(err)
								}
								f.Write(data)
								f.Close()
								break
							}
						}
					}
				}
			}

			if j%21 == 0 && j > 0 {
				first = uint8(data[j-21])
				second = uint8(data[j-20])
				third = uint8(data[j-19])
				fourth = uint8(data[j-18])
				fifth = uint8(data[j-17])
				sixth = uint8(data[j-16])
				seventh = uint8(data[j-15])
				eight = uint8(data[j-14])
				nine = uint8(data[j-13])
				ten = uint8(data[j-12])
				eleven = uint8(data[j-11])
				twelve = uint8(data[j-10])
				thirteen = uint8(data[j-9])
				fourtheen = uint8(data[j-8])
				fifthteen = uint8(data[j-7])
				sixteen = uint8(data[j-6])
				seventeen = uint8(data[j-5])
				eighteen = uint8(data[j-4])
				nineteen = uint8(data[j-3])
				twenty = uint8(data[j-2])
				twentyOne = uint8(data[j-1])
				addFirst += first
				addSecond += second
				addThird += third
				addFourth += fourth
				addFifth += fifth
				addSixth += sixth
				addSeventh += seventh
				addEight += eight
				addNinth += nine
				addTenth += ten
				addEleventh += eleven
				addTwelfth += twelve
				addThirteenth += thirteen
				addFourtheenth += fourtheen
				addFifteenth += fifthteen
				addSixteenth += sixteen
				addSeventeenth += seventeen
				addEighteenth += eighteen
				addNineteenth += nineteen
				addTwentieth += twenty
				addTwentyFirst += twentyOne

				if (seventh + eight + nine) == 0 {
					if first == 0 {
						if j/3-PreviousIndex == 56 {
							if twentyOne == 79 || twentyOne == 80 {
								fmt.Printf("j:%d First-->%d\n", j/3, addFirst)
								fmt.Printf("j:%d Second-->%d\n", j/3, addSecond)
								fmt.Printf("j:%d third-->%d\n", j/3, addThird)
								fmt.Printf("j:%d fourth-->%d\n", j/3, addFourth)
								fmt.Printf("j:%d fifth-->%d\n", j/3, addFifth)
								fmt.Printf("j:%d Sixth-->%d\n", j/3, addSixth)
								fmt.Printf("j:%d Seventh-->%d\n", j/3, addSeventh)
								fmt.Printf("j:%d Eight-->%d\n", j/3, addEight)
								fmt.Printf("j:%d Nineth-->%d\n", j/3, addNinth)
								fmt.Printf("j:%d Tenth-->%d\n", j/3, addTenth)
								fmt.Printf("j:%d Eleventh-->%d\n", j/3, addEleventh)
								fmt.Printf("j:%d Twelfth-->%d\n", j/3, addTwelfth)
								fmt.Printf("j:%d Thirtennth-->%d\n", j/3, addThirteenth)
								fmt.Printf("j:%d Fourtheenth-->%d\n", j/3, addFourtheenth)
								fmt.Printf("j:%d Fiftheenth-->%d\n", j/3, addFifteenth)
								fmt.Printf("j:%d Sixtheenth-->%d\n", j/3, addSixteenth)
								fmt.Printf("j:%d Seventeenth-->%d\n", j/3, addSeventeenth)
								fmt.Printf("j:%d Eighteenth-->%d\n", j/3, addEighteenth)
								fmt.Printf("j:%d NineTeenth-->%d\n", j/3, addNineteenth)
								fmt.Printf("j:%d Twentieth-->%d\n", j/3, addTwentieth)
								fmt.Printf("j:%d TwentyFirst-->%d\n", j/3, addTwentyFirst)
								fmt.Printf("TOTAL between 1 and 21 (exclusively) = %d\n\n", int(addSecond)+int(addThird)+int(addFourth)+int(addFifth)+int(addSixth)+int(addSeventh)+int(addEight)+int(addNinth)+int(addTenth)+int(addEleventh)+int(addTwelfth)+int(addThirteenth)+int(addFourtheenth)+int(addFifteenth)+int(addSixteenth)+int(addSeventeenth)+int(addEighteenth)+int(addNineteenth)+int(addTwentieth))
								AdjustementsDone = false
							}
						}
						PreviousIndex = j / 3
					}
				}
				addFirst = 0
				addSecond = 0
				addThird = 0
				addFourth = 0
				addFifth = 0
				addSixth = 0
				addSeventh = 0
				addEight = 0
				addNinth = 0
				addTenth = 0
				addEleventh = 0
				addTwelfth = 0
				addThirteenth = 0
				addFourtheenth = 0
				addFifteenth = 0
				addSixteenth = 0
				addSeventeenth = 0
				addEighteenth = 0
				addNineteenth = 0
				addTwentieth = 0
				addTwentyFirst = 0
				if (PreviousIndex*3+168 >= (int(OriginalQubitsFloatsSpan / 16))) && !AdjustementsDone {
					if (j+188 < len(data)) && (j > int(OriginalQubitsFloatsSpan/16)) {
						//Generating DNA Adjustments
						TenPoint4or5 := []int{79, 80}
						Choice := rand.Intn(2)
						if Choice == 0 {
							Choice++
						} else {
							Choice--
						}
						rand.NewSource(time.Now().UnixNano())
						data[PreviousIndex*3+168] = 0 //first
						data[PreviousIndex*3+175] = 0 //seventh
						data[PreviousIndex*3+176] = 0 //eight
						data[PreviousIndex*3+177] = 0 //nine
						data[PreviousIndex*3+188] = byte(TenPoint4or5[Choice])
						AdjustementsDone = true
					}
				}
			}
		}
		switch i {
		case 0:
			Realq1 = Float64sFromBytes(data)
		case 1:
			Realq2 = Float64sFromBytes(data)
		case 2:
			Imagq1 = Float64sFromBytes(data)
		case 3:
			Imagq2 = Float64sFromBytes(data)
		}
	}

	if len(Realq1)-len(Realq2) <= 0 {
		if len(Realq2)-len(Imagq1) <= 0 {
			if len(Imagq1)-len(Imagq2) <= 0 {
				fmt.Println("...Correct OK!")
			} else {
				chop1 := len(Imagq1) - len(Imagq2)
				Imagq1 = Imagq1[:len(Imagq1)-chop1]
				Realq2 = Realq2[:len(Realq1)-chop1]
				Realq1 = Realq1[:len(Realq1)-chop1]
			}
		} else {
			chop2 := len(Realq2) - len(Imagq1)
			Imagq1 = Imagq2[:len(Imagq1)-chop2]
			Realq2 = Realq1[:len(Realq2)-chop2]
			Realq1 = Realq2[:len(Realq1)-chop2]
		}
	} else {
		chop3 := len(Realq1) - len(Realq2)
		Imagq1 = Imagq2[:len(Imagq1)-chop3]
		Imagq2 = Realq1[:len(Imagq1)-chop3]
		Realq1 = Realq2[:len(Realq1)-chop3]
	}

	ProduceGenePools(Realq1, Realq2, Imagq1, Imagq2)

	WaveFmt := NewWaveFmt(1, 2, 48000, 16, nil)

	//WriteWaveFile(Float642Frame(cumSumQubits), WaveFmt, "Qubits_TimeFrame.wav")
	WriteWaveFile(Float642Frame(Realq1), WaveFmt, "Qubits_Realq1.wav")
	WriteWaveFile(Float642Frame(Realq2), WaveFmt, "Qubits_Realq2.wav")
	WriteWaveFile(Float642Frame(Imagq1), WaveFmt, "Qubits_Imagq1.wav")
	WriteWaveFile(Float642Frame(Imagq2), WaveFmt, "Qubits_Imagq2.wav")

	fmt.Println("Saved...Qubits_Realq1.wav")
	fmt.Println("Saved...Qubits_Realq1.wav")
	fmt.Println("Saved...Qubits_Imagq1.wav")
	fmt.Println("Saved...Qubits_Imagq2.wav")

	//fmt.Println("")

	//f.Write(float64ToBytes(Realq1))
	//f.Write(float64ToBytes(Realq2))
	//f.Write(float64ToBytes(Imagq1))
	//f.Write(float64ToBytes(Imagq2))

	//f.Close()
}

func ProduceGenePools(Realq1 []float64, Realq2 []float64, Imagq1 []float64, Imagq2 []float64) string {
	var Count int = 1
	filename := "Qubits"
	size := 70000
	// Create a FloatGenome from float slice
	SuperpositionGenome := NewFloatGenome(Realq1, Realq2, Imagq1, Imagq2)

	// Access information from FloatGenome
	fmt.Println("DNA Length:", len(SuperpositionGenome))
	//fmt.Println("Superposition DNA Sequence:", SuperpositionGenome) //DEBUG
	var j int = 0
	var SecondRun bool = false
	f, err := os.Create(filename + "_Split_" + strconv.Itoa(Count) + ".DNA")
	if err != nil {
		panic(err)
	}
	if len(SuperpositionGenome) < size {
		f.Write([]byte(SuperpositionGenome[:]))
		f.Close()
		fmt.Println("Saved Superposition DNA sequence to file", filename+"_Split_"+strconv.Itoa(Count)+".DNA")
	} else {
		for i := 0; i < len(SuperpositionGenome); i = i + size {
			if SecondRun {
				if i%size == 0 {
					f.Write([]byte(SuperpositionGenome[j:i]))
					f.Close()
					fmt.Println("Saved Superposition DNA sequence to file", filename+"_Split_"+strconv.Itoa(Count)+".DNA")
					Count++
					f, err = os.Create(filename + "_Split_" + strconv.Itoa(Count) + ".DNA")
					if err != nil {
						panic(err)
					}
				}
				if i+size >= len(SuperpositionGenome) {
					f.Write([]byte(SuperpositionGenome[i:]))
					f.Close()
					fmt.Println("Saved Superposition DNA sequence to file", filename+"_Split_"+strconv.Itoa(Count)+".DNA")
					break
				} else {
					j = i
				}
			} else {
				SecondRun = true
			}
		}
	}
	return SuperpositionGenome
}

// FloatToNucleotide encodes a float64 value to a nucleotide.
func FloatToNucleotides(combined float64) []alphabet.Letter {
	var converted []byte
	var nucleotides []alphabet.Letter
	converted = append(converted, float64ToBytes(combined)...)
	//converted = append(converted, byte(q2*math.Pow10(10)))
	nucleotides = append(nucleotides, alphabet.BytesToLetters(converted)...)
	//fmt.Println(nucleotides)
	return nucleotides
}

// NewFloatGenome creates a FloatGenome from a []float64.
func NewFloatGenome(Realq1 []float64, Realq2 []float64, Imagq1 []float64, Imagq2 []float64) string {
	var basePairs = map[string]string{
		"a": "t",
		"t": "a",
		"c": "g",
		"g": "c",
	}
	var bPrevious byte
	var bNow byte
	var bNext byte
	var bCatchPrevious []byte
	var bCatchNow []byte
	var bCatchNext []byte
	var bpWork [][]byte
	var bpFinal string
	var bpBuilder string
	var IncrementalScan int
	var Initialized bool = false
	var A bool = false
	var C bool = false
	var T bool = false
	var G bool = false
	var Regress bool = false
	var nucleotide1 []alphabet.Letter
	var nucleotide2 []alphabet.Letter

	for i := 0; i < len(Realq1); i++ {
		nucleotide1 = FloatToNucleotides(Realq1[i] + Realq2[i])
		nucleotide2 = FloatToNucleotides(Imagq1[i] + Imagq2[i])
		for _, nucleotide := range nucleotide1 {
			bpBuilder += string(alphabet.DNA.Letter((int(^nucleotide) % 4)))
		}
		for _, nucleotide := range nucleotide2 {
			bpBuilder += string(alphabet.DNA.Letter((int(^nucleotide) % 4)))
		}
	}
	fmt.Println(bpBuilder)
	for i := 0; i < len(bpBuilder); i++ {
		if i > 0 {
			bPrevious = bpBuilder[i-1]
			bNow = bpBuilder[i]
			if i+1 < len(bpBuilder) {
				bNext = bpBuilder[i+1]
				bCatchPrevious = append(bCatchPrevious, bPrevious)
				bCatchNow = append(bCatchNow, bNow)
				bCatchNext = append(bCatchNext, bNext)
				if !Initialized {
					if string(bCatchPrevious) == "a" {
						A = true
					}
					if string(bCatchPrevious) == "c" {
						C = true
					}
					if string(bCatchPrevious) == "t" {
						T = true
					}
					if string(bCatchPrevious) == "g" {
						G = true
					}
					if string(bCatchNow) == "a" {
						A = true
					}
					if string(bCatchNow) == "c" {
						C = true
					}
					if string(bCatchNow) == "t" {
						T = true
					}
					if string(bCatchNow) == "g" {
						G = true
					}
					if string(bCatchNext) == "a" {
						A = true
					}
					if string(bCatchNext) == "c" {
						C = true
					}
					if string(bCatchNext) == "t" {
						T = true
					}
					if string(bCatchNext) == "g" {
						G = true
					}
					bpWork = append(bpWork, bCatchPrevious)
					bpWork = append(bpWork, bCatchNow)
					bpWork = append(bpWork, bCatchNext)
					bCatchPrevious = nil
					bCatchNow = nil
					bCatchNext = nil
					if A && C && T && G {
						Initialized = true
					}
					bCatchPrevious = nil
					bCatchNow = nil
					bCatchNext = nil
				} else {
					bpWork = append(bpWork, bCatchPrevious)
					bpWork = append(bpWork, bCatchNow)
					bpWork = append(bpWork, bCatchNext)
					bCatchPrevious = nil
					bCatchNow = nil
					bCatchNext = nil
				}
				i = i + 2
			} else {
				bCatchPrevious = append(bCatchPrevious, bPrevious)
				bCatchNow = append(bCatchNow, bNow)
				bpWork = append(bpWork, bCatchPrevious)
				bpWork = append(bpWork, bCatchNow)
				//fmt.Println(bpWork) //TEST OK.
				var bCurrent []byte
				for j := 0; j < len(bpWork); j++ { //regression
					var Count int = 0
					for k := 0; k < len(bpWork); k++ {
						//fmt.Printf("%s<-->%s\n", basePairs[string(bCurrent)], string(bpWork[k]))
						if (basePairs[string(bCurrent)] == string(bpWork[k]) && Count == 0 && (j < k)) || (Regress && (basePairs[string(bCurrent)] == string(bpWork[k]) && Count == 0)) {
							//fmt.Printf("%s<-->%s\n", basePairs[string(bCurrent)], string(bpWork[k]))
							//os.Exit(1)
							bpFinal += string(bpWork[k]) + string(bpWork[j])
							bpWork[k] = nil
							bCurrent = bpWork[j]
							if bCurrent == nil {
								IncrementalScan = 0
								for {
									IncrementalScan++
									if bpWork[IncrementalScan] != nil {
										break
									}
								}
								bCurrent = bpWork[IncrementalScan]
							}
							if Regress {
								Regress = false
							}
							break
						} else {
							if string(bpWork[j]) == string(bpWork[k]) && (j < k) && bpWork[j] != nil {
								if j < k {
									bCurrent = bpWork[j]
									bpFinal += string(bpWork[k])
									bpWork[k] = nil
									Count++
									//fmt.Println(bpWork)
									//time.Sleep(time.Second)
								}
							} else {
								// we have a continuity in bases
								if basePairs[string(bCurrent)] == string(bpWork[k]) && bCurrent != nil && (j < k) {
									bpFinal += string(bpWork[k])
									bpWork[k] = nil
									Count--
									if Count == 0 {
										Regress = true
									}
								}
								if Count == 0 {
									if j == 0 {
										bpFinal += string(bpWork[j])
										bCurrent = bpWork[j]
										bpWork[j] = nil
										break
									}
									if Regress {
										//TODO:check condition
										IncrementalScan = 0
										for {
											IncrementalScan++
											if bpWork[IncrementalScan] != nil {
												break
											}
										}
										bCurrent = bpWork[IncrementalScan]
										k = IncrementalScan
										bpFinal += string(bpWork[k])
										bpWork[k] = nil
									}
								}
							}
						}
					}
					IncrementalScan = 0
					for {
						IncrementalScan++
						if IncrementalScan < len(bpWork) {
							if bpWork[IncrementalScan] != nil {
								break
							}
						} else {
							break
						}
					}
					if IncrementalScan < len(bpWork) {
						bCurrent = bpWork[IncrementalScan]
						bpFinal += string(bpWork[IncrementalScan])
						bpWork[IncrementalScan] = nil
					} else {
						break
					}
				}
			}
			//fmt.Println(bpWork)
		}
	}
	for _, base := range bpWork {
		if base != nil {
			bpFinal += string(base) //Get more Qubits from IBM Quantum.
		}
	}

	return bpFinal
}

func calculateMean(values []float64) float64 {
	var sum float64 = 0.0
	for _, value := range values {
		sum += value
	}
	return sum / float64(len(values))
}

func calculateSumOfSquareDeviations(values []float64, mean float64) float64 {
	var sum float64 = 0.0
	for _, value := range values {
		sum += (value - mean) * (value - mean)
	}
	return sum
}

func calculateCovariance(x []float64, y []float64, meanX float64, meanY float64) float64 {
	var sum float64 = 0.0
	for i := 0; i < len(x); i++ {
		sum += (x[i] - meanX) * (y[i] - meanY)
	}
	return sum / float64(len(x))
}

func calculateRegressionCoefficients(SSDX float64, SSDY float64, covariance float64, lenData float64, meanX float64, meanY float64) (float64, float64) {
	b1 := covariance / SSDX
	b0 := meanY - b1*meanX
	return b0, b1
}

func InterpretNextY(b0 float64, b1 float64, x float64) float64 {
	return b0 + (b1 * x)
}

func GenerateMoreValues(NumberOfNewValues int64, yValues []float64) []float64 {
	var xValues []float64
	var NewyValues []float64
	meanY := calculateMean(yValues)
	meanX := 1.0 //1 more cycle of throwing a Qubit. 1 more assumption.

	for i := 1; i <= len(yValues); i++ {
		xValues = append(xValues, float64(i))
	}

	SSDX := calculateSumOfSquareDeviations(xValues, meanX)
	SSDY := calculateSumOfSquareDeviations(yValues, meanY)

	covariance := calculateCovariance(xValues, yValues, meanX, meanY)

	b0, b1 := calculateRegressionCoefficients(SSDX, SSDY, covariance, float64(len(xValues)), meanX, meanY)

	for x := xValues[len(xValues)-1] + meanX; x <= float64(NumberOfNewValues); x++ {
		NewyValues = append(NewyValues, InterpretNextY(b0, b1, x))
	}
	return NewyValues
}
