package main

import (
	"bytes"
	"crypto/rand"
	"encoding/binary"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"image/color"
	"io"
	"math"
	"math/cmplx"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/shopspring/decimal"
	"gonum.org/v1/gonum/stat"
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
	Timestamp int64
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

func Tick(Q1A complex128, Q2A complex128, Q1O complex128, Q2O complex128, Q1B complex128, Q2B complex128, Q1T complex128, Q2T complex128, Moment int, magnitude chan chan bool, phase chan chan bool) {

	var mag bool = false
	var pha bool = false
	magnitude2 := make(chan bool, 1)
	phase2 := make(chan bool, 1)
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
		magnitude2 <- false
		phase2 <- false
	} else {
		if !mag && pha {
			magnitude2 <- false
			phase2 <- false
		}
		if mag && !pha {
			magnitude2 <- false
			phase2 <- false
		}
		if mag && pha {
			magnitude2 <- true
			phase2 <- true
		}
	}
	magnitude <- magnitude2
	phase <- phase2
}

func CollapsingQubits(wg *sync.WaitGroup, magnitude chan chan bool, phase chan chan bool, NumberCpu int, Originalq1 *[]Qubit, Originalq2 *[]Qubit) {

	for i := 0; i < NumberCpu; i++ { // /4 if using the 4 phases
		wg.Add(1)
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

		*Originalq1 = append(*Originalq1, Qubit{Alpha: cmplx.Sqrt(normalizeComplex(complex(A1, A1i))),
			Omega: cmplx.Sqrt(normalizeComplex(complex(O1, O1i))),
			Beta:  cmplx.Sqrt(normalizeComplex(complex(B1, B1i))),
			Theta: cmplx.Sqrt(normalizeComplex(complex(T1, T1i)))})

		*Originalq2 = append(*Originalq2, Qubit{Alpha: cmplx.Sqrt(normalizeComplex(complex(A2, A2i))),
			Omega: cmplx.Sqrt(normalizeComplex(complex(O2, O2i))),
			Beta:  cmplx.Sqrt(normalizeComplex(complex(B2, B2i))),
			Theta: cmplx.Sqrt(normalizeComplex(complex(T2, T2i)))})

		q1 := Qubit{Alpha: cmplx.Sqrt(normalizeComplex(complex(A1, A1i))),
			Omega: cmplx.Sqrt(normalizeComplex(complex(O1, O1i))),
			Beta:  cmplx.Sqrt(normalizeComplex(complex(B1, B1i))),
			Theta: cmplx.Sqrt(normalizeComplex(complex(T1, T1i)))}

		q2 := Qubit{Alpha: cmplx.Sqrt(normalizeComplex(complex(A2, A2i))),
			Omega: cmplx.Sqrt(normalizeComplex(complex(O2, O2i))),
			Beta:  cmplx.Sqrt(normalizeComplex(complex(B2, B2i))),
			Theta: cmplx.Sqrt(normalizeComplex(complex(T2, T2i)))}

		Q1A, Q2A := cmplx.Pow(q1.Alpha, 2), cmplx.Pow(q2.Alpha, 2)
		Q1O, Q2O := cmplx.Pow(q1.Omega, 2), cmplx.Pow(q2.Omega, 2)
		Q1B, Q2B := cmplx.Pow(q1.Beta, 2), cmplx.Pow(q2.Beta, 2)
		Q1T, Q2T := cmplx.Pow(q1.Theta, 2), cmplx.Pow(q2.Theta, 2)

		Q1A, Q2A = normalizeComplex(Q1A), normalizeComplex(Q2A)
		Q1O, Q2O = normalizeComplex(Q1O), normalizeComplex(Q2O)
		Q1B, Q2B = normalizeComplex(Q1B), normalizeComplex(Q2B)
		Q1T, Q2T = normalizeComplex(Q1T), normalizeComplex(Q2T)

		//go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 1, magnitude1, phase1)
		go Tick(Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 2, magnitude, phase)
		//go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 3, magnitude3, phase3)
		//go Tick(&wg, Q1A, Q2A, Q1O, Q2O, Q1B, Q2B, Q1T, Q2T, 4, magnitude4, phase4)

		wg.Done()
	}

}

func secureRandomFloat64() float64 {
	var randomBytes [8]byte
	var randomFloat float64
	for {
		_, err := rand.Read(randomBytes[:])
		if err != nil {
			panic(err)
		}

		randomFloat = math.Float64frombits(binary.BigEndian.Uint64(randomBytes[:]))
		if !math.IsNaN(randomFloat) {
			break
		}
	}

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

func WriteQubits(JSONfilename string, WriteBuffer int) [][]QubitRI {
	var wg sync.WaitGroup
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	var listQubit [][]QubitRI
	var q1 []Qubit
	var q2 []Qubit
	var Duration int64
	var QubitsRICombo []QubitRI
	var Then int64 = time.Now().UnixNano()
	var NumberCpu int = runtime.NumCPU()
	var CollectedBits int64 = 0
	var TotalCollectedQubits int64
	magnitude := make(chan chan bool, NumberCpu*WriteBuffer)
	phase := make(chan chan bool, NumberCpu*WriteBuffer)

	for {
		if m.Sys-m.Alloc < 512000 { // 500MB threshold
			break
		}

		for i := 0; i < WriteBuffer; i++ {
			wg.Add(1)
			CollapsingQubits(&wg, magnitude, phase, NumberCpu, &q1, &q2)
			wg.Done()
		}

		wg.Wait()

		for j := 0; j < NumberCpu*WriteBuffer; j++ {
			magnitude2 := <-magnitude
			phase2 := <-phase
			if <-magnitude2 && <-phase2 {
				Duration = time.Now().UnixNano() - Then //the time it takes to bring another collapse forward.
				q3 := RebuildFromTick(2, q1[j])

				q3RI := QubitRI{
					AlphaReal: real(q3.Alpha),
					AlphaImag: imag(q3.Alpha),
					OmegaReal: real(q3.Omega),
					OmegaImag: imag(q3.Omega),
					BetaReal:  real(q3.Beta),
					BetaImag:  imag(q3.Beta),
					ThetaReal: real(q3.Theta),
					ThetaImag: imag(q3.Theta),
					Tick:      2,
					Timestamp: Duration,
				}
				q2RI := QubitRI{
					AlphaReal: real(q2[j].Alpha),
					AlphaImag: imag(q2[j].Alpha),
					OmegaReal: real(q2[j].Omega),
					OmegaImag: imag(q2[j].Omega),
					BetaReal:  real(q2[j].Beta),
					BetaImag:  imag(q2[j].Beta),
					ThetaReal: real(q2[j].Theta),
					ThetaImag: imag(q2[j].Theta),
					Tick:      2,
					Timestamp: Duration,
				}

				QubitsRICombo = append(QubitsRICombo, q3RI, q2RI)
				//fmt.Println(QubitsRICombo) //DEBUG PRINT
				listQubit = append(listQubit, QubitsRICombo)
				QubitsRICombo = nil
				fmt.Printf("\rFound Qubits = %d", TotalCollectedQubits)
				TotalCollectedQubits++
				CollectedBits++

				if CollectedBits == int64(WriteBuffer) {
					jsonData, err := json.Marshal(listQubit)
					if err != nil {
						panic(err)
					} else {
						fmt.Printf("...Saving ... %d Qubits to file %s\n", TotalCollectedQubits, JSONfilename)
					}

					file, err := os.Create(JSONfilename)
					if err != nil {
						panic(err)
					}

					_, err = file.Write(jsonData)
					if err != nil {
						panic(err)
					}
					CollectedBits = 0
					file.Close()
				}
			}

		}
		q1 = nil
		q2 = nil
	}
	close(magnitude)
	close(phase)
	return listQubit
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

func cumulativeSum(numbers []float64) []float64 {
	var cumSum []float64
	var sum float64 = 0.0

	for _, num := range numbers {
		sum += num
		cumSum = append(cumSum, sum)
	}

	return cumSum
}

func Graph(cumSumQubits []float64, adjustedTicks []float64) {

	// Calculate the cumulative sum
	//cumSumTick := cumulativeSum(lstTick)
	//cumSumQubits := cumulativeSum(lstQubits)

	// Create a new plot
	p := plot.New()

	ptsLineQubits := make(plotter.XYs, len(cumSumQubits))
	for i, num := range cumSumQubits {
		ptsLineQubits[i].X = float64(i)
		ptsLineQubits[i].Y = num
	}

	ptsLineAdjustedTick := make(plotter.XYs, len(adjustedTicks))
	for i, num := range adjustedTicks {
		ptsLineAdjustedTick[i].X = float64(i)
		ptsLineAdjustedTick[i].Y = num
	}

	// Create lines
	lineQubits, err := plotter.NewLine(ptsLineQubits)
	if err != nil {
		panic(err)
	}

	lineTicks, err := plotter.NewLine(ptsLineAdjustedTick)
	if err != nil {
		panic(err)
	}

	// Set color for Line 1
	lineQubits.Color = color.RGBA{R: 255, G: 0, B: 0, A: 255} // Red

	// Set color for Line 2
	lineTicks.Color = color.RGBA{R: 0, G: 255, B: 0, A: 255} // Green

	// Add the scatter and line plots to the plot
	p.Add(lineQubits, lineTicks) //, lineTicks)

	// Set axis labels
	p.X.Label.Text = "Time"
	p.Y.Label.Text = "Qubits"

	// Save the plot to a PNG file
	fmt.Println("Building the plot data...Please wait!")
	if err := p.Save(12*vg.Inch, 8*vg.Inch, "Qubits_cumulative_sum_plot.png"); err != nil {
		panic(err)
	}
	fmt.Printf("Saved Qubits graph... Qubits_cumulative_sum_plot.png\n")
}

func logValues(values []float64) []float64 {
	loggedValues := make([]float64, len(values))
	for i, v := range values {
		loggedValues[i] = math.Log(v)
	}
	return loggedValues
}

func scaleValues(values []float64, scale float64) []float64 {
	scaledValues := make([]float64, len(values))
	for i, v := range values {
		scaledValues[i] = v * scale
	}
	return scaledValues
}

func int64ToFloat64Slice(input []int64) []float64 {
	result := make([]float64, len(input))
	for i, v := range input {
		result[i] = float64(v)
	}
	return result
}

func Regress(listQubits []float64, listTick []float64, SkipCumulativeSums bool) ([]float64, []float64) {

	var cumSumTick []float64 = listTick
	var cumSumQubits []float64 = listQubits
	var a1 float64
	var a2 float64
	//var b1 float64
	//var b2 float64

	if !SkipCumulativeSums {
		// Calculate the cumulative sum
		cumSumTick = cumulativeSum(listTick)
		cumSumQubits = cumulativeSum(listQubits)
	}

	var x []float64
	var yQubits []float64
	var yTicks []float64

	for i := 0; i < len(listQubits); i++ {
		x = append(x, float64(i))
	}
	for i := 0; i < len(listQubits); i++ {
		yQubits = append(yQubits, cumSumQubits[i])
	}
	for i := 0; i < len(listTick); i++ {
		yTicks = append(yTicks, cumSumTick[i])
	}

	xQubits := x
	xTicks := x

	// Perform linear regression for curve 1
	a1, _ = stat.LinearRegression(xQubits, logValues(yQubits), nil, false)

	// Perform linear regression for curve 2
	a2, _ = stat.LinearRegression(xTicks, logValues(yTicks), nil, false)

	// Calculate the difference between coefficients
	aDiff := math.Exp(a1 - a2)

	//fmt.Printf("Scale Factor Difference: %.2f\n", aDiff)
	adjustedTicks := scaleValues(yTicks, aDiff)

	return cumSumQubits, adjustedTicks
}

func Float642Frame(floats []float64) []Frame {
	Frames := []Frame{}
	for i := 0; i < len(floats); i++ {
		Frames = append(Frames, Frame(floats[i]))
	}
	return Frames
}

func findRepeatingPattern(nums []int64) []int64 {
	// Initialize the tortoise and hare pointers
	tortoise, hare := nums[0], nums[0]

	// Phase 1: Find the meeting point of tortoise and hare
	for {
		tortoise = nums[tortoise]
		hare = nums[nums[hare]]

		if tortoise == hare {
			break
		}
	}

	// Phase 2: Find the entrance to the cycle (repeating pattern)
	tortoise = nums[0]
	for tortoise != hare {
		tortoise = nums[tortoise]
		hare = nums[hare]
	}

	// The hare and tortoise are now at the entrance of the cycle
	return []int64{tortoise}
}

func findRepeatingSeries(lst []int64) []int64 {
	seen := make(map[int64]bool)
	repeatingSeries := make(map[int64]bool)

	for _, num := range lst {
		if seen[num] {
			repeatingSeries[num] = true
		} else {
			seen[num] = true
		}
	}

	var result []int64
	for num := range repeatingSeries {
		result = append(result, num)
	}

	return result
}

func QubitsToFloat64(listQubits [][]Qubit) ([]float64, []float64, []float64, []float64) {
	var listQ1Real []float64
	var listQ2Real []float64
	var listQ1Imag []float64
	var listQ2Imag []float64

	for _, Qubits := range listQubits {
		q1 := Qubits[0]
		q2 := Qubits[1]
		listQ1Real = append(listQ1Real, real(q1.Alpha), real(q1.Beta), real(q1.Omega), real(q1.Theta))
		listQ2Real = append(listQ2Real, real(q2.Alpha), real(q2.Beta), real(q2.Omega), real(q2.Theta))
		listQ1Imag = append(listQ1Imag, imag(q1.Alpha), imag(q1.Beta), imag(q1.Omega), imag(q1.Theta))
		listQ2Imag = append(listQ2Imag, imag(q2.Alpha), imag(q2.Beta), imag(q2.Omega), imag(q2.Theta))

	}

	return listQ1Real, listQ1Imag, listQ2Real, listQ2Imag
}

func main() {

	var WriteBufferThreshold int = 1000
	var lstQubit [][]QubitRI
	var listQubits [][]Qubit
	var RepeatQubits [][]Qubit
	var listNanoseconds []int64
	var listTick []int64
	var q1 Qubit
	var q2 Qubit

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
		lstQubit = WriteQubits(JSON_OUT_FILENAME, WriteBufferThreshold)
	} else {
		if r != "" {
			JSON_OUT_FILENAME := r
			lstQubit = ReadJSON(JSON_OUT_FILENAME)
		} else {
			fmt.Println("Usage: Qubit -w/-r <Qubits file>")
			os.Exit(1)
		}
	}

	var Foundq1 bool
	var Foundq2 bool
	for j, Qubits := range lstQubit {
		var Combination []Qubit
		var Repeat []Qubit
		q1 = Qubit{Alpha: complex(Qubits[0].AlphaReal, Qubits[0].AlphaImag),
			Beta:  complex(Qubits[0].BetaReal, Qubits[0].BetaImag),
			Omega: complex(Qubits[0].OmegaReal, Qubits[0].OmegaImag),
			Theta: complex(Qubits[0].ThetaReal, Qubits[0].ThetaImag),
		}

		q2 = Qubit{Alpha: complex(Qubits[1].AlphaReal, Qubits[1].AlphaImag),
			Beta:  complex(Qubits[1].BetaReal, Qubits[1].BetaImag),
			Omega: complex(Qubits[1].OmegaReal, Qubits[1].OmegaImag),
			Theta: complex(Qubits[1].ThetaReal, Qubits[1].ThetaImag),
		}

		listTick = append(listTick, int64(Qubits[0].Tick))
		listNanoseconds = append(listNanoseconds, Qubits[0].Timestamp)

		for i := j; i < len(lstQubit); i++ {
			if i+1 < len(lstQubit) {
				if q1 == (Qubit{Alpha: complex(lstQubit[i+1][0].AlphaReal, lstQubit[i+1][0].AlphaImag),
					Beta:  complex(lstQubit[i+1][0].BetaReal, lstQubit[i+1][0].BetaImag),
					Omega: complex(lstQubit[i+1][0].OmegaReal, lstQubit[i+1][0].OmegaImag),
					Theta: complex(lstQubit[i+1][0].ThetaReal, lstQubit[i+1][0].ThetaImag),
				}) {
					//fmt.Print("Found q1...")
					Foundq1 = true
					for k := i; k < len(lstQubit); k++ {
						if k+1 < len(lstQubit) {
							{
								if q2 == (Qubit{Alpha: complex(lstQubit[k+1][1].AlphaReal, lstQubit[k+1][1].AlphaImag),
									Beta:  complex(lstQubit[k+1][1].BetaReal, lstQubit[k+1][1].BetaImag),
									Omega: complex(lstQubit[k+1][1].OmegaReal, lstQubit[k+1][1].OmegaImag),
									Theta: complex(lstQubit[k+1][1].ThetaReal, lstQubit[k+1][1].ThetaImag),
								}) {
									//fmt.Println("Found q2!-->", i)
									Foundq2 = true
									break
								}
							}
						}
					}

					if Foundq1 && Foundq2 {
						Repeat = append(Repeat, q1, q2)
						RepeatQubits = append(RepeatQubits, Repeat)
						Foundq1 = false
						Foundq2 = false
					} else {
						if !Foundq1 && !Foundq2 {
							Combination = append(Combination, q1, q2)
							listQubits = append(listQubits, Combination)
						} else {
							Foundq1 = false
							Foundq2 = false
						}
					}
				}
			}
		}

	}
	fmt.Println("-==Repeated Qubits==-")
	for i, Repeat := range RepeatQubits {
		fmt.Printf("%d:%v+%v\n\n", i+1, Repeat[0], Repeat[1])
	}
	fmt.Println("-==Unique Qubits==-")
	for i, Combo := range listQubits {
		fmt.Printf("%d::%v+%v\n\n", i+1, Combo[0], Combo[1])
	}
	fmt.Println("Repeated Qubits harvested = ", len(RepeatQubits))
	fmt.Println("Unique Qubits harvested = ", len(listQubits))
	//repeatingTickPattern := findRepeatingPattern(listTick)
	//fmt.Printf("Repeating Tick pattern:%d\n", repeatingTickPattern)

	//repeatingTickSeries := findRepeatingSeries(listTick)
	//fmt.Printf("Repeating Tick series:%d\n", repeatingTickSeries)

	listTickFloat64 := int64ToFloat64Slice(listTick)
	listNanosecondsFloat64 := int64ToFloat64Slice(listNanoseconds)

	cumSumQubits, adjustedTicks := Regress(listNanosecondsFloat64, listTickFloat64, true)
	//CumSumQubits, adjustedTicks = Regress(CumSumQubits, adjustedTicks, false)

	/////TEST IT\\\\\\
	//var formula float64 = 1 * 0.05 + //math.Pow(math.Pi, math.Pi*math.Phi*2)

	//listQ1Real, listQ1Imag, listQ2Real, listQ2Imag := QubitsToFloat64(listQubits)

	//for i := 0; i < len(listQubits); i++ {
	//	fmt.Printf("%.7fi%.7f + %.7fi%.7f = %v\n", listQ1Real[i], listQ1Imag[i], listQ2Real[i], listQ2Imag[i], AddComplex(complex(listQ1Real[i], listQ1Imag[i]), complex(listQ2Real[i], listQ2Imag[i])))
	//}

	Graph(cumSumQubits, adjustedTicks) //Graph(CumSumQubits, adjustedTicks)

	WaveFmt := NewWaveFmt(1, 2, 44100, 16, nil)
	WriteWaveFile(Float642Frame(cumSumQubits), WaveFmt, "Qubits.wav")
	fmt.Println("Saved...Qubits.wav.")

	//for i := 1; i < len(CumSumQubits); i++ {
	//		fmt.Printf("Qubit-->%2.f<-->%2.f<-Tick\n", CumSumQubits[i], adjustedTicks[i])
	//	time.Sleep(time.Second)
	//}
	//if float64(cumSumQubits[1]) == (float64(cumSumTick[1])*interceptDiff + slopeDiff) {
	//	break
	//} else {
	//	fmt.Printf("Qubit-->%2.f<-->%2.f<-Tick\n", float64(cumSumQubits[1]), (float64(cumSumTick[1])))
	//}
}
