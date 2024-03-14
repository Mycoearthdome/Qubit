package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"slices"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
)

func Launch(file string, ZeroDNA []byte) {
	var DNASequence string
	var OutFile string
	var SpectrumKeys []float64
	fdata, err := os.ReadFile(file)
	if err != nil {
		panic(err)
	}
	data := strings.NewReader(string(fdata))
	template := linear.NewSeq("", nil, alphabet.DNAredundant)
	r := fasta.NewReader(data, template)

	sc := seqio.NewScanner(r)
	for sc.Next() {
		s := sc.Seq().(*linear.Seq)
		DNASequence += strings.ToLower(s.Seq.String())
	}

	OutFile = strings.Split(file, ".")[0] + "_Freq_Matched_DNA.Zero"

	if len(DNASequence) > len(ZeroDNA) {

		for i := range ZeroDNA { //Adjusting the loop.
			if ZeroDNA[i] == ZeroDNA[len(ZeroDNA)-1] {
				ZeroDNA = ZeroDNA[i : len(ZeroDNA)-1]
				break
			}
		}

		Aligned := false
		AlignedIndex := 0
		Count := 0
		Speed_of_Light := 299792458
		NucleotideSize := 0.33 / 1e9 //nm to meters
		var Frequencies []float64

		nucleotideIndex := 0

		f, err := os.Create(OutFile)
		if err != nil {
			panic(err)
		}

		for _, nucleotide := range DNASequence {
			if !Aligned {
				if string(nucleotide) != "n" {
					for i := 0; i < len(ZeroDNA); i++ {
						if byte(nucleotide) == ZeroDNA[i] { //aligned to DNA Zero
							Aligned = true
							AlignedIndex = i
						}
					}
				}
			} else {
				if ((nucleotideIndex + 1) < (len(DNASequence) - 1)) && (AlignedIndex < len(ZeroDNA)-1) {
					if byte(nucleotide) == ZeroDNA[AlignedIndex] && DNASequence[nucleotideIndex+1] == ZeroDNA[AlignedIndex+1] {
						WaveLength := float64(Count) * NucleotideSize
						Frequencies = append(Frequencies, float64(Speed_of_Light)/WaveLength)
						//f.write("%d=|%s-%s|" % (Count, nucleotide, DNAsequence[nucleotideIndex+1]))
						if AlignedIndex+1 == len(ZeroDNA)-1 { //loop it back
							AlignedIndex = 0
						} else {
							AlignedIndex = AlignedIndex + 1
						}
						Count = 0
					}
				} else {
					if AlignedIndex == len(ZeroDNA)-1 {
						AlignedIndex = 0
					}
				}
			}
			//f.write("%d=|%s|--------|%s|" % (Count, nucleotide, DNAsequence[nucleotideIndex+1]))
			Count = Count + 1
			nucleotideIndex = nucleotideIndex + 1
		}

		Spectrum := make(map[float64]int)
		for _, Frequency := range Frequencies {
			_, ok := Spectrum[Frequency]
			if !ok {
				Spectrum[Frequency] = 0
			} else {
				Spectrum[Frequency] = Spectrum[Frequency] + 1
			}
		}

		for keys := range Spectrum {
			SpectrumKeys = append(SpectrumKeys, keys)
		}

		slices.Sort(SpectrumKeys)

		for _, Frequency := range SpectrumKeys {
			outputString := fmt.Sprintf("Frequency=%f Count=%d\n", Frequency, Spectrum[Frequency])
			f.Write([]byte(outputString))
		}
		f.Close()
		fmt.Printf("Saved...%s\n", OutFile)
	}
}

func WalkMatch(root, pattern string) ([]string, error) {
	var matches []string
	err := filepath.WalkDir(root, func(path string, d os.DirEntry, err error) error { //info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() {
			return nil
		}
		if matched, err := filepath.Match(pattern, filepath.Base(path)); err != nil {
			return err
		} else if matched {
			matches = append(matches, path)
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	return matches, nil
}

func main() {
	var m runtime.MemStats
	DNAZeroRefFile := "Qubits_Split.DNA"
	Matches, err := WalkMatch("./", "*.fna")
	if err != nil {
		panic(err)
	}

	ZeroDNA, err := os.ReadFile(DNAZeroRefFile)
	if err != nil {
		panic(err)
	}

	ZeroDNA = ZeroDNA[22000 : len(ZeroDNA)-500000] // stripping leading and trailing leaving DNA Zero

	for _, fnaFile := range Matches {
		if len(strings.Split(fnaFile, ",")) > 1 {
			Rename := strings.Split(fnaFile, ",")[0] + strings.Split(fnaFile, ",")[1]
			os.Rename(fnaFile, Rename)
			fnaFile = Rename
		}
		fmt.Println(fnaFile)

		Launch(fnaFile, ZeroDNA)

		for {
			//time.Sleep(time.Second * 10)
			runtime.ReadMemStats(&m)
			if float64(m.Alloc) < float64(m.Sys)*0.8 { //80% RAM usage min
				break
			}
		}
	}
}
