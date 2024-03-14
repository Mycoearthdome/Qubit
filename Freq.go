package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

type Database struct {
	DBFrequencies map[float64][]string
	DBSpecies     map[string][]float64
	DBCount       map[string][]int
}
type float64Slice []float64

func (s float64Slice) Len() int           { return len(s) }
func (s float64Slice) Less(i, j int) bool { return s[i] < s[j] }
func (s float64Slice) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }

type stringSlice []string

func (s stringSlice) Len() int           { return len(s) }
func (s stringSlice) Less(i, j int) bool { return s[i] < s[j] }
func (s stringSlice) Swap(i, j int)      { s[i], s[j] = s[j], s[i] }

func convertFrequenciesMapToSlice(m map[float64][]string) []struct {
	Key   float64
	Value []string
} {
	var kvSlice []struct {
		Key   float64
		Value []string
	}

	keys := make(float64Slice, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}

	sort.Sort(keys)

	for _, k := range keys {
		kvSlice = append(kvSlice, struct {
			Key   float64
			Value []string
		}{Key: k, Value: m[k]})
	}

	return kvSlice
}

func convertSpeciesMapToSlice(m map[string][]float64) []struct {
	Key   string
	Value []float64
} {
	var kvSlice []struct {
		Key   string
		Value []float64
	}

	keys := make(stringSlice, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}

	sort.Sort(keys)

	for _, k := range keys {
		kvSlice = append(kvSlice, struct {
			Key   string
			Value []float64
		}{Key: k, Value: m[k]})
	}

	return kvSlice
}

func convertCountMapToSlice(m map[string][]int) []struct {
	Key   string
	Value []int
} {
	var kvSlice []struct {
		Key   string
		Value []int
	}

	keys := make(stringSlice, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}

	sort.Sort(keys)

	for _, k := range keys {
		kvSlice = append(kvSlice, struct {
			Key   string
			Value []int
		}{Key: k, Value: m[k]})
	}

	return kvSlice
}

func convertFrequencySliceToMap(s []struct {
	Key   float64
	Value []string
}) map[float64][]string {
	m := make(map[float64][]string)

	for _, kv := range s {
		m[kv.Key] = kv.Value
	}

	return m
}

func convertSpeciesSliceToMap(s []struct {
	Key   string
	Value []float64
}) map[string][]float64 {
	m := make(map[string][]float64)

	for _, kv := range s {
		m[kv.Key] = kv.Value
	}

	return m
}

func convertCountSliceToMap(s []struct {
	Key   string
	Value []int
}) map[string][]int {
	m := make(map[string][]int)

	for _, kv := range s {
		m[kv.Key] = kv.Value
	}

	return m
}

func ReadDatabase(Databasefilename string) Database {
	var database Database
	data, err := os.ReadFile(Databasefilename)
	if err != nil {
		panic(err)
	}

	// Unmarshal slice representation
	var tempDB struct {
		DBFrequencies []struct {
			Key   float64
			Value []string
		}
		DBSpecies []struct {
			Key   string
			Value []float64
		}
		DBCount []struct {
			Key   string
			Value []int
		}
	}

	err = json.Unmarshal(data, &tempDB)
	if err != nil {
		panic(err)
	}

	// Convert slice to map
	database.DBFrequencies = convertFrequencySliceToMap(tempDB.DBFrequencies)
	database.DBSpecies = convertSpeciesSliceToMap(tempDB.DBSpecies)
	database.DBCount = convertCountSliceToMap(tempDB.DBCount)

	return database
}

func SaveDatabase(Databasefilename string, database Database) {
	// Convert map to slice representation
	tempDB := struct {
		DBFrequencies []struct {
			Key   float64
			Value []string
		}
		DBSpecies []struct {
			Key   string
			Value []float64
		}
		DBCount []struct {
			Key   string
			Value []int
		}
	}{
		DBFrequencies: convertFrequenciesMapToSlice(database.DBFrequencies),
		DBSpecies:     convertSpeciesMapToSlice(database.DBSpecies),
		DBCount:       convertCountMapToSlice(database.DBCount),
	}

	jsonData, err := json.Marshal(tempDB)
	if err != nil {
		panic(err)
	}

	file, err := os.Create(Databasefilename)
	if err != nil {
		panic(err)
	}

	defer file.Close()

	_, err = file.Write(jsonData)
	if err != nil {
		panic(err)
	}
}

func WalkMatch(root, pattern string) ([]string, error) {
	var matches []string
	err := filepath.Walk(root, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		if info.IsDir() {
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

func HurtsHuman(dbSpecies map[string][]float64, Frequency float64) bool {
	for _, HumanFrequency := range dbSpecies["Human"] {
		if Frequency == HumanFrequency {
			return true
		}
	}
	return false
}

func HuntFrequenciesExcludeHumans(dbSpecies map[string][]float64, dbFrequencies map[float64][]string) map[float64][]string {
	dbCleared := make(map[float64][]string)

	for Specie, Frequencies := range dbSpecies {
		for _, Frequency := range Frequencies {
			if !HurtsHuman(dbSpecies, Frequency) {
				dbCleared[Frequency] = append(dbCleared[Frequency], Specie)
			}
		}
	}
	return dbCleared
}

func CheckNIH(Specie string) string {
	url := "https://www.ncbi.nlm.nih.gov/nuccore/" + Specie + "?report=fasta&log$=seqview&format=text"
	resp, err := http.Get(url)
	if err != nil {
		panic(err)
	}
	body, err := io.ReadAll(resp.Body)
	if err != nil {
		panic(err)
	}

	FrontMatter := strings.Split(string(body), "</title>")[0]

	return strings.Split(FrontMatter, "<title>")[1]
}

func main() {

	var DatabaseFilename = "Frequency.db"

	var inputFiles []string
	dbFrequencies := make(map[float64][]string) //the database of frequenciess to species
	dbSpecies := make(map[string][]float64)     //the database of species to frequencies
	dbCount := make(map[string][]int)           // the parallel database of the count per frequencies per species.
	var w string
	//var GigaHertzLOW = 1000000000.0
	//var GigahertzHIGH = 999999999999.0
	//var XRayLow = 3 * math.Pow10(16) //nothing found in X-RAY band
	//var XRayHigh = 3 * math.Pow10(19) //nothing found in X-RAY band
	//var UserInput string
	//var SpeciesInput bool = false
	//var FrequencyInput bool = false
	//var Exists bool = false
	//var Frequencies []float64
	//var Species []string
	var Count int

	flag.StringVar(&w, "w", "", "Write database")
	flag.Parse()

	if len(os.Args) < 1 {
		fmt.Println("Usage: Freq [OR] Freq -w <Zero file> or <'database'>[use of wildcards is accepted]")
		return
	}

	if w != "" && w != "database" {

		for i, arg := range os.Args {
			if i >= 2 {
				inputFiles = append(inputFiles, arg)
			}
		}

		for _, file := range inputFiles {

			data, err := os.ReadFile(file)
			if err != nil {
				panic(err)
			}

			Lines := strings.Split(string(data), "\n")

			dbID := strings.Split(file, "_Freq")[0]

			for _, line := range Lines {
				tempString := strings.Split(line, "=")
				if tempString[len(tempString)-1] != "" {
					Count, err = strconv.Atoi(tempString[len(tempString)-1])
					if err != nil {
						panic(err)
					}

					Frequency, err := strconv.ParseFloat((strings.Split(strings.Split(line, "=")[1], " ")[0]), 64)
					if err != nil {
						panic(err)
					}
					dbFrequencies[Frequency] = append(dbFrequencies[Frequency], dbID)
					dbSpecies[dbID] = append(dbSpecies[dbID], Frequency)
					dbCount[dbID] = append(dbCount[dbID], Count)
				}
			}
		}

		var database Database

		database.DBFrequencies = dbFrequencies
		database.DBSpecies = dbSpecies
		database.DBCount = dbCount

		SaveDatabase(DatabaseFilename, database)

	} else if w == "database" {
		inputFiles, err := WalkMatch("./", "*.Zero")
		if err != nil {
			panic(err)
		}

		for _, file := range inputFiles {

			data, err := os.ReadFile(file)
			if err != nil {
				panic(err)
			}

			Lines := strings.Split(string(data), "\n")

			dbID := strings.Split(file, "_Freq")[0]

			for _, line := range Lines {
				tempString := strings.Split(line, "=")
				if tempString[len(tempString)-1] != "" {
					Count, err = strconv.Atoi(tempString[len(tempString)-1])
					if err != nil {
						panic(err)
					}

					Frequency, err := strconv.ParseFloat((strings.Split(strings.Split(line, "=")[1], " ")[0]), 64)
					if err != nil {
						panic(err)
					}
					dbFrequencies[Frequency] = append(dbFrequencies[Frequency], dbID)
					dbSpecies[dbID] = append(dbSpecies[dbID], Frequency)
					dbCount[dbID] = append(dbCount[dbID], Count)
				}
			}
		}

		var database Database

		database.DBFrequencies = dbFrequencies
		database.DBSpecies = dbSpecies
		database.DBCount = dbCount

		SaveDatabase(DatabaseFilename, database)

	} else {
		database := ReadDatabase(DatabaseFilename)
		dbFrequencies = database.DBFrequencies
		dbSpecies = database.DBSpecies
		//dbCount = database.DBCount
	}

	dbCleared := HuntFrequenciesExcludeHumans(dbSpecies, dbFrequencies)

	f, err := os.Create("Frequencies.Cleared")
	if err != nil {
		panic(err)
	}
	for Frequency, Species := range dbCleared {
		//if Frequency >= GigaHertzLOW && Frequency <= GigahertzHIGH {
		writeme := fmt.Sprintf("Frequency=%f [%d]Species=", Frequency, len(Species))
		for _, Specie := range Species {
			if Specie != "Monkey" && Specie != "Mouse" {
				writeme += CheckNIH(Specie) + ","
			} else {
				writeme += Specie + ","
			}
		}
		writeme = writeme[:len(writeme)-1] + "\n\n"
		f.Write([]byte(writeme))
		//}
	}
	f.Close()
	fmt.Println("Saved...Frequencies.Cleared.")

	/*
		for {
			fmt.Print("[CTRL+C to EXIT] --> Enter the Frequency and/or species you want ot enquire about: ")
			_, err := fmt.Scan(&UserInput)
			if err != nil {
				panic(err)
			}

			Frequency, err := strconv.ParseFloat(UserInput, 64)
			if err != nil {
				SpeciesInput = true
				Frequencies, Exists = dbSpecies[UserInput]
			} else {
				FrequencyInput = true
				Species, Exists = dbFrequencies[Frequency]
			}

			if Exists {
				if SpeciesInput {
					for i, Frequency := range Frequencies {
						fmt.Printf("Freqency=%f Count=%d\n", Frequency, dbCount[UserInput][i])
					}
				}
				if FrequencyInput {
					i := 0
					for _, Specie := range Species {
						fmt.Printf("Specie #%d: %s\n", i, Specie)
						i++
					}
				}
			} else {
				//Dump the whole database for the user.
				for Frequency, Species := range dbFrequencies {
					fmt.Printf("Frequency=%f:", Frequency)
					for _, Specie := range Species {
						fmt.Printf("%s, ", Specie)
					}
					fmt.Printf("\n")
				}
			}
		}
	*/
}
