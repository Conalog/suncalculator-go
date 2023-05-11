package main

import (
	"fmt"
	"github.com/conalog/suncalculator-go/sun"
	"time"
)

func main() {
	lat, lng := 37.3358, 126.5839
	t := time.Now()
	times := sun.GetTimes(t, lng, lat)
	for name, timeValue := range times {
		fmt.Printf("%-20s:%s\n", name, timeValue.Format(time.RFC3339))
	}
}
