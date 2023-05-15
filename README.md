# Sun Calculator for Go

[Sun Calculator for Python](https://github.com/Conalog/sun_calculator)의 Go 버전입니다.

## Example

```go
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
```
