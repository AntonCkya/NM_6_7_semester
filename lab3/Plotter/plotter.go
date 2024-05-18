package plotter

import (
	"os"

	chart "github.com/wcharczuk/go-chart/v2"
)

func Plot2(x1, y1, x2, y2 []float64, filename string, plots []string) {
	graph := chart.Chart{
		Title:  filename,
		Width:  1280,
		Height: 720,
		Series: []chart.Series{
			chart.ContinuousSeries{
				Name:    plots[0],
				XValues: x1,
				YValues: y1,
			},
			chart.ContinuousSeries{
				Name:    plots[1],
				XValues: x2,
				YValues: y2,
			},
		},
	}

	graph.Elements = []chart.Renderable{
		chart.Legend(&graph),
	}

	f, _ := os.Create(filename + ".png")
	defer f.Close()

	graph.Render(chart.PNG, f)
}

func Plot1D(x1, y1, dx, dy []float64, filename string, plots []string) {
	graph := chart.Chart{
		Title:  filename,
		Width:  1280,
		Height: 720,
		Series: []chart.Series{
			chart.ContinuousSeries{
				Name:    plots[0],
				XValues: x1,
				YValues: y1,
			},
			chart.ContinuousSeries{
				Style: chart.Style{
					StrokeWidth: chart.Disabled,
					DotWidth:    5,
				},
				XValues: dx,
				YValues: dy,
			},
		},
	}

	graph.Elements = []chart.Renderable{
		chart.Legend(&graph),
	}

	f, _ := os.Create(filename + ".png")
	defer f.Close()

	graph.Render(chart.PNG, f)
}

func Plot2D(x1, y1, x2, y2, dx, dy []float64, filename string, plots []string) {
	graph := chart.Chart{
		Title:  filename,
		Width:  1280,
		Height: 720,
		Series: []chart.Series{
			chart.ContinuousSeries{
				Name:    plots[0],
				XValues: x1,
				YValues: y1,
			},
			chart.ContinuousSeries{
				Name:    plots[1],
				XValues: x2,
				YValues: y2,
			},
			chart.ContinuousSeries{
				Name:    plots[2],
				XValues: dx,
				YValues: dy,
			},
			chart.ContinuousSeries{
				Style: chart.Style{
					StrokeWidth: chart.Disabled,
					DotWidth:    5,
				},
				XValues: dx,
				YValues: dy,
			},
		},
	}

	graph.Elements = []chart.Renderable{
		chart.Legend(&graph),
	}

	f, _ := os.Create(filename + ".png")
	defer f.Close()

	graph.Render(chart.PNG, f)
}
