library(EARNMC)

# Generate the West Africa incidence summary graphs.
generateWestAfricaRawDataPlotPDF()

# Generate the West Africa EA-R comparison plots. 
makeR0ComparisonPlotPDFs()

# Generate the latest model prediction plot (at the time of manuscript 
# preparation).
makeWestAfricaPredictionPlotPDFs()

# Generate the summary of the spline basis functions used in estimation and
# prediction for West Africa
makeSplineBasisPlotPDF()

# Generate the data summary and EA-R comparison graphs for the Kikwit analyses
makeKikwitGraphPDFs()

