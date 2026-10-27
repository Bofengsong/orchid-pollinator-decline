# Run order

1. Run `code/Analysis.R` to regenerate `outputs/Results.csv` and `outputs/Results.rds`.
2. Run any `code/Figure*.R` to save that figure in `outputs/`.

Each plotting script works independently once results are available. To draw all statistical figures from the repository root:

```r
scripts <- list.files("code", pattern = "^Figure.*[.]R$", full.names = TRUE)
for (script in scripts) source(script, encoding = "UTF-8")
```

Figure 5 is the retained ecological illustration. The Chinese S1 review image is not a submission figure and is not included here.
