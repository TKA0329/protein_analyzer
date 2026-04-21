# Protein Properties Analyzer
Link to Streamlit Web App: https://proteinanalyzer-tfdxptasphzhddkxctm3ed.streamlit.app/

## How it works

Put your protein sequences in a CSV file with three columns:

1. sequence
2. mutation region (optional, e.g. `3-5, 6, 8, 10`)
3. number of random copies

```csv
sequence,mutation_region,num_random_copies
MPYEKHVEQTVVEKTEQGGSGGSYRHQTEAEAEKIRRELEKQGGGGSGGGGS,3-6,5
RQGGGAPAGGNIGGGQPQGGWGQPQQPQGGNQFSGGAQSRPQ,,1
MKTQRDGHSLGRWSLVLLLLGLVMPLAIIAQVLSYKEAVL,10,8
```

Upload it to the web app and click the mutation/analyze button. The app creates random substituted copies per row, then returns a table with properties for each generated sequence: molecular weight, isoelectric point (pI), GRAVY score, aromaticity, instability index, net charge at pH 7, and secondary structure fractions.

Useful for processing large numbers of sequences at once without having to submit them one by one to tools like ProtParam.

Run locally:

```bash
streamlit run analyzer.py
```
